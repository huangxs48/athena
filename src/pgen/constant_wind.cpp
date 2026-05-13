//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headersc
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../units/units.hpp"

//general variables to record mesh size and cell number
static int mesh_nx1, mesh_nx2, mesh_nx3;
static Real mesh_x1min, mesh_x1max, x1ratio;
static Real mesh_x2min, mesh_x2max;
static Real mesh_x3min, mesh_x3max;

//general variables converts c.g.s unit and unit-less variables in code
static Real kappa_es;
static Real temp_unit, l_unit, rho_unit, kappa_unit, vel_unit, time_unit;
static Real tfloor; //temperature floor used in radiation class
static Real dfloor, pfloor; //density, pressure floor used in hydro class

//initial background density and pressure
static Real rho_init, press_init;
static Real boundary_temp_lim; //optional, temperature uplimit at boundary

//prescribed wind base density, density profile index, total mdot
static Real rho_wind_base, rho_wind_index, mdot_wind, r_wind_in, vel_wind_base;
static Real lum_trapping; //assumed luminosity at trapping radius
static Real lum_base; //luminosity corresponding to the flux at boundary
static Real t_lum_base_ramp; //timescale for flux and mdot to ramp up

//opacity function
std::string  opacity_file;
std::string  opacity_type;
static int n_tem;
static int n_rho;

//frequency dependent free-free
Real kappa_ff_nu(Real nu, Real temp, Real rho);
void Multi_FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

//new combined opacity table
static AthenaArray<Real> combine_temp_grid;
static AthenaArray<Real> combine_rho_grid;
static AthenaArray<Real> combine_ross_table;
static AthenaArray<Real> combine_planck_table;
void combineopacity(const Real rho, const Real tgas, Real &kappa_ross, Real &kappa_planck);
void GetCombineOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

//the frequency grid
static AthenaArray<Real> fre_grid;

//frequency integrated opacity
Real kappa_ff_planck(Real temp, Real rho);
Real kappa_ff_ross(Real temp, Real rho);
void FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

// User-defined boundary conditions for hydro and radiation
void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void ConstMdotInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
		      Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void ConstFluxInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
		      const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
		      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//AMR condition
int RefinementCondition(MeshBlock *pmb);

//following are for reading in an initial csm profile
//buffer arrays for coordinate
std::vector<float> x1coord;
std::vector<float> x2coord;
std::vector<float> x3coord;

//recoding input density, temperature, velocity
static AthenaArray<Real> rho_init_buff;
static AthenaArray<Real> temp_init_buff;
static AthenaArray<Real> vel_init_buff;

//simple search for index of variable al in an array vec
int getindex(std::vector<float> vec, float val){
  std::vector<float>::iterator it = std::find(vec.begin(), vec.end(), val);
  int index = std::distance(vec.begin(), it);
  return index;
}


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================


void Mesh::InitUserMeshData(ParameterInput *pin) {
  int blocksizex1 = pin->GetOrAddInteger("meshblock", "nx1", 1);
  int blocksizex2 = pin->GetOrAddInteger("meshblock", "nx2", 1);
  int blocksizex3 = pin->GetOrAddInteger("meshblock", "nx3", 1);

  kappa_es = pin->GetReal("problem", "kappa_es");

  temp_unit = pin->GetReal("problem", "temp_unit");
  l_unit = pin->GetReal("problem", "l_unit");
  rho_unit = pin->GetReal("problem", "rho_unit");
  //kappa_unit: cm^2/g
  kappa_unit = 1.0/(rho_unit*l_unit);
  vel_unit = Constants::speed_of_light_cgs/pin->GetReal("radiation", "crat");

  tfloor = pin->GetOrAddReal("radiation", "tfloor", 1.0e-8);
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 1.0e-8);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 1.0e-8);

  rho_init = pin->GetOrAddReal("hydro", "rho_init", 1.0e-8);
  press_init = pin->GetOrAddReal("hydro", "press_init", 1.0e-8);

  boundary_temp_lim = pin->GetOrAddReal("problem", "boundary_temp_lim", 1.0e5);
  //mass_bottom_cgs = pin->GetReal("problem", "mass_bottom_cgs");

  //prescribed wind parameters
  rho_wind_base = pin->GetOrAddReal("problem", "rho_wind_base", 1.0e-1);
  rho_wind_index = pin->GetOrAddReal("problem", "rho_wind_index", -2.0);
  //r_wind_in = pin->GetOrAddReal("problem", "r_wind_in", 0.1);
  vel_wind_base = pin->GetOrAddReal("problem", "vel_wind_base", 1.0);
  mdot_wind = pin->GetOrAddReal("problem", "mdot_wind", 1.0);
  lum_trapping = pin->GetOrAddReal("problem", "lum_trapping", 0.0);
  lum_base = pin->GetOrAddReal("problem", "lum_base", 0.0);
  t_lum_base_ramp = pin->GetOrAddReal("problem", "t_lum_base_ramp", 1.0);

  //opacity
  opacity_file = pin->GetOrAddString("problem", "opacity_file", "None");
  opacity_type = pin->GetOrAddString("problem", "opacity_type", "table");
  n_rho = pin->GetOrAddInteger("problem" ,"n_rho", 70);
  n_tem = pin->GetOrAddInteger("problem" ,"n_tem", 140);

  // Enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, ConstMdotInnerX1); //HydroInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, HydroOuterX1);
  }

  // Enroll AMR condiation
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);


  if (NR_RADIATION_ENABLED){
    //Enroll rad boundaries
    if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, ConstFluxInnerX1); //RadInnerX1);
    }
    if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, RadOuterX1);
    }
  }

  // //do not use this "self-gravity" for now
  // EnrollUserExplicitSourceFunction(EnvGravity);
  
  // AllocateUserHistoryOutput(2);
  // EnrollUserHistoryOutput(0, einj_int, "einj_int");//
  // EnrollUserHistoryOutput(1, einj_dt, "einj_dt");//
  // EnrollUserHistoryOutput(3, massflux_Inj_x3, "massflux_Inj_x3");//mass flux outer
  // EnrollUserHistoryOutput(4, massfluxix1, "massfluxix1");//mass flux inner boundary
  // EnrollUserHistoryOutput(5, massfluxox1, "massfluxox1");//mass flux outer

  // read in mesh size and number of cells to prepare reading initial csm profile and graivty source term
  mesh_nx1 = pin->GetInteger("mesh", "nx1");
  mesh_x1min = pin->GetReal("mesh", "x1min");
  mesh_x1max = pin->GetReal("mesh", "x1max");
  mesh_nx2 = pin->GetInteger("mesh", "nx2");
  mesh_x2min = pin->GetReal("mesh", "x2min");
  mesh_x2max = pin->GetReal("mesh", "x2max");
  mesh_nx3 = pin->GetInteger("mesh", "nx3");
  mesh_x3min = pin->GetReal("mesh", "x3min");
  mesh_x3max = pin->GetReal("mesh", "x3max");
  x1ratio = pin->GetReal("mesh", "x1rat");

  //prepare three vectors for index finding of x1 x2 x3 coordinates

  //Real dx1 = (mesh_x1max - mesh_x1min)/mesh_nx1;
  Real dx2 = (mesh_x2max - mesh_x2min)/mesh_nx2;
  Real dx3 = (mesh_x3max - mesh_x3min)/mesh_nx3;
  //the vector are equivalent to pcoord->x1v, x2v, x3v
  // XS: NOTE here assumed logarithmic r_grid
  for(int i=0; i<mesh_nx1; i++){
    Real x1coord_now = (pow(x1ratio, i)-1.0)/(pow(x1ratio, mesh_nx1)-1.0) *
                       (mesh_x1max - mesh_x1min) + mesh_x1min;
    x1coord.push_back(x1coord_now);
  }
  for(int j=0; j<mesh_nx2; j++){
    x2coord.push_back(mesh_x2min+j*dx2);
  }
  for(int k=0; k<mesh_nx3; k++){
    x3coord.push_back(mesh_x3min+k*dx3);
  }

  //read in initial density, temperature, velocity
  AllocateRealUserMeshDataField(3);
  ruser_mesh_data[0].NewAthenaArray(mesh_nx1);
  ruser_mesh_data[1].NewAthenaArray(mesh_nx1);
  ruser_mesh_data[2].NewAthenaArray(mesh_nx1);
  //keep record of enclosed mass in each radius
  ruser_mesh_data[3].NewAthenaArray(mesh_nx1); //x1coordinate
  // ruser_mesh_data[4].NewAthenaArray(mesh_nx1); //mass in each shell
  // ruser_mesh_data[5].NewAthenaArray(mesh_nx1); //mass coordinate of each shell
  // ruser_mesh_data[6].NewAthenaArray(mesh_nx1); //enclosed mass in each radius
  // ruser_mesh_data[7].NewAthenaArray(mesh_nx1); //summbed b coefficient, not used,
  // ruser_mesh_data[8].NewAthenaArray(mesh_nx1); //added energy, not used

  for(int i=0; i<mesh_nx1; i++){
    ruser_mesh_data[3](i) = x1coord[i];
  }

  if (NR_RADIATION_ENABLED){

    if (opacity_file!="None"){
      //load combined opacity
      combine_temp_grid.NewAthenaArray(n_tem);
      combine_rho_grid.NewAthenaArray(n_rho);
      combine_ross_table.NewAthenaArray(n_tem, n_rho);
      combine_planck_table.NewAthenaArray(n_tem, n_rho);
      
      FILE *f_combineopacity;
      
      if ( (f_combineopacity=fopen(opacity_file.c_str(),"r"))==NULL )
	{
	  //printf("Open input file error combined opacity table %s\n", opacity_file.c_str());
	  //return;
	  std::stringstream msg;
	  msg << "FATAL ERROR: Could not open opacity file "<< opacity_file.c_str() << std::endl;
	  ATHENA_ERROR(msg);
	}
      
      //first two lines are n_temp, n_rho
      int buff;
      for(int i=0; i<2; i++){
	fscanf(f_combineopacity,"%d",&(buff));
      }
      
      //load temperature grid
      for(int i=0; i<n_tem; ++i){
	fscanf(f_combineopacity, "%lf", &(combine_temp_grid(i)));
      }
      
      //load density grid
      for(int i=0; i<n_rho; ++i){
	fscanf(f_combineopacity, "%lf", &(combine_rho_grid(i)));
      }
      
      //load grey Rosseland mean opacity
      for (int j=0; j<n_tem; ++j){
	for (int i=0; i<n_rho; ++i){
	  fscanf(f_combineopacity, "%lf", &(combine_ross_table(j, i)));
	}
      }
      
      //load grey Planck mean opacity
      for (int j=0; j<n_tem; ++j){
	for (int i=0; i<n_rho; ++i){
	  fscanf(f_combineopacity, "%lf", &(combine_planck_table(j, i)));
	}
      }
      
      fclose(f_combineopacity);
      
      // printf("testing opacity table\n");
      // Real kappa_p_test, kappa_r_test;
      // combineopacity(1.0e-15, 1.0e6, kappa_r_test, kappa_p_test);
      // printf("interpolated kappa_p:%g, kappa_r:%g\n", kappa_p_test, kappa_r_test);
    }

  }
    return;
}

//initialize user mesh block data
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  int blocksizex1 = pin->GetOrAddInteger("meshblock", "nx1", 1);
  int blocksizex2 = pin->GetOrAddInteger("meshblock", "nx2", 1);
  int blocksizex3 = pin->GetOrAddInteger("meshblock", "nx3", 1);

  blocksizex1 += 2*(NGHOST);
  if (blocksizex2 >1) blocksizex2 += 2*(NGHOST);
  if (blocksizex3 >1) blocksizex3 += 2*(NGHOST);
  
  AllocateRealUserMeshBlockDataField(4); //pre-allocate for diagnostic reason, probably ok to skip for now
  ruser_meshblock_data[0].NewAthenaArray(4,blocksizex3, blocksizex2, blocksizex1); //
  ruser_meshblock_data[1].NewAthenaArray(4,blocksizex3, blocksizex2, blocksizex1); //
  ruser_meshblock_data[2].NewAthenaArray(5,blocksizex3, blocksizex2, blocksizex1); //
  ruser_meshblock_data[3].NewAthenaArray(3,blocksizex3, blocksizex2, blocksizex1); //

  AllocateIntUserMeshBlockDataField(1);
  iuser_meshblock_data[0].NewAthenaArray(blocksizex1);//store the index of global coordinate array
  //iuser_meshblock_data[1].NewAthenaArray(blocksizex1);//flag to mass injection

  //enroll opacity function here
  if (NR_RADIATION_ENABLED){
    if (pnrrad->nfreq>1){
      pnrrad->EnrollOpacityFunction(Multi_FreeFreeOpacity);
    }else{
      if (opacity_type=="table"){
	pnrrad->EnrollOpacityFunction(GetCombineOpacity);//(FreeFreeOpacity);
      }else if (opacity_type=="freefree"){
	pnrrad->EnrollOpacityFunction(FreeFreeOpacity);
      }
      if (Globals::my_rank==0){
	std::cout<<"Using Opacity Type:"<<opacity_type<<std::endl;
      }
    }
  }

  //all for diagnostic, probably ok to skip 
  // AllocateUserOutputVariables(8);
  // SetUserOutputVariableName(0, "cellvol");
  // SetUserOutputVariableName(1, "dmass_r");
  // SetUserOutputVariableName(2, "mass_enclose_r");
  // SetUserOutputVariableName(3, "r_index");
  // SetUserOutputVariableName(4, "inject_flag");
  // SetUserOutputVariableName(5, "mass_coord_r");
  // SetUserOutputVariableName(6, "int_coefb");
  // SetUserOutputVariableName(7, "einj");
  return;
}

//an example refinement condition, find the density gradient maximum (simple shock finder)
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  int k=pmb->ks;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1)))/w(IDN,k,j,i);
      Real epsp= (std::abs(w(IPR,k,j,i+1)-2.0*w(IPR,k,j,i)+w(IPR,k,j,i-1)))/w(IPR,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  if (maxeps>1.0){
    printf("my_rank:%d, gid:%d, maxeps:%g\n", Globals::my_rank, pmb->gid, maxeps);
  }
  if(maxeps > 1.0) return 1;
  if(maxeps < 0.1) return -1;
  return 0;
}

// this block is mainly for calculating enclosed gravity, probably ok to skip it if not using EnvGravity, 
void Mesh::UserWorkInLoop(){

  
  MeshBlock *pmb = my_blocks(0);
  for(int nb=0; nb<nblocal; ++nb){ //loop over meshblocks on the same core
    pmb = my_blocks(nb);
    
    Hydro *phydro = pmb->phydro;
    Coordinates *pcoord = pmb->pcoord;
    int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;

    //NOTE only works for one-dimension problem right now
    Real mass_coord_now = 0.0;
    for(int i=is; i<=ie; i++){     
      
      Real r_now = pcoord->x1f(i);
      int index_rnow = pmb->iuser_meshblock_data[0](i);

      //then update mass in each shell
      Real dmass = 0.0; // mass_bottom_cgs/mass_unit;
        for(int k=ks; k<=ke; k++){
	  for(int j=js; j<=je; j++){
	    dmass += phydro->u(IDN,k,j,i) * pcoord->GetCellVolume(k,j,i);

	    //sanity check
	    for (int n=0; n<(NHYDRO);n++){
	      if (phydro->u(n,k,j,i) != phydro->u(n,k,j,i)){
		printf("block: %d, n: %d ,k: %d,j: %d,i: %d\n", pmb->gid,n,k,j,i);
		printf("x1v: %g, x2v:%g, x3v:%g\n",pmb->pcoord->x1v(i), pmb->pcoord->x2v(j),pmb->pcoord->x3v(k));
		//abort();
	      }   
	    }//end NHYDRO
	    
	  }//j
	}//k

    }//i
 
    
  }//loop over meshblocks
  
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){
  
  // for(int k=ks; k<=ke; k++){
  //   for(int j=js; j<=je; j++){
  //     for(int i=is; i<=ie; i++){
	
  // 	Real r_now = pcoord->x1f(i);
  // 	int index_rnow = iuser_meshblock_data[0](i);

  // 	//try calculate mass here
  // 	Real mass_rad = 0.0;
  // 	for (int ii=0; ii<index_rnow; ii++){
  // 	  mass_rad += mass_shell(ii);//pmy_mesh->ruser_mesh_data[4](ii);
  // 	}
  // 	//printf("index_now:%d, i:%d, gid:%d, my_rank:%d, mass_rad:%g\n", index_rnow, i, gid, Globals::my_rank, mass_rad);
      
  // 	// user_out_var(0,k,j,i) = pcoord->GetCellVolume(k,j,i);
  // 	// user_out_var(1,k,j,i) = pmy_mesh->ruser_mesh_data[4](index_rnow);
  // 	// //printf("output dmass:%g\n", pmy_mesh->ruser_mesh_data[4](index_rnow));s
  // 	// user_out_var(2,k,j,i) = pmy_mesh->ruser_mesh_data[6](index_rnow); //mass_rad;
  // 	// user_out_var(3,k,j,i) = index_rnow;
  // 	// user_out_var(4,k,j,i) = iuser_meshblock_data[1](i); //injection flag
  // 	// user_out_var(5,k,j,i) = pmy_mesh->ruser_mesh_data[5](index_rnow); //mass_rad;
  // 	// user_out_var(6,k,j,i) = pmy_mesh->ruser_mesh_data[7](index_rnow); //integrated coefficient b
  // 	// user_out_var(7,k,j,i) = ruser_meshblock_data[0](0,k,j,i);//injected energy
  
  //     }
  //    }
  //  }
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes shock profile as read-in data
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real gamma_gas = peos->GetGamma();

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

	// //load data, find current r index
        // Real x_now = pcoord->x1f(i);
        // int index_xnow = getindex(x1coord, x_now);
	// iuser_meshblock_data[0](i) = index_xnow;

	// Real rho_now = pmy_mesh->ruser_mesh_data[0](index_xnow);
	// Real temp_now = pmy_mesh->ruser_mesh_data[1](index_xnow);
	// Real vel_now = pmy_mesh->ruser_mesh_data[2](index_xnow);
        //printf("x_now:%g, index_x:%d, rho:%g, temp:%g, vel:%g\n", x_now, index_xnow, rho_now, temp_now, vel_now);

	// //find current radius, density and velocity
	// Real r_now = pcoord->x1v(i);
	// Real rho_now = rho_wind_base * pow(r_now/r_wind_in, rho_wind_index);
	// Real vel_now = mdot_wind / rho_now / (4.0*PI*r_now*r_now);
	// //apply floor
	// rho_now = std::max(rho_now, dfloor);
	// if (NR_RADIATION_ENABLED){
	//   vel_now = std::min(vel_now, 0.9*pnrrad->crat); //hard-coded for now
	// }
	// //set temperature by assuming a luminosity at trapping radius

	// Real kappa_es_code = kappa_es * rho_unit * l_unit; 
	// Real mass_load_wind = mdot_wind * (vel_now);
	// Real tgas4 =  kappa_es_code * lum_trapping / mass_load_wind / pow(r_now, 3) / pow(4.0*PI, 2) / (pnrrad->crat*pnrrad->prat);
	// Real temp_now = pow(tgas4, 0.25) ;
	// //printf("r_now:%g, mdot_wind:%g, mass_load:%g, vel_now:%g, tgas4:%g\n", r_now, mdot_wind, mass_load_wind, vel_now, tgas4);

	//assuming a low density background gas
	Real rho_now = rho_init;
	Real temp_now = press_init;

        phydro->u(IDN,k,j,i) = rho_now;
        phydro->u(IM1,k,j,i) = 0.0; //vel_now;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
	
	if (NR_RADIATION_ENABLED){

	  Real rho = rho_now;//phydro->w(IDN,k,j,i);
	  Real temp = temp_now;//phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);

	  Real rho_cgs = rho*rho_unit;
	  Real temp_cgs = temp*temp_unit;

	  Real kappa_s, kappa_ross, kappa_planck;

	  // electron scattering opacity, hard coded for now
	  kappa_s = 0.2 * (1.0 + 0.6);
	  Real T_ion = 1.0e4;//ionization temperature, below which assuming kappa_scatter=0
	  Real T_dust = 4.0e3;//where roughly opacity rises again due to dust, what to do for this?

	  if (pnrrad->nfreq>1){
	    for (int ifr=0; ifr<pnrrad->nfreq; ++ifr){

	      //get current frequency
	      //hard coded for now
	      Real evtohz = 2.41838e14;
	      Real nu_kev = 1.0; //fre_grid(ifr); //make this frequency grid 
	      Real nu_hz = nu_kev*1000*evtohz;
	      
	      Real kappa_ff_cgs = kappa_ff_nu(nu_hz, temp, rho);
	      //set rosseland mean and planck mean to be same for now, can be an issue
	      kappa_ross = kappa_ff_cgs;
	      kappa_planck = kappa_ff_cgs;
	      
	      pnrrad->sigma_s(k,j,i,ifr) = kappa_s * rho * rho_unit * l_unit; 
	      pnrrad->sigma_a(k,j,i,ifr) = kappa_ff_cgs * rho * rho_unit *l_unit; 
	      pnrrad->sigma_pe(k,j,i,ifr) = kappa_ff_cgs * rho * rho_unit *l_unit;
	      pnrrad->sigma_p(k,j,i,ifr) = kappa_ff_cgs * rho * rho_unit *l_unit;
	    }
	  }else{
	    
	    kappa_ross = kappa_ff_ross(temp, rho);
	    kappa_planck = kappa_ff_planck(temp, rho);
	    //one frequency, grey rhd for now
	    pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	    pnrrad->sigma_a(k,j,i,0) = kappa_ross * rho * rho_unit * l_unit; //rosseland mean
	    pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
	    pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean

	  }

	  //initialize intensity
	  for (int ifr=0; ifr<pnrrad->nfreq; ++ifr){
	    for(int n=0; n<pnrrad->nang; ++n){
	      int ang=ifr*pnrrad->nang+n;
	      pnrrad->ir(k,j,i,ang) = pow(temp, 4);//use temp_now^4 if assuming initial trad=tgas
	    }
	  }
     
	}//end rad

	if (NON_BAROTROPIC_EOS) {
	  phydro->u(IEN,k,j,i) = rho_now * temp_now /(gamma_gas - 1.0);
	  phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
	}//end non barotropic
      }//end i
      
      
    }//end j
  }//end k

  return;
}


void HydroOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=ks; k<=ke; ++k) {//phi
    Real phi_coord = pco->x3v(k);
    for (int j=js; j<=je; ++j) {//theta
      Real theta_coord = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {//R
        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
        prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
        prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
        prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
      if (NON_BAROTROPIC_EOS){
        prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie);
      }

      }//end R
    }//end theta
  }//end Phi
  
}

void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=ngh; ++i) {//R
        prim(IDN,k,j,is-i) = prim(IDN,k,j,is); 
        prim(IVX,k,j,is-i) = std::min(0.0, prim(IVX,k,j,is));//prim(IVX,k,j,is+i-1); //try only reflect velocity
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        if (NON_BAROTROPIC_EOS){
          prim(IPR,k,j,is-i) = std::max(prim(IPR,k,j,is), boundary_temp_lim/temp_unit * prim(IDN,k,j,is));
        }

      }//end R
    }//end theta
  }//end Phi

}


void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
  // copy radiation variables into ghost zones,
  // only allow outflow

  int &nang = pnrrad->nang; // angles per octant
  int &nfreq = pnrrad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=ngh; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir<0, inward: mu_dir>0
      Real mu_dir = pnrrad->mu(0,k,j,is,ang);
      if (mu_dir < 0.0){
        ir(k,j,is-i,ang) = ir(k,j,is,ang);
      }else{
        ir(k,j,is-i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}
}

void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  int &nang = pnrrad->nang; // angles per octant
  int &nfreq = pnrrad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=ngh; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir>0, inward: mu_dir<0
      Real mu_dir = pnrrad->mu(0,k,j,ie,ang);
      if (mu_dir > 0.0){
        ir(k,j,ie+i,ang) = ir(k,j,ie,ang);
      }else{
        ir(k,j,ie+i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}

  return;

}

void ConstFluxInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
		      const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
		      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
  
  //boundary condition to insert a constant flux corresponding to lum_base,
  //only for one band yet

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
  
	//try use cell-center values
	Real rho_local = w(IDN,k,j,is);//(w(IDN,k,j,is) + w(IDN,k,j,is-i))/2.0;
	Real sigma_local = 0.0;
	if (time>0.0)
	  sigma_local = pnrrad->sigma_a(k,j,is,0);
	
	Real dr = pco->dx1f(is);
	Real r_local = pco->x1f(is);

	//target flux
	Real frad_local = lum_base / (4.0*PI*r_local*r_local);
	//printf("frad_local:%g, lum_base:%g, r_local:%g\n", frad_local, lum_base, r_local);
	//add a grandually increasing factor
	if (time>0.0 && time<t_lum_base_ramp){
	    frad_local = frad_local * (1.0 - exp(-time/t_lum_base_ramp));
	    //printf("frad_now:%g\n", frad_local);
	}

	// //initialze moment array
	// if (time==0.0){
	//   pmb->pnrrad->CalculateMoment(pmb->pnrrad->ir);
	// }

	//get coefficient between Prr_rad/Erad, close to 1/3. when isotropic
	//energy density of first active cell
	Real er_is = 0.0;
	Real pr11_is = 0.0;
	//Real pr22_is = 0.0;
	//Real pr33_is = 0.0;
	for (int n=0; n<pnrrad->nang; ++n){//for single band
	  Real wmu = pnrrad->wmu(n);
	  Real mux = pnrrad->mu(0,k,j,is,n);
	  Real muy = pnrrad->mu(1,k,j,is,n);
	  Real muz = pnrrad->mu(2,k,j,is,n);
	  er_is += wmu * ir(k,j,is,n);
	  pr11_is += wmu * mux * mux * ir(k,j,is,n);
	  //pr22_is += wmu * muy * muy * ir(k,j,is,n);
	  //pr33_is += wmu * muz * muz * ir(k,j,is,n);
	}
	
        //Real er_is = pmb->pnrrad->rad_mom(IER,k,j,is);
	//Real pr11_is = pmb->pnrrad->rad_mom(IPR11,k,j,is);
	Real fedd_is = 3.0;
	
	// if (time > 0.0 and er_is >0.0 and pr11_is>0.0){
	//   fedd_is = (pr11_is/er_is); 
	//   //printf("fedd: %g, er_is/pr22_is:%g, er_is/pr33_is:%g\n", fedd_is, er_is/pr22_is, er_is/pr33_is);
	// }
	
	Real er_local = er_is + fedd_is * sigma_local * dr * frad_local;
	// if (i==1)
	//   printf("er_is=%g, pr11_is=%g, dr=%g, rho_local=%g, frad_local=%g, sigma_local=%g\n", er_is, pr11_is, dr, rho_local, frad_local, sigma_local);

	//get intensity coefficients
	Real coefa_u = 0.0, coefb_u = 0.0;
	Real coefa_d = 0.0, coefb_d = 0.0;

	for (int n=0; n<pnrrad->nang; ++n) {
	  Real mux = pnrrad->mu(0,k,j,is,n);
	  Real weight = pnrrad->wmu(n);
	  if (mux > 0.0){
	    coefa_u += weight;
	    coefb_u += mux * weight;
	  } else {
	    coefa_d += weight;
	    coefb_d += mux * weight;
	  }
	}//end angle

	for (int n=0; n<pnrrad->nang; ++n){
	  Real mux = pnrrad->mu(0,k,j,is-i,n);
	  if (mux > 0.0){
	    ir(k,j,is-i,n) = 0.5 * (er_local/coefa_u + frad_local/coefb_u);
	  }else{
	    ir(k,j,is-i,n) = 0.5 * (er_local/coefa_d + frad_local/coefb_d);
	  }
	  //printf("nang=%d, k=%d, j=%d, i=%d, ir=%g, er_local=%g, frad_local=%g, coefa_u=%g, coefb_u=%g, coefa_d=%g, coefb_d=%g\n", n, k, j, i, ir(k,j,is-i, n), er_local, frad_local, coefa_u, coefb_u, coefa_d, coefb_d);
	  // if (i==1 && (n==0||n==7))
	  //   printf("nang=%d, i=%d, ir=%g, er_local=%g, frad_local=%g,  er_is=%g, pr11_is=%g, er/pr=%g\n", n, i, ir(k,j,is-i, n), er_local, frad_local,  er_is, pr11_is, er_is/pr11_is);
	}//end angle

	//debug
	Real fr_now = 0.0;
	for (int n=0; n<pnrrad->nang; ++n){//for single band
	  Real wmu = pnrrad->wmu(n);
	  Real mux = pnrrad->mu(0,k,j,is-i,n);
	  Real muy = pnrrad->mu(1,k,j,is-i,n);
	  Real muz = pnrrad->mu(2,k,j,is-i,n);
	  fr_now += wmu * mux * ir(k,j,is-i,n);
	}
	//printf("i=%d, fr(is-i):%g\n", i, fr_now);
	
      }//i
    }//j
  }//k

}

void ConstMdotInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=ngh; ++i) {//R

	//find current radius, density and velocity
	//use is value to be consistent with flux
	Real r_now = pco->x1f(is);
	//Real rho_now = rho_wind_base;
	Real vel_now = vel_wind_base;
	Real mdot_wind_now = mdot_wind;
	Real lum_base_now = lum_base;
	if (time>0.0 && time<t_lum_base_ramp){
	   mdot_wind_now = mdot_wind * (1.0 - exp(-time/t_lum_base_ramp));
	   lum_base_now = lum_base * (1.0 - exp(-time/t_lum_base_ramp));
	}
	//Real vel_now = mdot_wind_now / rho_now / (4.0*PI*r_now*r_now);
	
	Real rho_now = mdot_wind_now / vel_now / (4.0*PI*r_now*r_now);
	//printf("mdot_now:%g, rho_now:%g, vel_now:%g, r_now:%g\n", mdot_wind_now, rho_now, vel_now, r_now);

	//estimate gas temperature
	Real mass_load_wind = mdot_wind_now / vel_now;
	Real temp_now = press_init/rho_init;
	
	if (NR_RADIATION_ENABLED){
	  Real kappa_es_code = kappa_es * rho_unit * l_unit; 
	  Real tgas4 =  kappa_es_code * lum_base / mass_load_wind / pow(r_now, 3) / pow(4.0*PI, 2) / (pmb->pnrrad->crat*pmb->pnrrad->prat);
	  temp_now = pow(tgas4, 0.25) ;
	}
	
        prim(IDN,k,j,is-i) = rho_now; 
        prim(IVX,k,j,is-i) = vel_now; //prim(IVX,k,j,is);
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);

	if (NON_BAROTROPIC_EOS){
          prim(IPR,k,j,is-i) =  prim(IPR,k,j,is);//std::max(prim(IPR,k,j,is), boundary_temp_lim/temp_unit * prim(IDN,k,j,is)); //std::max(rho_now*temp_now, boundary_temp_lim/temp_unit * rho_now);
        }

      }//end R
    }//end theta
  }//end Phi

}

//input code unit, output cgs
Real kappa_ff_nu(Real nu, Real temp, Real rho){

  Real h_planck = 6.626196e-27 ;
  Real evtohz = 2.41838e14;
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real m_p = 1.6726e-24;
  Real k_B = 1.3807e-16;
  
  Real  gff = 1.0;
  Real  z = 1.0;

  Real  he_adbund = 0.04;
  Real  nh = rho_cgs/m_p/(1.0 + 4.0*he_adbund);
  Real  nhe = nh*he_adbund;
  Real  ne = nh + 2.0*nhe;
  Real  n_rho = rho_cgs/m_p/0.62;

  Real  e_ff = 3.7e8 * pow(temp_cgs, -0.5) * pow(z, 2) * pow(n_rho, 2) * pow(nu, -3) * (1.0 - exp(-h_planck*nu/k_B/temp_cgs)) * gff;

  return e_ff/rho_cgs;
}

void Multi_FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  NRRadiation *prad = pmb->pnrrad;
  int il = pmb->is; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie; int ju = pmb->je; int ku = pmb->ke;
  il -= NGHOST;
  iu += NGHOST;
  if(ju > jl){
    jl -= NGHOST;
    ju += NGHOST;
  }
  if(ku > kl){
    kl -= NGHOST;
    ku += NGHOST;
  }

  // electron scattering opacity
  Real kappas = 0.2 * (1.0 + 0.6);
  Real kappaa = 0.0;
  Real T_ion = 1.0e4;//ionization temperature, below which assuming kappa_scatter=0
  Real T_llim = 1.0e4;//lower lim of TOPs data, temperature at which switch Combined opacity grey opacity including dust
  Real T_dust = 4.0e3;//where roughly opacity rises again due to dust, what to do for this?
  
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
  for (int ifr=0; ifr<prad->nfreq; ++ifr){
    Real rho  = prim(IDN,k,j,i);
    Real tgas = std::max(prim(IEN,k,j,i)/rho, tfloor);

    Real kappa_ff_cgs;

    Real rho_cgs = rho * rho_unit;
    Real tgas_cgs = tgas * temp_unit;

    //hard coded for now
    Real evtohz = 2.41838e14;
    Real nu_kev = 1.0; //fre_grid(ifr); //make this frequency grid 
    Real nu_hz = nu_kev*1000*evtohz;

    kappa_ff_cgs = kappa_ff_nu(nu_hz, tgas, rho);

    prad->sigma_s(k,j,i,ifr) = kappa_es * rho * rho_unit * l_unit; 
    //assuming planck mean and rossland mean are same, make change to adapt your problem
    prad->sigma_a(k,j,i,ifr) = kappa_ff_cgs * rho * rho_unit *l_unit; 
    prad->sigma_pe(k,j,i,ifr) = kappa_ff_cgs * rho * rho_unit *l_unit;
    prad->sigma_p(k,j,i,ifr) = kappa_ff_cgs * rho * rho_unit *l_unit;
  }    

 }}}

}

void FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim){

  NRRadiation *pnrrad=pmb->pnrrad;
  int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
  //int kl=pmb->kl, ku=pmb->ku, js=pmb->jl, je=pmb->ju, is=pmb->il, ie=pmb->iu;
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl=js, ju=je;
  int kl=ks, ku=ke;
  if (pmb->pmy_mesh->f2){
    jl = js - NGHOST;
    ju = je + NGHOST;
  }
  if (pmb->pmy_mesh->f3){
    kl = ks - NGHOST;
    ku = ke + NGHOST;
  }

  //printf("kappaes:%g\n", kappa_es);
  for (int k=kl; k<=ku; k++){
    for (int j=jl; j<=ju; j++){
      for (int i=il; i<=iu; i++){
	
	Real rho = prim(IDN,k,j,i);
	Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i); //std::max(prim(IPR,k,j,i)/prim(IDN,k,j,i), tfloor);
	Real rho_cgs = rho*rho_unit;
	Real temp_cgs = temp*temp_unit;
	
	Real kappa_s, kappa_ross, kappa_planck;
	kappa_s = kappa_es;
	kappa_ross = kappa_ff_ross(temp, rho);
	kappa_planck = kappa_ff_planck(temp, rho);
	
	//one frequency
	pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	pnrrad->sigma_a(k,j,i,0) = kappa_ross * rho * rho_unit * l_unit; //rosseland mean
	pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
        pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean
      
      }//end i
    }//end j
  }//end k

}

 //planck mean free free absorption
Real kappa_ff_planck(Real temp, Real rho){
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real kappa_cgs = 2.86e-5*(rho_cgs/1.0e-8)*pow(temp_cgs/1.0e6, -3.5);

  return kappa_cgs;
}

//rosseland mean free free absorption
Real kappa_ff_ross(Real temp, Real rho){
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real kappa_cgs = 7.73e-7*(rho_cgs/1.0e-8)*pow(temp_cgs/1.0e6, -3.5);

  return kappa_cgs;
}

void GetCombineOpacity(MeshBlock *pmb, AthenaArray<Real> &prim){

  NRRadiation *pnrrad=pmb->pnrrad;
  int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
  //int kl=pmb->kl, ku=pmb->ku, js=pmb->jl, je=pmb->ju, is=pmb->il, ie=pmb->iu;
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl=js, ju=je;
  int kl=ks, ku=ke;
  //ONE-DIMENSIONAL ONLY!
  if (pmb->pmy_mesh->f2){
    jl = js - NGHOST;
    ju = je + NGHOST;
  }
  if (pmb->pmy_mesh->f3){
    kl = ks - NGHOST;
    ku = ke + NGHOST;
  }

  // electron scattering opacity
  Real kappas = 0.2 * (1.0 + 0.6);
  Real kappa_sct_ross = 0.0;
  Real kappa_ross = 0.0;
  Real T_ion = 1.0e4;//ionization temperature, below which assuming kappa_scatter=0
  Real T_llim = 1.0e4;//lower lim of TOPs data, temperature at which switch Combined opacity grey opacity including dust
  Real T_dust = 4.0e3;//where roughly opacity rises again due to dust, what to do for this?
  
  for (int k=kl; k<=ku; k++){
    for (int j=jl; j<=ju; j++){
      for (int i=il; i<=iu; i++){
	Real rho = prim(IDN,k,j,i);
	Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i); //std::max(prim(IPR,k,j,i)/prim(IDN,k,j,i), tfloor);
	Real rho_cgs = rho*rho_unit;
	Real temp_cgs = temp*temp_unit;
  
	Real kappa_s, kappa_ross, kappa_planck;
	combineopacity(rho_cgs, temp_cgs, kappa_sct_ross, kappa_planck);
  
	if(kappa_sct_ross < kappa_es){
	  if(temp < T_ion/temp_unit){
	    kappa_ross = kappa_sct_ross;
	    kappa_s = 0.0;
	  }else{
	    kappa_ross = 0.0;
	    kappa_s = kappa_sct_ross;
	  }
	}else{
	  kappa_ross = kappa_sct_ross - kappa_es;
	  kappa_s = kappa_es;
	}
	// //assume a simple constant opacity for test
	// kappa_s = 0.32;
	// kappa_ross = 1.0;
	// kappa_planck = 1.0;

	//one frequency
	pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	pnrrad->sigma_a(k,j,i,0) = kappa_ross * rho * rho_unit * l_unit; //rosseland mean
	pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
	pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean
      
      }//end i
    }//end j
  }//end k

}

void combineopacity(const Real rho, const Real tgas, Real &kappa_ross, Real &kappa_planck){

  //STEP1: find index of temperature and density range
  //index searching segment in rho grid
  int nrho1 = 0;
  int nrho2 = 0;

  while(( rho > combine_rho_grid(nrho2)) && (nrho2 < n_rho-1)){
    nrho1 = nrho2;
    nrho2++;
  }
  //if hits the end of table, set two index equal
  if(nrho2==n_rho-1 && (rho > combine_rho_grid(nrho2))){
    nrho1=nrho2;
  }


  //index searching segments in temperature grid
  int nt1 = 0;
  int nt2 = 0;
  while((tgas > combine_temp_grid(nt2)) && (nt2 < n_tem-1)){
    nt1 = nt2;
    nt2++;
  }
  //if hits the end of table, set two index equal
  if(nt2==n_tem-1 && (tgas > combine_temp_grid(nt2))){
    nt1=nt2;
  }

  //STEP2: read the templated opacities, get ready for interpolation
  
  Real kappa_t1_rho1_gray=combine_ross_table(nt1,nrho1);
  Real kappa_t1_rho2_gray=combine_ross_table(nt1,nrho2);
  Real kappa_t2_rho1_gray=combine_ross_table(nt2,nrho1);
  Real kappa_t2_rho2_gray=combine_ross_table(nt2,nrho2);

  Real planck_t1_rho1_gray=combine_planck_table(nt1,nrho1);
  Real planck_t1_rho2_gray=combine_planck_table(nt1,nrho2);
  Real planck_t2_rho1_gray=combine_planck_table(nt2,nrho1);
  Real planck_t2_rho2_gray=combine_planck_table(nt2,nrho2);

  //in the case the temperature is out of range, extrapolate planck mean opacity by T^-3.5
  Real logt = log10(tgas);
  Real logtlim_table = log10(combine_temp_grid(n_tem-1));
  if(nt2 == n_tem-1 && (logt > logtlim_table)){
    Real scaling = pow(10.0, -3.5*(logt - logtlim_table));
    planck_t1_rho1_gray *= scaling;
    planck_t1_rho2_gray *= scaling;
    planck_t2_rho1_gray *= scaling;
    planck_t2_rho2_gray *= scaling;
  }

  //Note that if density is below the tabulated value, will use the lowest temperature in table

  Real rho_1 = combine_rho_grid(nrho1);
  Real rho_2 = combine_rho_grid(nrho2);

  Real t_1 = combine_temp_grid(nt1);
  Real t_2 = combine_temp_grid(nt2);

  //printf("rho1:%g, rho2:%g, t1:%g, t2:%g\n", rho_1, rho_2, t_1, t_2);

  //SPEP 3: Rossland opacity interpolation
  if (nrho1 == nrho2){ //if density both on lower or upper end of table 
    if (nt1 == nt2){ //if temperature also on lower or upper end of table
      kappa_ross = kappa_t1_rho1_gray; //use the only value, don't interpolate
    }else{ //interpolate only on temperature
      kappa_ross = kappa_t1_rho1_gray + (kappa_t2_rho1_gray - kappa_t1_rho1_gray) 
            * (tgas - t_1)/(t_2 - t_1);
    }
  }else{ //if two densitites are different
    if(nt1 == nt2){ //if temperature index are the same, only interpolate density
      kappa_ross = kappa_t1_rho1_gray + (kappa_t1_rho2_gray - kappa_t1_rho1_gray) 
                                * (rho - rho_1)/(rho_2 - rho_1);
    }else{ //interpolate both density and temperature

      kappa_ross = kappa_t1_rho1_gray * (t_2 - tgas) * (rho_2 - rho)  
                           /((t_2 - t_1) * (rho_2 - rho_1))
           + kappa_t2_rho1_gray * (tgas - t_1) * (rho_2 - rho)
                                /((t_2 - t_1) * (rho_2 - rho_1))
           + kappa_t1_rho2_gray * (t_2 - tgas) * (rho - rho_1)
                                /((t_2 - t_1) * (rho_2 - rho_1))
           + kappa_t2_rho2_gray * (tgas - t_1) * (rho - rho_1)
                     /((t_2 - t_1) * (rho_2 - rho_1));
    }
  }

  //STEP4: Planck opacity interpolation
    if (nrho1 == nrho2){ //if density both on lower or upper end of table 
      if (nt1 == nt2){ //if temperature also on lower or upper end of table
        kappa_planck = planck_t1_rho1_gray;
      }else{ //interpolate only on temperature
        kappa_planck = planck_t1_rho1_gray + (planck_t2_rho1_gray - planck_t1_rho1_gray)
               *(tgas - t_1)/(t_2 - t_1);
      }
    }else{//if two densitites are different
      if (nt1 == nt2){
        kappa_planck = planck_t1_rho1_gray + (planck_t1_rho2_gray - planck_t1_rho1_gray)
               *(rho - rho_1)/(rho_2 - rho_1);
      }else{ //interpolate both density and temperature
        kappa_planck = planck_t1_rho1_gray * (t_2 - tgas) * (rho_2 - rho)
                                /((t_2 - t_1) * (rho_2 - rho_1))
              + planck_t2_rho1_gray * (tgas - t_1) * (rho_2 - rho)
                                /((t_2 - t_1) * (rho_2 - rho_1))
              + planck_t1_rho2_gray * (t_2 - tgas) * (rho - rho_1)
                                /((t_2 - t_1) * (rho_2 - rho_1))
              + planck_t2_rho2_gray * (tgas - t_1) * (rho - rho_1)
                                /((t_2 - t_1) * (rho_2 - rho_1));
      }
    }
}

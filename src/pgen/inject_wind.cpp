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

#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()

//general variables to record mesh size and cell number
static int mesh_nx1, mesh_nx2, mesh_nx3;
static Real mesh_x1min, mesh_x1max, x1ratio;
static Real mesh_x2min, mesh_x2max;
static Real mesh_x3min, mesh_x3max;

//general variables converts c.g.s unit and unit-less variables in code
static Real kappa_es;
static Real temp_unit, l_unit, rho_unit, kappa_unit, vel_unit;
static Real tfloor; //temperature floor used in radiation class
static Real dfloor, pfloor; //density, pressure floor used in hydro class

//initial background density and pressure
static Real rho_init, press_init;
static Real boundary_temp_lim; //optional, temperature uplimit at boundary
//injected wind velocity and density
static Real rho_wind, vx_wind;

//gravity and mass, these are used to calculate mass within the shell, and then added as source terms
//static AthenaArray<Real> totmass;
static AthenaArray<Real> mass_shell, mass_coord;
static Real mass_bottom_cgs;
void EnvGravity(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

//opacity function
//frequency dependent free-free
Real kappa_ff_nu(Real nu, Real temp, Real rho);
void Multi_FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

//the frequency grid
static AthenaArray<Real> fre_grid;

//frequency integrated opacity
Real kappa_ff_planck(Real temp, Real rho);
Real kappa_ff_ross(Real temp, Real rho);
void FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

// User-defined boundary conditions for hydro and radiation
void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void WindInjInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


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
  vel_unit = 2.99792458e10/pin->GetReal("radiation", "crat");

  tfloor = pin->GetOrAddReal("radiation", "tfloor", 1.0e-8);
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 1.0e-8);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 1.0e-8);

  rho_init = pin->GetOrAddReal("hydro", "rho_init", 1.0e-8);
  press_init = pin->GetOrAddReal("hydro", "press_init", 1.0e-8);

  boundary_temp_lim = pin->GetOrAddReal("problem", "boundary_temp_lim", 1.0e5);
  mass_bottom_cgs = pin->GetReal("problem", "mass_bottom_cgs");

  rho_wind = pin->GetOrAddReal("problem", "rho_wind", 1.0e-1);
  vx_wind = pin->GetOrAddReal("problem", "vx_wind", 1.0);

  // Enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    //EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HydroInnerX1);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, WindInjInnerX1);
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
      EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadInnerX1);
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

  //prepare mass array
  mass_shell.NewAthenaArray(mesh_nx1);
  mass_coord.NewAthenaArray(mesh_nx1);

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
  AllocateRealUserMeshDataField(9);
  ruser_mesh_data[0].NewAthenaArray(mesh_nx1);
  ruser_mesh_data[1].NewAthenaArray(mesh_nx1);
  ruser_mesh_data[2].NewAthenaArray(mesh_nx1);
  //keep record of enclosed mass in each radius
  ruser_mesh_data[3].NewAthenaArray(mesh_nx1); //x1coordinate
  ruser_mesh_data[4].NewAthenaArray(mesh_nx1); //mass in each shell
  ruser_mesh_data[5].NewAthenaArray(mesh_nx1); //mass coordinate of each shell
  ruser_mesh_data[6].NewAthenaArray(mesh_nx1); //enclosed mass in each radius
  ruser_mesh_data[7].NewAthenaArray(mesh_nx1); //summbed b coefficient, not used,
  ruser_mesh_data[8].NewAthenaArray(mesh_nx1); //added energy, not used

  for(int i=0; i<mesh_nx1; i++){
    ruser_mesh_data[3](i) = x1coord[i];
  }
  
  FILE *f_init_rho;
  if ( (f_init_rho=fopen("./init_rho.txt","r"))==NULL )
    {
      printf("Open input file error: initial density, init_rho.txt");
      return;
    }

  FILE *f_init_temp;
  if ( (f_init_temp=fopen("./init_temp.txt","r"))==NULL )
    {
      printf("Open input file error: initial temperature, init_temp.txt");
      return;
    }

  FILE *f_init_vel;
  if ( (f_init_vel=fopen("./init_vel.txt","r"))==NULL )
    {
      printf("Open input file error: initial velocity, init_vel.txt");
      return;
    }

  rho_init_buff.NewAthenaArray(mesh_nx1);
  temp_init_buff.NewAthenaArray(mesh_nx1);
  vel_init_buff.NewAthenaArray(mesh_nx1);
  
  //load density, temperature, velocity
  for(int i=0; i<mesh_nx1; ++i){
    fscanf(f_init_rho, "%lf", &(rho_init_buff(i)));
    fscanf(f_init_temp, "%lf", &(temp_init_buff(i)));
    fscanf(f_init_vel, "%lf", &(vel_init_buff(i)));   
  }

  for(int i=0; i<mesh_nx1; ++i){
    ruser_mesh_data[0](i) = rho_init_buff(i) / rho_unit;
    ruser_mesh_data[1](i) = temp_init_buff(i) / temp_unit;
    ruser_mesh_data[2](i)= vel_init_buff(i) / vel_unit;   
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
      pnrrad->EnrollOpacityFunction(FreeFreeOpacity);
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

  for(int i=0; i<mesh_nx1; ++i){
    mass_shell(i) = 0.0;
    mass_coord(i) = 0.0;
  }
  
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
	
	//mass_coord_now += dmass;
	//mass_coord(index_rnow) = mass_coord_now;
	//pmb->pmy_mesh->ruser_mesh_data[5](index_rnow) = mass_coord_now;
      
      pmb->pmy_mesh->ruser_mesh_data[3](index_rnow) = pcoord->x1v(i);
      //printf("rad:%g, dmass:%g\n", pcoord->x1v(i), dmass);
      //pmb->pmy_mesh->ruser_mesh_data[4](index_rnow) = dmass;
      mass_shell(index_rnow) = dmass;

      // Real mass_rad = 0.0;
      // for (int ii=0; ii<index_rnow; ii++){
      // 	mass_rad += mass_shell(ii);//pmy_mesh->ruser_mesh_data[4](ii);
      // }
      // mass_coord(index_rnow) = mass_rad;
      // pmb->pmy_mesh->ruser_mesh_data[5](index_rnow) = mass_rad;
      
    }//i
 
    
  }//loop over meshblocks

  //broadcast all nodes with mass in each shell
#ifdef MPI_PARALLEL
    AthenaArray<Real> cachearray0, cachearray1;
    cachearray0.NewAthenaArray(mesh_nx1);
    //cachearray1.NewAthenaArray(mesh_nx1);

    //MPI_Bcast(&(mass_shell(0)), mesh_nx1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&(mass_shell(0)), &(cachearray0(0)), mesh_nx1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&(mass_coord(0)), &(cachearray1(0)), mesh_nx1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i=0; i<mesh_nx1; ++i){
      mass_shell(i) = cachearray0(i);
      //mass_coord(i) = cachearray1(i);
    }

    cachearray0.DeleteAthenaArray();
    //cachearray1.DeleteAthenaArray();
#endif  

  
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

	//load data, find current r index
        Real x_now = pcoord->x1f(i);
        int index_xnow = getindex(x1coord, x_now);
	iuser_meshblock_data[0](i) = index_xnow;

	Real rho_now = pmy_mesh->ruser_mesh_data[0](index_xnow);
	Real temp_now = pmy_mesh->ruser_mesh_data[1](index_xnow);
	Real vel_now = pmy_mesh->ruser_mesh_data[2](index_xnow);
        //printf("x_now:%g, index_x:%d, rho:%g, temp:%g, vel:%g\n", x_now, index_xnow, rho_now, temp_now, vel_now);

        phydro->u(IDN,k,j,i) = rho_now;
        phydro->u(IM1,k,j,i) = rho_now * vel_now;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
	
	if (NR_RADIATION_ENABLED){

	  Real rho = rho_now;//phydro->w(IDN,k,j,i);
	  Real temp = temp_now;//phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);

	  Real rho_cgs = rho*rho_unit;
	  Real temp_cgs = temp*temp_unit;
  
	  Real kappaa = 0.0;
	  Real kappa_s, kappa_ross, kappa_planck;
	  combineopacity(rho_cgs, temp_cgs, kappa_ross, kappa_planck);
	  Real t_ion = 1.0e4;
  
	  if(kappa_ross < kappa_es){
	    if(temp < t_ion/temp_unit){
	      kappaa = kappa_ross;
	      kappa_s = 0.0;
	    }else{
	      kappaa = 0.0;
	      kappa_s = kappa_ross;
	    }
	  }else{
	    kappaa = kappa_ross - kappa_es;
	    kappa_s = kappa_es;
	  }
	  //printf("gid:%d, i:%d, rho_cgs:%g, temp_cgs:%g, rho:%g, temp:%g, rho_unit:%g, temp_unit:%g, kappas:%g, kappar:%g, kappap:%g\n", gid, i, rho_cgs, temp_cgs, rho, temp, rho_unit, temp_unit, kappa_s, kappaa, kappa_planck);
  
	  //one frequency, grey rhd for now
	  pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	  pnrrad->sigma_a(k,j,i,0) = kappaa * rho * rho_unit * l_unit; //rosseland mean
	  pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
	  pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean

	  //initialize intensity
	  for (int ifr=0; ifr<pnrrad->nfreq; ++ifr){
	    for(int n=0; n<pnrrad->nang; ++n){
	      int ang=ifr*pnrrad->nang+n;
	      pnrrad->ir(k,j,i,ang) = 0.0;//use temp_now^4 if assuming initial trad=tgas
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
        prim(IVX,k,j,is-i) = prim(IVX,k,j,is+i-1); //try only reflect velocity
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        if (NON_BAROTROPIC_EOS){
          prim(IPR,k,j,is-i) = prim(IPR,k,j,is);
        }

      }//end R
    }//end theta
  }//end Phi

}

void WindInjInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=ngh; ++i) {//R
        prim(IDN,k,j,is-i) = rho_wind; //prim(IDN,k,j,is); 
        prim(IVX,k,j,is-i) = vx_wind; //prim(IVX,k,j,is); 
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        if (NON_BAROTROPIC_EOS){
          prim(IPR,k,j,is-i) = prim(IPR,k,j,is); //fmin(prim(IPR,k,j,is), (boundary_temp_lim/temp_unit)*rho_wind);
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

// // Temporary function to copy intensity
// void CopyIntensity_(Real *iri, Real *iro, int li, int lo, int n_ang) {
//   // here ir is only intensity for each cell and each frequency band
//   for (int n=0; n<n_ang; ++n) {
//     int angi = li * n_ang + n;
//     int ango = lo * n_ang + n;
//     iro[angi] = iri[ango];
//     iro[ango] = iri[angi];
//   }
// }

// //use reflecting radiation condition
// void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
//                 const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
//                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
//   // copy radiation variables into ghost zones,
//   // reflect rays along angles with opposite nx
//   const int noct = pmb->pnrrad->noct;
//   int n_ang = pmb->pnrrad->nang/noct; // angles per octant
//   const int nfreq = pmb->pnrrad->nfreq; // number of frequency bands

//   for (int k=ks; k<=ke; ++k) {
//     for (int j=js; j<=je; ++j) {
//       for (int i=1; i<=ngh; ++i) {
//         for (int ifr=0; ifr<nfreq; ++ifr) {
//           //AthenaArray<Real> &var = *ir;
//           Real *iri = &ir(k,j,(is+i-1),ifr*pmb->pnrrad->nang);
//           Real *iro = &ir(k,j, is-i, ifr*pmb->pnrrad->nang);
//           CopyIntensity_(iri, iro, 0, 1, n_ang);
//           if (noct > 2) {
//             CopyIntensity_(iri, iro, 2, 3, n_ang);
//           }
//           if (noct > 3) {
//             CopyIntensity_(iri, iro, 4, 5, n_ang);
//             CopyIntensity_(iri, iro, 6, 7, n_ang);
//           }
//         }
//       }
//     }
//   }
//   return;
// }


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


void EnvGravity(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
    const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar){

  //similar to default point mass, just use enclosed mass in each radius
  Hydro *phydro=pmb->phydro;
  Coordinates *pcoord=pmb->pcoord;

  //hard-coded parameters, Eq17 SNEC_NOTE
  //Real Rt = 100.0;
  //Real Rm = 100.0;
  
  Real mass_unit = rho_unit * pow(l_unit, 3);
  Real time_unit = l_unit / vel_unit;
  Real G_code = 6.67259e-8/(pow(l_unit,3)*pow(mass_unit,-1)*pow(time_unit, -2)); //hard coded for now
  Real mass_bottom_code = mass_bottom_cgs/mass_unit;

  ///* Try Inject Energy As Explicit Source */
  
  // Real inj_mass_spread_code = inj_mass_spread/mass_unit;
  // //getting distribution of mass injection
  // Real a_coef = std::log(Rm) / inj_mass_spread_code; // a' coefficient in (Eq18) SNEC_NOTE

  // //for each meshblock, find the global coordinate where mass is within mass injection
  // Real mass_inj_accumulate = 0.0;
  // Real mass_inj_tot = 0.0;
  // Real coef_b_sum_factor = 0.0; //denominator of b' coefficient in (Eq18) SNEC_NOTE
  // Real mass_coord = 0.0;
  // int count = 0;
  // for (int ii=0; ii<mesh_nx1; ii++){    
  //   if (mass_inj_accumulate <= inj_mass_spread_code ){
  //     count += 1;
  //     mass_inj_tot += mass_shell(ii);

  //     //ge the mass coordinate for injection
  //     mass_coord += mass_shell(ii);
  //     pmb->pmy_mesh->ruser_mesh_data[5](ii) = mass_coord;
  //     coef_b_sum_factor += std::exp(-a_coef * mass_coord) * mass_shell(ii);
  //     pmb->pmy_mesh->ruser_mesh_data[7](ii) = coef_b_sum_factor;
  //     // //debug lines
  //     // if (mass_shell(ii) < 100 && pmb->pmy_mesh->ncycle>1){
  //     // 	printf("rank:%d, gid:%d, mass shell is small in ii=%d, mass_shell:%g, count:%d\n", Globals::my_rank, pmb->gid, ii, mass_shell(ii), count);
  //     // }
      
  //   }
  //   mass_inj_accumulate += mass_shell(ii);
  // }
  // //printf("my_rank:%d, coef_b_sum_factor:%g, acoef:%g\n", Globals::my_rank, coef_b_sum_factor, a_coef);
  // Real x1_coord_inj_lim = pmb->pmy_mesh->ruser_mesh_data[3](count); //x1coord[count];
  // // printf("my_rank:%d, count:%d, mass_x1coord:%g, mass_inj_tot:%g\n", Globals::my_rank, count, x1_coord_inj_lim, mass_inj_tot);

  
  // //first, calculate total energy injected at this time

  // Real tb_start = inj_tstart / time_unit;
  // Real tb_end = inj_tend / time_unit;
  // Real energy_unit = rho_unit * pow(vel_unit, 2) * pow(l_unit, 3);
  // Real Etot = final_energy / energy_unit;
  // Real c_coef = std::log(Rt) / (tb_end -  tb_start);
  // Real d_coef = c_coef * Etot / (std::exp(-c_coef * tb_start) - std::exp(-c_coef * tb_end));
  // Real Pb_now = d_coef * std::exp(-c_coef * time); //in unit of energy per unit time
  // Real Pb_per_mass = Pb_now / mass_inj_tot; //total energy per mass per time

  // if (Globals::my_rank==0){
  //   printf("time_now:%g, pb_now:%g, smooth_mass:%d\n", time, Pb_now, MASS_SMOOTH_FLAG);
  // }

  // calculate accumulated mass at each radius 
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = phydro->w(IDN,k,j,i);
	
	//get the enclosed mass
	Real r_now = pcoord->x1f(i);
	int index_rnow = pmb->iuser_meshblock_data[0](i);
	pmb->pmy_mesh->ruser_mesh_data[4](index_rnow) = mass_shell(index_rnow);

	Real mass_rad = 0.0;			 
	for (int ii=0; ii<index_rnow; ii++){
	  mass_rad += mass_shell(ii);//pmy_mesh->ruser_mesh_data[4](ii);
	  
	}//end ii
	
	Real gm_en = G_code * (mass_rad+mass_bottom_code);//add the mass where it get truncated

	//record enclosed mass
	pmb->pmy_mesh->ruser_mesh_data[6](index_rnow) = mass_rad;
	//(TODO) hard-coded changed the protected variables in Coordinates class for now
        Real src = dt*den*pmb->pcoord->coord_src1_i_(i)*gm_en/pmb->pcoord->x1v(i);
        cons(IM1,k,j,i) -= src; //momentum
	//keep record of enclosed graivty momentum source term
	pmb->ruser_meshblock_data[0](1,k,j,i) = -den*pmb->pcoord->coord_src1_i_(i)*gm_en/pmb->pcoord->x1v(i);
        if (NON_BAROTROPIC_EOS) {
          cons(IEN,k,j,i) -=
              dt*0.5*(pmb->pcoord->phy_src1_i_(i)*phydro->flux[X1DIR](IDN,k,j,i)*gm_en
                      +pmb->pcoord->phy_src2_i_(i)*phydro->flux[X1DIR](IDN,k,j,i+1)*gm_en);

        }
	//keep record of enclosed gravity KE source term 
	pmb->ruser_meshblock_data[0](2,k,j,i) = 0.5*(pmb->pcoord->phy_src1_i_(i)*phydro->flux[X1DIR](IDN,k,j,i)*gm_en
                      +pmb->pcoord->phy_src2_i_(i)*phydro->flux[X1DIR](IDN,k,j,i+1)*gm_en);

	// //inject energy as SNEC did
	// if (pmb->pmy_mesh->ncycle>1){
	//   pmb->iuser_meshblock_data[1](i) = 0; //reset user mesh data for the flag
	//   pmb->ruser_meshblock_data[3](1,k,j,i) = 0.0; //reset user meshblock data for the energy injection record
	//   pmb->ruser_meshblock_data[0](0,k,j,i) = 0.0;
	//   if ((r_now < x1_coord_inj_lim)){
	//     Real einj_now = 0.0;
	//     if (MASS_SMOOTH_FLAG==0){
	//       einj_now = Pb_per_mass * den * pcoord->GetCellVolume(k,j,i);
	//     }else if (MASS_SMOOTH_FLAG==1){
	//     //mass distribution
	//       Real b_coef = std::exp(-a_coef * mass_rad);
	//       einj_now = Pb_now * b_coef * den * pcoord->GetCellVolume(k,j,i)/coef_b_sum_factor;
	//     }
	    
	//     if (NON_BAROTROPIC_EOS) {
	//       cons(IEN,k,j,i) += einj_now * dt;
	//     }
	//     pmb->ruser_meshblock_data[0](0,k,j,i) = einj_now;
	//     //record the flag
	//     pmb->iuser_meshblock_data[1](i) = 1;
	//     //record time-integrated energy injection
	//     pmb->ruser_meshblock_data[3](0,k,j,i) += einj_now * dt;
	//     pmb->ruser_meshblock_data[3](1,k,j,i) = einj_now;
	    
	//     if (einj_now == 0){
	//       Real b_coef = std::exp(-a_coef * mass_rad);
	//       printf("zero energy injection x1:%g, Pb_per_mass:%g, den:%g, cellvol:%g, Pb_now:%g, b_coef:%g, coef_b_sum:%g\n", pcoord->x1v(i), Pb_per_mass, den, pcoord->GetCellVolume(k,j,i), Pb_now, b_coef, coef_b_sum_factor);
	//     }
	    
	//   }
	// }
	
      }
    }
  }

  
  return;
}

// Real einj_int(MeshBlock *pmb, int iout)
// {
//   Real einj=0;
//   int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

//   for(int k=ks; k<=ke; k++) {
//     for(int j=js; j<=je; j++) {
//       for(int i=is; i<=ie; i++) {
// 	einj += pmb->ruser_meshblock_data[3](0,k,j,i);
//       }
//     }
//   }

//   return einj;
// }

// Real einj_dt(MeshBlock *pmb, int iout)
// {
//   Real einj=0;
//   int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

//   for(int k=ks; k<=ke; k++) {
//     for(int j=js; j<=je; j++) {
//       for(int i=is; i<=ie; i++) {
// 	einj += pmb->ruser_meshblock_data[3](1,k,j,i);
//       }
//     }
//   }

//   return einj;
// }

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

  for (int k=kl; k<=ku; k++){
    for (int j=jl; j<=ju; j++){
      for (int i=il; i<=iu; i++){
	
	Real rho = prim(IDN,k,j,i);
	Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i); //std::max(prim(IPR,k,j,i)/prim(IDN,k,j,i), tfloor);
	Real rho_cgs = rho*rho_unit;
	Real temp_cgs = temp*temp_unit;
	
	Real kappa_s, kappa_ross, kappa_planck;
	kappa_s = kappa_es;
	kappa_ross = kappa_ff_ross(temp_cgs, rho_cgs);
	kappa_planck = kappa_ff_planck(temp_cgs, rho_cgs);
	
	//one frequency
	pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	pnrrad->sigma_a(k,j,i,0) = kappa_ross * rho * rho_unit * l_unit; //rosseland mean
	pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
        pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean
      
      }//end i
    }//end j
  }//end k

}

 //input code unit, output code unit, planck mean free free absorption
Real kappa_ff_planck(Real temp, Real rho){
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real kappa_cgs = 2.86e-5*(rho_cgs/1.0e-8)*pow(temp_cgs/1.0e6, -3.5);

  return kappa_cgs/kappa_unit;
}

//input code unit, output code unit, rosseland mean free free absorption
Real kappa_ff_ross(Real temp, Real rho){
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real kappa_cgs = 7.73e-7*(rho_cgs/1.0e-8)*pow(temp_cgs/1.0e6, -3.5);

  return kappa_cgs/kappa_unit;
}

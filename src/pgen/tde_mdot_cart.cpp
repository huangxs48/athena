//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// python configure.py --prob tde_mdot_cart -nr_radiation -mpi -hdf5 --hdf5_path=/usr/local/hdf5-mpi
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

// Athena++ headers
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

static Real GM1, gm;
static Real kappa_es;
static Real temp_unit,l_unit,rho_unit;
static Real kappa_unit, kappa_es_code;
static Real tfloor;
static Real dfloor, pfloor;
static Real user_dt;
static Real r_isco;

//Injection point
static Real x_inj, y_inj, z_inj, local_dens, local_vx, local_vy, local_vz, local_press;
//Initialize domain
static Real rho_init, press_init;
static Real temp_stream;

//opacity tables
static AthenaArray<Real> opacitytable;
static AthenaArray<Real> planckopacity;
static AthenaArray<Real> logttable;
static AthenaArray<Real> logrhottable;
static AthenaArray<Real> logttable_planck;
static AthenaArray<Real> logrhottable_planck;

//new combined opacity table
static AthenaArray<Real> combine_temp_grid;
static AthenaArray<Real> combine_rho_grid;
static AthenaArray<Real> combine_ross_table;
static AthenaArray<Real> combine_planck_table;
static int n_rho = 140;
static int n_tem = 70;
void combineopacity(const Real rho, const Real tgas, Real &kappa_ross, Real &kappa_planck);
void GetCombineOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);
int replace_low_dens_ff;
Real rho_cut_ff; 

//refinement
static Real x_thresh, y_thresh, z_thresh;
static int blocksizex1_base, blocksizex2_base, blocksizex3_base;
static Real delta_theta_base;
static Real x_coarse_thresh, y_coarse_thresh, z_coarse_thresh;

//stream injection boundary
static Real r0, inj_thresh, b0;
static Real x_inject2; 
int rho1_flag;

// User-defined boundary conditions for disk simulations
void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StreamInjectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void HydroOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadInnerX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadOuterX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadInnerX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadOuterX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void GeneralNewtonianPotentialCart(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                              const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
			      AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

Real kappa_ff_planck(Real temp, Real rho);
Real kappa_ff_ross(Real temp, Real rho);
Real kappa_ff_nu(Real nu, Real temp, Real rho);

//user history 
Real massflux_AInj_x1(MeshBlock *pmb, int iout);
Real massflux_AInj_x3(MeshBlock *pmb, int iout);
Real massflux_Inj_x1(MeshBlock *pmb, int iout);
Real massflux_Inj_x3(MeshBlock *pmb, int iout);

Real massfluxix1(MeshBlock *pmb, int iout);
Real massfluxox1(MeshBlock *pmb, int iout);

//AMR condition
int RefinementCondition3(MeshBlock *pmb);
Real orbit_fit_r(Real phi);
Real orbit_fit_phi(Real r);

//try small time step like rad
Real RadTimeStep(MeshBlock *pmb);

//time dependent mdot
static int mdot_table_ng; //number of segment groups from spline fitting
static Real mdot_rho1; //when setting rho=1, the total mdot through phi and theta
static AthenaArray<Real> mdot_table;
static AthenaArray<Real> mdot_time_table;
Real GetMdot(MeshBlock *pmb, Real tnow);


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

Real RadTimeStep(MeshBlock *pmb){
  //Real dt_rad = 2.2039989460532532e-06;
  return user_dt;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  int blocksizex1 = pin->GetOrAddInteger("meshblock", "nx1", 1);
  int blocksizex2 = pin->GetOrAddInteger("meshblock", "nx2", 1);
  int blocksizex3 = pin->GetOrAddInteger("meshblock", "nx3", 1);

  blocksizex1_base = pin->GetOrAddInteger("meshblock", "nx1", 1);
  blocksizex2_base = pin->GetOrAddInteger("meshblock", "nx2", 1);
  blocksizex3_base = pin->GetOrAddInteger("meshblock", "nx3", 1);
  // int theta_max = pin->GetReal("mesh","x2max");
  // int theta_min = pin->GetReal("mesh","x2min");
  // int theta_base = pin->GetInteger("mesh","nx2");
  // delta_theta_base = (theta_max - theta_min)/theta_base * blocksizex2_base;

  GM1 = pin->GetOrAddReal("problem","GM1",1.0);
  kappa_es = pin->GetReal("problem", "kappa_es");

  temp_unit = pin->GetReal("problem", "temp_unit");
  l_unit = pin->GetReal("problem", "l_unit");
  rho_unit = pin->GetReal("problem", "rho_unit");

  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.001);
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 0.001);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 0.001);
  r_isco = pin->GetOrAddReal("problem", "r_isco", 3.0);

  //opacity flags
  replace_low_dens_ff = pin->GetOrAddInteger("problem", "replace_low_dens_ff", 0);
  rho_cut_ff = pin->GetOrAddReal("problem", "rho_cut_ff", 1.0e-7);  

  user_dt = pin->GetOrAddReal("problem", "user_dt", 1.0e-6);
  //Initialize the injection point
  x_inj = pin->GetReal("problem", "x_inj");
  y_inj = pin->GetReal("problem", "y_inj");
  z_inj = pin->GetReal("problem", "z_inj");
  local_dens = pin->GetReal("problem", "local_dens");
  local_vx = pin->GetReal("problem", "local_vx");
  local_vy = pin->GetReal("problem", "local_vy");
  local_vz = pin->GetReal("problem", "local_vz");
  //local_press = pin->GetReal("problem", "local_press");
  rho1_flag = pin->GetOrAddInteger("problem", "rho1_flag", 0);

  //x_thresh = pin->GetOrAddReal("problem", "x_thresh", 0.1);  
  //y_thresh = pin->GetOrAddReal("problem", "y_thresh", 0.01);
  //z_thresh = pin->GetOrAddReal("problem", "z_thresh", 0.01);
  x_coarse_thresh = pin->GetOrAddReal("problem", "x_coarse_thresh", 0.1);
  y_coarse_thresh = pin->GetOrAddReal("problem", "y_coarse_thresh", 0.1);
  z_coarse_thresh = pin->GetOrAddReal("problem", "z_coarse_thresh", 0.1);

  //stream injection boundary
  b0 = pin->GetOrAddReal("problem", "b0", 0.1);
  r0 = pin->GetOrAddReal("problem", "r0", 7.53);
  inj_thresh = pin->GetOrAddReal("problem", "inj_thresh", 7.53);
  x_inject2 = pin->GetOrAddReal("problem", "x_inject2", 407.886);

  //Initialize domain
  rho_init = pin->GetReal("problem", "rho_init");
  press_init = pin->GetReal("problem", "press_init");
  temp_stream = pin->GetReal("problem", "temp_stream");  
  local_press = local_dens*temp_stream/temp_unit;

  //kappa_unit: cm^2/g
  kappa_unit = 1.0/(rho_unit*l_unit);
  kappa_es_code = kappa_es/kappa_unit;

  //mdot table
  mdot_table_ng = pin->GetInteger("problem", "mdot_table_ng");
  mdot_rho1 = pin->GetReal("problem", "mdot_rho1");

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HydroInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, StreamInjectOuterX1);//HydroOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, HydroInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, HydroOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, HydroInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, HydroOuterX3);
  }

  //EnrollUserTimeStepFunction(RadTimeStep);

  // Enroll AMR
  if(adaptive==true)
   EnrollUserRefinementCondition(RefinementCondition3);

  //load mdot table
  mdot_table.NewAthenaArray(mdot_table_ng);//segment bounds for cubic spline fitting
  mdot_time_table.NewAthenaArray(mdot_table_ng);//coeficients for cubic spline fitting

  FILE *fmdot_time_table;
  if ( (fmdot_time_table=fopen("./mdot_tday_table.txt","r"))==NULL ){
    printf("Open input file error mdot time table");
    return;
  }
  //read mdot
  // time in days
  for(int i=0; i<mdot_table_ng; ++i){
    fscanf(fmdot_time_table,"%lf",&(mdot_time_table(i)));
  }
  FILE *fmdot_table;
  if ( (fmdot_table=fopen("./mdot_table.txt","r"))==NULL ){
    printf("Open input file error mdot table");
    return;
  }
  // mdot in cgs
  for(int i=0; i<mdot_table_ng; ++i){
    fscanf(fmdot_table,"%lf",&(mdot_table(i)));
  }

  //printf("testing spline txt, mdot_table_ng(7):%g, mdot_coef:%g %g %g %g\n", mdot_x_breaks(7), mdot_coef(7, 0), mdot_coef(7, 1), mdot_coef(7, 2), mdot_coef(7, 3));

  if (NR_RADIATION_ENABLED){
    //Enroll rad boundaries
    if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadInnerX1);
    }
    if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, RadOuterX1);
    }
    if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::inner_x2, RadInnerX2);
    }
    if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::outer_x2, RadOuterX2);
    }
    if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::inner_x3, RadInnerX3);
    }
    if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::outer_x3, RadOuterX3);
    }

    //load combined opacity
    combine_temp_grid.NewAthenaArray(n_tem);
    combine_rho_grid.NewAthenaArray(n_rho);
    combine_ross_table.NewAthenaArray(n_tem, n_rho);
    combine_planck_table.NewAthenaArray(n_tem, n_rho);

    FILE *f_combineopacity;
    
    if ( (f_combineopacity=fopen("./output_grey_combined.txt","r"))==NULL )
    {
      printf("Open input file error combined opacity table");
      return;
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

    //printf("testing frequency table\n");
    //printf("grey planck:%g, ross:%g\n", combine_planck_table(30, 50), combine_ross_table(30, 50));
    //Real kappa_p_test, kappa_r_test;
    //combineopacity(1.0e-15, 1.0e6, kappa_r_test, kappa_p_test);
    //printf("interpolated kappa_p:%g, kappa_r:%g\n", kappa_p_test, kappa_r_test);
    
 
  }
  

  EnrollUserExplicitSourceFunction(GeneralNewtonianPotentialCart);//(PointMassPotential);

  AllocateUserHistoryOutput(6);
  EnrollUserHistoryOutput(0, massflux_AInj_x1, "massflux_AInj_x1");//mass flux outer
  EnrollUserHistoryOutput(1, massflux_AInj_x3, "massflux_AInj_x3");//mass flux outer
  EnrollUserHistoryOutput(2, massflux_Inj_x1, "massflux_Inj_x1");//mass flux outer
  EnrollUserHistoryOutput(3, massflux_Inj_x3, "massflux_Inj_x3");//mass flux outer
  EnrollUserHistoryOutput(4, massfluxix1, "massfluxix1");//mass flux inner boundary
  EnrollUserHistoryOutput(5, massfluxox1, "massfluxox1");//mass flux outer
  
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  int blocksizex1 = pin->GetOrAddInteger("meshblock", "nx1", 1);
  int blocksizex2 = pin->GetOrAddInteger("meshblock", "nx2", 1);
  int blocksizex3 = pin->GetOrAddInteger("meshblock", "nx3", 1);

  blocksizex1 += 2*(NGHOST);
  if (blocksizex2 >1) blocksizex2 += 2*(NGHOST);
  if (blocksizex3 >1) blocksizex3 += 2*(NGHOST);

  gm = 5.0e5;
  if (NR_RADIATION_ENABLED) gm = 0.5 * pnrrad->crat * pnrrad->crat;
  
  AllocateRealUserMeshBlockDataField(3); 
  ruser_meshblock_data[0].NewAthenaArray(4,blocksizex3, blocksizex2, blocksizex1); //store external source
  ruser_meshblock_data[1].NewAthenaArray(4,blocksizex3, blocksizex2, blocksizex1); //injection face area and flux
  ruser_meshblock_data[2].NewAthenaArray(5,blocksizex3, blocksizex2, blocksizex1); //diagnostic energy source terms, dt

  //printf("user real data\n");
  AllocateIntUserMeshBlockDataField(1);
  iuser_meshblock_data[0].NewAthenaArray(3); //labels special meshblocks
  iuser_meshblock_data[0](0) = 0; //inejection point meshblocks
  iuser_meshblock_data[0](1) = 0; //near polar meshblocks
  iuser_meshblock_data[0](2) = 0; //orbit refinement, fix grid
 // printf("user_int_data\n");  

  for(int k=0; k<blocksizex3; k++){
   for(int j=0; j<blocksizex2; j++){
     for(int i=0; i<blocksizex1; i++){

       //initialize the storing flux and area
       ruser_meshblock_data[1](0,k,j,i) = 0.0;
       ruser_meshblock_data[1](1,k,j,i) = 0.0;
       ruser_meshblock_data[1](2,k,j,i) = 0.0;
       ruser_meshblock_data[1](3,k,j,i) = 0.0;
       
       Real z_coord = pcoord->x3v(k);
       Real y_coord = pcoord->x2v(j);
       Real x_coord = pcoord->x1v(i);
       Real rad = sqrt(x_coord*x_coord + y_coord*y_coord + z_coord*z_coord);

       //Real dr = (pcoord->dx1f(i-1)+pcoord->dx1f(i+1))/2.0;
       //Real dphi = (pcoord->dx3f(k-1)+pcoord->dx3f(k+1))/2.0;
       //Real dth = (pcoord->dx2f(j-1)+pcoord->dx2f(j+1))/2.0;	  

       //injection point
       if (fabs(x_coord-x_inj) <= x_coarse_thresh){
       	 if (fabs(y_coord-y_inj) <= y_coarse_thresh){// if within L1 point
       	   if (fabs(z_coord - 0.0)<= z_coarse_thresh){
       	     //printf("coarse\n");
	     //printf("r:%g, ph:%g, th:%g, dis_r:%g, dis_phi:%g, dis_th:%g, i:%d, gid:%d, dr:%g, dph:%g, dth:%g\n", rad, phi_coord, theta_coord, fabs(rad-x_inject), fabs(phi_coord - y_inject), fabs(theta_coord - PI/2.0), i, gid, pcoord->dx1f(i), pcoord->dx3f(k), pcoord->dx2f(j));
       	     iuser_meshblock_data[0](0) = 1;
       	   }
       	 }
       }//end injection point


     }//k
   }//j
  }//i

  AllocateUserOutputVariables(9);
  SetUserOutputVariableName(0, "GravSrc_IM1");
  SetUserOutputVariableName(1, "GravSrc_IM2");
  SetUserOutputVariableName(2, "GravSrc_IM3");
  SetUserOutputVariableName(3, "GravSrc_IEN");
  SetUserOutputVariableName(4, "GravSrc_etot");
  SetUserOutputVariableName(5, "GravSrc_ein");
  SetUserOutputVariableName(6, "dt");
  SetUserOutputVariableName(7, "RadSrc_ein");
  SetUserOutputVariableName(8, "RadSrc_ein_new");

  if (NR_RADIATION_ENABLED){
    pnrrad->EnrollOpacityFunction(GetCombineOpacity);//(opalopacity);
  }

  return;
}

//density curvature refinement test
int RefinementCondition3(MeshBlock *pmb){
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  Real min_r = 9999.0;
  Real max_dens = 0.0;

  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
	Real epsr = (std::abs(w(IDN,k,j,i+1) - 2.0*w(IDN,k,j,i) + w(IDN,k,j,i-1))
		     + std::abs(w(IDN,k,j+1,i) - 2.0*w(IDN,k,j,i) + w(IDN,k,j-1,i))
		     + std::abs(w(IDN,k+1,j,i) - 2.0*w(IDN,k,j,i) + w(IDN,k-1,j,i)))/w(IDN,k,j,i);
	Real dist_r = sqrt(SQR(pmb->pcoord->x1v(i)) + SQR(pmb->pcoord->x2v(j)) + SQR(pmb->pcoord->x3v(k)));
	min_r = std::min(dist_r, min_r);
	// Real epsp = (std::abs(w(IPR,k,j,i+1) - 2.0*w(IPR,k,j,i) + w(IPR,k,j,i-1))
	//             + std::abs(w(IPR,k,j+1,i) - 2.0*w(IPR,k,j,i) + w(IPR,k,j-1,i)))
	//             /w(IPR,k,j,i);
	// Real eps = std::max(epsr, epsp);
	maxeps = std::max(maxeps, epsr);
	max_dens = std::max(max_dens, w(IDN,k,j,i));
      }//i
    }//j
  }//k

  //in the initialization step, the flag is based on coarse cell size
  int injection_flag = pmb->iuser_meshblock_data[0](0);
  //after initialization, the injection grids are refined, refinement criterial should 
  //follow finer cell size, so switch to m_data(1)
  if (pmb->pmy_mesh->ncycle>1) injection_flag = pmb->iuser_meshblock_data[0](1);

  //if (injection_flag==1)
  //printf("gid:%d, min_r:%g, maxdens:%g\n", pmb->gid, min_r, max_dens);

  if ((maxeps>1.0 && max_dens>1.0e-3)||injection_flag ==1||min_r<3.0) return 1; //run to t=0.3
  if (maxeps<1.0e-3) return -1;//run to t=0.3
  
  Real dis_z_l = std::fabs(pmb->pcoord->x2f(pmb->js) - PI/2.0);
  Real dis_z_r = std::fabs(pmb->pcoord->x2f(pmb->je) - PI/2.0);
  Real dis_z_min = std::min(dis_z_l, dis_z_r);

  //run to t=2.9
 // if ((maxeps>1.0 && max_dens>1.0e-5) && r_min>4.0 && dis_th_min<0.088) return 1;
  //if (maxeps<1.0e-4 || r_min<4.0) return -1; //run to t=3.0 for level=4 run

  //run to t=3.3
  //if ((maxeps>5.0 && max_dens>1.0e-3) && r_min>4.5 && dis_th_min<0.088) return 1;
  //if (maxeps<0.5 || r_min<4.0) return -1; //run to t=5.

  //run to t=3.5
  //if ((maxeps>10.0 && max_dens>1.0e-3) && r_min>4.5 && dis_th_min<0.088) return 1;
  //if (maxeps<1.0 || r_min<4.0 || max_dens<5.0e-4) return -1; //run to t=1.9

  //run to t=4.02
  //if ((maxeps>1.0 && max_dens>5.0e-3) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.1 || r_min<4.8 || max_dens<1.0e-4) return -1; //run to t=5.

  return 0;

}

void Mesh::UserWorkInLoop(){
  MeshBlock *pmb = my_blocks(0);
  for(int nb=0; nb<nblocal; ++nb){ //loop over meshblocks on the same core
    pmb = my_blocks(nb);
    
    Hydro *phydro = pmb->phydro;
    Coordinates *pcoord = pmb->pcoord;
    int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
    
    //get current loacla dens
    Real t_current = pmb->pmy_mesh->time;
    Real local_dens_now_ = GetMdot(pmb, t_current);

    //hard-code local density to 1.0 in runs to calculate mass flux normalization.
    if (rho1_flag==1){
      local_dens_now_ = 1.0;
    }
    for(int k=ks; k<=ke; k++){
      for(int i=is; i<=ie; i++){
	for(int j=js; j<=je; j++){
	  for (int n=0; n<(NHYDRO);n++){
	    
	    if (phydro->u(n,k,j,i) != phydro->u(n,k,j,i)){
	      printf("block: %d, n: %d ,k: %d,j: %d,i: %d\n", pmb->gid,n,k,j,i);
	      printf("x1v: %g, x2v:%g, x3v:%g, idn:%g, im1:%g, im2:%g, im3:%g, ien:%g\n",pmb->pcoord->x1v(i), pmb->pcoord->x2v(j),pmb->pcoord->x3v(k), pmb->phydro->u(IDN,k,j,i), pmb->phydro->u(IM1,k,j,i), pmb->phydro->u(IM2,k,j,i), pmb->phydro->u(IM3,k,j,i), pmb->phydro->u(IEN,k,j,i));
	      //abort();
	    }	  
	  }//end NHYDRO

	  //apply floor
	  Real x1_now = pmb->pcoord->x1v(i);
	  Real x2_now = pmb->pcoord->x2v(j);
	  Real x3_now = pmb->pcoord->x3v(k);
	  Real r_now = sqrt(SQR(x1_now) + SQR(x2_now) + SQR(x3_now));

	  if (r_now < r_isco){
	    pmb->phydro->u(IDN,k,j,i) = dfloor;
	    pmb->phydro->u(IM1,k,j,i) = 0.0;
	    pmb->phydro->u(IM2,k,j,i) = 0.0;
	    pmb->phydro->u(IM3,k,j,i) = 0.0;
	    pmb->phydro->u(IEN,k,j,i) = pfloor /(pmb->peos->GetGamma()-1.0);
	    
	    pmb->phydro->w(IDN,k,j,i) = dfloor;
	    pmb->phydro->w(IVX,k,j,i) = 0.0;
	    pmb->phydro->w(IVY,k,j,i) = 0.0;
	    pmb->phydro->w(IVZ,k,j,i) = 0.0;
	    pmb->phydro->w(IPR,k,j,i) = pfloor;

	    if (NR_RADIATION_ENABLED){
	      for (int n=0; n<=pmb->pnrrad->n_fre_ang; n++){
		pmb->pnrrad->ir(k,j,i,n) = 0.0;
	      }
	    }
	    
	  }

	}//j
      }//i
    }//k
  }//loop over meshblocks
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){
  for(int k=ks; k<=ke; k++){
    for(int j=js; j<=je; j++){
      for(int i=is; i<=ie; i++){
  	//user_out_var(0,k,j,i) = pcoord->GetCellVolume(k,j,i);
  	user_out_var(0,k,j,i) = ruser_meshblock_data[0](0,k,j,i);
   	user_out_var(1,k,j,i) = ruser_meshblock_data[0](1,k,j,i);
   	user_out_var(2,k,j,i) = ruser_meshblock_data[0](2,k,j,i);
   	user_out_var(3,k,j,i) = ruser_meshblock_data[0](3,k,j,i);
        user_out_var(4,k,j,i) = ruser_meshblock_data[2](0,k,j,i);
        user_out_var(5,k,j,i) = ruser_meshblock_data[2](1,k,j,i);
        user_out_var(6,k,j,i) = ruser_meshblock_data[2](2,k,j,i);
        user_out_var(7,k,j,i) = ruser_meshblock_data[2](3,k,j,i);
        user_out_var(8,k,j,i) = ruser_meshblock_data[2](4,k,j,i);     
      }
     }
   }
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  //Real rho_floor = 1.527944e-05;
  //Real press_init = 2.806747e-08;//2.656888e-5;
  Real gamma_gas = peos->GetGamma();

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
       
        phydro->u(IDN,k,j,i) = rho_init;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

	if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = press_init/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }//end non barotropic

	if (NR_RADIATION_ENABLED){

	  Real rho = rho_init;
	  Real temp = press_init/rho_init;

	  for (int n=0; n<pnrrad->n_fre_ang; n++){
	    pnrrad->ir(k,j,i,n) = temp * temp * temp * temp;
	  }

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
	
	  //one frequency
	  pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	  pnrrad->sigma_a(k,j,i,0) = kappaa * rho * rho_unit * l_unit; //rosseland mean
	  pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
	  pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean
	
     
	}//end rad

      }//end i
      
      
    }//end j
  }//end k

  return;
}

Real GetMdot(MeshBlock *pmb, Real tnow){
  //tnow in code unit, convert to days
  //dm/dt in cgs
  Real vel_unit = 2.99792458e10/pmb->pnrrad->crat;
  Real time_unit = l_unit / vel_unit;
  Real t_days  =  tnow * time_unit / 24 / 3600;
  Real mdot_unit = rho_unit * pow(l_unit, 3) / time_unit;
  Real msun = 1.9891e33;
  Real year = 3.155815e7;

  int ngroup = mdot_table_ng;
  int ig1 = 0;
  int ig2 = 0;
  while ((t_days>mdot_time_table(ig2)) && (ig2<ngroup)){
    ig1 = ig2;
    ig2 += 1;
  }

  //if hits end of table
  if (ig2 == ngroup) {
    ig1 = ig2;
  }

  //read templated time and mdot
  Real mdot1 = mdot_table(ig1);
  Real mdot2 = mdot_table(ig2);

  if ((ig2 == ngroup) && (t_days > mdot_time_table(ig2))){
    mdot1 = mdot_table(ngroup);
    mdot2 = mdot_table(ngroup);
  }

  if ((ig2 == 0) && (t_days > mdot_time_table(0))){
    mdot1 = mdot_table(0);
    mdot2 = mdot_table(0);
  }

  Real mdot_interp = 0.0;
  //linear interpolate
  Real mdot_t1 = mdot_time_table(ig1);
  Real mdot_t2 = mdot_time_table(ig2);

  if ((ig1 == ig2) && (ig2 == ngroup)){
    mdot_interp = mdot_table(ngroup);
  }else if ((ig1 == ig2) && (ig2==0)){
    mdot_interp = mdot_table(0);
  }else{
    mdot_interp = mdot1 + (t_days - mdot_t1) * (mdot2 - mdot1)/(mdot_t2 - mdot_t1);
  }

  Real mdot_cgs = mdot_interp * (msun/year);
  Real mdot_code = mdot_cgs / mdot_unit;
  Real local_dens_now = mdot_code / mdot_rho1;
  //printf("t now:%g, t days:%g, mdot now:%g, local_dens:%g, mdot_code:%g, mdot_rho1:%g\n", tnow, t_days, mdot_interp, local_dens_now, mdot_code, mdot_rho1);

  
  return local_dens_now;
  
}

void GeneralNewtonianPotentialCart(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar){
  //printf("user source\n");
  //potential as Tejeda & Rosswog 13
  
  // Gravitational acceleration from orbital motion
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real rho = prim(IDN,k,j,i);
   
	Real x = pmb->pcoord->x1v(i);
	Real y = pmb->pcoord->x2v(j);
	Real z = pmb->pcoord->x3v(k);

	Real r = sqrt(x*x + y*y + z*z);

	if (r > r_isco){

	  Real dx = prim(IVX, k,j,i);
	  Real dy = prim(IVY, k,j,i);
	  Real dz = prim(IVZ, k,j,i);
	  
	  //coefficients of acceleration
	  Real coef1 = -gm * pow(1.0 - 1.0/r, 2) / pow(r, 3);
	  Real coef2 = (x*dx + y*dy + z*dz)/(pow(r, 2)*(r - 1.0));
	  Real coef3 = - (3.0/2.0) * (1.0/pow(r, 5)) 
	               * (pow(z*dx - x*dz, 2) + pow(y*dx - x*dy, 2) +pow(y*dz - z*dy, 2));
	  //accelerations
	  Real ax = coef1 * x + coef2 * dx + coef3 * x;
	  Real ay = coef1 * y + coef2 * dy + coef3 * y;
	  Real az = coef1 * z + coef2 * dz + coef3 * z;
	  
	  
	  cons(IM1,k,j,i) += dt * rho * ax;
	  cons(IM2,k,j,i) += dt * rho * ay;
	  cons(IM3,k,j,i) += dt * rho * az;
	  
	  cons(IEN,k,j,i) += dt * ax * rho * prim(IVX,k,j,i);
	  cons(IEN,k,j,i) += dt * ay * rho * prim(IVY,k,j,i) + dt * az * rho * prim(IVZ,k,j,i);
	}

      }
    }
  }
  
}


void StreamInjectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,int is, int ie, int js, int je, int ks, int ke, int ngh){

  Real t_current = pmb->pmy_mesh->time;
  Real local_dens_now_;
  
  if (rho1_flag==1){
    local_dens_now_ = 1.0; //for rho1 run
  }else{
    local_dens_now_ = GetMdot(pmb, t_current);
  }

  Real rinj_thresh=inj_thresh;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
	
	//get currtent Cartesian coordinates
	Real x_now = pco->x1v(ie+i);
	Real y_now = pco->x2v(j);
	Real z_now = pco->x3v(k);

	//Then what to do?
	//we can try saying that if the distance to injection point is within a range rinj_thresh
	Real d_inj = sqrt((x_now - x_inj)*(x_now - x_inj)+
			  (y_now - y_inj)*(y_now - y_inj)+
			  (z_now - z_inj)*(z_now - z_inj));

	//printf("in boundary, d_inj:%g, rinj_thresh:%g\n", d_inj, rinj_thresh);
	if (d_inj < rinj_thresh){ //injection cells
          
          //Real dist = sqrt(r_coord*r_coord + r_inj*r_inj - 2.0*r_coord*r_inj*(sin(th_coord)*sin(th_inj)*cos(ph_coord-ph_inj) + 
	  //								    cos(th_coord)*cos(th_inj))); //real distance from injection point
                                                                                                        // (fc) to injection cell (vc)
	  Real dist = d_inj;
	  //printf("boundary: x:%g, y:%g, z:%g, d_inj:%g, x_inj:%g, y_inj:%g, z_inj:%g\n", x_now, y_now, z_now, d_inj, x_inj, y_inj, z_inj);

	   prim(IDN,k,j,ie+i) = local_dens_now_*exp(-(dist*dist)/(r0*r0));
	   prim(IVX,k,j,ie+i) = local_vx; 
	   prim(IVY,k,j,ie+i) = local_vy;
	   prim(IVZ,k,j,ie+i) = local_vz;
	   // prim(IDN,k,j,ie+i+1) = local_dens_now_*exp(-(dist*dist)/(r0*r0));
	   // prim(IVX,k,j,ie+i+1) = local_vx; 
	   // prim(IVY,k,j,ie+i+1) = local_vy;
	   // prim(IVZ,k,j,ie+i+1) = local_vz;
	   if (NON_BAROTROPIC_EOS){
	     prim(IPR,k,j,ie+i) = prim(IDN,k,j,ie+i)*(temp_stream/temp_unit);
	     //prim(IPR,k,j,ie+i+1) = prim(IDN,k,j,ie+i+1)*(temp_stream/temp_unit);
	   }
		   
	 
	}else{//one-direction outflow
	  prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
	  prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
	  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);//pco->x1v(ie+i)*sqrt(GM1/pow(pco->x1v(ie+i),3));
	  prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
	  if (NON_BAROTROPIC_EOS)
	    prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie); //is a temperature floor necessary?
	  
	}
	//maybe label these cells so it doesn't have to select them again in b fields? 


      }//end R
    }//end theta
  }//end Phi
  
  if (MAGNETIC_FIELDS_ENABLED) {

    //(TODO) didn't check MHD version
    Real ph_inj = y_inj;
    Real th_inj = PI/2.0;
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=ngh; ++i){
	Real r_inj;
	Real radiusac;
	if (i==1){
	  r_inj = x_inj;
	  //r_inj=pco->x1f(ie+2);
	  radiusac=pco->x1f(ie+2);
	} else {
	  r_inj = x_inj;
	  //r_inj=pco->x1f(ie+3);
	  radiusac=pco->x1f(ie+3);
	}
	Real x_inj = r_inj * cos(ph_inj) * sin(th_inj);
	Real y_inj = r_inj * sin(ph_inj) * sin(th_inj);
	Real z_inj = r_inj * cos(th_inj);
        Real radius=pco->x1v(ie+i);
        Real theta=pco->x2v(j);
        Real phi=pco->x3v(k);
        Real x_now=radiusac*sin(theta)*cos(phi);
        Real y_now=radiusac*sin(theta)*sin(phi);
        Real z_now=radiusac*cos(theta);
        Real dinj= sqrt((x_now - x_inj)*(x_now - x_inj)+
			  (y_now- y_inj)*(y_now - y_inj)+
			  (z_now - z_inj)*(z_now - z_inj));
	if (dinj < rinj_thresh){ //injection cells
          Real dist = sqrt(radius*radius + r_inj*r_inj - 2.0*radius*r_inj*(sin(theta)*sin(th_inj)*cos(phi-ph_inj) + 
									   cos(theta)*cos(th_inj)));
          Real atheta = b0*r0*exp(-(dist*dist)/(r0*r0));
	  Real ddist2dphi = 2.0*radius*r_inj*sin(theta)*sin(th_inj)*sin(phi-ph_inj);
	  b.x1f(k,j,ie+i+1) = 1.0/(radius*sin(theta)) * atheta * (1.0/(r0*r0)) * ddist2dphi;
	}
	else{
	  b.x1f(k,j,ie+i+1) = b.x1f(k,j,ie+1); //outflow BC for non-injection cells
	}
      }}}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma simd
      for(int i=1; i<=ngh; ++i){
        Real radius=pco->x1v(ie+i);
        Real theta=pco->x2v(j); 
        Real phi=pco->x3v(k);
	Real r_inj;
	Real radiusac;
	if (i==1){
	  r_inj = x_inj;
	  //r_inj=pco->x1f(ie+2);
	  radiusac=pco->x1f(ie+2);
	} else {
	  r_inj = x_inject2;
	  //r_inj=pco->x1f(ie+3);
	  radiusac=pco->x1f(ie+3);
	}
        Real x_now=radiusac*sin(theta)*cos(phi);
        Real y_now=radiusac*sin(theta)*sin(phi);
        Real z_now=radiusac*cos(theta);
	Real x_inj = r_inj * cos(ph_inj) * sin(th_inj);
	Real y_inj = r_inj * sin(ph_inj) * sin(th_inj);
	Real z_inj = r_inj * cos(th_inj);
        Real dinj= sqrt((x_now - x_inj)*(x_now - x_inj)+
			  (y_now- y_inj)*(y_now - y_inj)+
			  (z_now - z_inj)*(z_now - z_inj));
 
	if (dinj < rinj_thresh){ //injection cells
	  b.x2f(k,j,ie+i) = 0.0;
	}else{
	  b.x2f(k,j,ie+i) = b.x2f(k,j,ie); //outflow 
	}
      }}}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=ngh; ++i){
        Real radius=pco->x1v(ie+i); 
        Real theta=pco->x2v(j); 
        Real phi=pco->x3v(k);
	Real r_inj;
	Real radiusac;
	if (i==1){
	  //r_inj = x_inject;
	  r_inj=pco->x1f(ie+2);
	  radiusac=pco->x1f(ie+2);
	} else {
	  //r_inj = x_inject2;
	  r_inj=pco->x1f(ie+3);
	  radiusac=pco->x1f(ie+3);
	}
        Real x_now=radiusac*sin(theta)*cos(phi);
        Real y_now=radiusac*sin(theta)*sin(phi);
        Real z_now=radiusac*cos(theta);
	Real x_inj = r_inj * cos(ph_inj) * sin(th_inj);
	Real y_inj = r_inj * sin(ph_inj) * sin(th_inj);
	Real z_inj = r_inj * cos(th_inj);
        Real dinj= sqrt((x_now - x_inj)*(x_now - x_inj)+
			  (y_now- y_inj)*(y_now - y_inj)+
			  (z_now - z_inj)*(z_now - z_inj));
	if (dinj < rinj_thresh){ //injection cells
	  Real dist = sqrt(radius*radius + r_inj*r_inj - 2.0*radius*r_inj*(sin(theta)*sin(th_inj)*cos(phi-ph_inj) + 
									   cos(theta)*cos(th_inj)));
	  Real atheta = b0*r0*exp(-(dist*dist)/(r0*r0));
	  Real ddist2dr = 2.0*radius - 2*r_inj*(sin(theta)*sin(th_inj)*cos(phi-ph_inj) + cos(theta)*cos(th_inj));
	  b.x3f(k,j,ie+i) = atheta/radius - atheta*(1.0/(r0*r0))*ddist2dr;
	}else{
	  b.x3f(k,j,ie+i) = b.x3f(k,j,ie); //outflow
	}
      }}}
  }

}

void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  Real local_dens = 0.1;
  Real local_vr = 0.0;
  Real local_vphi = 0.036;
  Real local_press = 1.0e-4;
  //printf("hydro ix1\n");
  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=(NGHOST); ++i) {//R
	  prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
	  prim(IVX,k,j,is-i) = std::min(0.0, prim(IVX,k,j,is));
	  prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
	  prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,j,is-i) = prim(IPR,k,j,is);
	  }

      }//end R
    }//end theta
  }//end Phi
  
}

void HydroOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  Real local_dens = 0.1;
  Real local_vr = 0.0;
  Real local_vphi = 0.036;
  Real local_press = 1.0e-4;
  //printf("hydro ix1\n");
  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=(NGHOST); ++i) {//R
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

void HydroInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=1; j<=(NGHOST); ++j) {//theta
      for (int i=is; i<=ie; ++i) {//R
	prim(IDN,k,js-j,i) = prim(IDN,k,js,i); 
	prim(IVX,k,js-j,i) = prim(IVX,k,js,i); 
	prim(IVZ,k,js-j,i) = prim(IVZ,k,js,i);
	prim(IVY,k,js-j,i) = std::min(0.0, prim(IVY,k,js,i));
	if (NON_BAROTROPIC_EOS){
	  prim(IPR,k,js-j,i) = prim(IPR,k,js,i);
	}

      }//end R
    }//end theta
  }//end Phi

}

void HydroOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=NGHOST; ++j) {
      for (int i=is; i<=ie; ++i) {
	  prim(IDN,k,je+j,i) = prim(IDN,k,je,i);
	  prim(IVX,k,je+j,i) = prim(IVX,k,je,i);
	  prim(IVZ,k,je+j,i) = prim(IVZ,k,je,i);
	  prim(IVY,k,je+j,i) = std::max(0.0, prim(IVY,k,je,i));
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,je+j,i) = prim(IPR,k,je,i);
	  }

      }
    }
  }
  
}

void HydroInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
	prim(IDN,ks-k,j,i) = prim(IDN,ks,j,i); 
	prim(IVX,ks-k,j,i) = prim(IVX,ks,j,i); 
	prim(IVZ,ks-k,j,i) = std::min(0.0, prim(IVZ,ks,j,i));
	prim(IVY,ks-k,j,i) = prim(IVY,ks,j,i);
	if (NON_BAROTROPIC_EOS){
	  prim(IPR,ks-k,j,i) = prim(IPR,ks,j,i);
	}

      }
    }
  }

}

void HydroOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
	prim(IDN,ke+k,j,i) = prim(IDN,ke,j,i); 
	prim(IVX,ke+k,j,i) = prim(IVX,ke,j,i); 
	prim(IVZ,ke+k,j,i) = std::min(0.0, prim(IVZ,ke,j,i));
	prim(IVY,ke+k,j,i) = prim(IVY,ke,j,i);
	if (NON_BAROTROPIC_EOS){
	  prim(IPR,ke+k,j,i) = prim(IPR,ke,j,i);
	}

      }
    }
  }

}

void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
  // copy radiation variables into ghost zones,
  // only allow outflow
  int &nang = prad->nang; // angles per octant
  int &nfreq = prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=ngh; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir<0, inward: mu_dir>0
      Real mu_dir = prad->mu(0,k,j,is-i,ang);
      if (mu_dir < 0.0){
	ir(k,j,is-i,ang) = ir(k,j,is,ang);
      }else{
	ir(k,j,is-i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}
}

void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  int &nang = prad->nang; // angles per octant
  int &nfreq = prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=ngh; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir>0, inward: mu_dir<0
      Real mu_dir = prad->mu(0,k,j,ie+i,ang);
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

void RadInnerX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
  // copy radiation variables into ghost zones,
  // only allow outflow
  //printf("rad ix2\n");
  int &nang = prad->nang; // angles per octant
  int &nfreq = prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=ngh; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir<0, inward: mu_dir>0
      Real mu_dir = prad->mu(1,k,js-j,i,ang);
      if (mu_dir < 0.0){
	ir(k,js-j,i,ang) = ir(k,js,i,ang);
      }else{
	ir(k,js-j,i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}
}

void RadOuterX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  int &nang = prad->nang; // angles per octant
  int &nfreq = prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=ngh; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir>0, inward: mu_dir<0
      Real mu_dir = prad->mu(1,k,je+j,i,ang);
      if (mu_dir > 0.0){
	ir(k,je+j,i,ang) = ir(k,je,i,ang);
      }else{
       	ir(k,je+j,i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}

  return;

}



void RadInnerX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){
  int &nang = prad->nang; // angles per octant
  int &nfreq = prad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir<0, inward: mu_dir>0
      Real mu_dir = prad->mu(2,ks-k,j,i,ang);
      if (mu_dir < 0.0){
	ir(ks-k,j,i,ang) = ir(ks,j,i,ang);
      }else{
	ir(ks-k,j,i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}
}

void RadOuterX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  int &nang = prad->nang; // angles per octant
  int &nfreq = prad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      //if directs outwards: mu_dir>0, inward: mu_dir<0
      Real mu_dir = prad->mu(2,ke+k,j,i,ang);
      if (mu_dir > 0.0){
	ir(ke+k,j,i,ang) = ir(ke,j,i,ang);
      }else{
       	ir(ke+k,j,i,ang) = 0.0;
      }
    }// end n
  }// end ifr
  }}}

  return;

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

    //if replace_low_dens_ff is true, replace the low density opacity with free-free
    if (replace_low_dens_ff==1){
      if (rho <= rho_cut_ff){
	kappa_ross = kappa_es + kappa_ff_ross(tgas, rho);
	kappa_planck = kappa_ff_planck(tgas, rho);
      }
    }
    
}


void GetCombineOpacity(MeshBlock *pmb, AthenaArray<Real> &prim){

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

  Real kappaa = 0.0;

  for (int k=kl; k<=ku; k++){
    for (int j=jl; j<=ju; j++){
      for (int i=il; i<=iu; i++){
	Real rho = prim(IDN,k,j,i);
	Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i); //std::max(prim(IPR,k,j,i)/prim(IDN,k,j,i), tfloor);
	Real rho_cgs = rho*rho_unit;
	Real temp_cgs = temp*temp_unit;
	
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
	
	//if (temp_cgs > 5.0e6){
	//  printf("temp:%g, rho:%g, kappa_p:%g, kappa_r+kappaes:%g, r:%g, ph:%g\n", temp_cgs, rho_cgs, kappa_planck, kappa_s+kappaa, pmb->pcoord->x1v(i), pmb->pcoord->x3v(k));
	//}

	//one frequency
	pnrrad->sigma_s(k,j,i,0) = kappa_s * rho * rho_unit * l_unit; //scatter
	pnrrad->sigma_a(k,j,i,0) = kappaa * rho * rho_unit * l_unit; //rosseland mean
	pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //planck mean
        pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;//planck mean
      
      }//end i
    }//end j
  }//end k

}



Real massflux_AInj_x1(MeshBlock *pmb, int iout){
  Real area = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      for (int i=is; i<=ie; i++){
	area += pmb->ruser_meshblock_data[1](0,k,j,i);
      }
    }
  }
    
  return area;
}

Real massflux_Inj_x1(MeshBlock *pmb, int iout){
  Real flux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      for (int i=is; i<=ie; i++){
        flux += pmb->ruser_meshblock_data[1](1,k,j,i);
      }
    }
  }
    
  return flux;
}

Real massflux_AInj_x3(MeshBlock *pmb, int iout){
  Real area = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
    
  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      for (int i=is; i<=ie; i++){
	area += pmb->ruser_meshblock_data[1](2,k,j,i);
      }
    }
  }
    
  return area;
}

Real massflux_Inj_x3(MeshBlock *pmb, int iout){
  Real flux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      for (int i=is; i<=ie; i++){
        flux += pmb->ruser_meshblock_data[1](3,k,j,i);
      }
    }
  }
    
  return flux;
}

Real massfluxox1(MeshBlock *pmb, int iout){
  Real massflux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1;
  //face1.NewAthenaArray((ie-is)+2*NGHOST+1);

  AthenaArray<Real> x1flux = pmb->phydro->flux[X1DIR];
  //bool boundaryflag;
  //boundaryflag = pmb->pbval->apply_bndry_fn_[BoundaryFace::outer_x1];
  //printf("outer boundary flag:%d\n", boundaryflag);
  if (pmb->pbval->get_apply_bndry_fn_flag(BoundaryFace::outer_x1)){
    //printf("here\n");
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
	//pmb->pcoord->Face1Area(k , j, is, ie, face1);
	for (int i=1; i<=1; i++){
      massflux += pmb->pcoord->GetFace1Area(k,j,ie+i)*x1flux(IDN,k,j,ie+i); //x1flux(IDN) is the density flux, multiply by volume to get mass
          //printf("face:%g, x1flux:%g\n", face1(ie), x1flux(0, k, j, ie));
	}

      }
    }
  }

  //printf("outter mass flux:%g\n", massflux);

  //face1.DeleteAthenaArray();
  x1flux.DeleteAthenaArray();

  return massflux;


}

Real massfluxix1(MeshBlock *pmb, int iout){
  Real massflux = 0.0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1;
  face1.NewAthenaArray((ie-is)+2*NGHOST+2);

  AthenaArray<Real> x1flux = pmb->phydro->flux[X1DIR];
  bool boundaryflag = pmb->pbval->get_apply_bndry_fn_flag(BoundaryFace::inner_x1);
  //printf("inner boundary flag:%d\n", boundaryflag);
  if (pmb->pbval->get_apply_bndry_fn_flag(BoundaryFace::inner_x1)){
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
	pmb->pcoord->Face1Area(k , j, is, ie, face1);
	for (int i=is; i<=is; i++){
	  massflux += face1(is)*x1flux(0,k,j,is);//x1flux(0) is the density flux, multiply by volume to get mass
        }

      }
    }
  }

  //printf("inner mass flux:%g\n", massflux);
  face1.DeleteAthenaArray();
  x1flux.DeleteAthenaArray();

  return massflux;


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


//input code unit, output cgs
Real kappa_ff_nu(Real nu, Real temp, Real rho){

  Real h_planck = 6.626196e-27 ;
  Real evtohz = 2.41838e14;
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real m_p = m_p = 1.6726e-24;
  Real k_B = 1.3807e-16;
  
  Real  gff = 1.0;
  Real  z = 1.0;

  Real  he_adbund = 0.04;
  Real  nh = rho_cgs/m_p/(1.0 + 4.0*he_adbund);
  Real  nhe = nh*he_adbund;
  Real  ne = nh + 2.0*nhe;
  Real  n_rho = rho_cgs/m_p/0.62;

  Real  e_ff = 3.7e8 * pow(temp_cgs, -0.5) * pow(z, 2) * pow(n_rho, 2) * pow(nu, -3) * (1.0 - exp(-h_planck*nu/k_B/temp_cgs)) * gff;

  return std::max(e_ff/rho_cgs, 1.0e-10); //add a floor
}

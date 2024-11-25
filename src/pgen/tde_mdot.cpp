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

//Injection point
static Real rad_inject, phi_inject, local_dens, local_vr, local_vphi, local_press;
static Real f_dr, f_dth, f_dph;
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

//refinement
static Real rad_thresh, ph_thresh, th_thresh;
static int blocksizex1_base, blocksizex2_base, blocksizex3_base;
static Real delta_theta_base;
static Real rad_coarse_thresh, ph_coarse_thresh, th_coarse_thresh;

//stream injection boundary
static Real r0, inj_thresh, b0;
static Real rad_inject2; //rad_inject,  phi_inject;

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
void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void RadOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void PointMassPotential(MeshBlock *pmb, const Real time, const Real dt,
			const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void GeneralNewtonianPotential(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                              const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
			      AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

void opacity(MeshBlock *pmb, AthenaArray<Real> &prim);
void opalopacity(MeshBlock *pmb, AthenaArray<Real> &prim);
Real kappa_ff(Real temp, Real rho);
void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck);


//user history 
Real massflux_AInj_x1(MeshBlock *pmb, int iout);
Real massflux_AInj_x3(MeshBlock *pmb, int iout);
Real massflux_Inj_x1(MeshBlock *pmb, int iout);
Real massflux_Inj_x3(MeshBlock *pmb, int iout);

Real massfluxix1(MeshBlock *pmb, int iout);
Real massfluxox1(MeshBlock *pmb, int iout);

//AMR condition
int RefinementCondition(MeshBlock *pmb);
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
  int theta_max = pin->GetReal("mesh","x2max");
  int theta_min = pin->GetReal("mesh","x2min");
  int theta_base = pin->GetInteger("mesh","nx2");
  delta_theta_base = (theta_max - theta_min)/theta_base * blocksizex2_base;

  GM1 = pin->GetOrAddReal("problem","GM1",1.0);
  kappa_es = pin->GetReal("problem", "kappa_es");

  temp_unit = pin->GetReal("problem", "temp_unit");
  l_unit = pin->GetReal("problem", "l_unit");
  rho_unit = pin->GetReal("problem", "rho_unit");

  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.001);
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 0.001);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 0.001);

  user_dt = pin->GetOrAddReal("problem", "user_dt", 1.0e-6);
  //Initialize the injection point
  rad_inject = pin->GetReal("problem", "rad_inject");
  phi_inject = pin->GetReal("problem", "phi_inject");
  local_dens = pin->GetReal("problem", "local_dens");
  local_vr = pin->GetReal("problem", "local_vr");
  local_vphi = pin->GetReal("problem", "local_vphi");
  //local_press = pin->GetReal("problem", "local_press");

  f_dr = pin->GetOrAddInteger("problem", "f_dr", 1);  
  f_dth = pin->GetOrAddInteger("problem", "f_dth", 1);
  f_dph = pin->GetOrAddInteger("problem", "f_dph", 1);  

  rad_thresh = pin->GetOrAddReal("problem", "rad_thresh", 0.1);  
  ph_thresh = pin->GetOrAddReal("problem", "ph_thresh", 0.01);
  th_thresh = pin->GetOrAddReal("problem", "th_thresh", 0.01);
  rad_coarse_thresh = pin->GetOrAddReal("problem", "rad_coarse_thresh", 0.1);
  ph_coarse_thresh = pin->GetOrAddReal("problem", "ph_coarse_thresh", 0.01);
  th_coarse_thresh = pin->GetOrAddReal("problem", "th_coarse_thresh", 0.01);

  //stream injection boundary
  b0 = pin->GetOrAddReal("problem", "b0", 0.1);
  r0 = pin->GetOrAddReal("problem", "r0", 7.53);
  inj_thresh = pin->GetOrAddReal("problem", "inj_thresh", 7.53);
  //rad_inject = pin->GetOrAddReal("problem", "rad_inject", 400.0);
  rad_inject2 = pin->GetOrAddReal("problem", "rad_inject2", 407.886);
  //phi_inject = pin->GetOrAddReal("problem", "phi_inj", 3.14);

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

    //opacity tables
    opacitytable.NewAthenaArray(212,46);
    planckopacity.NewAthenaArray(138,37);
    logttable.NewAthenaArray(212);
    logrhottable.NewAthenaArray(46);
    logttable_planck.NewAthenaArray(138);
    logrhottable_planck.NewAthenaArray(37);
    
    // read in the opacity table
    FILE *fkappa, *flogt, *flogrhot, *fplanck, *flogt_planck, *flogrhot_planck;
      
    if ( (fkappa=fopen("./aveopacity_combined.txt","r"))==NULL )
    {
      printf("Open input file error aveopacity_combined");
      return;
    }

    if ( (fplanck=fopen("./PlanckOpacity.txt","r"))==NULL )
    {
      printf("Open input file error PlanckOpacity");
      return;
    }

    if ( (flogt=fopen("./logT.txt","r"))==NULL )
    {
      printf("Open input file error logT");
      return;
    }

    if ( (flogrhot=fopen("./logRhoT.txt","r"))==NULL )
    {
      printf("Open input file error logRhoT");
      return;
    }

    if ( (flogt_planck=fopen("./logT_planck.txt","r"))==NULL )
    {
      printf("Open input file error logT_planck");
      return;
    }

    if ( (flogrhot_planck=fopen("./logRhoT_planck.txt","r"))==NULL )
    {
      printf("Open input file error logRhoT_planck");
      return;
    }

    for(int j=0; j<212; j++){
      for(int i=0; i<46; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
      }
    }

    for(int j=0; j<138; j++){
      for(int i=0; i<37; i++){
          fscanf(fplanck,"%lf",&(planckopacity(j,i)));
      }
     }


    for(int i=0; i<46; i++){
      fscanf(flogrhot,"%lf",&(logrhottable(i)));
    }

    for(int i=0; i<212; i++){
      fscanf(flogt,"%lf",&(logttable(i)));
    }

    for(int i=0; i<37; i++){
      fscanf(flogrhot_planck,"%lf",&(logrhottable_planck(i)));
    }

    for(int i=0; i<138; i++){
      fscanf(flogt_planck,"%lf",&(logttable_planck(i)));
    }

    fclose(fkappa);
    fclose(flogt);
    fclose(flogrhot);
    fclose(fplanck);
    fclose(flogt_planck);
    fclose(flogrhot_planck);


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
  

  EnrollUserExplicitSourceFunction(GeneralNewtonianPotential);//(PointMassPotential);

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

  // Real r_l= pcoord->x1f(is); 
  // Real r_u = pcoord->x1f(ie);
  // Real th_l = pcoord->x2f(js);
  // Real th_u = pcoord->x2f(je);
  // Real ph_l = pcoord->x3f(ks);
  // Real ph_u = pcoord->x3f(ke);

  // Real th_edge_dis = std::max(fabs(th_l-PI/2.0), fabs(th_u-PI/2.0));
  // int ph_peri_half = 0;
  // if (ph_l>=PI/2.0 &&  ph_u<1.5*PI){
  //   ph_peri_half = 1;
  // }

  // Real max_dis_th = 0.0;
  // Real max_dis_phi = 0.0;
  // Real max_dis_phi_p = 0.0;
  // Real max_dis_phi_m = 0.0;
  // Real max_dth, max_dphi = 0.0;

  for(int k=0; k<blocksizex3; k++){
   for(int j=0; j<blocksizex2; j++){
     for(int i=0; i<blocksizex1; i++){

       //initialize the storing flux and area
       ruser_meshblock_data[1](0,k,j,i) = 0.0;
       ruser_meshblock_data[1](1,k,j,i) = 0.0;
       ruser_meshblock_data[1](2,k,j,i) = 0.0;
       ruser_meshblock_data[1](3,k,j,i) = 0.0;
       
       Real phi_coord = pcoord->x3v(k);
       Real theta_coord = pcoord->x2v(j);
       Real rad = pcoord->x1v(i);

       Real dr = (pcoord->dx1f(i-1)+pcoord->dx1f(i+1))/2.0;
       Real dphi = (pcoord->dx3f(k-1)+pcoord->dx3f(k+1))/2.0;
       Real dth = (pcoord->dx2f(j-1)+pcoord->dx2f(j+1))/2.0;	  

       //injection point
       if (fabs(rad-rad_inject) <= rad_coarse_thresh){
       	 if (fabs(phi_coord-phi_inject) <= ph_coarse_thresh){// if within L1 point
       	   if (fabs(theta_coord - PI/2.0)<= th_coarse_thresh){
       	     //printf("coarse\n");
		//printf("r:%g, ph:%g, th:%g, dis_r:%g, dis_phi:%g, dis_th:%g, i:%d, gid:%d, dr:%g, dph:%g, dth:%g\n", rad, phi_coord, theta_coord, fabs(rad-rad_inject), fabs(phi_coord - phi_inject), fabs(theta_coord - PI/2.0), i, gid, pcoord->dx1f(i), pcoord->dx3f(k), pcoord->dx2f(j));
       	     iuser_meshblock_data[0](0) = 1;
       	   }
       	 }
       }//end injection point

       // Real phi_now = pcoord->x3v(k);
       // Real th_now = pcoord->x2v(j);
       // Real r_now = pcoord->x1v(i);
       // Real r_l = pcoord->x1v(is);
       // Real r_u = pcoord->x1v(ie);

       // Real phi_fit = orbit_fit_phi(r_now);
       // Real dis_phi = fabs(phi_now - phi_fit);
       // Real dis_th = fabs(th_now - PI/2.0);

       // //the positive and negative distantce
       // Real dis_phi_p = 0.0;
       // if (phi_now > phi_fit){
       // 	 dis_phi_p = fabs(phi_now - phi_fit);
       // }
       // Real dis_phi_m = 0.0;
       // if (phi_now <= phi_fit){
       // 	 dis_phi_m = fabs(phi_now - phi_fit);
       // }


       // 	//find the max dr/dphi in this meshblock to r
       // 	if (dphi > max_dphi){
       // 	  max_dphi = dphi;
       // 	}
       // 	if (dth > max_dth){
       // 	  max_dth = dth;
       // 	}

       // 	//find the max dr/dth in this meshblock
       // 	if (dis_phi > max_dis_phi){
       // 	  max_dis_phi = dis_phi;
       // 	}
       // 	if (dis_phi_m > max_dis_phi_m){
       // 	  max_dis_phi_m = dis_phi_m;
       // 	}
       // 	if (dis_phi_p > max_dis_phi_p){
       // 	  max_dis_phi_p = dis_phi_p;
       // 	}
       // 	if (dis_th > max_dis_th){
       // 	  max_dis_th = dis_th;
       // 	}

     }//k
   }//j
  }//i

  // //level 3
  // if (max_dis_phi <= 0.68 && r_u<50.0 && r_l>22.0 && max_dis_th<0.54){ //inner part of stream
  //   iuser_meshblock_data[0](2) = 1;
  //   //printf("refine 1, r_u:%g, r_l:%g, max_dis_phi:%g\n", r_u, r_l, max_dis_phi);
  // }else if (max_dis_phi <= 1.2 && r_u<22.0 && r_l>11.0 && max_dis_th<0.54){
  //   iuser_meshblock_data[0](2) = 1;
  //   printf("refine 2, r_u:%g, r_l:%g, max_dis_phi:%g\n", r_u, r_l, max_dis_phi);
  // }else if (r_u<14.0 && r_l>2.0 && max_dis_th<0.42){//outer part of stream
  //   iuser_meshblock_data[0](2) = 1;
  // }

 //level 3, different upper phi and lower phi to cover collision region
  // if (max_dis_phi_p <= 0.68 && max_dis_phi_p > 0.0 && r_u<50.0 && r_l>11.0 && max_dis_th<0.54){ //inner part of stream
  //   iuser_meshblock_data[0](2) = 1;
  // }else if (max_dis_phi_m <= 0.68 && max_dis_phi_m >0.0 && r_u<50.0 && r_l>22.0 && max_dis_th<0.54){
  //   iuser_meshblock_data[0](2) = 1;
  // }else if (max_dis_phi_m <= 1.4 && max_dis_phi_m >0.0 && r_u<22.0 && r_l>11.0 && max_dis_th<0.54){
  //   iuser_meshblock_data[0](2) = 1;
  // }else if (r_u<14.0 && r_l>2.0 && max_dis_th<0.42){//outer part of stream
  //   iuser_meshblock_data[0](2) = 1;
  // }
  
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
  Real max_r = 0.0;
  Real max_th = 0.0;
  Real max_ph = 0.0;
  Real max_dens = 0.0;
  
  int current_level = pmb->pmy_mesh->get_current_level();//pmb->loc.level - pmb->pmy_mesh->root_level;

  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
	Real epsr = (std::abs(w(IDN,k,j,i+1) - 2.0*w(IDN,k,j,i) + w(IDN,k,j,i-1))
		     + std::abs(w(IDN,k,j+1,i) - 2.0*w(IDN,k,j,i) + w(IDN,k,j-1,i))
		     + std::abs(w(IDN,k+1,j,i) - 2.0*w(IDN,k,j,i) + w(IDN,k-1,j,i)))/w(IDN,k,j,i);
	// Real epsp = (std::abs(w(IPR,k,j,i+1) - 2.0*w(IPR,k,j,i) + w(IPR,k,j,i-1))
	//             + std::abs(w(IPR,k,j+1,i) - 2.0*w(IPR,k,j,i) + w(IPR,k,j-1,i)))
	//             /w(IPR,k,j,i);
	// Real eps = std::max(epsr, epsp);
	if (epsr >= maxeps){
	  max_r = pmb->pcoord->x1v(i);
	  max_th = pmb->pcoord->x2v(j);
	  max_ph = pmb->pcoord->x3v(k);
	}
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
  //  printf("gid:%d, maxeps:%g, maxdens:%g, current level:%d, pmb->root:%d, pmb->cuurent:%d\n", pmb->gid, maxeps, max_dens, current_level, pmb->pmy_mesh->root_level, pmb->pmy_mesh->current_level);

  //if ((maxeps>1.0 && max_dens>1.0e-5)||injection_flag ==1) return 1; //run to t=0.3
  //if (maxeps<1.0e-3) return -1;//run to t=0.3
  
  Real r_min = pmb->pcoord->x1f(pmb->is);
  Real dis_th_l = std::fabs(pmb->pcoord->x2f(pmb->js) - PI/2.0);
  Real dis_th_r = std::fabs(pmb->pcoord->x2f(pmb->je) - PI/2.0);
  Real dis_th_min = std::min(dis_th_l, dis_th_r);

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

  //run to t=4.33
  //if ((maxeps>1.0 && max_dens>5.0e-2) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.1 || r_min<4.8 || max_dens<1.0e-4) return -1; //run to t=5.

  //run to t=4.7
  //if ((maxeps>1.0 && max_dens>5.0e-2) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.1 || r_min<4.8 || max_dens<1.0e-3) return -1; //run to t=5.

  //run to t=5.0
  //if ((maxeps>0.1 && max_dens>1.0e-1) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.01 || r_min<4.8 || max_dens<1.0e-3) return -1; //run to t=5.

  //run to t-5.39
  //if ((maxeps>0.1 && max_dens>2.0e-1) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.01 || r_min<4.8 || max_dens<2.0e-3) return -1; //run to t=5.

  //run to t=5.4
  //if ((maxeps>0.1 && max_dens>3.0e-1) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.01 || r_min<4.8 || max_dens<3.0e-3) return -1; //run to t=5

  //run to t=5.87
  //if ((maxeps>0.1 && max_dens>5.0e-1) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.01 || r_min<4.8 || max_dens<4.0e-3) return -1; //run to t=5

  //run to t=6.07
  //if ((maxeps>0.2 && max_dens>6.0e-1) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.02 || r_min<4.8 || max_dens<5.0e-3) return -1; //run to t=5

  //run to t=6.42
  //if ((maxeps>0.5 && max_dens>8.0e-1) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.05 || r_min<4.8 || max_dens<8.0e-3) return -1; 

  //run to t=10.05
  //if ((maxeps>1.0 && max_dens>1.0) && r_min>6.0 && dis_th_min<0.088) return 1;
  //if (maxeps<0.05 || r_min<4.8 || max_dens<1.0e-2) return -1;

  if ((maxeps>5.0 && max_dens>3.0) && r_min>6.0 && dis_th_min<0.088) return 1;
  if (maxeps<0.05 || r_min<4.8 || max_dens<5.0e-2) return -1;

  //only used in trying restart a .rst file with numlevel=6 -> numlevel=5
  //if (pmb->pmy_mesh->current_level>8 && (maxeps>1.0 && max_dens>5.0e-4) && r_min>3.0 && dis_th_min<0.088) return -1;
  //if (pmb->pmy_mesh->current_level>8) return -1;
  
   //return derefine to polar region?
  //distance to theta=0.0 for max/min theta of current meshblock
  Real dis_topole1_l = std::fabs(pmb->pcoord->x2f(pmb->js) - 0.0);
  Real dis_topole1_r = std::fabs(pmb->pcoord->x2f(pmb->je) - 0.0); 
  Real dis_topole1_min = std::min(dis_topole1_l, dis_topole1_r);

  //distance to theta=pi for max/min theta of current meshblock
  Real dis_topole2_l = std::fabs(pmb->pcoord->x2f(pmb->js) - PI);
  Real dis_topole2_r = std::fabs(pmb->pcoord->x2f(pmb->je) - PI);
  Real dis_topole2_min = std::min(dis_topole2_l, dis_topole2_r);
  
  //PI/16.0 is a randomly choosed threshold, should test different values
  //if ((dis_topole1_min<PI/16.0) || (dis_topole2_min<PI/16.0)){  //stay off the poles (>11deg)
  if ((dis_topole1_min<delta_theta_base) || (dis_topole2_min<delta_theta_base)){  //stay off the poles (>11deg)
    // Real current_level = pmb->current_level(gid);
    // if current_level != root_level; return -1;
    return 0; //-1
  }

  return 0;

}

void Mesh::UserWorkInLoop(){
  MeshBlock *pmb = my_blocks(0);
  for(int nb=0; nb<nblocal; ++nb){ //loop over meshblocks on the same core
    pmb = my_blocks(nb);
    
    Hydro *phydro = pmb->phydro;
    Coordinates *pcoord = pmb->pcoord;
    int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
    
    if (pmb->iuser_meshblock_data[0](0) !=0 ){
    //get current loacla dens
    Real t_current = pmb->pmy_mesh->time;
    Real local_dens_now_ = GetMdot(pmb, t_current);
    //for rho=1
    //local_dens_now_ = 1.0;
    for(int k=ks; k<=ke; k++){
      for(int i=is; i<=ie; i++){
	for(int j=js; j<=je; j++){
	  for (int n=0; n<(NHYDRO);n++){
	    if (phydro->u(n,k,j,i) != phydro->u(n,k,j,i)){
	      printf("block: %d, n: %d ,k: %d,j: %d,i: %d\n", pmb->gid,n,k,j,i);
	      printf("x1v: %g, x2v:%g, x3v:%g\n",pmb->pcoord->x1v(i), pmb->pcoord->x2v(j),pmb->pcoord->x3v(k));
	      abort();
	    }	  
	  }//end NHYDRO
	  //inject stream
	  Real phi_coord = pmb->pcoord->x3v(k);
	  Real theta_coord = pmb->pcoord->x2v(j);
	  Real rad = pmb->pcoord->x1v(i);

	  Real dr = (pcoord->dx1f(i-1)+pcoord->dx1f(i+1))/2.0;
	  Real dphi = (pcoord->dx3f(k-1)+pcoord->dx3f(k+1))/2.0;
	  Real dth = (pcoord->dx2f(j-1)+pcoord->dx2f(j+1))/2.0;

	  //get current loacla dens
	  //Real t_current = pmb->pmy_mesh->time;
	  //Real local_dens_now_ = GetMdot(pmb, t_current);
 
	  // store user meshblock datas for diagnostic
	  AthenaArray<Real> &x1flux = phydro->flux[X1DIR];
	  AthenaArray<Real> &x2flux = phydro->flux[X2DIR];
	  AthenaArray<Real> &x3flux = phydro->flux[X3DIR];

	  /* change to boundary injection */
	  // if (fabs(rad-rad_inject) <= rad_thresh){
	  //   if (fabs(phi_coord - phi_inject) <= ph_thresh){
	  //     if (fabs(theta_coord - PI/2.0)<= th_thresh){
	  // 	pmb->iuser_meshblock_data[0](1) = 1;
	  //       printf("r:%g, ph:%g, th:%g, dis_r:%g, dis_phi:%g, dis_th:%g, i:%d, gid:%d, dr:%g, dph:%g, dth:%g\n", rad, phi_coord, theta_coord, fabs(rad-rad_inject), fabs(phi_coord - phi_inject), fabs(theta_coord - PI/2.0), i, pmb->gid, pmb->pcoord->dx1f(i), pmb->pcoord->dx3f(k), pmb->pcoord->dx2f(j));
	  // 	pmb->phydro->u(IDN,k,j,i) += local_dens_now_*pmb->pmy_mesh->dt;
	  // 	pmb->phydro->u(IM1,k,j,i) += local_dens_now_*local_vr*pmb->pmy_mesh->dt;
	  // 	pmb->phydro->u(IM3,k,j,i) += local_dens_now_*local_vphi*pmb->pmy_mesh->dt;
	  // 	if (NON_BAROTROPIC_EOS){
	  // 	  Real ke_new = 0.5*(1.0/pmb->phydro->u(IDN,k,j,i))*(SQR(pmb->phydro->u(IM1,k,j,i)) + SQR(pmb->phydro->u(IM2,k,j,i)) + SQR(pmb->phydro->u(IM3,k,j,i)));
	  // 	  Real tgas_stream = temp_stream/temp_unit;
	  // 	  Real etot_new = ke_new + tgas_stream * phydro->w(IDN,k,j,i)/(pmb->peos->GetGamma()-1.0);
	  // 	  phydro->u(IEN,k,j,i) = etot_new;

	  // 	}
	  //       //store the injection area and flux
  	  // 	pmb->ruser_meshblock_data[1](0,k,j,i) = pcoord->GetFace1Area(k,j,i);
  	  // 	pmb->ruser_meshblock_data[1](1,k,j,i) = x1flux(IDN,k,j,i) * pcoord->GetFace1Area(k,j,i) - x1flux(IDN,k,j,i+1) * pcoord->GetFace1Area(k,j,i+1); 
  	  // 	pmb->ruser_meshblock_data[1](2,k,j,i) = pcoord->GetFace3Area(k,j,i);
  	  // 	pmb->ruser_meshblock_data[1](3,k,j,i) = x3flux(IDN,k,j,i) * pcoord->GetFace3Area(k,j,i) - x3flux(IDN,k+1,j,i) * pcoord->GetFace3Area(k+1,j,i); 
		
	  //     }//in radius range
	  //   }//in phi range
	  // }//in thetarange


	}//j
      }//i
    }//k
    }//injection flag
    // pmb = pmb->next;
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

	if (NR_RADIATION_ENABLED){
	  // // electron scattering opacity
	  // Real kappa_es = 0.2 * (1.0 + 0.6);
	  // Real kappaa = 0.0;
	  // Real rho = phydro->w(IDN,k,j,i);
	  // Real temp = std::max(phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i), tfloor);
	  // Real kappa, kappa_planck;
	  // //rossopacity(rho, temp, kappa, kappa_planck);
	
	  // if(kappa < kappa_es){
	  //   if(temp < 0.14){
	  //     kappaa = kappa;
	  //     kappa = 0.0;
	  //   }else{
	  //     kappaa = 0.0;
	  //   }
	  // }else{
	  //   kappaa = kappa - kappa_es;
	  //   kappa = kappa_es;
	  // }
	  // //one frequency
	  // pnrrad->sigma_s(k,j,i,0) = kappa * rho * rho_unit * l_unit; //scatter
	  // pnrrad->sigma_a(k,j,i,0) = kappaa * rho * rho_unit * l_unit; //Rossland mean
	  // pnrrad->sigma_ae(k,j,i,0) = pnrrad->sigma_a(k,j,i,0); //Rossland mean
	  // //Planck mean - Rossland mean
	  // if(kappaa < kappa_planck){
	  //   pnrrad->sigma_planck(k,j,i,0) = (kappa_planck-kappaa)*rho*rho_unit*l_unit;
	  // }else{
	  //   pnrrad->sigma_planck(k,j,i,0) = 0.0;
	  // }

	  Real rho = phydro->w(IDN,k,j,i);
	  Real temp = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);

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

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = press_init/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }//end non barotropic
      }//end i
      
      
      if (NR_RADIATION_ENABLED){
      	for (int i=ie; i>=is; i--){
      	  pnrrad->ir(k,j,i,0) = 0.0;
      	}
      }// end rad i loop
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
  while ((t_days>mdot_time_table(ig2)) & (ig2<ngroup)){
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

  if ((ig2 == ngroup) & (t_days > mdot_time_table(ig2))){
    mdot1 = mdot_table(ngroup);
    mdot2 = mdot_table(ngroup);
  }

  if ((ig2 == 0) & (t_days > mdot_time_table(0))){
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

void GeneralNewtonianPotential(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar){
  //printf("user source\n");
  //potential as Tejeda & Rosswog 13
  
  // Gravitational acceleration from orbital motion
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {

	Real r = pmb->pcoord->x1v(i);
	Real th = pmb->pcoord->x2v(j);
	Real ph = pmb->pcoord->x3v(k);

	Real sinth = sin(th);
	Real costh = cos(th);
	Real tanth = sinth/costh;

	Real vr = prim(IVX, k,j,i);//pmb->phydro->u(IM1,k,j,i)/pmb->phydro->u(IDN,k,j,i);
	Real vth = prim(IVY, k,j,i);//pmb->phydro->u(IM2,k,j,i)/pmb->phydro->u(IDN,k,j,i);
	Real vph = prim(IVZ, k,j,i);//pmb->phydro->u(IM3,k,j,i)/pmb->phydro->u(IDN,k,j,i);

        Real gm1=GM1;
	if (NR_RADIATION_ENABLED){	
            gm1 = pmb->pnrrad->crat * pmb->pnrrad->crat / 2.0;
	}
	Real rho = prim(IDN,k,j,i);

	//get r_dd, th_dd, ph_dd as in Eq 2.14, 2.15, 2.16
	//v_th = th_d * r, v_ph = ph_d * r * sinth
	Real r_dd = - gm1 * (1.0 - 1.0/r) * (1.0 - 1.0/r)/(r * r);
	r_dd += vr * vr / (r*(r-1.0));
	r_dd += (vth * vth + vph * vph) * (r - 3.0/2.0)/(r * r);

	Real th_dd = -(2.0 * vr * vth/(r * r)) * (r-3.0/2.0)/(r-1.0) + vph*vph/(r * r)*(costh/sinth);
	Real ph_dd = -(2.0 * vr * vph/(r * r)/sinth) * (r-3.0/2.0)/(r-1.0) - (2.0/tanth)*vph*vth/(r * r)/sinth;
	
	//Real acc1 = gm1/((r-1.0)*(r-1.0));
	//acceleration in all directions
	Real acc1 = -gm1 * pow(1.0/r, 2) * pow(1.0 - 1.0/r, 2) + pow(vr, 2)/r/(r-1.0) - (3.0/2.0)*pow(1.0/r, 2)*(pow(vth, 2) + pow(vph, 2));//r_dd - r*pow(vth/r, 2) - r*pow(vph/r, 2);
	Real acc2 = vr * vth /r/(r - 1.0); //r*th_dd + 2.0*vr*(vth/r) - r*sinth*costh*pow(vph/r/sinth, 2);
	Real acc3 = vr * vph /r/(r - 1.0);//r*sinth*ph_dd + 2.0*vr*(vph/r) + 2.0*r*costh*(vph/r/sinth)*(vth/r);

	//current internal energy
	Real etot_now = cons(IEN,k,j,i);
        Real ke_now = 0.5*(1.0/cons(IDN,k,j,i))*(SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
	Real ie_now = etot_now - ke_now;

	//check total energy
	Real etot_update = cons(IEN,k,j,i) + dt * acc1 * rho * prim(IVX,k,j,i) +  acc2 * rho * prim(IVY,k,j,i) + dt * acc3 * rho * prim(IVZ,k,j,i);
	Real im1_update = cons(IM1,k,j,i) + dt * rho * acc1;
	Real im2_update = cons(IM2,k,j,i) + dt * rho * acc2;
	Real im3_update = cons(IM3,k,j,i) + dt * rho * acc3;
        
        Real ek_update = 0.5*(1.0/cons(IDN,k,j,i))*(SQR(im1_update) + SQR(im2_update) + SQR(im3_update));
	Real update_press = (pmb->peos->GetGamma()-1.0) * (etot_update-ek_update);

	Real vel1_update = im1_update/cons(IDN,k,j,i);
	Real vel2_update = im2_update/cons(IDN,k,j,i);
	Real vel3_update = im3_update/cons(IDN,k,j,i);
	Real vel_update = sqrt(SQR(vel1_update) + SQR(vel2_update) + SQR(vel3_update));
	

        //do dE/dt, or add dt
	//if (update_press >= pfloor){

	//only update when vel_tot < c_light
	Real vel_ceiling = 0.95*pmb->pnrrad->crat;
	if (vel_update <= vel_ceiling){
	
	  cons(IM1,k,j,i) += dt * rho * acc1;
	  cons(IM2,k,j,i) += dt * rho * acc2;
	  cons(IM3,k,j,i) += dt * rho * acc3;
          cons(IEN,k,j,i) = ie_now + ek_update;

	  //cons(IEN,k,j,i) += dt * acc1 * rho * prim(IVX,k,j,i);//0.5 * (pmb->phydro->flux[X1DIR](IDN,k,j,i) + pmb->phydro->flux[X1DIR](IDN,k,j,i+1));
	  //cons(IEN,k,j,i) += dt * acc2 * rho * prim(IVY,k,j,i) + dt * acc3 * rho * prim(IVZ,k,j,i);
	//}else{
	//  cons(IM1,k,j,i) += dt * rho * acc1;
	//  cons(IM2,k,j,i) += dt * rho * acc2;
	//  cons(IM3,k,j,i) += dt * rho * acc3;

	 // cons(IEN,k,j,i) = pfloor/(pmb->peos->GetGamma()-1.0) + ek_update;
	//}

	  Real etot_new = cons(IEN,k,j,i);
	  Real ke_new = 0.5*(1.0/cons(IDN,k,j,i))*(SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
	  Real ie_new = etot_new - ke_new;

	  pmb->ruser_meshblock_data[0](0,k,j,i) = rho * acc1;
	  pmb->ruser_meshblock_data[0](1,k,j,i) = rho * acc2;
	  pmb->ruser_meshblock_data[0](2,k,j,i) = rho * acc3;
	  pmb->ruser_meshblock_data[0](3,k,j,i) = acc1 * rho * prim(IVX,k,j,i) + acc2 * rho * prim(IVY,k,j,i) + acc3 * rho * prim(IVZ,k,j,i);

	  //store energy change due to gravity for diagnostic
	  pmb->ruser_meshblock_data[2](0,k,j,i) = etot_new - etot_now;
	  pmb->ruser_meshblock_data[2](1,k,j,i) = ie_new - ie_now;
	  //pmb->ruser_meshblock_data[2](2,k,j,i) = dt;
	}//vel ceiling
      }
    }
  }
  
}



void PointMassPotential(MeshBlock *pmb, const Real time, const Real dt,const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons){
  //Psuedo Newtonian potential: Phi = -GM/(r/rs-1.)
  //acc at r direction : acc = dPhi/dr = (Phi_l - Phi_r)/dr, forward
  //Momentum_r = -dt * rho * acc
  
  //Energy (dt * rhov * acc/2.0) = 0.5*(dt * rhov * divphi) is calculated by finite volume method:
  //rhov * divPhi = div(rhov*Phi) - Phi*div(rhov) 
  //Phi*div(rhov) = Phi(i) * (Area(i+1)*x1Flux(rhov, i+1) - Area(i)*x1Flux(rhov, i)) 
  //div(rhov*Phi) = Phi(i+1) * x1Flux(rhov, i+1) * Area(i+1) - Phi(i) * x1Flux(rhov, i) * Area(i)
  //Energy = -0.5 * GM * dt * (1/vol(i)) * ( div(rhov*Phi) - Phi*div(rhov) )

  //get the x1flux(rho), rhov
  AthenaArray<Real> &x1flux = pmb->phydro->flux[X1DIR];

  //for each cell, add energy and momentum
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
	// Real rho = prim(IDN,k,j,i);
	// Real Area_m = pmb->pcoord->GetFace1Area(k,j,i);
	// Real Area_p = pmb->pcoord->GetFace1Area(k,j,i+1);
	// Real vol = pmb->pcoord->GetCellVolume(k,j,i);
        // Real r_c = pmb->pcoord->x1v(i);
	// Real r_m = pmb->pcoord->x1f(i);
	// Real r_p = pmb->pcoord->x1f(i+1);
	// Real dr_c = pmb->pcoord->dx1f(i);
	
	// Real Phi_c = -gm/(r_c - 1.0);
	// Real Phi_m = -gm/(r_m - 1.0);
	// Real Phi_p = -gm/(r_p - 1.0);

	// update momentum
	// Real acc = (Phi_p - Phi_m)/dr_c;
	// cons(IM1,k,j,i) += -dt * rho * acc;
	// store external momentum source
	// pmb->ruser_meshblock_data[0](0,k,j,i) = -rho*acc;

	// update energy
	// Real PhiDivRhoV = Phi_c * (Area_p * x1flux(IDN,k,j,i+1) - Area_m * x1flux(IDN,k,j,i)) / vol;
	// Real DivRhoVPhi =(Phi_p * x1flux(IDN,k,j,i+1) * Area_p - Phi_m * x1flux(IDN,k,j,i) * Area_m) / vol;
	
	// initialize user output data
	// Real src = (PhiDivRhoV - DivRhoVPhi);
	// cons(IEN,k,j,i) += dt * src;
	// pmb->ruser_meshblock_data[0](1,k,j,i) = src;

        Real rho = prim(IDN,k,j,i);
        Real phic = -gm/(pmb->pcoord->x1v(i)-1.0);
        Real phil = -gm/(pmb->pcoord->x1f(i)-1.0);
        Real phir = -gm/(pmb->pcoord->x1f(i+1)-1.0);
        Real rr = pmb->pcoord->x1f(i+1);
        Real rl = pmb->pcoord->x1f(i);
        
        Real areal = rl * rl;
        Real arear = rr * rr;
        Real vol = (rr*rr*rr-rl*rl*rl)/3.0;
        Real src = - dt * rho * (phir - phil)/pmb->pcoord->dx1f(i);
        cons(IM1,k,j,i) += src; //make it src/dt/rho
	pmb->ruser_meshblock_data[0](0,k,j,i) = src;

        Real phidivrhov = (arear*x1flux(IDN,k,j,i+1) -
                           areal*x1flux(IDN,k,j,i))*phic/vol;
        Real divrhovphi = (arear*x1flux(IDN,k,j,i+1)*phir -
                           areal*x1flux(IDN,k,j,i)*phil)/vol;
        cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
	pmb->ruser_meshblock_data[0](1,k,j,i) = (phidivrhov - divrhovphi);

	
	//compare old implementation
	Real r_ = pmb->pcoord->x1v(i);
	Real th_ = pmb->pcoord->x2v(j);
	Real ph_ = pmb->pcoord->x3v(k);

	Real rho_ = prim(IDN,k,j,i);
	Real acc1_ = gm/((r_-1.0)*(r_-1.0));

	pmb->ruser_meshblock_data[0](2,k,j,i) = -rho*acc1_;
	pmb->ruser_meshblock_data[0](3,k,j,i) = -acc1_*0.5*(pmb->phydro->flux[X1DIR](IDN,k,j,i) + pmb->phydro->flux[X1DIR](IDN,k,j,i+1));
	
      }//end i
    }//end j
  }//end k

						   

}


void HydroOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh){

  Real local_dens = 0.1;
  Real local_vr = 0.0;
  Real local_vphi = 0.036;
  Real local_press = 1.0e-4;

  for (int k=ks; k<=ke; ++k) {//phi
    Real phi_coord = pco->x3v(k);
    for (int j=js; j<=je; ++j) {//theta
      Real theta_coord = pco->x2v(j);
      for (int i=1; i<=(NGHOST); ++i) {//R
	//printf("phi: %g, theta: %g\n", phi_coord, theta_coord);
	//if (phi_coord<=0.05 || (2*PI-phi_coord<=0.05)){// if within L1 point
	//  if ((theta_coord <= PI/2.+0.05) && (theta_coord >= PI/2.-0.05)){
	    //printf("phi: %g, theta: %g, k: %d, j:%d, i:%d\n", phi_coord, theta_coord,k,j,i);

	//  prim(IDN,k,j,ie+i) = local_dens;
	//  prim(IVX,k,j,ie+i) = local_vr;
	//  prim(IVY,k,j,ie+i) = 0.0;
	//  prim(IVZ,k,j,ie+i) = local_vphi; 
	  
	//  if (NON_BAROTROPIC_EOS) {
	//    prim(IPR,k,j,ie+i) = local_press;
	//  }
	//}else{
	//  prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
	//  prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
	//  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);//pco->x1v(ie+i)*sqrt(GM1/pow(pco->x1v(ie+i),3));
	//  prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
	//  if (NON_BAROTROPIC_EOS){
	//    prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie);
	//  }

	//}
	//}else{//one-direction outflow
	  prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie);
	  prim(IVX,k,j,ie+i) = std::max(0.0, prim(IVX,k,j,ie));
	  prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);//pco->x1v(ie+i)*sqrt(GM1/pow(pco->x1v(ie+i),3));
	  prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie);
	  }
	  //}
	


      }//end R
    }//end theta
  }//end Phi
  
}

void StreamInjectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b, Real time, Real dt,int is, int ie, int js, int je, int ks, int ke, int ngh){

  Real t_current = pmb->pmy_mesh->time;
  Real local_dens_now_ = GetMdot(pmb, t_current);
  //local_dens_now_ = 1.0; //for rho1 run

  Real rinj_thresh=inj_thresh;//2x2x3 (r,th,phi) ghost cells for a wider stream than it is thick
  for (int k=ks; k<=ke; ++k) {//phi
    Real ph_coord = pco->x3v(k);
    for (int j=js; j<=je; ++j) {//theta
      Real th_coord = pco->x2v(j);
      for (int i=1; i<=NGHOST; ++i) {//R

	//try constraining distance to injection point first
	//only inject at theta, phi direction cell for radial boundary
 
	Real r_inj; 
	Real rad; //used to select injection cells only, r_inj and rad are very close (1e-2-1e-1) for i=1 and exactly the same for i=2
	// //XS: my understanding is because boundary cell size can be very different for dx1f(ie+1) and dx1f(ie+2), so here we're setting different r_inj for two ghost cells
	// if (i==1) {
	//   r_inj = rad_inject; //this is very close to rad
	//   //r_inj=pco->x1f(ie+2);
	//   rad=pco->x1f(ie+2);
	// } else {
	//   r_inj = rad_inject2;
	//   //r_inj=pco->x1f(ie+3);
	//   rad=pco->x1f(ie+3);
	// }

	//XS: alternatively, only find closest one cell. If it's i=1, then inject at both (k,j,ie+i) and (k,j,ie+i+1); otherwise, if i=2, inject at (k,j,ie+i) nad (k,j,ie+i-1)
	r_inj = rad_inject; 
	rad=pco->x1f(ie+2);

	
        Real r_coord = pco->x1v(ie+i); //ghost zone radius for density
	//get currtent Cartesian coordinates
	Real x_now = rad * cos(ph_coord) * sin(th_coord);
	Real y_now = rad * sin(ph_coord) * sin(th_coord);
	Real z_now = rad * cos(th_coord);
	//get injection cell coordinates
	
	Real ph_inj = phi_inject;
	Real th_inj = PI/2.0;
	//get injection cell Cartesian coordinates
	Real x_inj = r_inj * cos(ph_inj) * sin(th_inj);
	Real y_inj = r_inj * sin(ph_inj) * sin(th_inj);
	Real z_inj = r_inj * cos(th_inj);

	//Then what to do?
	//we can try saying that if the distance to injection point is within a range rinj_thresh
	Real d_inj = sqrt((x_now - x_inj)*(x_now - x_inj)+
			  (y_now - y_inj)*(y_now - y_inj)+
			  (z_now - z_inj)*(z_now - z_inj));

	//printf("in boundary, d_inj:%g, rinj_thresh:%g\n", d_inj, rinj_thresh);
	if (d_inj < rinj_thresh){ //injection cells
          
          Real dist = sqrt(r_coord*r_coord + r_inj*r_inj - 2.0*r_coord*r_inj*(sin(th_coord)*sin(th_inj)*cos(ph_coord-ph_inj) + 
									    cos(th_coord)*cos(th_inj))); //real distance from injection point
                                                                                                        // (fc) to injection cell (vc)

	  //printf("boundary: r:%g, phi:%g, theta:%g, x:%g, y:%g, z:%g, d_inj:%g, x_inj:%g, y_inj:%g, z_inj:%g\n", r_coord, ph_coord, th_coord, x_now, y_now, z_now, d_inj, x_inj, y_inj, z_inj);
	  //HERE CHOOSE rinj_thresh THAT ONLY ALLOW SINGLE INDEX OF i 
	 //printf("boundary, gid:%d, k:%d, j:%d, i:%d, x1v:%g, x1v_l:%g, x1v_r:%g, x2v:%g, x3v:%g, d_inj:%g\n", pmb->gid, k, j, i, pmb->pcoord->x1v(ie+i), pmb->pcoord->x1v(ie+i-1), pmb->pcoord->x1v(ie+i+1), pmb->pcoord->x2v(j), pmb->pcoord->x3v(k), d_inj);

	  //similar thing here
	 if (i==1){
	   prim(IDN,k,j,ie+i) = local_dens_now_*exp(-(dist*dist)/(r0*r0));
	   prim(IVX,k,j,ie+i) = local_vr; 
	   prim(IVY,k,j,ie+i) = 0.0;
	   prim(IVZ,k,j,ie+i) = local_vphi;
	   prim(IDN,k,j,ie+i+1) = local_dens_now_*exp(-(dist*dist)/(r0*r0));
	   prim(IVX,k,j,ie+i+1) = local_vr; 
	   prim(IVY,k,j,ie+i+1) = 0.0;
	   prim(IVZ,k,j,ie+i+1) = local_vphi;
	   if (NON_BAROTROPIC_EOS){
	     prim(IPR,k,j,ie+i) = prim(IDN,k,j,ie+i)*(temp_stream/temp_unit);
	     prim(IPR,k,j,ie+i+1) = prim(IDN,k,j,ie+i+1)*(temp_stream/temp_unit);
	   }
		   
	 }else if (i==2){
	   prim(IDN,k,j,ie+i) = local_dens_now_*exp(-(dist*dist)/(r0*r0));
	   prim(IVX,k,j,ie+i) = local_vr; 
	   prim(IVY,k,j,ie+i) = 0.0;
	   prim(IVZ,k,j,ie+i) = local_vphi;
	   prim(IDN,k,j,ie+i-1) = local_dens_now_*exp(-(dist*dist)/(r0*r0));
	   prim(IVX,k,j,ie+i-1) = local_vr; 
	   prim(IVY,k,j,ie+i-1) = 0.0;
	   prim(IVZ,k,j,ie+i-1) = local_vphi;
           if (NON_BAROTROPIC_EOS){
             prim(IPR,k,j,ie+i) = prim(IDN,k,j,ie+i)*(temp_stream/temp_unit);
             prim(IPR,k,j,ie+i-1) = prim(IDN,k,j,ie+i-1)*(temp_stream/temp_unit);
           }
	 }else{
	   printf("Do not find X1 outer boundary cell to inject stream\n");
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
    Real ph_inj = phi_inject;
    Real th_inj = PI/2.0;
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=ngh; ++i){
	Real r_inj;
	Real radiusac;
	if (i==1){
	  r_inj = rad_inject;
	  //r_inj=pco->x1f(ie+2);
	  radiusac=pco->x1f(ie+2);
	} else {
	  r_inj = rad_inject2;
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
	  r_inj = rad_inject;
	  //r_inj=pco->x1f(ie+2);
	  radiusac=pco->x1f(ie+2);
	} else {
	  r_inj = rad_inject2;
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
	  //r_inj = rad_inject;
	  r_inj=pco->x1f(ie+2);
	  radiusac=pco->x1f(ie+2);
	} else {
	  //r_inj = rad_inject2;
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

  for (int k=ks; k<=ke; ++k) {//phi
    for (int j=js; j<=je; ++j) {//theta
      for (int i=1; i<=(NGHOST); ++i) {//R
	prim(IDN,k,j,is-i) = prim(IDN,k,j,is); 
	prim(IVX,k,j,is-i) = std::min(0.0, prim(IVX,k,j,is));
	prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);//pco->x1v(is-i)*sqrt(GM1/pow(pco->x1v(is-i),3));
	prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
	if (NON_BAROTROPIC_EOS){
	  prim(IPR,k,j,is-i) = prim(IPR,k,j,is);
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
      Real mu_dir = pnrrad->mu(0,k,j,is-i,ang);
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
      Real mu_dir = pnrrad->mu(0,k,j,ie+i,ang);
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


//input code unit, output code unit
Real kappa_ff(Real temp, Real rho){
  Real rho_cgs = rho*rho_unit;
  Real temp_cgs =  temp*temp_unit;
  Real kappa_cgs = 2.86e-5*(rho_cgs/1.0e-8)*pow(temp_cgs/1.0e6, -3.5);

  return kappa_cgs/kappa_unit;
}

void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck)
{
  
    
    Real logt = log10(tgas * temp_unit);
    Real logrhot = log10(rho* rho_unit) - 3.0* logt + 18.0;
    int nrhot1_planck = 0;
    int nrhot2_planck = 0;
    
    int nrhot1 = 0;
    int nrhot2 = 0;

    while((logrhot > logrhottable_planck(nrhot2_planck)) && (nrhot2_planck < 36)){
      nrhot1_planck = nrhot2_planck;
      nrhot2_planck++;
    }
    if(nrhot2_planck==36 && (logrhot > logrhottable_planck(nrhot2_planck)))
      nrhot1_planck=nrhot2_planck;

    while((logrhot > logrhottable(nrhot2)) && (nrhot2 < 45)){
      nrhot1 = nrhot2;
      nrhot2++;
    }
    if(nrhot2==45 && (logrhot > logrhottable(nrhot2)))
      nrhot1=nrhot2;
  
  /* The data point should between NrhoT1 and NrhoT2 */
    int nt1_planck = 0;
    int nt2_planck = 0;
    int nt1 = 0;
    int nt2 = 0;
    while((logt > logttable_planck(nt2_planck)) && (nt2_planck < 137)){
      nt1_planck = nt2_planck;
      nt2_planck++;
    }
    if(nt2_planck==137 && (logt > logttable_planck(nt2_planck)))
      nt1_planck=nt2_planck;

    while((logt > logttable(nt2)) && (nt2 < 211)){
      nt1 = nt2;
      nt2++;
    }
    if(nt2==211 && (logt > logttable(nt2)))
      nt1=nt2;

  

    Real kappa_t1_rho1=opacitytable(nt1,nrhot1);
    Real kappa_t1_rho2=opacitytable(nt1,nrhot2);
    Real kappa_t2_rho1=opacitytable(nt2,nrhot1);
    Real kappa_t2_rho2=opacitytable(nt2,nrhot2);

    Real planck_t1_rho1=planckopacity(nt1_planck,nrhot1_planck);
    Real planck_t1_rho2=planckopacity(nt1_planck,nrhot2_planck);
    Real planck_t2_rho1=planckopacity(nt2_planck,nrhot1_planck);
    Real planck_t2_rho2=planckopacity(nt2_planck,nrhot2_planck);


    // in the case the temperature is out of range
    // the planck opacity should be smaller by the 
    // ratio T^-3.5
    if(nt2_planck == 137 && (logt > logttable_planck(nt2_planck))){
       Real scaling = pow(10.0, -3.5*(logt - logttable_planck(137)));
       planck_t1_rho1 *= scaling;
       planck_t1_rho2 *= scaling;
       planck_t2_rho1 *= scaling;
       planck_t2_rho2 *= scaling;
    }


    Real rho_1 = logrhottable(nrhot1);
    Real rho_2 = logrhottable(nrhot2);
    Real t_1 = logttable(nt1);
    Real t_2 = logttable(nt2);

    
    if(nrhot1 == nrhot2){
      if(nt1 == nt2){
        kappa = kappa_t1_rho1;
      }else{
        kappa = kappa_t1_rho1 + (kappa_t2_rho1 - kappa_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1 == nt2){
        kappa = kappa_t1_rho1 + (kappa_t1_rho2 - kappa_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);
      }else{
        kappa = kappa_t1_rho1 * (t_2 - logt) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t2_rho1 * (logt - t_1) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t1_rho2 * (t_2 - logt) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t2_rho2 * (logt - t_1) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1));
      }
    }/* end same rhoT */

    rho_1 = logrhottable_planck(nrhot1_planck);
    rho_2 = logrhottable_planck(nrhot2_planck);
    t_1 = logttable_planck(nt1_planck);
    t_2 = logttable_planck(nt2_planck);
 
  /* Now do the same thing for Planck mean opacity */
    if(nrhot1_planck == nrhot2_planck){
      if(nt1_planck == nt2_planck){
        kappa_planck = planck_t1_rho1;
      }else{
        kappa_planck = planck_t1_rho1 + (planck_t2_rho1 - planck_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1_planck == nt2_planck){
        kappa_planck = planck_t1_rho1 + (planck_t1_rho2 - planck_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);

      }else{        
        kappa_planck = planck_t1_rho1 * (t_2 - logt) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t2_rho1 * (logt - t_1) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t1_rho2 * (t_2 - logt) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t2_rho2 * (logt - t_1) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1));
      }
    }/* end same rhoT */

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
}

void opacity(MeshBlock *pmb, AthenaArray<Real> &prim){
  NRRadiation *pnrrad=pmb->pnrrad;
  int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
  //int kl=pmb->kl, ku=pmb->ku, js=pmb->jl, je=pmb->ju, is=pmb->il, ie=pmb->iu;
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js - NGHOST;
  int ju = je + NGHOST;
  int kl = ks - NGHOST;
  int ku = ke + NGHOST;

  for (int k=kl; k<=ku; k++){
    for (int j=jl; j<=ju; j++){
      for (int i=il; i<=iu; i++){
	Real rho = prim(IDN,k,j,i);
	Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i);
	
	//one frequency
	pnrrad->sigma_s(k,j,i,0) = rho*kappa_es_code; //scatter
	pnrrad->sigma_a(k,j,i,0) = 0.0; //flux mean
	pnrrad->sigma_pe(k,j,i,0) = rho*kappa_ff(temp, rho); //energy mean
	pnrrad->sigma_p(k,j,i,0) = rho*kappa_ff(temp, rho); //Planck mean 
	
      }//end i
    }//end j
  }//end k

}


void opalopacity(MeshBlock *pmb, AthenaArray<Real> &prim){
  NRRadiation *pnrrad=pmb->pnrrad;
  int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
  //int kl=pmb->kl, ku=pmb->ku, js=pmb->jl, je=pmb->ju, is=pmb->il, ie=pmb->iu;
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js - NGHOST;
  int ju = je + NGHOST;
  int kl = ks - NGHOST;
  int ku = ke + NGHOST;

  // electron scattering opacity
  //Real kappa_es = 0.2 * (1.0 + 0.7);
  Real kappaa = 0.0;

  for (int k=kl; k<=ku; k++){
    for (int j=jl; j<=ju; j++){
      for (int i=il; i<=iu; i++){
	Real rho = prim(IDN,k,j,i);
	Real temp = std::max(prim(IPR,k,j,i)/prim(IDN,k,j,i), tfloor);
	Real kappa, kappa_planck;
	rossopacity(rho, temp, kappa, kappa_planck);
        Real t_ion = 1.0e4;
	
	if(kappa < kappa_es){
	  if(temp < t_ion/temp_unit){
	    kappaa = kappa;
	    kappa = 0.0;
	  }else{
	    kappaa = 0.0;
	  }
	}else{
	  kappaa = kappa - kappa_es;
	  kappa = kappa_es;
	}
	
	//one frequency
	pnrrad->sigma_s(k,j,i,0) = kappa * rho * rho_unit * l_unit; //scatter
	pnrrad->sigma_a(k,j,i,0) = kappaa * rho * rho_unit * l_unit; //Rossland mean
	pnrrad->sigma_pe(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit; //Rossland mean
        pnrrad->sigma_p(k,j,i,0) = kappa_planck * rho * rho_unit * l_unit;
	// //Planck mean - Rossland mean
	// if(kappaa < kappa_planck){
	//   pnrrad->sigma_planck(k,j,i,0) = (kappa_planck-kappaa)*rho*rho_unit*l_unit;
	// }else{
	//   pnrrad->sigma_planck(k,j,i,0) = 0.0;
        // }
      
      }//end i
    }//end j
  }//end k

}

void GetCombineOpacity(MeshBlock *pmb, AthenaArray<Real> &prim){

  NRRadiation *pnrrad=pmb->pnrrad;
  int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
  //int kl=pmb->kl, ku=pmb->ku, js=pmb->jl, je=pmb->ju, is=pmb->il, ie=pmb->iu;
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl=is, ju=je;
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
  face1.NewAthenaArray((ie-is)+2*NGHOST+2);

  AthenaArray<Real> x1flux = pmb->phydro->flux[X1DIR];
  //bool boundaryflag;
  //boundaryflag = pmb->pbval->apply_bndry_fn_[BoundaryFace::outer_x1];
  //printf("outer boundary flag:%d\n", boundaryflag);
  if (pmb->pbval->get_apply_bndry_fn_flag(BoundaryFace::outer_x1)){
    //printf("here\n");
    for (int k=ks; k<=ke; k++){
      for (int j=js; j<=je; j++){
	pmb->pcoord->Face1Area(k , j, is, ie, face1);
	for (int i=ie; i<=ie; i++){
	  massflux += face1(ie)*x1flux(0,k,j,ie);//x1flux(0) is the density flux, multiply by volume to get mass
          //printf("face:%g, x1flux:%g\n", face1(ie), x1flux(0, k, j, ie));
	}

      }
    }
  }

  //printf("outter mass flux:%g\n", massflux);

  face1.DeleteAthenaArray();
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

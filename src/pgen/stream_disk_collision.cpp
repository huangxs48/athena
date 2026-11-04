// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details

// A stream colliding with a section of accretion disk
// --prob stream_disk_collision

//========================================================================================
// C header

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
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../units/units.hpp"

static Real temp_unit, l_unit, rho_unit;
static Real vel_unit, time_unit; //derived units
static Real tfloor;
static Real dfloor, pfloor, d_amb, p_amb;

static int bh_grav_flag;
static Real GM_bh;
static Real x_bh, y_bh, z_bh; //BH location

static int uniform_stream;
static Real r_inj_thresh;
static Real vx_inj, vy_inj, vz_inj;
static Real rho_inj, temp_inj, h_inj, x_inj, y_inj, z_inj;
static Real rho_norm;

//time step related 
static int force_rhd_timestep; //force radiation timestep (not implemented yet)
Real user_dt;

//numericla setup related
//disk set up
static int include_disk;
static Real h_disk;
static Real temp_disk_mid;
static Real rho_disk_mid;
static Real disk_scalar_frac;
static int disk_grav_flag;
static Real grav_disk; //disk gravity
static Real z_disk_mid;

void Phi_BH_DiskGrav(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

// User-defined boundary conditions for disk simulations
void StreamInjectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
			 Real time, Real dt,
			 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void RadOuterX2(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
	        const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
        	Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);


Real RadTimeStep(MeshBlock *pmb){
  // //leave the option to force radiation timestep
  // if (force_rhd_timestep==1) {
  //   Real crat = Constants::speed_of_light_cgs/vel_unit; 
  //   Real cfl = pmb->pmy_mesh->cfl_number;
  // //loop over cells to find minimal timestep
  // }
  
  if (user_dt>0.0){
    return user_dt;
  }else{
    return pmb->pmy_mesh->dt;
  }
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
  Real x1length = mesh_size.nx1+1+2*NGHOST;

  temp_unit = pin->GetReal("problem", "temp_unit");
  l_unit = pin->GetReal("problem", "l_unit");
  rho_unit = pin->GetReal("problem", "rho_unit");
  vel_unit = pin->GetReal("problem", "vel_unit");
  time_unit = l_unit/vel_unit;

  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.001); 
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 0.001);
  d_amb = pin->GetOrAddReal("problem", "d_amb", 0.0001);
  p_amb = pin->GetOrAddReal("problem", "p_amb", 0.001);

  bh_grav_flag = pin->GetOrAddInteger("problem", "bh_grav_flag", 0);
  GM_bh = pin->GetOrAddReal("problem", "GM_bh", 1.0);
  x_bh = pin->GetOrAddReal("problem", "x_bh", 100.0);
  y_bh = pin->GetOrAddReal("problem", "y_bh", 100.0);
  z_bh = pin->GetOrAddReal("problem", "z_bh", 100.0);

  include_disk = pin->GetOrAddInteger("problem", "include_disk", 0);
  h_disk = pin->GetOrAddReal("problem", "h_disk", 1.0);
  temp_disk_mid = pin->GetOrAddReal("problem", "temp_disk_mid", 0.0);
  rho_disk_mid = pin->GetOrAddReal("problem", "rho_disk_mid", 0.0);
  z_disk_mid = pin->GetOrAddReal("problem", "z_disk_mid", 0.0);
  disk_grav_flag = pin->GetOrAddInteger("problem", "disk_grav_flag", 0);
  grav_disk = pin->GetOrAddReal("problem", "grav_disk", 0.0);
  disk_scalar_frac = pin->GetOrAddReal("problem", "grav_disk", 4.0); //number of scale heigh to paint disk scalar

  uniform_stream = pin->GetOrAddReal("problem", "uniform_stream", 0.0);
  r_inj_thresh = pin->GetOrAddReal("problem", "r_inj_thresh", 0.0);
  h_inj = pin->GetOrAddReal("problem", "h_inj", 1.0);
  x_inj = pin->GetOrAddReal("problem", "x_inj", 0.0);
  y_inj = pin->GetOrAddReal("problem", "y_inj", 0.0);
  z_inj = pin->GetOrAddReal("problem", "z_inj", 0.0);
  rho_inj = pin->GetOrAddReal("problem", "rho_inj", 0.0);
  rho_norm = pin->GetOrAddReal("problem", "rho_norm", 1.0);
  temp_inj = pin->GetOrAddReal("problem", "temp_inj", 0.0);
  vx_inj = pin->GetOrAddReal("problem", "vx_inj", 0.0);
  vy_inj = pin->GetOrAddReal("problem", "vy_inj", 0.0);
  vz_inj = pin->GetOrAddReal("problem", "vz_inj", 0.0);

  user_dt = pin->GetOrAddReal("problem", "user_dt", -1.0);
  
  EnrollUserExplicitSourceFunction(Phi_BH_DiskGrav);

  EnrollUserTimeStepFunction(RadTimeStep);
  //if (adaptive==true)
  //  EnrollUserRefinementCondition(RefinementCondition);


  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, StreamInjectOuterX2);
  }

  if (NR_RADIATION_ENABLED){
    //Enroll rad boundaries
    if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
      EnrollUserRadBoundaryFunction(BoundaryFace::outer_x2, RadOuterX2);
    }
  }
  
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

  
  // AllocateRealUserMeshBlockDataField(1); 
  //ruser_meshblock_data[0].NewAthenaArray(6,blocksizex3, blocksizex2, blocksizex1); //store external source

  //store fofc flag   //fofc?
  //AllocateIntUserMeshBlockDataField(1); 
  //iuser_meshblock_data[0].NewAthenaArray(1,blocksizex3, blocksizex2, blocksizex1); //store external source

  Real x1_l= pcoord->x1f(is); 
  Real x1_u = pcoord->x1f(ie);
  Real x2_l = pcoord->x2f(js);
  Real x2_u = pcoord->x2f(je);
  Real x3_l = pcoord->x3f(ks);
  Real x3_u = pcoord->x3f(ke);

  // AllocateUserOutputVariables(7);
  // SetUserOutputVariableName(2, "GravSrc_m2_IM1");
  // SetUserOutputVariableName(3, "GravSrc_m2_IM2");
  // SetUserOutputVariableName(4, "GravSrc_m2_IM3");
  // SetUserOutputVariableName(5, "GravSrc_m2_IEN");

  return;
}

int RefinementCondition(MeshBlock *pmb){ //not used yet, not sure if it works with self-gravity

  AthenaArray<Real> &w = pmb->phydro->w;
  Real max_d2 = 0.0;//maximum length to m2 in a meshblock
  for(int k=pmb->ks; k<=pmb->ke; k++) {
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {

	//current coordinate
	Real x_now = pmb->pcoord->x1v(i);
	Real y_now = pmb->pcoord->x2v(j);
	Real z_now = pmb->pcoord->x3v(k);

	Real r_now = sqrt(pow(x_now, 2) + pow(y_now, 2) + pow(z_now, 2));

	max_d2 = std::max(max_d2, r_now);
      }
    }
  }
  //printf("max_d2:%g\n", max_d2);
  if(max_d2 < 0.8) return 1;
  if(max_d2 > 0.8) return -1;
  return 0;
}


void Mesh::UserWorkInLoop(){
  //find pressure increase region
  for(int nb=0; nb<nblocal; ++nb){ //loop over meshblocks on the same core
    MeshBlock *pmb = my_blocks(nb);
    ParameterInput *pin;

    Hydro *phydro = pmb->phydro;
    Coordinates *pcoord = pmb->pcoord;
    int ks=pmb->ks, ke=pmb->ke, js=pmb->js, je=pmb->je, is=pmb->is, ie=pmb->ie;
         
    for(int k=pmb->ks; k<=pmb->ke; k++){
      for(int i=pmb->is; i<=pmb->ie; i++){
    	for(int j=pmb->js; j<=pmb->je; j++){
	  Real ycoord = pcoord->x2v(j);
	  Real rho = phydro->u(IDN,k,j,i);
	  Real press = phydro->w(IPR,k,j,i);
	  Real temp = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);

	  // // ////debug
	  for (int n=0; n<(NHYDRO); n++){
	    if (std::isnan(phydro->u(n,k,j,i))){
	      printf("x1v:%g, x2v:%g, n:%d eint:%g, rho:%g, gid:%d, k:%d, j:%d, i:%d\n", pcoord->x1v(i), pcoord->x2v(j), n, phydro->u(IEN,k,j,i), phydro->u(IDN,k,j,i), pmb->gid, k, j, i);
		abort();
	      }
	  }

	  //inject energy
	  

	  Real vel1 = pmb->phydro->w(IVX,k,j,i);
	  Real vel2 = pmb->phydro->w(IVY,k,j,i);
	  Real vel3 = pmb->phydro->w(IVZ,k,j,i);
	  Real c_light = Constants::speed_of_light_cgs/vel_unit; 
	  Real vel_mag = sqrt(vel1*vel1 + vel2*vel2 + vel3*vel3);

	  if (vel_mag >= c_light){
	    //printf("bad cell, rho:%g, press:%g, ie:%g, temp:%g, vel1:%g, vel2:%g, x1v:%g, x2v:%g\n", rho, press, phydro->u(IEN,k,j,i), temp, phydro->w(IVX,k,j,i), phydro->w(IVY,k,j,i), pmb->pcoord->x1v(i), pmb->pcoord->x2v(j));
	    pmb->phydro->w(IVX,k,j,i) = (vel1 >=0) ? std::min(vel1, 0.99*c_light) : std::max(vel1, -0.99*c_light);
	    pmb->phydro->w(IVY,k,j,i) = (vel2 >=0) ? std::min(vel2, 0.99*c_light) : std::max(vel2, -0.99*c_light);
	    pmb->phydro->w(IVZ,k,j,i) = vel3 = (vel3 >=0) ? std::min(vel3, 0.99*c_light) : std::max(vel3, -0.99*c_light);
	  }



	}//j
      }//i
    }//k


  }//loop over meshblocks

}

// void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){
//   for(int k=ks; k<=ke; k++){
//     for(int j=js; j<=je; j++){
//       for(int i=is; i<=ie; i++){
//   	user_out_var(0,k,j,i) = ruser_meshblock_data[0](0,k,j,i);
//    	user_out_var(1,k,j,i) = ruser_meshblock_data[0](1,k,j,i);
//    	user_out_var(2,k,j,i) = ruser_meshblock_data[0](2,k,j,i);
//    	user_out_var(3,k,j,i) = ruser_meshblock_data[0](3,k,j,i);
//         user_out_var(4,k,j,i) = ruser_meshblock_data[0](4,k,j,i);
//         user_out_var(5,k,j,i) = ruser_meshblock_data[0](5,k,j,i);
//         user_out_var(6,k,j,i) = iuser_meshblock_data[0](0,k,j,i);
//       }
//      }
//    }


// }

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real gamma_gas = peos->GetGamma();

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

	//current coordinate
	Real x_now = pcoord->x1v(i);
	Real y_now = pcoord->x2v(j);
	Real z_now = pcoord->x3v(k);
	Real dis_now = (y_now - z_disk_mid);

	Real dens_now = pin->GetReal("hydro", "dfloor");
	Real pgas_now = pin->GetReal("hydro", "pfloor");

	//add the disk if a disk is included
	if (include_disk==1 && fabs(dis_now)<(disk_scalar_frac * h_disk)){
	  dens_now = rho_disk_mid * exp(-SQR(dis_now)/SQR(h_disk)/2.0);
	  pgas_now = dens_now * temp_disk_mid;
	}
	
	Real d0 = 0.0; //paint the disk gas
        d0 = fabs(dis_now) > disk_scalar_frac * h_disk ? d0 : 1.0;

	dens_now = std::max(dens_now, d_amb);
	pgas_now = std::max(pgas_now, p_amb);

        phydro->u(IDN,k,j,i) = dens_now;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
	
	if (NON_BAROTROPIC_EOS) {
          //phydro->u(IEN,k,j,i) = (rho_now*temp_iso)/(gamma_gas - 1.0);
	  phydro->u(IEN,k,j,i) = pgas_now/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
	}//end non barotropic

        if (NSCALARS>0){ //scalar to mark the disk and stream, not used if not compiled with --nscalar=1
	  pscalars->s(0,k,j,i) = phydro->u(IDN,k,j,i) * d0;
	}

      }//end i
      

    }//end j
  }//end k

  return;
}

/* utility functions */

/* end utility functions */

void Phi_BH_DiskGrav(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar){

  //Newtonian graivity
  if (bh_grav_flag==0)
    return;
  
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real rho = prim(IDN,k,j,i);

	Real x = pmb->pcoord->x1v(i);
	Real y = pmb->pcoord->x2v(j);
	Real z = pmb->pcoord->x3v(k);

	Real r = sqrt((x - x_bh) * (x - x_bh)
		      + (y - y_bh) * (y - y_bh)
		      + (z - z_bh) * (z - z_bh));

	Real dx = prim(IVX, k,j,i);
	Real dy = prim(IVY, k,j,i);
	Real dz = prim(IVZ, k,j,i);
	  
	//acceleration components
	Real ax = - GM_bh / (r * r) * (x - x_bh)/r;
	Real ay = - GM_bh / (r * r) * (y - y_bh)/r;
	Real az = - GM_bh / (r * r) * (z - z_bh)/r;
	
	//force conserving internal energy
	//current internal energy , not adapt to generalized EoS
	Real etot_now = cons(IEN,k,j,i);
	Real ke_now = 0.5*(1.0/cons(IDN,k,j,i))*(SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
	Real ie_now = etot_now - ke_now;
	
	Real im1_update = cons(IM1,k,j,i) + dt * rho * ax;
	Real im2_update = cons(IM2,k,j,i) + dt * rho * ay;
	Real im3_update = cons(IM3,k,j,i) + dt * rho * az;
	
	Real ek_update = 0.5*(1.0/cons(IDN,k,j,i))*(SQR(im1_update) + SQR(im2_update) + SQR(im3_update));
	
	cons(IM1,k,j,i) += dt * rho * ax;
	cons(IM2,k,j,i) += dt * rho * ay;
	cons(IM3,k,j,i) += dt * rho * az;
	cons(IEN,k,j,i) = ie_now + ek_update;
	
	//}

      }
    }
  }

}



/* boundary functions */

void StreamInjectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
			 Real time, Real dt,
			 int il, int iu, int jl, int ju, int kl, int ku, int ngh){

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=NGHOST; ++j) {
      for (int i=il; i<=iu; ++i) {
	
	//get currtent Cartesian coordinates
	Real x_now = pco->x1v(i);
	Real y_now = pco->x2v(ju+j);
	Real z_now = pco->x3v(k);

	//distance to injection point
	Real dis_to_inj = sqrt((x_now - x_inj)*(x_now - x_inj)
			       + (y_now - y_inj)*(y_now - y_inj)
			       + (z_now - z_inj)*(z_now - z_inj)
			       );
	if (dis_to_inj < r_inj_thresh){
	  //printf("xnow:%g, ynow:%g, znow:%g, dis_to_inj:%g, r_inj_thresh:%g\n", x_now, y_now, z_now,  dis_to_inj, r_inj_thresh);
	  Real local_dens = rho_inj * rho_norm;

	  //stream density, if not uniform
	  if (uniform_stream==0){
	    local_dens = local_dens * std::exp( -std::pow(dis_to_inj/h_inj, 2)/2.0 );
	  }

	  prim(IDN,k,ju+j,i) = local_dens;
	  prim(IVX,k,ju+j,i) = vx_inj;
	  prim(IVY,k,ju+j,i) = vy_inj;
	  prim(IVZ,k,ju+j,i) = vz_inj;

	  //printf("local_dens:%g, vx_inj:%g, vy_inj:%g, vz_inj:%g\n", local_dens, vx_inj, vy_inj, vz_inj);
	  
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,ju+j,i) = local_dens * temp_inj;
	  }
	}else{
	  prim(IDN,k,ju+j,i) = prim(IDN,k,ju,i);
	  prim(IVX,k,ju+j,i) = prim(IVX,k,ju,i);
	  prim(IVZ,k,ju+j,i) = std::max(0.0, prim(IVZ,k,ju,i));
	  prim(IVY,k,ju+j,i) = prim(IVY,k,ju,i);
	  if (NON_BAROTROPIC_EOS){
	    prim(IPR,k,ju+j,i) = prim(IPR,k,ju,i); 
	  }
	  
	}//r_inj_thresh
	
      }//i
    }//j
  }//k

}

/*
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
      //Real mu_dir = prad->mu(1,k,je+j,i,n);
      
      ir(k,je+j,i,ang) = ir(k,je,i,ang);

    }// end n
  }// end ifr
  }}}

  return;

}
*/

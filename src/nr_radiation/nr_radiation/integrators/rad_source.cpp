//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================


// C headers

// C++ headers
#include <algorithm>
#include <string>
#include <vector>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../implicit/radiation_implicit.hpp"
#include "../radiation.hpp"

// class header
#include "./rad_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void RadIntegrator::CalSourceTerms(MeshBlock *pmb, const Real dt,
        const int k, const int j, const int i, AthenaArray<Real> &u,
        AthenaArray<Real> &ir_ini, AthenaArray<Real> &ir) {
  NRRadiation *prad=pmb->pnrrad;
  Real invcrat = 1.0/prad->crat;

  Real *sigma_at, *sigma_s, *sigma_p, *sigma_pe;
  Real *lab_ir;

  const int &nang =prad->nang;
  const int &nfreq=prad->nfreq;

  //xiaoshan: for debugging
  Real invredc = 1.0/prad->reduced_c;
  Real invredfactor = invredc/invcrat;

  // Get the temporary arrays
  AthenaArray<Real> &wmu_cm = wmu_cm_;
  AthenaArray<Real> &tran_coef = tran_coef_;
  AthenaArray<Real> &ir_cm = ir_cm_;
  AthenaArray<Real> &cm_to_lab = cm_to_lab_;

  AthenaArray<Real> split_ratio;
  AthenaArray<int> map_start, map_end;

  // for implicit update, using the quantities from the partially
  // updated u, not from w
  Real rho = u(IDN,k,j,i);
  Real vx = vel_source_(k,j,i,0);
  Real vy = vel_source_(k,j,i,1);
  Real vz = vel_source_(k,j,i,2);
  Real vsq = vx*vx + vy*vy + vz*vz;

  Real lorzsq = 1.0/(1.0 - vsq  * invcrat * invcrat);
  Real lorz = std::sqrt(lorzsq);

  sigma_at = &(prad->sigma_a(k,j,i,0));
  sigma_s = &(prad->sigma_s(k,j,i,0));
  sigma_p = &(prad->sigma_p(k,j,i,0));
  sigma_pe =&(prad->sigma_pe(k,j,i,0));

  //xiaoshan: add a cell-wise doppler flag for debugging
  //doppler_flag_cell_(k,j,i) = 1;
  pmb->ruser_meshblock_data[1](0,k,j,i) = 1.0;
  //only do the Comoving<->Lab step for certain density and temperature range, don't change intensity in the funnel region
  //first constrain by density and radial direction velocity, then density and temperature, using aggresvie condition
  Real tgas_now = tgas_(k,j,i) * prad->tunit;
  Real rho_now = rho * prad->rhounit;
  Real vx_now = vx / prad->crat;//use the velocity goes into Doppler shift	
  if ((rho_now<=1.0e-14 && vx_now>=0.15) || (rho_now<=1.0e-14 && tgas_now>1.0e7)){
    //doppler_flag_cell_(k,j,i) = 0;
    pmb->ruser_meshblock_data[1](0,k,j,i) = 0.0;
  }
  pmb->user_out_var(29,k,j,i) += pmb->ruser_meshblock_data[1](0,k,j,i); //doppler_flag_cell_(k,j,i);

  // Prepare the transformation coefficients
  Real numsum = 0.0;

  for (int n=0; n<nang; ++n) {
    Real vdotn = vx * prad->mu(0,k,j,i,n) + vy * prad->mu(1,k,j,i,n)
                + vz * prad->mu(2,k,j,i,n);
    Real vnc = 1.0 - vdotn * invcrat;
    tran_coef(n) = lorz * vnc;
    wmu_cm(n) = prad->wmu(n)/(tran_coef(n) * tran_coef(n));
    numsum += wmu_cm(n);
    cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
  }
           // Normalize weight in co-moving frame to make sure the sum is one
  numsum = 1.0/numsum;

  for (int n=0; n<nang; ++n) {
    wmu_cm(n) *= numsum;
  }

  for (int ifr=0; ifr<nfreq; ++ifr) {
    lab_ir=&(ir_ini(k,j,i,ifr*nang));
    for (int n=0; n<nang; ++n)
      ir_cm(n+ifr*nang) = lab_ir[n];
    if (IM_RADIATION_ENABLED) {
      // add the explicit flux divergence term
      Real *divflux = &(divflx_(k,j,i,ifr*nang));
      for (int n=0; n<nang; ++n) {
        ir_cm(n+ifr*nang) += divflux[n];
      }
      // add angular flux
      if (prad->angle_flag == 1) {
        Real *p_angflx  = &(ang_flx_(k,j,i,ifr*nang));
        for (int n=0; n<nang; ++n) {
          ir_cm(n+ifr*nang) += p_angflx[n];
        }
      }
    }
    for (int n=0; n<nang; ++n) {
      ir_cm(n+ifr*nang) *= cm_to_lab(n);
    }
    // apply floor to ir_cm
    for (int n=0; n<nang; ++n) {
      ir_cm(n+ifr*nang) = std::max(ir_cm(n+ifr*nang),static_cast<Real>(TINY_NUMBER));
    }

    for (int n=0; n<nang; n++) {
      implicit_coef_(n+ifr*nang) = 1.0;
      if (IM_RADIATION_ENABLED) {
        implicit_coef_(n+ifr*nang) += const_coef_(k,j,i,n+ifr*nang);
        if (prad->angle_flag == 1) {
          implicit_coef_(n+ifr*nang) += imp_ang_coef_(k,j,i,n);
        }
      }
    }
  }
  if (nfreq == 1) {
  // Add absorption and scattering opacity source
    tgas_new_(k,j,i) = AbsorptionScattering(wmu_cm,tran_coef, sigma_at, sigma_p,
                                            sigma_pe, sigma_s, dt, lorz, rho,
                                            tgas_(k,j,i),
                                            implicit_coef_, ir_cm);
    // Add compton scattering
    if (compton_flag_ > 0) {
      Compton(wmu_cm,tran_coef, sigma_s, dt, lorz, rho, tgas_new_(k,j,i), ir_cm);
    }
  } else {
    // map frequency grid
    if (doppler_flag_ > 0) {

      //record ifr=5 and ifr=18 group total intensity and r direction flux increment
      Real ir_5_before = 0.0;
      Real ir_18_before = 0.0;
      Real fr_5_before = 0.0;
      Real fr_18_before = 0.0;

      Real ir_5_after = 0.0;
      Real ir_18_after = 0.0;
      Real fr_5_after = 0.0;
      Real fr_18_after = 0.0;
      
      for (int n=0; n<nang; ++n) {
	//get the angular weight for debugging flux
	//weight of domega*cos(mu)_x
	Real wmu = prad->wmu(n);
	Real wmu_cosx = prad->mu(0,k,j,i,n) * prad->wmu(n);
        for (int ifr=0; ifr<nfreq; ++ifr) {
          ir_ori_(ifr) = ir_cm(ifr*nang+n);
	  if (ifr==5) ir_5_before += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==5) fr_5_before += ir_cm(ifr*nang+n)*wmu_cosx;
	  if (ifr==18) ir_18_before += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==18) fr_18_before += ir_cm(ifr*nang+n)*wmu_cosx;
        }
        split_ratio.InitWithShallowSlice(split_ratio_,3,n,1);
        map_start.InitWithShallowSlice(map_bin_start_,2,n,1);
        map_end.InitWithShallowSlice(map_bin_end_,2,n,1);

        MapLabToCmFrequency(tran_coef(n), split_ratio, map_start, map_end,
                                                     ir_ori_, ir_done_);

        for (int ifr=0; ifr<nfreq; ++ifr) {
	  //only do the step for certain density and temperature range, don't change intensity in the funnel region
	  //first constrain by density and radial direction velocity, then density and temperature, using aggresvie condition
	  Real tgas_now = tgas_(k,j,i) * prad->tunit;
	  Real rho_now = u(IDN,k,j,i) * prad->rhounit;
          Real vx_now = vel_source_(k,j,i,0) / prad->crat;//use the velocity goes into Doppler shift

	  if (pmb->ruser_meshblock_data[1](0,k,j,i)==0.0){
	    //std::cout<<"doppler flag zero in gid:"<<pmb->gid<<" k: "<<k<<" j: "<<j<<" i: "<<i<<" rho: "<<pmb->phydro->w(IDN,k,j,i)<<" temp: "<<pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i)<<std::endl;
	    ir_done_(ifr) = ir_cm(ifr*nang+n); //ir_ori_(ifr);
	  }//zero doppler flag 
          ir_cm(ifr*nang+n) = ir_done_(ifr);
	  if (ifr==5) ir_5_after += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==5) fr_5_after += ir_cm(ifr*nang+n)*wmu_cosx;
	  if (ifr==18) ir_18_after += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==18) fr_18_after += ir_cm(ifr*nang+n)*wmu_cosx;
        }
      }
      //write diff into user mesh block data, make them the same unit as delta_source terms
      pmb->user_out_var(25,k,j,i) += (ir_5_after - ir_5_before) * invredfactor; //energy moment
      pmb->user_out_var(26,k,j,i) += (ir_18_after - ir_18_before) * invredfactor; 
      pmb->user_out_var(27,k,j,i) += (fr_5_after - fr_5_before) * invredc; //flux moment
      pmb->user_out_var(28,k,j,i) += (fr_18_after - fr_18_before) * invredc;
      
    }

    Real ir_5_before = 0.0;
    Real ir_18_before = 0.0;
    Real fr_5_before = 0.0;
    Real fr_18_before = 0.0;

    Real ir_5_after = 0.0;
    Real ir_18_after = 0.0;
    Real fr_5_after = 0.0;
    Real fr_18_after = 0.0;

    for (int n=0; n<nang; ++n) {
      //get the angular weight for debugging flux
      //weight of domega*cos(mu)_x
      Real wmu = prad->wmu(n);
      Real wmu_cosx = prad->mu(0,k,j,i,n) * prad->wmu(n);
      for (int ifr=0; ifr<nfreq; ++ifr) {
	if (ifr==5) ir_5_before += ir_cm(ifr*nang+n)*wmu;
	if (ifr==5) fr_5_before += ir_cm(ifr*nang+n)*wmu_cosx;
	if (ifr==18) ir_18_before += ir_cm(ifr*nang+n)*wmu;
	if (ifr==18) fr_18_before += ir_cm(ifr*nang+n)*wmu_cosx;
      }
    }
    
    // calculate the source term
    tgas_new_(k,j,i) = MultiGroupAbsScat(wmu_cm,tran_coef, sigma_at, sigma_p,
                               sigma_pe, sigma_s, dt, lorz, rho, tgas_(k,j,i),
                                                     implicit_coef_,ir_cm);

    for (int n=0; n<nang; ++n){
      //get the angular weight for debugging flux
      //weight of domega*cos(mu)_x
      Real wmu = prad->wmu(n);
      Real wmu_cosx = prad->mu(0,k,j,i,n) * prad->wmu(n);
      for (int ifr=0; ifr<nfreq; ++ifr) {
	if (ifr==5) ir_5_after += ir_cm(ifr*nang+n)*wmu;
	if (ifr==5) fr_5_after += ir_cm(ifr*nang+n)*wmu_cosx;
	if (ifr==18) ir_18_after += ir_cm(ifr*nang+n)*wmu;
	if (ifr==18) fr_18_after += ir_cm(ifr*nang+n)*wmu_cosx;
      }
    }

    //write diff into user mesh block data, make them the same unit as delta_source terms
    pmb->user_out_var(21,k,j,i) += (ir_5_after - ir_5_before) * invredfactor; //energy moment
    pmb->user_out_var(22,k,j,i) += (ir_18_after - ir_18_before) * invredfactor; //energy moment
    pmb->user_out_var(23,k,j,i) += (fr_5_after - fr_5_before) * invredc; //flux moment
    pmb->user_out_var(24,k,j,i) += (fr_18_after - fr_18_before) * invredc; //flux moment
    
    // Add compton scattering
    // Compton scattering for implicit scheme is added separately
    if ((compton_flag_ > 0) && (split_compton_ == 0)) {
      Real t_ini = tgas_new_(k,j,i);
      Real t_old;
      int count=0;
      Real relative_error = 1;
      while((count < iteration_compton_) && (relative_error > compton_error_)) {
        ir_buff_ = ir_cm;
        t_old = tgas_new_(k,j,i);
        MultiGroupCompton(wmu_cm,tran_coef,dt,lorz,rho,t_ini,tgas_new_(k,j,i),
                                                                    ir_buff_);
        count++;
        relative_error = std::abs(t_old-tgas_new_(k,j,i))/tgas_new_(k,j,i);
      }
      ir_cm = ir_buff_;
      // inverseshift
    }
    // map frequency grid
    if (doppler_flag_ > 0) {
      //record ifr=5 and ifr=18 group total intensity and r direction flux increment
      Real ir_5_before = 0.0;
      Real ir_18_before = 0.0;
      Real fr_5_before = 0.0;
      Real fr_18_before = 0.0;

      Real ir_5_after = 0.0;
      Real ir_18_after = 0.0;
      Real fr_5_after = 0.0;
      Real fr_18_after = 0.0;
      
      for (int n=0; n<nang; ++n) {
	//get the angular weight for debugging flux
	//weight of domega*cos(mu)_x
	Real wmu = prad->wmu(n);
	Real wmu_cosx = prad->mu(0,k,j,i,n) * prad->wmu(n);
        for (int ifr=0; ifr<nfreq; ++ifr) {
          ir_ori_(ifr) = ir_cm(ifr*nang+n);
	  if (ifr==5) ir_5_before += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==5) fr_5_before += ir_cm(ifr*nang+n)*wmu_cosx;
	  if (ifr==18) ir_18_before += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==18) fr_18_before += ir_cm(ifr*nang+n)*wmu_cosx;
        }
	
        split_ratio.InitWithShallowSlice(split_ratio_,3,n,1);
        map_start.InitWithShallowSlice(map_bin_start_,2,n,1);
        map_end.InitWithShallowSlice(map_bin_end_,2,n,1);

        bool invertible = FreMapMatrix(split_ratio, tran_coef(n), map_start,
                                        map_end, map_count_, fre_map_matrix_);
        if (invertible) {
          bool success = InverseMapFrequency(tran_coef(n), map_count_,
                                          fre_map_matrix_, ir_ori_, ir_done_);

          if (!success)
            MapCmToLabFrequency(tran_coef(n),split_ratio, map_start, map_end,
                                                        ir_ori_,ir_done_);

        } else {
          MapCmToLabFrequency(tran_coef(n),split_ratio, map_start, map_end,
                                                        ir_ori_,ir_done_);
        }

        for (int ifr=0; ifr<nfreq; ++ifr) {
          //only do the step for certain density and temperature range, don't change intensity in the funnel region
          //first constrain by density and radial direction velocity, then density and temperature, using aggresvie condition
          Real tgas_now = tgas_(k,j,i) * prad->tunit;
          Real rho_now = u(IDN,k,j,i) * prad->rhounit;
          Real vx_now = vel_source_(k,j,i,0) / prad->crat;//use the velocity goes into Doppler shift

          if (pmb->ruser_meshblock_data[1](0,k,j,i)==0.0){
	    ir_done_(ifr) = ir_cm(ifr*nang+n); //ir_ori_(ifr);
	  }//no doppler flag
          ir_cm(ifr*nang+n) = ir_done_(ifr);
	  if (ifr==5) ir_5_after += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==5) fr_5_after += ir_cm(ifr*nang+n)*wmu_cosx;
	  if (ifr==18) ir_18_after += ir_cm(ifr*nang+n)*wmu;
	  if (ifr==18) fr_18_after += ir_cm(ifr*nang+n)*wmu_cosx;
        }
      }
    
      //write diff into user mesh block data, make them the same unit as delta_source terms
      pmb->user_out_var(13,k,j,i) += (ir_5_after - ir_5_before) * invredfactor; //energy moment
      pmb->user_out_var(14,k,j,i) += (ir_18_after - ir_18_before) * invredfactor; 
      pmb->user_out_var(15,k,j,i) += (fr_5_after - fr_5_before) * invredc; //flux moment
      pmb->user_out_var(16,k,j,i) += (fr_18_after - fr_18_before) * invredc;
    }
  }

  // update specific intensity in the lab frame
  // do not modify ir_ini
  Real omega = 1.0;
  if (IM_RADIATION_ENABLED) {
    omega = pmb->pmy_mesh->pimrad->omega;
  }
  Real omega_1 = 1.0 - omega;

  if (std::abs(omega_1) < TINY_NUMBER) {
    Real ir_5_before = 0.0;
    Real ir_18_before = 0.0;
    Real fr_5_before = 0.0;
    Real fr_18_before = 0.0;

    Real ir_5_after = 0.0;
    Real ir_18_after = 0.0;
    Real fr_5_after = 0.0;
    Real fr_18_after = 0.0;
    for (int ifr=0; ifr<nfreq; ++ifr) {
      lab_ir = &(ir(k,j,i,nang*ifr));
      for (int n=0; n<nang; ++n) {
	Real wmu_cosx = prad->mu(0,k,j,i,n) * prad->wmu(n);
	if (ifr==5) ir_5_before += ir_cm(n+ifr*nang)/cm_to_lab(n);
	if (ifr==5) fr_5_before += wmu_cosx*ir_cm(n+ifr*nang)/cm_to_lab(n);
	if (ifr==18) ir_18_before += ir_cm(n+ifr*nang)/cm_to_lab(n);
	if (ifr==18) fr_18_before += wmu_cosx*ir_cm(n+ifr*nang)/cm_to_lab(n);	       
        lab_ir[n] = std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),
                             static_cast<Real>(TINY_NUMBER));
	if (ifr==5) ir_5_after += std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),static_cast<Real>(TINY_NUMBER));
	if (ifr==5) fr_5_after += wmu_cosx*std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),static_cast<Real>(TINY_NUMBER));
	if (ifr==18) ir_18_after += std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),static_cast<Real>(TINY_NUMBER));
	if (ifr==18) fr_18_after += wmu_cosx*std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),static_cast<Real>(TINY_NUMBER));	 
      }
    }
    // //write diff into user mesh block data
    // pmb->user_out_var(17,k,j,i) += ir_5_after - ir_5_before;
    // pmb->user_out_var(18,k,j,i) += ir_18_after - ir_18_before;
    // pmb->user_out_var(19,k,j,i) += fr_5_after - fr_5_before;
    // pmb->user_out_var(20,k,j,i) += fr_18_after - fr_18_before;
  } else {
    for (int ifr=0; ifr<nfreq; ++ifr) {
      lab_ir = &(ir(k,j,i,nang*ifr));
      for (int n=0; n<nang; ++n) {
        Real ir_old = lab_ir[n];
        Real ir_new = ir_cm(n+ifr*nang)/cm_to_lab(n);
        ir_new = omega_1 * ir_old + omega * ir_new;
        lab_ir[n] = std::max(ir_new, static_cast<Real>(TINY_NUMBER));
      }
    }
    tgas_new_(k,j,i) = omega_1*tgas_(k,j,i) + omega*tgas_new_(k,j,i);
  }
}



void RadIntegrator::AddMultiGroupCompt(MeshBlock *pmb, const Real dt,
        AthenaArray<Real> &u, AthenaArray<Real> &ir) {
  // need to transform lab frame ir to co-moving frame
  NRRadiation *prad=pmb->pnrrad;

  const Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;
  Real invredfactor = prad->crat/prad->reduced_c;

  AthenaArray<Real> split_ratio;
  AthenaArray<int> map_start, map_end;

  Real *lab_ir;

  const int& nang =prad->nang;
  const int& nfreq=prad->nfreq;

  // only apply for multi-grou case
  if ((nfreq > 1) && (compton_flag_ > 0) && (split_compton_ > 0)) {
    // Get the temporary arrays
    AthenaArray<Real> &wmu_cm = wmu_cm_;
    AthenaArray<Real> &tran_coef = tran_coef_;
    AthenaArray<Real> &ir_cm = ir_cm_;
    AthenaArray<Real> &cm_to_lab = cm_to_lab_;

    int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
    int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real rho = u(IDN,k,j,i);
          Real vx = vel_source_(k,j,i,0);
          Real vy = vel_source_(k,j,i,1);
          Real vz = vel_source_(k,j,i,2);
          Real vsq = vx*vx + vy*vy + vz*vz;

          Real lorzsq = 1.0/(1.0 - vsq  * invcrat * invcrat);
          Real lorz = std::sqrt(lorzsq);

          // Prepare the transformation coefficients
          Real numsum = 0.0;

          for (int n=0; n<nang; ++n) {
             Real vdotn = vx * prad->mu(0,k,j,i,n) + vy * prad->mu(1,k,j,i,n)
                        + vz * prad->mu(2,k,j,i,n);
             Real vnc = 1.0 - vdotn * invcrat;
             tran_coef(n) = lorz * vnc;
             wmu_cm(n) = prad->wmu(n)/(tran_coef(n) * tran_coef(n));
             numsum += wmu_cm(n);
             cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
          }
           // Normalize weight in co-moving frame to make sure the sum is one
          numsum = 1.0/numsum;
#pragma omp simd
          for (int n=0; n<nang; ++n) {
            wmu_cm(n) *= numsum;
          }
          for (int ifr=0; ifr<nfreq; ++ifr) {
            lab_ir=&(ir(k,j,i,ifr*nang));
            for (int n=0; n<nang; ++n)
              ir_cm(n+ifr*nang) = std::max(lab_ir[n] * cm_to_lab(n),
                                           static_cast<Real>(TINY_NUMBER));
            // store the moments
            Real er_fr = 0.0;
            for (int n=0; n<nang; ++n) {
              Real ir_weight = lab_ir[n]*prad->wmu(n);
              er_fr += ir_weight;
            }
            delta_source(0,ifr) = er_fr;
          }
          // map frequency grid
          for (int n=0; n<nang; ++n) {
            for (int ifr=0; ifr<nfreq; ++ifr)
              ir_ori_(ifr) = ir_cm(ifr*nang+n);

            split_ratio.InitWithShallowSlice(split_ratio_,3,n,1);
            map_start.InitWithShallowSlice(map_bin_start_,2,n,1);
            map_end.InitWithShallowSlice(map_bin_end_,2,n,1);

            MapLabToCmFrequency(tran_coef(n), split_ratio, map_start, map_end,
                                                            ir_ori_, ir_done_);


            for (int ifr=0; ifr<nfreq; ++ifr) {
              ir_cm(ifr*nang+n) = ir_done_(ifr);
            }
          }

          // Add compton scattering
          Real t_ini = tgas_new_(k,j,i);
          Real t_old;
          int count=0;
          Real relative_error = 1;
          while ((count < iteration_compton_) && (relative_error > compton_error_)) {
            ir_buff_ = ir_cm;
            t_old = tgas_new_(k,j,i);
            MultiGroupCompton(wmu_cm,tran_coef,dt,lorz,rho,t_ini,tgas_new_(k,j,i),
                                                                        ir_buff_);
            count++;
            relative_error = std::abs(t_old-tgas_new_(k,j,i))/tgas_new_(k,j,i);
          }
          ir_cm = ir_buff_;

          // map frequency grid
          for (int n=0; n<nang; ++n) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              ir_ori_(ifr) = ir_cm(ifr*nang+n);
            }

            split_ratio.InitWithShallowSlice(split_ratio_,3,n,1);
            map_start.InitWithShallowSlice(map_bin_start_,2,n,1);
            map_end.InitWithShallowSlice(map_bin_end_,2,n,1);

            bool invertible = FreMapMatrix(split_ratio, tran_coef(n), map_start,
                                            map_end, map_count_, fre_map_matrix_);
            if (invertible) {
              bool success = InverseMapFrequency(tran_coef(n), map_count_,
                                          fre_map_matrix_, ir_ori_, ir_done_);
              if (!success) {
                  MapCmToLabFrequency(tran_coef(n),split_ratio, map_start, map_end,
                                                        ir_ori_,ir_done_);
              }
            } else {
              MapCmToLabFrequency(tran_coef(n),split_ratio, map_start, map_end,
                                                             ir_ori_,ir_done_);
            }

            for (int ifr=0; ifr<nfreq; ++ifr) {
              ir_cm(ifr*nang+n) = ir_done_(ifr);
            }
          }
          for (int ifr=0; ifr<nfreq; ++ifr) {
            lab_ir = &(ir(k,j,i,nang*ifr));
            for (int n=0; n<nang; ++n) {
              lab_ir[n] = std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),
                                   static_cast<Real>(TINY_NUMBER));
            }
          }

          if (prad->set_source_flag > 0) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              Real *p_ir =  &(ir(k,j,i,ifr*nang));
              Real er_fr = 0.0;
              for (int n=0; n<nang; ++n) {
                Real ir_weight = p_ir[n] * prad->wmu(n);
                er_fr  += ir_weight;
              }
              rad_source(0,k,j,i) += (-prat*(er_fr-delta_source(0,ifr))*invredfactor);
            }
          }
        }
      }
    }
  }
}

// ir_ini and ir only differ by the source term for explicit scheme
// for implicit scheme, ir also includes the flux divergence term

void RadIntegrator::GetHydroSourceTerms(MeshBlock *pmb,
                       AthenaArray<Real> &ir_ini, AthenaArray<Real> &ir) {
  NRRadiation *prad=pmb->pnrrad;
  const Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;
  Real invredc = 1.0/prad->reduced_c;
  Real invredfactor = invredc/invcrat;

  const int& nang =prad->nang;
  const int& nfreq=prad->nfreq;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // first, calculate Er and Fr in lab frame before the step
        for (int ifr=0; ifr<nfreq; ++ifr) {
          Real *p_ir0 =  &(ir_ini(k,j,i,ifr*nang));
          Real er_fr = 0.0;
          Real frx_fr = 0.0;
          Real fry_fr = 0.0;
          Real frz_fr = 0.0;
          if (IM_RADIATION_ENABLED) {
            for (int n=0; n<nang; ++n) {
              Real ir_weight = p_ir0[n];
              ir_weight += divflx_(k,j,i,ifr*nang+n);
              if (prad->angle_flag == 1) {
                ir_weight += ang_flx_(k,j,i,ifr*nang+n);
                ir_weight -= ((imp_ang_coef_(k,j,i,n))
                           * ir(k,j,i,ifr*nang+n));
              }
              ir_weight -= (const_coef_(k,j,i,ifr*nang+n)
                           * ir(k,j,i,ifr*nang+n));
              ir_weight *= prad->wmu(n);
              er_fr  += ir_weight;
              frx_fr += ir_weight * prad->mu(0,k,j,i,n);
              fry_fr += ir_weight * prad->mu(1,k,j,i,n);
              frz_fr += ir_weight * prad->mu(2,k,j,i,n);
            }
          } else { // end IM_Radiation
#pragma omp simd reduction (+:er_fr,frx_fr,fry_fr,frz_fr)
            for (int n=0; n<nang; ++n) {
              Real ir_weight = p_ir0[n];
              ir_weight *= prad->wmu(n);
              er_fr  += ir_weight;
              frx_fr += ir_weight * prad->mu(0,k,j,i,n);
              fry_fr += ir_weight * prad->mu(1,k,j,i,n);
              frz_fr += ir_weight * prad->mu(2,k,j,i,n);
            }
          }
          delta_source(0,ifr) = er_fr;
          delta_source(1,ifr) = frx_fr;
          delta_source(2,ifr) = fry_fr;
          delta_source(3,ifr) = frz_fr;
        }
        for (int ifr=0; ifr<nfreq; ++ifr) {
          Real *p_ir =  &(ir(k,j,i,ifr*nang));
          Real er_fr = 0.0;
          Real frx_fr = 0.0;
          Real fry_fr = 0.0;
          Real frz_fr = 0.0;
#pragma omp simd reduction (+:er_fr,frx_fr,fry_fr,frz_fr)
          for (int n=0; n<nang; ++n) {
            Real ir_weight = p_ir[n] * prad->wmu(n);
            er_fr  += ir_weight;
            frx_fr += ir_weight * prad->mu(0,k,j,i,n);
            fry_fr += ir_weight * prad->mu(1,k,j,i,n);
            frz_fr += ir_weight * prad->mu(2,k,j,i,n);
          }
          delta_source(0,ifr) = er_fr  - delta_source(0,ifr);
          delta_source(1,ifr) = frx_fr - delta_source(1,ifr);
          delta_source(2,ifr) = fry_fr - delta_source(2,ifr);
          delta_source(3,ifr) = frz_fr - delta_source(3,ifr);
        }

        Real delta_er = 0.0, delta_frx=0.0, delta_fry = 0.0, delta_frz = 0.0;

        for (int ifr=0; ifr<nfreq; ++ifr) {
          delta_er  += delta_source(0,ifr);
          delta_frx += delta_source(1,ifr);
          delta_fry += delta_source(2,ifr);
          delta_frz += delta_source(3,ifr);
	  if (ifr==5){
            pmb->user_out_var(17,k,j,i) += (prat*delta_source(0,ifr) * invredfactor);
            pmb->user_out_var(18,k,j,i) += (prat*delta_source(1,ifr) * invredc);
          }
          if (ifr==18){
            pmb->user_out_var(19,k,j,i) += (prat*delta_source(0,ifr) * invredfactor);
            pmb->user_out_var(20,k,j,i) += (prat*delta_source(1,ifr) * invredc);
          }
        }

        // Now apply the radiation source terms to gas with energy and
        // momentum conservation
        rad_source(0,k,j,i) = (-prat*delta_er  * invredfactor);
        rad_source(1,k,j,i) = (-prat*delta_frx * invredc);
        rad_source(2,k,j,i) = (-prat*delta_fry * invredc);
        rad_source(3,k,j,i) = (-prat*delta_frz * invredc);
	
      }
    }
  }
}


void RadIntegrator::AddSourceTerms(MeshBlock *pmb, AthenaArray<Real> &u) {
  NRRadiation *prad=pmb->pnrrad;
  Real invcrat = 1.0/prad->crat;
  Field *pfield = pmb->pfield;

  Real gm1 = pmb->peos->GetGamma() - 1.0;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        u(IM1,k,j,i) += rad_source(1,k,j,i);
        u(IM2,k,j,i) += rad_source(2,k,j,i);
        u(IM3,k,j,i) += rad_source(3,k,j,i);

        // limit the velocity by speed of light
        Real vx = u(IM1,k,j,i)/u(IDN,k,j,i);
        Real vy = u(IM2,k,j,i)/u(IDN,k,j,i);
        Real vz = u(IM3,k,j,i)/u(IDN,k,j,i);
        Real vsq = vx*vx+vy*vy+vz*vz;
        Real vel = std::sqrt(vsq);
        Real ratio = vel * invcrat;
        if (ratio > prad->vmax) {
          Real factor = prad->vmax/ratio;
          u(IM1,k,j,i) *= factor;
          u(IM2,k,j,i) *= factor;
          u(IM3,k,j,i) *= factor;
        }

        Real ekin = 0.5 *(SQR(u(IM1,k,j,i))+SQR(u(IM2,k,j,i))
                          + SQR(u(IM3,k,j,i)))/u(IDN,k,j,i);
        Real pb = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                    + SQR(pfield->bcc(IB3,k,j,i)));
        }

        if (prad->set_source_flag == 2) {
          Real eint = tgas_new_(k,j,i) * u(IDN,k,j,i)/gm1;
          u(IEN,k,j,i) = eint + pb + ekin;
        } else {
          Real e_source = rad_source(0,k,j,i);
          // first check that gas internal energy will not become negative
          Real eint = u(IEN,k,j,i) + e_source - ekin - pb;
          Real tgas = eint * gm1/u(IDN,k,j,i);
          if (eint < 0.0) {
            eint = tgas_new_(k,j,i) * u(IDN,k,j,i)/gm1;
            u(IEN,k,j,i) = eint + pb + ekin;
          } else if (tgas > prad->t_ceiling_(k,j,i)) {
            eint = prad->t_ceiling_(k,j,i) * u(IDN,k,j,i)/gm1;
            u(IEN,k,j,i) = ekin + pb + eint;
          } else {
            u(IEN,k,j,i) += e_source;
          }
        }
      }
    }
  }
  
}

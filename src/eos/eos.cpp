//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos.cpp
//! \brief FOFC-support helpers: conservative->primitive test (writes fofc_ flags)
//!        and reusable SingleConservativeToPrimitive{Hydro,MHD} building blocks.

// C headers

// C++ headers
#include <cmath>   // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ConservedToPrimitiveTest
//! \brief Flag cells where the ConsToPrim would have hit the density or pressure floor.
//!        cons(IEN) must be e_int + e_k (no e_mag contribution) even for MHD.

void EquationOfState::ConservedToPrimitiveTest(
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bcc,
    int il, int iu, int jl, int ju, int kl, int ku) {
  if (!fofc_.IsAllocated()) {
    fofc_.NewAthenaArray(pmy_block_->ncells3, pmy_block_->ncells2, pmy_block_->ncells1);
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real u_d  = cons(IDN,k,j,i);
        Real u_m1 = cons(IM1,k,j,i);
        Real u_m2 = cons(IM2,k,j,i);
        Real u_m3 = cons(IM3,k,j,i);
        Real u_e  = cons(IEN,k,j,i);
        Real e_mag = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          Real bcc1 = bcc(IB1,k,j,i);
          Real bcc2 = bcc(IB2,k,j,i);
          Real bcc3 = bcc(IB3,k,j,i);
          e_mag = 0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
        }

        Real w_d, w_vx, w_vy, w_vz, w_p, dp;
        bool dfloor_used = false, pfloor_used = false;
        SingleConservativeToPrimitiveMHD(u_d, u_m1, u_m2, u_m3, u_e,
                                         w_d, w_vx, w_vy, w_vz, w_p,
                                         dp, dfloor_used, pfloor_used, e_mag);
        fofc_(k,j,i) = dfloor_used || pfloor_used;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::SingleConservativeToPrimitiveHydro
//! \brief Pure-hydro conservative->primitive for a single cell; sets floor-use flags.
void EquationOfState::SingleConservativeToPrimitiveHydro(
  Real &u_d, Real &u_m1, Real &u_m2, Real &u_m3, Real &u_e,
  Real &w_d, Real &w_vx, Real &w_vy, Real &w_vz, Real &w_p,
  Real &dp, bool &dfloor_used, bool &pfloor_used) {
  if (u_d < density_floor_) {
    u_d = density_floor_;
    dfloor_used = true;
  }
  w_d = u_d;

  Real di = 1.0/u_d;
  w_vx = u_m1*di;
  w_vy = u_m2*di;
  w_vz = u_m3*di;

  if (NON_BAROTROPIC_EOS) {
    Real gm1 = gamma_ - 1.0;
    Real e_k = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
    w_p = gm1*(u_e - e_k);

    if (w_p < pressure_floor_) {
      dp = pressure_floor_ - w_p;
      w_p = pressure_floor_;
      u_e = w_p/gm1 + e_k;
      pfloor_used = true;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::SingleConservativeToPrimitiveMHD
//! \brief MHD conservative->primitive for a single cell; subtracts the supplied e_mag
//!        from the total energy before applying the pressure floor.
void EquationOfState::SingleConservativeToPrimitiveMHD(
  Real &u_d, Real &u_m1, Real &u_m2, Real &u_m3, Real &u_e,
  Real &w_d, Real &w_vx, Real &w_vy, Real &w_vz, Real &w_p,
  Real &dp, bool &dfloor_used, bool &pfloor_used, const Real e_mag) {
  if (u_d < density_floor_) {
    u_d = density_floor_;
    dfloor_used = true;
  }
  w_d = u_d;

  Real di = 1.0/u_d;
  w_vx = u_m1*di;
  w_vy = u_m2*di;
  w_vz = u_m3*di;

  if (NON_BAROTROPIC_EOS) {
    Real gm1 = gamma_ - 1.0;
    Real e_k = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
    w_p = gm1*(u_e - e_k - e_mag);

    if (w_p < pressure_floor_) {
      dp = pressure_floor_ - w_p;
      w_p = pressure_floor_;
      u_e = w_p/gm1 + e_k + e_mag;
      pfloor_used = true;
    }
  }
}

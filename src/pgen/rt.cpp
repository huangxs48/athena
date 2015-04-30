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

// Primary header
#include "../mesh.hpp"

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // GetGamma
#include "../field/field.hpp"      // magnetic field
#include "../bvals/bvals.hpp"      // EnrollFluidBoundaryFunction
#include "../fluid/srcterms/srcterms.hpp"  // GetG2, GetG3

double ran2(long int *idum); // random number generator from NR
void reflect_ix2(MeshBlock *pmb, AthenaArray<Real> &a);
void reflect_ox2(MeshBlock *pmb, AthenaArray<Real> &a);
void reflect_ix3(MeshBlock *pmb, AthenaArray<Real> &a);
void reflect_ox3(MeshBlock *pmb, AthenaArray<Real> &a);

// made global to share with BC functions
static Real gm1;
static Real grav_acc;

//======================================================================================
//! \file rt.c
//  \brief Problem generator for RT instabilty.
//
// Note the gravitational acceleration is hardwired to be 0.1. Density difference is
// hardwired to be 2.0 in 2D, and is set by the input parameter <problem>/rhoh in 3D
// (default value is 3.0). This reproduces 2D results of Liska & Wendroff, 3D results of
// Dimonte et al.
// 
// FOR 2D HYDRO:
// Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to match Liska
// & Wendroff. Interface is at y=0; perturbation added to Vy. Gravity acts in y-dirn.
// Special reflecting boundary conditions added in x2 to improve hydrostatic eqm
// (prevents launching of weak waves) Atwood number A=(d2-d1)/(d2+d1)=1/3. Options:
//     iprob = 1  -- Perturb V2 using single mode
//     iprob != 1 -- Perturb V2 using multiple mode
//
// FOR 3D:
// Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1, gamma=5/3 to
// match Dimonte et al.  Interface is at z=0; perturbation added to Vz. Gravity acts in
// z-dirn. Special reflecting boundary conditions added in x3.  A=1/2.  Options:
//     iprob = 1 -- Perturb V3 using single mode
//     iprob = 2 -- Perturb V3 using multiple mode
//     iprob = 3 -- B rotated by "angle" at interface, multimode perturbation
//
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  long int iseed = -1;
  Real gamma = pfl->pf_eos->GetGamma();
  gm1 = gamma - 1.0;
  
  Real kx = 2.0*(PI)/(pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min);
  Real ky = 2.0*(PI)/(pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min);
  Real kz = 2.0*(PI)/(pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min);

  // Read perturbation amplitude, problem switch, density ratio
  Real amp = pin->GetReal("problem","amp");
  int iprob = pin->GetInteger("problem","iprob");
  Real drat = pin->GetOrAddReal("problem","drat",3.0);


// 2D PROBLEM ---------------------------------------------------------------

  if (pmb->block_size.nx3 == 1) {
    grav_acc = pfl->pf_srcterms->GetG2();
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real den=1.0;
        if (pmb->x2v(j) > 0.0) den *= drat;

        if (iprob == 1) {
          pfl->u(IM2,k,j,i) = (1.0+cos(kx*pmb->x1v(i)))*(1.0+cos(ky*pmb->x2v(j)))/4.0;
        } else {
          pfl->u(IM2,k,j,i) = (ran2(&iseed) - 0.5)*(1.0+cos(ky*pmb->x2v(j)));
        }

        pfl->u(IDN,k,j,i) = den;
        pfl->u(IM1,k,j,i) = 0.0;
        pfl->u(IM2,k,j,i) *= (den*amp);
        pfl->u(IM3,k,j,i) = 0.0;
        pfl->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pmb->x2v(j)))/gm1;
        pfl->u(IEN,k,j,i) += 0.5*SQR(pfl->u(IM2,k,j,i))/den;
      }
    }}

    // initialize interface B, same for all iprob
    if (MAGNETIC_FIELDS_ENABLED) {
      // Read magnetic field strength, angle [in degrees, 0 is along +ve X-axis]
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfd->b.x1f(k,j,i) = b0;
      }}}
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        pfd->b.x2f(k,j,i) = 0.0;
      }}}
      for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfd->b.x3f(k,j,i) = 0.0;
      }}}
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfl->u(IEN,k,j,i) += 0.5*b0*b0;
        }}}
      }
    }

    // Enroll special BCs
    pmb->pbval->EnrollFluidBoundaryFunction(inner_x2, reflect_ix2);
    pmb->pbval->EnrollFluidBoundaryFunction(outer_x2, reflect_ox2);

// 3D PROBLEM ----------------------------------------------------------------

  } else {
    grav_acc = pfl->pf_srcterms->GetG3();
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real den=1.0;
        if (pmb->x3v(k) > 0.0) den *= drat;

        if (iprob == 1) {
          pfl->u(IM3,k,j,i) = (1.0+cos(kx*(pmb->x1v(i))))/8.0
                     *(1.0+cos(ky*pmb->x2v(j)))*(1.0+cos(kz*pmb->x3v(k)));
        } else {
          pfl->u(IM3,k,j,i) = amp*(ran2(&iseed) - 0.5)*(1.0+cos(kz*pmb->x3v(k)));
        }

        pfl->u(IDN,k,j,i) = den;
        pfl->u(IM1,k,j,i) = 0.0;
        pfl->u(IM2,k,j,i) = 0.0;
        pfl->u(IM3,k,j,i) *= (den*amp);
        pfl->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pmb->x3v(k)))/gm1;
        pfl->u(IEN,k,j,i) += 0.5*SQR(pfl->u(IM3,k,j,i))/den;
      }
    }}

    // initialize interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      // Read magnetic field strength, angle [in degrees, 0 is along +ve X-axis]
      Real b0 = pin->GetReal("problem","b0");
      Real angle = pin->GetReal("problem","angle");
      angle = (angle/180.)*PI;
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        if (pmb->x3v(k) > 0.0) {
          pfd->b.x1f(k,j,i) = b0;
        } else {
          pfd->b.x1f(k,j,i) = b0*cos(angle);
        }
      }}}
      for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        if (pmb->x3v(k) > 0.0) {
          pfd->b.x2f(k,j,i) = 0.0;
        } else {
          pfd->b.x2f(k,j,i) = b0*sin(angle);
        }
      }}}
      for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfd->b.x3f(k,j,i) = 0.0;
      }}}
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfl->u(IEN,k,j,i) += 0.5*b0*b0;
        }}}
      }
    }

    // Enroll special BCs
    pmb->pbval->EnrollFluidBoundaryFunction(inner_x3, reflect_ix3);
    pmb->pbval->EnrollFluidBoundaryFunction(outer_x3, reflect_ox3);

  } /* end of 3D initialization */

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void reflect_ix2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void reflect_ix2(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM2)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM2,k,js-j,i) = -a(IM2,k,js+j-1,i);  // reflect 2-mom
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IEN,k,js-j,i) = a(IEN,k,js+j-1,i) 
             - a(IDN,k,js+j-1,i)*grav_acc*(2*j-1)*pmb->dx2f(j)/gm1;
        }
      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,k,js-j,i) = a(n,k,js+j-1,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void reflect_ox2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void reflect_ox2(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM2)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM2,k,je+j,i) = -a(IM2,k,je-j+1,i);  // reflect 2-mom
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IEN,k,je+j,i) = a(IEN,k,je-j+1,i) 
             + a(IDN,k,je-j+1,i)*grav_acc*(2*j-1)*pmb->dx2f(j)/gm1;
        }
      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,k,je+j,i) = a(n,k,je-j+1,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void reflect_ix3()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void reflect_ix3(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM3)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM3,ks-k,j,i) = -a(IM3,ks+k-1,j,i);  // reflect 3-mom
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IEN,ks-k,j,i) = a(IEN,ks+k-1,j,i) 
             - a(IDN,ks+k-1,j,i)*grav_acc*(2*k-1)*pmb->dx3f(k)/gm1;
        }
      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,ks-k,j,i) = a(n,ks+k-1,j,i);
        }
      }

    }
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void reflect_ox3()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void reflect_ox3(MeshBlock *pmb, AthenaArray<Real> &a)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ke = pmb->ke;

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int n=0; n<(NFLUID); ++n) {

      if (n==(IM3)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IM3,ke+k,j,i) = -a(IM3,ke-k+1,j,i);  // reflect 3-mom
        }
      } else if (n==(IEN)) {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(IEN,ke+k,j,i) = a(IEN,ke-k+1,j,i)
             + a(IDN,ke-k+1,j,i)*grav_acc*(2*k-1)*pmb->dx3f(k)/gm1;
        }
      } else {
#pragma simd
        for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
          a(n,ke+k,j,i) = a(n,ke-k+1,j,i);
        }
      }

    }
  }}

  return;
}
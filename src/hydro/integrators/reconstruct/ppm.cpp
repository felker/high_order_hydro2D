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
//! \file ppm.cpp
//  \brief piecewise parabolic reconstruction, adapted from remap in CMHOG code
//======================================================================================

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../hydro.hpp"
#include "../../../mesh.hpp"
#include "../../../coordinates/coordinates.hpp"

// this class header
#include "../hydro_integrator.hpp"

//limiter options for all directions (define only one)
#define SEKORA_6
#undef MIGNONE
#undef VL_LIMITER

//--------------------------------------------------------------------------------------
//! \fn HydroIntegrator::ReconstructionFuncX1()
//  \brief 

void HydroIntegrator::PiecewiseParabolicX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  AthenaArray<Real> dd_, dph_, c1_, c2_, c3_, c4_, c5_, c6_; 
  Real qa,qb,qc,qd, qe; // do I want to keep these variables only inside the loop scope?

  double lim_const = 1.25; //Colella/Sekora constant used in second derivative limiter, >1, independent of dx

  int ncells1 = iu - il +2*NGHOST; //=264 //(iu+2) - (il-1) +1; 
  //allocate scratch arrays for parabolic interpolation
  //continue using underscores for private variable??
  dd_.NewAthenaArray(ncells1);
  dph_.NewAthenaArray(ncells1);
  c1_.NewAthenaArray(ncells1);
  c2_.NewAthenaArray(ncells1);
  c3_.NewAthenaArray(ncells1);
  c4_.NewAthenaArray(ncells1);
  c5_.NewAthenaArray(ncells1);
  c6_.NewAthenaArray(ncells1);
  //Temporary array of L/R parabolic interpolation coefficients that follows CMHOG and Colella notation
  AthenaArray<Real> qplus_, qminus_; 
  qplus_.NewAthenaArray(ncells1);
  qminus_.NewAthenaArray(ncells1);
  
  //Mesh spacing information
  Coordinates *pco = pmy_hydro->pmy_block->pcoord;
 
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il-2; i<=iu+1; ++i){
      //Compute coefficients used in interpolation formulae
      Real& dx_im2 = pco->dx1v(i-2); 
      Real& dx_im1 = pco->dx1v(i-1);
      Real& dx_i   = pco->dx1v(i);
      Real& dx_ip1   = pco->dx1v(i+1);
      qa = dx_i/(dx_im1 + dx_i + dx_ip1); //first factor in eq 1.7
      c1_(i) = qa*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i); 
      c2_(i) = qa*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i); 
    }
    //The q* are used to construct the interpolation coefficients c*_
    //All of these variables only depend on the mesh spacing (and can be reused for all prim)
    for (int i=il-1; i<=iu+1; ++i){
      Real& dx_im2 = pco->dx1v(i-2); 
      Real& dx_im1 = pco->dx1v(i-1);
      Real& dx_i   = pco->dx1v(i);
      Real& dx_ip1   = pco->dx1v(i+1);
      qa = dx_im2 + dx_im1 + dx_i + dx_ip1;
      qb = dx_im1/(dx_im1 + dx_i); 
      qc = (dx_im2 + dx_im1)/(2.0*dx_im1 + dx_i); 
      qd = (dx_ip1 + dx_i)/(2.0*dx_i + dx_im1); 
      qb = qb + 2.0*dx_i*qb/qa*(qc-qd); 
      c3_(i) = 1.0 - qb; 
      c4_(i) = qb; 
      c5_(i) = dx_i/qa*qd; 
      c6_(i) = -dx_im1/qa*qc; 
    }

    for (int i=il-2; i<=iu+1; ++i){
      // Compute average linear slopes (eqn 1.7)
      // Only used in this loop to compute and monotonize average zone slope, dd_(i)
      Real dplus = q(n,k,j,i+1) - q(n,k,j,i); 
      Real dmnus = q(n,k,j,i) - q(n,k,j,i-1);
      dd_(i) = c1_(i)*dplus + c2_(i)*dmnus;
      
      // Monotonize, van Leer limiting (eqn 1.8)
#ifdef VL_LIMITER
      qa = std::min( fabs(dd_(i)), 2.0*std::min(fabs(dmnus),fabs(dplus)) ); 
      dd_(i) = 0.0;
      if (dplus*dmnus > 0.0) 
      dd_(i) = qa*SIGN(dd_(i));  
#endif
    }
    
    // Initialize interface values (eqn 1.6)
    // a_j-1/2 value that is assigned to most L/R Riemann states
    for (int i=il-1; i<=(iu+1); ++i){ //need the first two ghost cells at far end of domain to do reconstruction on A_l / R Riemann state at last boundary
      Real& dx_i   = pco->dx1v(i);
      //Interpolation function for nonuniform mesh spacing
#if defined(VL_LIMITER) || defined(MIGNONE)
      dph_(i)= c3_(i)*q(n,k,j,i-1) + c4_(i)*q(n,k,j,i) + c5_(i)*dd_(i-1) + c6_(i)*dd_(i); 
      //Simplified expression for uniform mesh, no van Leer monotonization of average zone slopes
      //      dph_(i)= 7.0/12.0*(q(n,k,j,i-1) + q(n,k,j,i)) - 1.0/12.0*(q(n,k,j,i+1) + q(n,k,j,i-2));
#endif
      
      //6th order approximation to interface state, uniform mesh spacing
      //Colella Sekora 2008 equation 17 (i -=1)
#ifdef SEKORA_6
      dph_(i) = 37.0/60.0*(q(n,k,j,i-1) + q(n,k,j,i)) - 2.0/15.0*(q(n,k,j,i-2)+ q(n,k,j,i+1)) +
      1.0/60.0*(q(n,k,j,i-3)+ q(n,k,j,i+2));
      // Colella Sekora eqn 13: Check monotonization constraint 
      if (std::min(q(n,k,j,i-1), q(n,k,j,i)) > dph_(i) || std::max(q(n,k,j,i-1), q(n,k,j,i)) < dph_(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_i*dx_i)*(q(n,k,j,i-1) - 2.0*dph_(i) + q(n,k,j,i)); 
	//left
	qb = 1.0/(dx_i*dx_i)*(q(n,k,j,i-2) - 2.0*q(n,k,j,i-1) + q(n,k,j,i)); 
	//right
	qc = 1.0/(dx_i*dx_i)*(q(n,k,j,i-1) - 2.0*q(n,k,j,i) + q(n,k,j,i+1)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_(i) = 0.5*(q(n,k,j,i-1) + q(n,k,j,i)) - dx_i*dx_i/3.0 * qd; 
      } //end sekora limiter for extrema 
#endif      
      //To print out the initial a_{j-1/2}
      /*      if (j==25 && n == 0){
	std::cout << i << " , " << dph_(i) << "\n"; 
	} */
    } //end loop for interpolation of interface states

    //Initialize parabolic interpolation coefficients a_{L,j} = a_{j-1/2}, a_{R,j} = a_{j+3/2}
    for (int i=il-1; i<=iu; ++i){ 
      //For use in Sekora and standard PPM limiters
      qminus_(i) = dph_(i  );
      qplus_(i) = dph_(i+1 );  
      //Alternative monotonization: Mignone 2014, eq 45. Do not need to check above monotonization conditions
      //since minmax assignments automatically adjust violations while preserving normal cases
#ifdef MIGNONE
      qminus_(i) = std::min(qminus_(i), std::max(q(n,k,j,i),q(n,k,j,i-1)));
      qplus_(i) = std::min(qplus_(i), std::max(q(n,k,j,i),q(n,k,j,i+1))); 

      qminus_(i) = std::max(qminus_(i), std::min(q(n,k,j,i),q(n,k,j,i-1)));
      qplus_(i) = std::max(qplus_(i), std::min(q(n,k,j,i),q(n,k,j,i+1)));  
#endif
    } 

    // Apply original PPM limiters to jth cell reconstruction 
#ifndef SEKORA_6
    for (int i=il-1; i<=iu; ++i){ 
      qa = (qplus_(i) - q(n,k,j,i))*(q(n,k,j,i) - qminus_(i));
      qd = qplus_(i)-qminus_(i);
      qe = 6.0*(q(n,k,j,i) - 0.5*(qplus_(i) + qminus_(i))); //this is a6 parabolic coefficient

      //Parabolic overshoot monotonization (eqn 1.10). They are mutually exclusive
      //Alternative formulation of condition commented out. See Colella Sekora
      if ( (qd*(qd - qe)) < 0.0) { 	//	  if (fabs(qminus_(i) - q(n,k,j,i)) >= 2*fabs(qplus_(i) - q(n,k,j,i))){	
	qminus_(i) = 3.0*q(n,k,j,i) - 2.0*qplus_(i); 
      }      
      if ( (qd*(qd + qe)) < 0.0) {  //    if (fabs(qplus_(i) - q(n,k,j,i)) >= 2*fabs(qminus_(i) - q(n,k,j,i))){
	qplus_(i) = 3.0*q(n,k,j,i) - 2.0*qminus_(i);
      }
      // Local extrema monotonization (eqn 1.10)
      if (qa <= 0.0){
	qminus_(i) = q(n,k,j,i);
	qplus_(i) = q(n,k,j,i);
      }      
    } //end loop over i, original PPM limiters
#endif    

    //Apply Colella and Sekora limiters to parabolic interpolant in jth cell
#ifdef SEKORA_6
    for (int i=il-1; i<=iu; ++i){ 
      Real& dx_i   = pco->dx1v(i); 
      qa = (qplus_(i) - q(n,k,j,i))*(q(n,k,j,i) - qminus_(i));
      qd = qplus_(i)-qminus_(i);
      qe = 6.0*(q(n,k,j,i) - 0.5*(qplus_(i) + qminus_(i))); //this is a6 parabolic coefficient          
      //eq 20 additional local extrema check
      qb = (q(n,k,j,i-1) - q(n,k,j,i))*(q(n,k,j,i) - q(n,k,j,i+1));

      // Local extrema monotonization, eq 20
      if (qa <= 0.0 || qb <= 0.0 ){              
	//Compute various approximations to the second derivative at the cell center
	//Reuse the temporary scratch variables from above
	//centered
	qa = 1.0/(dx_i*dx_i)*(q(n,k,j,i-1) - 2.0*q(n,k,j,i) + q(n,k,j,i+1));
	//left
	qb = 1.0/(dx_i*dx_i)*(q(n,k,j,i-2) - 2.0*q(n,k,j,i-1) + q(n,k,j,i));
	//right
	qc = 1.0/(dx_i*dx_i)*(q(n,k,j,i) - 2.0*q(n,k,j,i+1) + q(n,k,j,i+2));      
	//other centered approx. qe=a6 parabolic coefficient gives measure of second derivative
	qd = -2.0*qe/(dx_i*dx_i);

	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)){//transitive property--> all have same sign
	  //limited second derivative as a nonlinear combination of the 4 approximations
	  qe = SIGN(qd)* std::min(std::min(lim_const*fabs(qa),lim_const*fabs(qb)),std::min(lim_const*fabs(qc),fabs(qd)));
	}
	else {
	  qe =0.0;
	}
	if (qd != 0.0){ // Do not divide by zero
	  qminus_(i) = q(n,k,j,i) + (qminus_(i) - q(n,k,j,i))*qe/qd; 
	  qplus_(i) = q(n,k,j,i) + (qplus_(i) - q(n,k,j,i))*qe/qd; 
	}
	else{
	  qminus_(i) = q(n,k,j,i); 
	  qplus_(i) = q(n,k,j,i); 
	}
      }//end local extrema loop

      //Construct the parabolic interpolant using a standard monotonicity preserving limiter
      else { //Eq 26, less restrictive condition than van Leer limiting
	//alpha_-,+, interpolant coefficients renormalized to the cell avearge
	qa = qminus_(i) - q(n,k,j,i);
	qb = qplus_(i) - q(n,k,j,i);
	// Overshoot L/- state 
	if (fabs(qa) >= 2*fabs(qb)){
	  //eq 25, \delta I_{ext}
	  qc = -qa*qa/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qminus_(i) - q(n,k,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j,i+1) - q(n,k,j,i-1)); 
	  if (qc*qe >= qe*qd)
	    qminus_(i) = q(n,k,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qb));
	}
	// Overshoot R/+ state 
	if (fabs(qb) >= 2*fabs(qa)){
	  //eq 25, \delta I_{ext}
	  qc = -qb*qb/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qplus_(i) - q(n,k,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j,i+1) - q(n,k,j,i-1)); 
	  if (qc*qe >= qe*qd)
	    qplus_(i) = q(n,k,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qa));
	}
      } //end else statement (not a local extrema)      
    }//end loop over i direction, Sekora limiters (setting interpolant values a_+,-)
#endif    

    //Convert a_L, a_R interpolation coefficients to L/R Riemann states
    for (int i=il; i<=iu; ++i){ //same bounds as PLM interpolation final loops
      ql(n,i) = qplus_(i-1); 
      qr(n,i) = qminus_(i);
    } 
  } //end loop over primitive variables
  
  //free 1D arrays
  qplus_.DeleteAthenaArray(); 
  qminus_.DeleteAthenaArray(); 
  dd_.DeleteAthenaArray(); 
  dph_.DeleteAthenaArray(); 
  c1_.DeleteAthenaArray(); 
  c2_.DeleteAthenaArray(); 
  c3_.DeleteAthenaArray(); 
  c4_.DeleteAthenaArray(); 
  c5_.DeleteAthenaArray(); 
  c6_.DeleteAthenaArray(); 
  return;
}

//--------------------------------------------------------------------------------------
//! \fn HydroIntegrator::ReconstructionFuncX2()
//  \brief 

void HydroIntegrator::PiecewiseParabolicX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  AthenaArray<Real> dd_, dph_, c1_, c2_, c3_, c4_, c5_, c6_;
  AthenaArray<Real> dd_jp1, dph_jp1, c1_jp1, c2_jp1, c3_jp1, c4_jp1, c5_jp1, c6_jp1;
  AthenaArray<Real> dd_jm1, dph_jm1, c1_jm1, c2_jm1, c3_jm1, c4_jm1, c5_jm1, c6_jm1;
  AthenaArray<Real> dd_jm2, c1_jm2, c2_jm2; //need average linear slopes to compute nonuniform interpolation at j-3/2 
  Real qa,qb,qc,qd, qe; 

  double lim_const = 1.25; //Colella/Sekora constant used in second derivative limiter, >1, independent of dx

  int ncells1 = iu - il +2*NGHOST; //=264 //(iu+2) - (il-1) +1; 
  //allocate scratch arrays for parabolic interpolation
  dd_.NewAthenaArray(ncells1);
  dd_jp1.NewAthenaArray(ncells1);
  dd_jm1.NewAthenaArray(ncells1);
  dd_jm2.NewAthenaArray(ncells1);
  dph_.NewAthenaArray(ncells1);
  dph_jp1.NewAthenaArray(ncells1);
  dph_jm1.NewAthenaArray(ncells1);
  //
  c1_.NewAthenaArray(ncells1);
  c2_.NewAthenaArray(ncells1);
  c3_.NewAthenaArray(ncells1);
  c4_.NewAthenaArray(ncells1);
  c5_.NewAthenaArray(ncells1);
  c6_.NewAthenaArray(ncells1);
  //
  c1_jp1.NewAthenaArray(ncells1);
  c2_jp1.NewAthenaArray(ncells1);
  c3_jp1.NewAthenaArray(ncells1);
  c4_jp1.NewAthenaArray(ncells1);
  c5_jp1.NewAthenaArray(ncells1);
  c6_jp1.NewAthenaArray(ncells1);
  //
  c1_jm1.NewAthenaArray(ncells1);
  c2_jm1.NewAthenaArray(ncells1);
  c3_jm1.NewAthenaArray(ncells1);
  c4_jm1.NewAthenaArray(ncells1);
  c5_jm1.NewAthenaArray(ncells1);
  c6_jm1.NewAthenaArray(ncells1);
  //
  c1_jm2.NewAthenaArray(ncells1);
  c2_jm2.NewAthenaArray(ncells1);
  //Temporary array of L/R parabolic interpolation coefficients that follows CMHOG and Colella notation
  AthenaArray<Real> qplus_, qminus_; 
  qplus_.NewAthenaArray(ncells1);
  qminus_.NewAthenaArray(ncells1);
  AthenaArray<Real> qplus_jm1, qminus_jm1; 
  qplus_jm1.NewAthenaArray(ncells1);
  qminus_jm1.NewAthenaArray(ncells1);
  
  //Mesh spacing information
  Coordinates *pco = pmy_hydro->pmy_block->pcoord;
 
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){ //this is all redundant along x1 array, if you assume dx2 does not change in i, only j
      //This loop computes the coefficients for nonuniform van Leer limited linear slopes
      Real& dx_jm2 = pco->dx2v(j-2); 
      Real& dx_jm1 = pco->dx2v(j-1);
      Real& dx_j   = pco->dx2v(j);
      Real& dx_jp1   = pco->dx2v(j+1);

      //The following are needed to compute coefficients for j-1, j+1 van Leer
      Real& dx_jm3 = pco->dx2v(j-3); 
      Real& dx_jp2 = pco->dx2v(j+2); 

      //j-2
      qa = dx_jm2/(dx_jm3 + dx_jm2 + dx_jm1); //first factor in eq 1.7
      c1_jm2(i) = qa*(2.0*dx_jm3+dx_jm2)/(dx_jm1 + dx_jm2); 
      c2_jm2(i) = qa*(2.0*dx_jm1+dx_jm2)/(dx_jm3 + dx_jm2); 
      //j-1
      qa = dx_jm1/(dx_jm2 + dx_jm1 + dx_j); //first factor in eq 1.7
      c1_jm1(i) = qa*(2.0*dx_jm2+dx_jm1)/(dx_j + dx_jm1); 
      c2_jm1(i) = qa*(2.0*dx_j+dx_jm1)/(dx_jm2 + dx_jm1); 
      //j
      qa = dx_j/(dx_jm1 + dx_j + dx_jp1); //first factor in eq 1.7
      c1_(i) = qa*(2.0*dx_jm1+dx_j)/(dx_jp1 + dx_j); 
      c2_(i) = qa*(2.0*dx_jp1+dx_j)/(dx_jm1 + dx_j); 
      //j+1
      qa = dx_jp1/(dx_j + dx_jp1 + dx_jp2); //first factor in eq 1.7
      c1_jp1(i) = qa*(2.0*dx_j+dx_jp1)/(dx_jp2 + dx_jp1); 
      c2_jp1(i) = qa*(2.0*dx_jp2+dx_jp1)/(dx_j + dx_jp1); 
    }
    //The q* are used to construct the interpolation coefficients c*_
    //All of these variables only depend on the mesh spacing (and can be reused for all prim)
    for (int i=il; i<=iu; ++i){
      Real& dx_jm2 = pco->dx2v(j-2); 
      Real& dx_jm1 = pco->dx2v(j-1);
      Real& dx_j   = pco->dx2v(j);
      Real& dx_jp1   = pco->dx2v(j+1);
      //The following are needed to compute coefficients for j-3/2, j+1/2
      Real& dx_jm3 = pco->dx2v(j-3); 
      Real& dx_jp2 = pco->dx2v(j+2); 
      //j-3/2
      qa = dx_jm3 + dx_jm2 + dx_jm1 + dx_j;
      qb = dx_jm2/(dx_jm2 + dx_jm1); 
      qc = (dx_jm3 + dx_jm2)/(2.0*dx_jm2 + dx_jm1); 
      qd = (dx_j + dx_jm1)/(2.0*dx_jm1 + dx_jm2); 
      qb = qb + 2.0*dx_jm1*qb/qa*(qc-qd); 
      c3_jm1(i) = 1.0 - qb; 
      c4_jm1(i) = qb; 
      c5_jm1(i) = dx_jm1/qa*qd; 
      c6_jm1(i) = -dx_jm2/qa*qc; 
      //j-1/2
      qa = dx_jm2 + dx_jm1 + dx_j + dx_jp1;
      qb = dx_jm1/(dx_jm1 + dx_j); 
      qc = (dx_jm2 + dx_jm1)/(2.0*dx_jm1 + dx_j); 
      qd = (dx_jp1 + dx_j)/(2.0*dx_j + dx_jm1); 
      qb = qb + 2.0*dx_j*qb/qa*(qc-qd); 
      c3_(i) = 1.0 - qb; 
      c4_(i) = qb; 
      c5_(i) = dx_j/qa*qd; 
      c6_(i) = -dx_jm1/qa*qc; 
      //j+1/2
      qa = dx_jm1 + dx_j + dx_jp1 + dx_jp2;
      qb = dx_j/(dx_j + dx_jp1); 
      qc = (dx_jm1 + dx_j)/(2.0*dx_j + dx_jp1); 
      qd = (dx_jp2 + dx_jp1)/(2.0*dx_jp1 + dx_j); 
      qb = qb + 2.0*dx_jp1*qb/qa*(qc-qd); 
      c3_jp1(i) = 1.0 - qb; 
      c4_jp1(i) = qb; 
      c5_jp1(i) = dx_jp1/qa*qd; 
      c6_jp1(i) = -dx_j/qa*qc; 
    }

    for (int i=il; i<=iu; ++i){
      // Compute average linear slopes in cell (eqn 1.7)
      // Only used in this loop to compute and monotonize average zone slope, dd_(i)
      //j-2
      Real dplus_jm2 = q(n,k,j-1,i) - q(n,k,j-2,i); 
      Real dmnus_jm2 = q(n,k,j-2,i) - q(n,k,j-3,i);
      dd_jm2(i) = c1_jm2(i)*dplus_jm2 + c2_jm2(i)*dmnus_jm2;
      //j-1
      Real dplus_jm1 = q(n,k,j,i) - q(n,k,j-1,i); 
      Real dmnus_jm1 = q(n,k,j-1,i) - q(n,k,j-2,i);
      dd_jm1(i) = c1_jm1(i)*dplus_jm1 + c2_jm1(i)*dmnus_jm1;
      //j
      Real dplus = q(n,k,j+1,i) - q(n,k,j,i); 
      Real dmnus = q(n,k,j,i) - q(n,k,j-1,i);
      dd_(i) = c1_(i)*dplus + c2_(i)*dmnus;
      //j+1
      Real dplus_jp1 = q(n,k,j+2,i) - q(n,k,j+1,i); 
      Real dmnus_jp1 = q(n,k,j+1,i) - q(n,k,j,i);
      dd_jp1(i) = c1_jp1(i)*dplus_jp1 + c2_jp1(i)*dmnus_jp1;
      
      // Monotonize, van Leer limiting (eqn 1.8)
#ifdef VL_LIMITER
      //j-2
      qa = std::min( fabs(dd_jm2(i)), 2.0*std::min(fabs(dmnus_jm2),fabs(dplus_jm2)) ); 
      dd_jm2(i) = 0.0;
      if (dplus_jm2*dmnus_jm2 > 0.0) 
	dd_jm2(i) = qa*SIGN(dd_jm2(i));  
      //j-1
      qa = std::min( fabs(dd_jm1(i)), 2.0*std::min(fabs(dmnus_jm1),fabs(dplus_jm1)) ); 
      dd_jm1(i) = 0.0;
      if (dplus_jm1*dmnus_jm1 > 0.0) 
	dd_jm1(i) = qa*SIGN(dd_jm1(i));  
      //j
      qa = std::min( fabs(dd_(i)), 2.0*std::min(fabs(dmnus),fabs(dplus)) ); 
      dd_(i) = 0.0;
      if (dplus*dmnus > 0.0) 
	dd_(i) = qa*SIGN(dd_(i));  
      //j+1
      qa = std::min( fabs(dd_jp1(i)), 2.0*std::min(fabs(dmnus_jp1),fabs(dplus_jp1)) ); 
      dd_jp1(i) = 0.0;
      if (dplus_jp1*dmnus_jp1 > 0.0) 
      dd_jp1(i) = qa*SIGN(dd_jp1(i));   
#endif
    }
    
    // Initialize interface values (eqn 1.6)
    // a_j-1/2 value that is assigned to most L/R Riemann states
    for (int i=il; i<=iu; ++i){ //need the first two ghost cells at far end of domain to do reconstruction on A_l / R Riemann state at last boundary
      Real& dx_j   = pco->dx2v(j);
      //Interpolation function for nonuniform mesh spacing
#if defined(VL_LIMITER) || defined(MIGNONE)
      //Nonuniform mesh option disabled for x2, x3 directions
      dph_jm1(i)= c3_jm1(i)*q(n,k,j-2,i) + c4_jm1(i)*q(n,k,j-1,i) + c5_jm1(i)*dd_jm2(i) + c6_jm1(i)*dd_jm1(i); 
      dph_(i)= c3_(i)*q(n,k,j-1,i) + c4_(i)*q(n,k,j,i) + c5_(i)*dd_jm1(i) + c6_(i)*dd_(i); 
      dph_jp1(i)= c3_jp1(i)*q(n,k,j,i) + c4_jp1(i)*q(n,k,j+1,i) + c5_jp1(i)*dd_(i) + c6_jp1(i)*dd_jp1(i); 
      //Simplified expression for uniform mesh, no van Leer monotonization of average zone slopes
      /*      dph_jm1(i)= 7.0/12.0*(q(n,k,j-2,i) + q(n,k,j-1,i)) - 1.0/12.0*(q(n,k,j,i) + q(n,k,j-3,i));
      dph_(i)= 7.0/12.0*(q(n,k,j-1,i) + q(n,k,j,i)) - 1.0/12.0*(q(n,k,j+1,i) + q(n,k,j-2,i));
      dph_jp1(i)= 7.0/12.0*(q(n,k,j,i) + q(n,k,j+1,i)) - 1.0/12.0*(q(n,k,j+2,i) + q(n,k,j-1,i));  */
#endif

      //6th order approximation to interface state, uniform mesh spacing
      //Colella Sekora 2008 equation 17 
#ifdef SEKORA_6
      // eq 17, j -=2
      dph_jm1(i) = 37.0/60.0*(q(n,k,j-2,i) + q(n,k,j-1,i)) - 2.0/15.0*(q(n,k,j-3,i)+ q(n,k,j,i)) +
	1.0/60.0*(q(n,k,j-4,i)+ q(n,k,j+1,i));
      // eq 17, j -=1
      dph_(i) = 37.0/60.0*(q(n,k,j-1,i) + q(n,k,j,i)) - 2.0/15.0*(q(n,k,j-2,i)+ q(n,k,j+1,i)) +
	1.0/60.0*(q(n,k,j-3,i)+ q(n,k,j+2,i));
      // eq 17
      dph_jp1(i) = 37.0/60.0*(q(n,k,j,i) + q(n,k,j+1,i)) - 2.0/15.0*(q(n,k,j-1,i)+ q(n,k,j+2,i)) +
	1.0/60.0*(q(n,k,j-2,i)+ q(n,k,j+3,i));
      // Colella Sekora eqn 13: Check monotonization constraint for j-3/2 boundary
      if (std::min(q(n,k,j-2,i), q(n,k,j-1,i)) > dph_jm1(i) || std::max(q(n,k,j-2,i), q(n,k,j-1,i)) < dph_jm1(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_j*dx_j)*(q(n,k,j-2,i) - 2.0*dph_jm1(i) + q(n,k,j-1,i)); 
	//left
	qb = 1.0/(dx_j*dx_j)*(q(n,k,j-3,i) - 2.0*q(n,k,j-2,i) + q(n,k,j-1,i)); 
	//right
	qc = 1.0/(dx_j*dx_j)*(q(n,k,j-2,i) - 2.0*q(n,k,j-1,i) + q(n,k,j,i)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_jm1(i) = 0.5*(q(n,k,j-2,i) + q(n,k,j-1,i)) - dx_j*dx_j/3.0 * qd; 
      } //end sekora limiter for extrema at j-3/2 interface

      // Colella Sekora eqn 13: Check monotonization constraint for j-1/2 boundary
      if (std::min(q(n,k,j-1,i), q(n,k,j,i)) > dph_(i) || std::max(q(n,k,j-1,i), q(n,k,j,i)) < dph_(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_j*dx_j)*(q(n,k,j-1,i) - 2.0*dph_(i) + q(n,k,j,i)); 
	//left
	qb = 1.0/(dx_j*dx_j)*(q(n,k,j-2,i) - 2.0*q(n,k,j-1,i) + q(n,k,j,i)); 
	//right
	qc = 1.0/(dx_j*dx_j)*(q(n,k,j-1,i) - 2.0*q(n,k,j,i) + q(n,k,j+1,i)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_(i) = 0.5*(q(n,k,j-1,i) + q(n,k,j,i)) - dx_j*dx_j/3.0 * qd; 
      } //end sekora limiter for extrema at j-1/2 interface

      // Colella Sekora eqn 13: Check monotonization constraint for j+3/2 boundary
      if (std::min(q(n,k,j,i), q(n,k,j+1,i)) > dph_jp1(i) || std::max(q(n,k,j,i), q(n,k,j+1,i)) < dph_jp1(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_j*dx_j)*(q(n,k,j,i) - 2.0*dph_jp1(i) + q(n,k,j+1,i)); 
	//left
	qb = 1.0/(dx_j*dx_j)*(q(n,k,j-1,i) - 2.0*q(n,k,j,i) + q(n,k,j+1,i)); 
	//right
	qc = 1.0/(dx_j*dx_j)*(q(n,k,j,i) - 2.0*q(n,k,j+1,i) + q(n,k,j+2,i)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_jp1(i) = 0.5*(q(n,k,j,i) + q(n,k,j+1,i)) - dx_j*dx_j/3.0 * qd; 
      } //end sekora limiter for extrema at j-1/2 interface
#endif      
    } //end loop for interpolation of interface states

    //Initialize parabolic interpolation coefficients a_{L,j} = a_{j-1/2}, a_{R,j} = a_{j+1/2}
    for (int i=il; i<=iu; ++i){ 
      qminus_(i) = dph_(i);
      qplus_(i) = dph_jp1(i);

    /* This is where things get tricky in x2 reconstruction along x1 slices.
       Need to interpolate j+1/2 boundary for reconstruction in j. 
       Need to do entire reconstruction (interpolation and limiting) for j-1 cell.
       Currently, these are wasted at end of routine and recomputed in subsequent
       call to PPMx2() in van_leer2.cpp */

      //Initialize parabolic interpolation coefficients a_{L,j-1} = a_{j-3/2}, a_{R,j-1} = a_{j-1/2}
      qminus_jm1(i) = dph_jm1(i);
      qplus_jm1(i) = dph_(i);
      
      //Alternative monotonization: Mignone 2014, eq 45. Do not need to check above monotonization conditions
      //since minmax assignments automatically adjust violations while preserving normal cases
#ifdef MIGNONE
      //j-3/2 boundary
      qminus_jm1(i) = std::min(qminus_jm1(i), std::max(q(n,k,j-2,i),q(n,k,j-1,i)));
      qplus_jm1(i) = std::min(qplus_jm1(i), std::max(q(n,k,j-1,i),q(n,k,j,i))); 

      qminus_jm1(i) = std::max(qminus_jm1(i), std::min(q(n,k,j-2,i),q(n,k,j-1,i)));
      qplus_jm1(i) = std::max(qplus_jm1(i), std::min(q(n,k,j-1,i),q(n,k,j,i)));  

      //j-1/2 boundary
      qminus_(i) = std::min(qminus_(i), std::max(q(n,k,j-1,i),q(n,k,j,i)));
      qplus_(i) = std::min(qplus_(i), std::max(q(n,k,j+1,i),q(n,k,j,i))); 

      qminus_(i) = std::max(qminus_(i), std::min(q(n,k,j-1,i),q(n,k,j,i)));
      qplus_(i) = std::max(qplus_(i), std::min(q(n,k,j+1,i),q(n,k,j,i)));  
#endif
    } 
    
#ifndef SEKORA_6
    // Apply original PPM limiters to jth cell reconstruction 
    for (int i=il; i<=iu; ++i){ 
      qa = (qplus_(i) - q(n,k,j,i))*(q(n,k,j,i) - qminus_(i));
      qd = qplus_(i)-qminus_(i);
      qe = 6.0*(q(n,k,j,i) - 0.5*(qplus_(i) + qminus_(i))); //this is a6 parabolic coefficient

      //Parabolic overshoot monotonization (eqn 1.10). They are mutually exclusive
      //Alternative formulation of condition commented out. See Colella Sekora
      if ( (qd*(qd - qe)) < 0.0) { 	//	  if (fabs(qminus_(i) - q(n,k,j,i)) >= 2*fabs(qplus_(i) - q(n,k,j,i))){	
	qminus_(i) = 3.0*q(n,k,j,i) - 2.0*qplus_(i); 
      }      
      if ( (qd*(qd + qe)) < 0.0) {  //    if (fabs(qplus_(i) - q(n,k,j,i)) >= 2*fabs(qminus_(i) - q(n,k,j,i))){
	qplus_(i) = 3.0*q(n,k,j,i) - 2.0*qminus_(i);
      }
      // Local extrema monotonization (eqn 1.10)
      if (qa <= 0.0){
	qminus_(i) = q(n,k,j,i);
	qplus_(i) = q(n,k,j,i);
      }      

      // Apply original PPM limiters to j-1 cell reconstruction 
      //bug here somewhere
      qa = (qplus_jm1(i) - q(n,k,j-1,i))*(q(n,k,j-1,i) - qminus_jm1(i));
      qd = qplus_jm1(i)-qminus_jm1(i);
      qe = 6.0*(q(n,k,j-1,i) - 0.5*(qplus_jm1(i) + qminus_jm1(i))); //this is a6 parabolic coefficient
      //Parabolic overshoot monotonization (eqn 1.10). They are mutually exclusive
      //Alternative formulation of condition commented out. See Colella Sekora
      if ( (qd*(qd - qe)) < 0.0) { 	//	  if (fabs(qminus_jm1(i) - q(n,k,j-1,i)) >= 2*fabs(qplus_jm1(i) - q(n,k,j-1,i))){	
	qminus_jm1(i) = 3.0*q(n,k,j-1,i) - 2.0*qplus_jm1(i); 
      }      
      if ( (qd*(qd + qe)) < 0.0) {  //    if (fabs(qplus_jm1(i) - q(n,k,j-1,i)) >= 2*fabs(qminus_jm1(i) - q(n,k,j-1,i))){
	qplus_jm1(i) = 3.0*q(n,k,j-1,i) - 2.0*qminus_jm1(i);
      }
      // Local extrema monotonization (eqn 1.10)
      if (qa <= 0.0){
	qminus_jm1(i) = q(n,k,j-1,i);
	qplus_jm1(i) = q(n,k,j-1,i);
      }      

      } //end loop over i, original PPM limiters
#endif    


#ifdef SEKORA_6
    //Apply Colella and Sekora limiters to parabolic interpolant in jth cell
    for (int i=il; i<=iu; ++i){ 
      Real& dx_j   = pco->dx2v(j);
      qa = (qplus_(i) - q(n,k,j,i))*(q(n,k,j,i) - qminus_(i));
      qd = qplus_(i)-qminus_(i);
      qe = 6.0*(q(n,k,j,i) - 0.5*(qplus_(i) + qminus_(i))); //this is a6 parabolic coefficient          
      //eq 20 additional local extrema check
      qb = (q(n,k,j-1,i) - q(n,k,j,i))*(q(n,k,j,i) - q(n,k,j+1,i));

      // Local extrema monotonization, eq 20
      if (qa <= 0.0 || qb <= 0.0 ){              
	//Compute various approximations to the second derivative at the cell center
	//Reuse the temporary scratch variables from above
	//centered
	qa = 1.0/(dx_j*dx_j)*(q(n,k,j-1,i) - 2.0*q(n,k,j,i) + q(n,k,j+1,i));
	//left
	qb = 1.0/(dx_j*dx_j)*(q(n,k,j-2,i) - 2.0*q(n,k,j-1,i) + q(n,k,j,i));
	//right
	qc = 1.0/(dx_j*dx_j)*(q(n,k,j,i) - 2.0*q(n,k,j+1,i) + q(n,k,j+2,i));      
	//other centered approx. qe=a6 parabolic coefficient gives measure of second derivative
	qd = -2.0*qe/(dx_j*dx_j);

	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)){//transitive property--> all have same sign
	  //limited second derivative as a nonlinear combination of the 4 approximations
	  qe = SIGN(qd)* std::min(std::min(lim_const*fabs(qa),lim_const*fabs(qb)),std::min(lim_const*fabs(qc),fabs(qd)));
	}
	else {
	  qe =0.0;
	}
	if (qd != 0.0){ // Do not divide by zero
	  qminus_(i) = q(n,k,j,i) + (qminus_(i) - q(n,k,j,i))*qe/qd; 
	  qplus_(i) = q(n,k,j,i) + (qplus_(i) - q(n,k,j,i))*qe/qd; 
	}
	else{
	  qminus_(i) = q(n,k,j,i); 
	  qplus_(i) = q(n,k,j,i); 
	}
      }//end local extrema loop

      //Construct the parabolic interpolant using a standard monotonicity preserving limiter
      else { //Eq 26, less restrictive condition than van Leer limiting
	//alpha_-,+, interpolant coefficients renormalized to the cell avearge
	qa = qminus_(i) - q(n,k,j,i);
	qb = qplus_(i) - q(n,k,j,i);
	// Overshoot L/- state 
	if (fabs(qa) >= 2*fabs(qb)){
	  //eq 25, \delta I_{ext}
	  qc = -qa*qa/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qminus_(i) - q(n,k,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j+1,i) - q(n,k,j-1,i)); 
	  if (qc*qe >= qe*qd)
	    qminus_(i) = q(n,k,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qb));
	}
	// Overshoot R/+ state 
	if (fabs(qb) >= 2*fabs(qa)){
	  //eq 25, \delta I_{ext}
	  qc = -qb*qb/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qplus_(i) - q(n,k,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j+1,i) - q(n,k,j-1,i)); 
	  if (qc*qe >= qe*qd)
	    qplus_(i) = q(n,k,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qa));
	}
      } //end else statement (not a local extrema, j cell)      

      //Apply Colella and Sekora limiters to parabolic interpolant in j-1 cell
      qa = (qplus_jm1(i) - q(n,k,j-1,i))*(q(n,k,j-1,i) - qminus_jm1(i));
      qd = qplus_jm1(i)-qminus_jm1(i);
      qe = 6.0*(q(n,k,j-1,i) - 0.5*(qplus_jm1(i) + qminus_jm1(i))); //this is a6 parabolic coefficient          
      //eq 20 additional local extrema check
      qb = (q(n,k,j-2,i) - q(n,k,j-1,i))*(q(n,k,j-1,i) - q(n,k,j,i));

      // Local extrema monotonization, eq 20
      if (qa <= 0.0 || qb <= 0.0 ){              
	//Compute various approximations to the second derivative at the cell center
	//Reuse the temporary scratch variables from above
	//centered
	qa = 1.0/(dx_j*dx_j)*(q(n,k,j-2,i) - 2.0*q(n,k,j-2,i) + q(n,k,j,i));
	//left
	qb = 1.0/(dx_j*dx_j)*(q(n,k,j-3,i) - 2.0*q(n,k,j-1,i) + q(n,k,j-1,i));
	//right
	qc = 1.0/(dx_j*dx_j)*(q(n,k,j-1,i) - 2.0*q(n,k,j,i) + q(n,k,j+1,i));      
	//other centered approx. qe=a6 parabolic coefficient gives measure of second derivative
	qd = -2.0*qe/(dx_j*dx_j);

	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)){//transitive property--> all have same sign
	  //limited second derivative as a nonlinear combination of the 4 approximations
	  qe = SIGN(qd)* std::min(std::min(lim_const*fabs(qa),lim_const*fabs(qb)),std::min(lim_const*fabs(qc),fabs(qd)));
	}
	else {
	  qe =0.0;
	}
	if (qd != 0.0){ // Do not divide by zero
	  qminus_jm1(i) = q(n,k,j-1,i) + (qminus_jm1(i) - q(n,k,j-1,i))*qe/qd; 
	  qplus_jm1(i) = q(n,k,j-1,i) + (qplus_jm1(i) - q(n,k,j-1,i))*qe/qd; 
	}
	else{
	  qminus_jm1(i) = q(n,k,j-1,i); 
	  qplus_jm1(i) = q(n,k,j-1,i); 
	}
      }//end local extrema loop

      //Construct the parabolic interpolant using a standard monotonicity preserving limiter
      else { //Eq 26, less restrictive condition than van Leer limiting
	//alpha_-,+, interpolant coefficients renormalized to the cell avearge
	qa = qminus_jm1(i) - q(n,k,j-1,i);
	qb = qplus_jm1(i) - q(n,k,j-1,i);
	// Overshoot L/- state 
	if (fabs(qa) >= 2*fabs(qb)){
	  //eq 25, \delta I_{ext}
	  qc = -qa*qa/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qminus_jm1(i) - q(n,k,j-1,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j,i) - q(n,k,j-2,i)); 
	  if (qc*qe >= qe*qd)
	    qminus_jm1(i) = q(n,k,j-1,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qb));
	}
	// Overshoot R/+ state 
	if (fabs(qb) >= 2*fabs(qa)){
	  //eq 25, \delta I_{ext}
	  qc = -qb*qb/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qplus_jm1(i) - q(n,k,j-1,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j,i) - q(n,k,j-2,i)); 
	  if (qc*qe >= qe*qd)
	    qplus_jm1(i) = q(n,k,j-1,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qa));
	}
      } //end else statement (not a local extrema, j-1)      
    }//end loop over i direction, Sekora limiters (setting interpolant values a_+,-)
#endif    

    //Convert a_L, a_R interpolation coefficients to L/R Riemann states
    for (int i=il; i<=iu; ++i){ //same bounds as PLM interpolation final loops
      ql(n,i) = qplus_jm1(i); 
      qr(n,i) = qminus_(i);
    } 
  } //end loop over primitive variables
  
  //free 1D arrays
  qplus_.DeleteAthenaArray(); 
  qminus_.DeleteAthenaArray(); 
  qplus_jm1.DeleteAthenaArray(); 
  qminus_jm1.DeleteAthenaArray(); 
  //
  dd_.DeleteAthenaArray(); 
  dd_jp1.DeleteAthenaArray(); 
  dd_jm1.DeleteAthenaArray(); 
  dd_jm2.DeleteAthenaArray(); 
  dph_jm1.DeleteAthenaArray(); 
  dph_.DeleteAthenaArray(); 
  dph_jp1.DeleteAthenaArray(); 
  c1_.DeleteAthenaArray(); 
  c2_.DeleteAthenaArray(); 
  c3_.DeleteAthenaArray(); 
  c4_.DeleteAthenaArray(); 
  c5_.DeleteAthenaArray(); 
  c6_.DeleteAthenaArray(); 
  //
  c1_jp1.DeleteAthenaArray(); 
  c2_jp1.DeleteAthenaArray(); 
  c3_jp1.DeleteAthenaArray(); 
  c4_jp1.DeleteAthenaArray(); 
  c5_jp1.DeleteAthenaArray(); 
  c6_jp1.DeleteAthenaArray(); 
  //
  c1_jm1.DeleteAthenaArray(); 
  c2_jm1.DeleteAthenaArray(); 
  c3_jm1.DeleteAthenaArray(); 
  c4_jm1.DeleteAthenaArray(); 
  c5_jm1.DeleteAthenaArray(); 
  c6_jm1.DeleteAthenaArray(); 
  //
  c1_jm2.DeleteAthenaArray(); 
  c2_jm2.DeleteAthenaArray(); 


  return;
}

//--------------------------------------------------------------------------------------
//! \fn HydroIntegrator::ReconstructionFuncX3()
//  \brief 

void HydroIntegrator::PiecewiseParabolicX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  AthenaArray<Real> dd_, dph_, c1_, c2_, c3_, c4_, c5_, c6_;
  AthenaArray<Real> dd_kp1, dph_kp1, c1_kp1, c2_kp1, c3_kp1, c4_kp1, c5_kp1, c6_kp1;
  AthenaArray<Real> dd_km1, dph_km1, c1_km1, c2_km1, c3_km1, c4_km1, c5_km1, c6_km1;
  AthenaArray<Real> dd_km2, c1_km2, c2_km2; //need average linear slopes to compute nonuniform interpolation at j-3/2 
  Real qa,qb,qc,qd, qe; 

  double lim_const = 1.25; //Colella/Sekora constant used in second derivative limiter, >1, independent of dx

  int ncells1 = iu - il +2*NGHOST; //=264 //(iu+2) - (il-1) +1; 
  //allocate scratch arrays for parabolic interpolation
  dd_.NewAthenaArray(ncells1);
  dd_kp1.NewAthenaArray(ncells1);
  dd_km1.NewAthenaArray(ncells1);
  dd_km2.NewAthenaArray(ncells1);
  dph_.NewAthenaArray(ncells1);
  dph_kp1.NewAthenaArray(ncells1);
  dph_km1.NewAthenaArray(ncells1);
  //
  c1_.NewAthenaArray(ncells1);
  c2_.NewAthenaArray(ncells1);
  c3_.NewAthenaArray(ncells1);
  c4_.NewAthenaArray(ncells1);
  c5_.NewAthenaArray(ncells1);
  c6_.NewAthenaArray(ncells1);
  //
  c1_kp1.NewAthenaArray(ncells1);
  c2_kp1.NewAthenaArray(ncells1);
  c3_kp1.NewAthenaArray(ncells1);
  c4_kp1.NewAthenaArray(ncells1);
  c5_kp1.NewAthenaArray(ncells1);
  c6_kp1.NewAthenaArray(ncells1);
  //
  c1_km1.NewAthenaArray(ncells1);
  c2_km1.NewAthenaArray(ncells1);
  c3_km1.NewAthenaArray(ncells1);
  c4_km1.NewAthenaArray(ncells1);
  c5_km1.NewAthenaArray(ncells1);
  c6_km1.NewAthenaArray(ncells1);
  //
  c1_km2.NewAthenaArray(ncells1);
  c2_km2.NewAthenaArray(ncells1);
  //Temporary array of L/R parabolic interpolation coefficients that follows CMHOG and Colella notation
  AthenaArray<Real> qplus_, qminus_; 
  qplus_.NewAthenaArray(ncells1);
  qminus_.NewAthenaArray(ncells1);
  AthenaArray<Real> qplus_km1, qminus_km1; 
  qplus_km1.NewAthenaArray(ncells1);
  qminus_km1.NewAthenaArray(ncells1);
  
  //Mesh spacing information
  Coordinates *pco = pmy_hydro->pmy_block->pcoord;
 
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){ //this is all redundant along x1 array, if you assume dx2 does not change in i, only j
      //This loop computes the coefficients for nonuniform van Leer limited linear slopes
      Real& dx_km2 = pco->dx3v(k-2); 
      Real& dx_km1 = pco->dx3v(k-1);
      Real& dx_k   = pco->dx3v(k);
      Real& dx_kp1   = pco->dx3v(k+1);

      //The following are needed to compute coefficients for k-1, k+1 van Leer
      Real& dx_km3 = pco->dx3v(k-3); 
      Real& dx_kp2 = pco->dx3v(k+2); 

      //k-2
      qa = dx_km2/(dx_km3 + dx_km2 + dx_km1); //first factor in eq 1.7
      c1_km2(i) = qa*(2.0*dx_km3+dx_km2)/(dx_km1 + dx_km2); 
      c2_km2(i) = qa*(2.0*dx_km1+dx_km2)/(dx_km3 + dx_km2); 
      //k-1
      qa = dx_km1/(dx_km2 + dx_km1 + dx_k); //first factor in eq 1.7
      c1_km1(i) = qa*(2.0*dx_km2+dx_km1)/(dx_k + dx_km1); 
      c2_km1(i) = qa*(2.0*dx_k+dx_km1)/(dx_km2 + dx_km1); 
      //k
      qa = dx_k/(dx_km1 + dx_k + dx_kp1); //first factor in eq 1.7
      c1_(i) = qa*(2.0*dx_km1+dx_k)/(dx_kp1 + dx_k); 
      c2_(i) = qa*(2.0*dx_kp1+dx_k)/(dx_km1 + dx_k); 
      //k+1
      qa = dx_kp1/(dx_k + dx_kp1 + dx_kp2); //first factor in eq 1.7
      c1_kp1(i) = qa*(2.0*dx_k+dx_kp1)/(dx_kp2 + dx_kp1); 
      c2_kp1(i) = qa*(2.0*dx_kp2+dx_kp1)/(dx_k + dx_kp1); 
    }
    //The q* are used to construct the interpolation coefficients c*_
    //All of these variables only depend on the mesh spacing (and can be reused for all prim)
    for (int i=il; i<=iu; ++i){
      Real& dx_km2 = pco->dx3v(k-2); 
      Real& dx_km1 = pco->dx3v(k-1);
      Real& dx_k   = pco->dx3v(k);
      Real& dx_kp1   = pco->dx3v(k+1);

      //The following are needed to compute coefficients for k-1, k+1 van Leer
      Real& dx_km3 = pco->dx3v(k-3); 
      Real& dx_kp2 = pco->dx3v(k+2); 
      //k-3/2
      qa = dx_km3 + dx_km2 + dx_km1 + dx_k;
      qb = dx_km2/(dx_km2 + dx_km1); 
      qc = (dx_km3 + dx_km2)/(2.0*dx_km2 + dx_km1); 
      qd = (dx_k + dx_km1)/(2.0*dx_km1 + dx_km2); 
      qb = qb + 2.0*dx_km1*qb/qa*(qc-qd); 
      c3_km1(i) = 1.0 - qb; 
      c4_km1(i) = qb; 
      c5_km1(i) = dx_km1/qa*qd; 
      c6_km1(i) = -dx_km2/qa*qc; 
      //k-1/2
      qa = dx_km2 + dx_km1 + dx_k + dx_kp1;
      qb = dx_km1/(dx_km1 + dx_k); 
      qc = (dx_km2 + dx_km1)/(2.0*dx_km1 + dx_k); 
      qd = (dx_kp1 + dx_k)/(2.0*dx_k + dx_km1); 
      qb = qb + 2.0*dx_k*qb/qa*(qc-qd); 
      c3_(i) = 1.0 - qb; 
      c4_(i) = qb; 
      c5_(i) = dx_k/qa*qd; 
      c6_(i) = -dx_km1/qa*qc; 
      //k+1/2
      qa = dx_km1 + dx_k + dx_kp1 + dx_kp2;
      qb = dx_k/(dx_k + dx_kp1); 
      qc = (dx_km1 + dx_k)/(2.0*dx_k + dx_kp1); 
      qd = (dx_kp2 + dx_kp1)/(2.0*dx_kp1 + dx_k); 
      qb = qb + 2.0*dx_kp1*qb/qa*(qc-qd); 
      c3_kp1(i) = 1.0 - qb; 
      c4_kp1(i) = qb; 
      c5_kp1(i) = dx_kp1/qa*qd; 
      c6_kp1(i) = -dx_k/qa*qc; 
    }

    for (int i=il; i<=iu; ++i){
      // Compute average linear slopes in cell (eqn 1.7)
      // Only used in this loop to compute and monotonize average zone slope, dd_(i)
      //k-2
      Real dplus_km2 = q(n,k-1,j,i) - q(n,k-2,j,i); 
      Real dmnus_km2 = q(n,k-2,j,i) - q(n,k-3,j,i);
      dd_km2(i) = c1_km2(i)*dplus_km2 + c2_km2(i)*dmnus_km2;
      //k-1
      Real dplus_km1 = q(n,k,j,i) - q(n,k-1,j,i); 
      Real dmnus_km1 = q(n,k-1,j,i) - q(n,k-2,j,i);
      dd_km1(i) = c1_km1(i)*dplus_km1 + c2_km1(i)*dmnus_km1;
      //k
      Real dplus = q(n,k+1,j,i) - q(n,k,j,i); 
      Real dmnus = q(n,k,j,i) - q(n,k-1,j,i);
      dd_(i) = c1_(i)*dplus + c2_(i)*dmnus;
      //k+1
      Real dplus_kp1 = q(n,k+2,j,i) - q(n,k+1,j,i); 
      Real dmnus_kp1 = q(n,k+1,j,i) - q(n,k,j,i);
      dd_kp1(i) = c1_kp1(i)*dplus_kp1 + c2_kp1(i)*dmnus_kp1;
      
      // Monotonize, van Leer limiting (eqn 1.8)
#ifdef VL_LIMITER
      //k-2
      qa = std::min( fabs(dd_km2(i)), 2.0*std::min(fabs(dmnus_km2),fabs(dplus_km2)) ); 
      dd_km2(i) = 0.0;
      if (dplus_km2*dmnus_km2 > 0.0) 
	dd_km2(i) = qa*SIGN(dd_km2(i));  
      //k-1
      qa = std::min( fabs(dd_km1(i)), 2.0*std::min(fabs(dmnus_km1),fabs(dplus_km1)) ); 
      dd_km1(i) = 0.0;
      if (dplus_km1*dmnus_km1 > 0.0) 
	dd_km1(i) = qa*SIGN(dd_km1(i));  
      //k
      qa = std::min( fabs(dd_(i)), 2.0*std::min(fabs(dmnus),fabs(dplus)) ); 
      dd_(i) = 0.0;
      if (dplus*dmnus > 0.0) 
	dd_(i) = qa*SIGN(dd_(i));  
      //k+1
      qa = std::min( fabs(dd_kp1(i)), 2.0*std::min(fabs(dmnus_kp1),fabs(dplus_kp1)) ); 
      dd_kp1(i) = 0.0;
      if (dplus_kp1*dmnus_kp1 > 0.0) 
      dd_kp1(i) = qa*SIGN(dd_kp1(i));   
#endif
    }
    
    // Initialize interface values (eqn 1.6)
    // a_k-1/2 value that is assigned to most L/R Riemann states
    for (int i=il; i<=iu; ++i){ //need the first two ghost cells at far end of domain to do reconstruction on A_l / R Riemann state at last boundary
      Real& dx_k   = pco->dx3v(k);
      //Interpolation function for nonuniform mesh spacing
#if defined(VL_LIMITER) || defined(MIGNONE)
      //Nonuniform mesh option disabled for x2, x3 directions
      dph_km1(i)= c3_km1(i)*q(n,k-2,j,i) + c4_km1(i)*q(n,k-1,j,i) + c5_km1(i)*dd_km2(i) + c6_km1(i)*dd_km1(i); 
      dph_(i)= c3_(i)*q(n,k-1,j,i) + c4_(i)*q(n,k,j,i) + c5_(i)*dd_km1(i) + c6_(i)*dd_(i); 
      dph_kp1(i)= c3_kp1(i)*q(n,k,j,i) + c4_kp1(i)*q(n,k+1,j,i) + c5_kp1(i)*dd_(i) + c6_kp1(i)*dd_kp1(i); 
      //Simplified expression for uniform mesh, no van Leer monotonization of average zone slopes
      /*      dph_km1(i)= 7.0/12.0*(q(n,k-2,j,i) + q(n,k-1,j,i)) - 1.0/12.0*(q(n,k,j,i) + q(n,k-3,j,i));
      dph_(i)= 7.0/12.0*(q(n,k-1,j,i) + q(n,k,j,i)) - 1.0/12.0*(q(n,k+1,j,i) + q(n,k-2,j,i));
      dph_kp1(i)= 7.0/12.0*(q(n,k,j,i) + q(n,k+1,j,i)) - 1.0/12.0*(q(n,k+2,j,i) + q(n,k-1,j,i));  */
#endif

      //6th order approximation to interface state, uniform mesh spacing
      //Colella Sekora 2008 equation 17 
#ifdef SEKORA_6
      // eq 17, k -=2
      dph_km1(i) = 37.0/60.0*(q(n,k-2,j,i) + q(n,k-1,j,i)) - 2.0/15.0*(q(n,k-3,j,i)+ q(n,k,j,i)) +
	1.0/60.0*(q(n,k,j-4,i)+ q(n,k+1,j,i));
      // eq 17, k -=1
      dph_(i) = 37.0/60.0*(q(n,k-1,j,i) + q(n,k,j,i)) - 2.0/15.0*(q(n,k-2,j,i)+ q(n,k+1,j,i)) +
	1.0/60.0*(q(n,k-3,j,i)+ q(n,k+2,j,i));
      // eq 17
      dph_kp1(i) = 37.0/60.0*(q(n,k,j,i) + q(n,k+1,j,i)) - 2.0/15.0*(q(n,k-1,j,i)+ q(n,k+2,j,i)) +
	1.0/60.0*(q(n,k-2,j,i)+ q(n,k+3,j,i));
      // Colella Sekora eqn 13: Check monotonization constraint for k-3/2 boundary
      if (std::min(q(n,k-2,j,i), q(n,k-1,j,i)) > dph_km1(i) || std::max(q(n,k-2,j,i), q(n,k-1,j,i)) < dph_km1(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_k*dx_k)*(q(n,k-2,j,i) - 2.0*dph_km1(i) + q(n,k-1,j,i)); 
	//left
	qb = 1.0/(dx_k*dx_k)*(q(n,k-3,j,i) - 2.0*q(n,k-2,j,i) + q(n,k-1,j,i)); 
	//right
	qc = 1.0/(dx_k*dx_k)*(q(n,k-2,j,i) - 2.0*q(n,k-1,j,i) + q(n,k,j,i)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_km1(i) = 0.5*(q(n,k-2,j,i) + q(n,k-1,j,i)) - dx_k*dx_k/3.0 * qd; 
      } //end sekora limiter for extrema at k-3/2 interface

      // Colella Sekora eqn 13: Check monotonization constraint for k-1/2 boundary
      if (std::min(q(n,k-1,j,i), q(n,k,j,i)) > dph_(i) || std::max(q(n,k-1,j,i), q(n,k,j,i)) < dph_(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_k*dx_k)*(q(n,k-1,j,i) - 2.0*dph_(i) + q(n,k,j,i)); 
	//left
	qb = 1.0/(dx_k*dx_k)*(q(n,k-2,j,i) - 2.0*q(n,k-1,j,i) + q(n,k,j,i)); 
	//right
	qc = 1.0/(dx_k*dx_k)*(q(n,k-1,j,i) - 2.0*q(n,k,j,i) + q(n,k+1,j,i)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_(i) = 0.5*(q(n,k-1,j,i) + q(n,k,j,i)) - dx_k*dx_k/3.0 * qd; 
      } //end sekora limiter for extrema at k-1/2 interface

      // Colella Sekora eqn 13: Check monotonization constraint for j+3/2 boundary
      if (std::min(q(n,k,j,i), q(n,k+1,j,i)) > dph_kp1(i) || std::max(q(n,k,j,i), q(n,k+1,j,i)) < dph_kp1(i)){
	//Compute various approximations to the second derivative at the interface (reusing temp variables)       
	//centered
	qa = 3.0/(dx_k*dx_k)*(q(n,k,j,i) - 2.0*dph_kp1(i) + q(n,k+1,j,i)); 
	//left
	qb = 1.0/(dx_k*dx_k)*(q(n,k-1,j,i) - 2.0*q(n,k,j,i) + q(n,k+1,j,i)); 
	//right
	qc = 1.0/(dx_k*dx_k)*(q(n,k,j,i) - 2.0*q(n,k+1,j,i) + q(n,k+2,j,i)); 
	//built in macro function in defs.hpp
	//#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 ). It returns -1 for 0
	// This is not a problem for the if statement, but may be a problem for the assignment of qd
	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
	  //limited second derivative as a nonlinear combination of the 3 approximations
	  //the left and right approximations should be multiplied by some constant C>1, but not sure what it is
	  qd = SIGN(qa)* std::min(lim_const*fabs(qb),std::min(lim_const*fabs(qc),fabs(qa)));// if all are 0, then qd= -1*0.0
	}
	else {
	  qd = 0.0;
	}
	dph_kp1(i) = 0.5*(q(n,k,j,i) + q(n,k+1,j,i)) - dx_k*dx_k/3.0 * qd; 
      } //end sekora limiter for extrema at k-1/2 interface
#endif      
    } //end loop for interpolation of interface states

    //Initialize parabolic interpolation coefficients a_{L,j} = a_{k-1/2}, a_{R,j} = a_{k+1/2}
    for (int i=il; i<=iu; ++i){ 
      qminus_(i) = dph_(i);
      qplus_(i) = dph_kp1(i);

    /* This is where things get tricky in x2 reconstruction along x1 slices.
       Need to interpolate k+1/2 boundary for reconstruction in j. 
       Need to do entire reconstruction (interpolation and limiting) for j-1 cell.
       Currently, these are wasted at end of routine and recomputed in subsequent
       call to PPMx2() in van_leer2.cpp */

      //Initialize parabolic interpolation coefficients a_{L,j-1} = a_{k-3/2}, a_{R,j-1} = a_{k-1/2}
      qminus_km1(i) = dph_km1(i);
      qplus_km1(i) = dph_(i);
      
      //Alternative monotonization: Mignone 2014, eq 45. Do not need to check above monotonization conditions
      //since minmax assignments automatically adjust violations while preserving normal cases
#ifdef MIGNONE
      //k-3/2 boundary
      qminus_km1(i) = std::min(qminus_km1(i), std::max(q(n,k-2,j,i),q(n,k-1,j,i)));
      qplus_km1(i) = std::min(qplus_km1(i), std::max(q(n,k-1,j,i),q(n,k,j,i))); 

      qminus_km1(i) = std::max(qminus_km1(i), std::min(q(n,k-2,j,i),q(n,k-1,j,i)));
      qplus_km1(i) = std::max(qplus_km1(i), std::min(q(n,k-1,j,i),q(n,k,j,i)));  

      //k-1/2 boundary
      qminus_(i) = std::min(qminus_(i), std::max(q(n,k-1,j,i),q(n,k,j,i)));
      qplus_(i) = std::min(qplus_(i), std::max(q(n,k+1,j,i),q(n,k,j,i))); 

      qminus_(i) = std::max(qminus_(i), std::min(q(n,k-1,j,i),q(n,k,j,i)));
      qplus_(i) = std::max(qplus_(i), std::min(q(n,k+1,j,i),q(n,k,j,i)));  
#endif
    } 
    
#ifndef SEKORA_6
    // Apply original PPM limiters to jth cell reconstruction 
    for (int i=il; i<=iu; ++i){ 
      qa = (qplus_(i) - q(n,k,j,i))*(q(n,k,j,i) - qminus_(i));
      qd = qplus_(i)-qminus_(i);
      qe = 6.0*(q(n,k,j,i) - 0.5*(qplus_(i) + qminus_(i))); //this is a6 parabolic coefficient

      //Parabolic overshoot monotonization (eqn 1.10). They are mutually exclusive
      //Alternative formulation of condition commented out. See Colella Sekora
      if ( (qd*(qd - qe)) < 0.0) { 	//	  if (fabs(qminus_(i) - q(n,k,j,i)) >= 2*fabs(qplus_(i) - q(n,k,j,i))){	
	qminus_(i) = 3.0*q(n,k,j,i) - 2.0*qplus_(i); 
      }      
      if ( (qd*(qd + qe)) < 0.0) {  //    if (fabs(qplus_(i) - q(n,k,j,i)) >= 2*fabs(qminus_(i) - q(n,k,j,i))){
	qplus_(i) = 3.0*q(n,k,j,i) - 2.0*qminus_(i);
      }
      // Local extrema monotonization (eqn 1.10)
      if (qa <= 0.0){
	qminus_(i) = q(n,k,j,i);
	qplus_(i) = q(n,k,j,i);
      }      

      // Apply original PPM limiters to k-1 cell reconstruction 
      //bug here somewhere
      qa = (qplus_km1(i) - q(n,k-1,j,i))*(q(n,k-1,j,i) - qminus_km1(i));
      qd = qplus_km1(i)-qminus_km1(i);
      qe = 6.0*(q(n,k-1,j,i) - 0.5*(qplus_km1(i) + qminus_km1(i))); //this is a6 parabolic coefficient
      //Parabolic overshoot monotonization (eqn 1.10). They are mutually exclusive
      //Alternative formulation of condition commented out. See Colella Sekora
      if ( (qd*(qd - qe)) < 0.0) { 	//	  if (fabs(qminus_km1(i) - q(n,k-1,j,i)) >= 2*fabs(qplus_km1(i) - q(n,k-1,j,i))){	
	qminus_km1(i) = 3.0*q(n,k-1,j,i) - 2.0*qplus_km1(i); 
      }      
      if ( (qd*(qd + qe)) < 0.0) {  //    if (fabs(qplus_km1(i) - q(n,k-1,j,i)) >= 2*fabs(qminus_km1(i) - q(n,k-1,j,i))){
	qplus_km1(i) = 3.0*q(n,k-1,j,i) - 2.0*qminus_km1(i);
      }
      // Local extrema monotonization (eqn 1.10)
      if (qa <= 0.0){
	qminus_km1(i) = q(n,k-1,j,i);
	qplus_km1(i) = q(n,k-1,j,i);
      }      

      } //end loop over i, original PPM limiters
#endif    


#ifdef SEKORA_6
    //Apply Colella and Sekora limiters to parabolic interpolant in jth cell
    for (int i=il; i<=iu; ++i){ 
      Real& dx_k   = pco->dx3v(k);
      qa = (qplus_(i) - q(n,k,j,i))*(q(n,k,j,i) - qminus_(i));
      qd = qplus_(i)-qminus_(i);
      qe = 6.0*(q(n,k,j,i) - 0.5*(qplus_(i) + qminus_(i))); //this is a6 parabolic coefficient          
      //eq 20 additional local extrema check
      qb = (q(n,k-1,j,i) - q(n,k,j,i))*(q(n,k,j,i) - q(n,k+1,j,i));

      // Local extrema monotonization, eq 20
      if (qa <= 0.0 || qb <= 0.0 ){              
	//Compute various approximations to the second derivative at the cell center
	//Reuse the temporary scratch variables from above
	//centered
	qa = 1.0/(dx_k*dx_k)*(q(n,k-1,j,i) - 2.0*q(n,k,j,i) + q(n,k+1,j,i));
	//left
	qb = 1.0/(dx_k*dx_k)*(q(n,k-2,j,i) - 2.0*q(n,k-1,j,i) + q(n,k,j,i));
	//right
	qc = 1.0/(dx_k*dx_k)*(q(n,k,j,i) - 2.0*q(n,k+1,j,i) + q(n,k+2,j,i));      
	//other centered approx. qe=a6 parabolic coefficient gives measure of second derivative
	qd = -2.0*qe/(dx_k*dx_k);

	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)){//transitive property--> all have same sign
	  //limited second derivative as a nonlinear combination of the 4 approximations
	  qe = SIGN(qd)* std::min(std::min(lim_const*fabs(qa),lim_const*fabs(qb)),std::min(lim_const*fabs(qc),fabs(qd)));
	}
	else {
	  qe =0.0;
	}
	if (qd != 0.0){ // Do not divide by zero
	  qminus_(i) = q(n,k,j,i) + (qminus_(i) - q(n,k,j,i))*qe/qd; 
	  qplus_(i) = q(n,k,j,i) + (qplus_(i) - q(n,k,j,i))*qe/qd; 
	}
	else{
	  qminus_(i) = q(n,k,j,i); 
	  qplus_(i) = q(n,k,j,i); 
	}
      }//end local extrema loop

      //Construct the parabolic interpolant using a standard monotonicity preserving limiter
      else { //Eq 26, less restrictive condition than van Leer limiting
	//alpha_-,+, interpolant coefficients renormalized to the cell avearge
	qa = qminus_(i) - q(n,k,j,i);
	qb = qplus_(i) - q(n,k,j,i);
	// Overshoot L/- state 
	if (fabs(qa) >= 2*fabs(qb)){
	  //eq 25, \delta I_{ext}
	  qc = -qa*qa/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qminus_(i) - q(n,k,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k+1,j,i) - q(n,k-1,j,i)); 
	  if (qc*qe >= qe*qd)
	    qminus_(i) = q(n,k,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qb));
	}
	// Overshoot R/+ state 
	if (fabs(qb) >= 2*fabs(qa)){
	  //eq 25, \delta I_{ext}
	  qc = -qb*qb/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qplus_(i) - q(n,k,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k+1,j,i) - q(n,k-1,j,i)); 
	  if (qc*qe >= qe*qd)
	    qplus_(i) = q(n,k,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qa));
	}
      } //end else statement (not a local extrema, j cell)      

      //Apply Colella and Sekora limiters to parabolic interpolant in k-1 cell
      qa = (qplus_km1(i) - q(n,k-1,j,i))*(q(n,k-1,j,i) - qminus_km1(i));
      qd = qplus_km1(i)-qminus_km1(i);
      qe = 6.0*(q(n,k-1,j,i) - 0.5*(qplus_km1(i) + qminus_km1(i))); //this is a6 parabolic coefficient          
      //eq 20 additional local extrema check
      qb = (q(n,k-2,j,i) - q(n,k-1,j,i))*(q(n,k-1,j,i) - q(n,k,j,i));

      // Local extrema monotonization, eq 20
      if (qa <= 0.0 || qb <= 0.0 ){              
	//Compute various approximations to the second derivative at the cell center
	//Reuse the temporary scratch variables from above
	//centered
	qa = 1.0/(dx_k*dx_k)*(q(n,k-2,j,i) - 2.0*q(n,k-2,j,i) + q(n,k,j,i));
	//left
	qb = 1.0/(dx_k*dx_k)*(q(n,k-3,j,i) - 2.0*q(n,k-1,j,i) + q(n,k-1,j,i));
	//right
	qc = 1.0/(dx_k*dx_k)*(q(n,k-1,j,i) - 2.0*q(n,k,j,i) + q(n,k+1,j,i));      
	//other centered approx. qe=a6 parabolic coefficient gives measure of second derivative
	qd = -2.0*qe/(dx_k*dx_k);

	if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)){//transitive property--> all have same sign
	  //limited second derivative as a nonlinear combination of the 4 approximations
	  qe = SIGN(qd)* std::min(std::min(lim_const*fabs(qa),lim_const*fabs(qb)),std::min(lim_const*fabs(qc),fabs(qd)));
	}
	else {
	  qe =0.0;
	}
	if (qd != 0.0){ // Do not divide by zero
	  qminus_km1(i) = q(n,k-1,j,i) + (qminus_km1(i) - q(n,k-1,j,i))*qe/qd; 
	  qplus_km1(i) = q(n,k-1,j,i) + (qplus_km1(i) - q(n,k-1,j,i))*qe/qd; 
	}
	else{
	  qminus_km1(i) = q(n,k-1,j,i); 
	  qplus_km1(i) = q(n,k-1,j,i); 
	}
      }//end local extrema loop

      //Construct the parabolic interpolant using a standard monotonicity preserving limiter
      else { //Eq 26, less restrictive condition than van Leer limiting
	//alpha_-,+, interpolant coefficients renormalized to the cell avearge
	qa = qminus_km1(i) - q(n,k-1,j,i);
	qb = qplus_km1(i) - q(n,k-1,j,i);
	// Overshoot L/- state 
	if (fabs(qa) >= 2*fabs(qb)){
	  //eq 25, \delta I_{ext}
	  qc = -qa*qa/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qminus_km1(i) - q(n,k-1,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j,i) - q(n,k-2,j,i)); 
	  if (qc*qe >= qe*qd)
	    qminus_km1(i) = q(n,k-1,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qb));
	}
	// Overshoot R/+ state 
	if (fabs(qb) >= 2*fabs(qa)){
	  //eq 25, \delta I_{ext}
	  qc = -qb*qb/(4.0*(qa+qb)); 
	  //eq 25, \delta a
	  qd = qplus_km1(i) - q(n,k-1,j,i);	  
	  //eq 24, s
	  qe = SIGN(q(n,k,j,i) - q(n,k-2,j,i)); 
	  if (qc*qe >= qe*qd)
	    qplus_km1(i) = q(n,k-1,j,i) - (2.0*qd + 2.0*qe*sqrt(qd*qd - qd*qa));
	}
      } //end else statement (not a local extrema, k-1)      
    }//end loop over i direction, Sekora limiters (setting interpolant values a_+,-)
#endif    

    //Convert a_L, a_R interpolation coefficients to L/R Riemann states
    for (int i=il; i<=iu; ++i){ //same bounds as PLM interpolation final loops
      ql(n,i) = qplus_km1(i); 
      qr(n,i) = qminus_(i);
    } 
  } //end loop over primitive variables
  
  //free 1D arrays
  qplus_.DeleteAthenaArray(); 
  qminus_.DeleteAthenaArray(); 
  qplus_km1.DeleteAthenaArray(); 
  qminus_km1.DeleteAthenaArray(); 
  //
  dd_.DeleteAthenaArray(); 
  dd_kp1.DeleteAthenaArray(); 
  dd_km1.DeleteAthenaArray(); 
  dd_km2.DeleteAthenaArray(); 
  dph_km1.DeleteAthenaArray(); 
  dph_.DeleteAthenaArray(); 
  dph_kp1.DeleteAthenaArray(); 
  c1_.DeleteAthenaArray(); 
  c2_.DeleteAthenaArray(); 
  c3_.DeleteAthenaArray(); 
  c4_.DeleteAthenaArray(); 
  c5_.DeleteAthenaArray(); 
  c6_.DeleteAthenaArray(); 
  //
  c1_kp1.DeleteAthenaArray(); 
  c2_kp1.DeleteAthenaArray(); 
  c3_kp1.DeleteAthenaArray(); 
  c4_kp1.DeleteAthenaArray(); 
  c5_kp1.DeleteAthenaArray(); 
  c6_kp1.DeleteAthenaArray(); 
  //
  c1_km1.DeleteAthenaArray(); 
  c2_km1.DeleteAthenaArray(); 
  c3_km1.DeleteAthenaArray(); 
  c4_km1.DeleteAthenaArray(); 
  c5_km1.DeleteAthenaArray(); 
  c6_km1.DeleteAthenaArray(); 
  //
  c1_km2.DeleteAthenaArray(); 
  c2_km2.DeleteAthenaArray(); 


  return;
}


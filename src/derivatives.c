/*
##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include "probFun.h"


/*
// To Delete

// Function for computing binomial coefficients
// taken  from this website "https://rosettacode.org/wiki/Evaluate_binomial_coefficients#C"
static unsigned long gcd_ui(unsigned long x, unsigned long y) {
unsigned long t;
if (y < x) { t = x; x = y; y = t; }
while (y > 0) {
t = y;  y = x % y;  x = t;  // y1 <- x0 % y0 ; x1 <- y0
}
return x;
}

unsigned long binomial(unsigned long a, unsigned long b) {
unsigned long n, d, g, r = 1;
n = a + b;
if (a == 0) return 1;
if (a == 1) return n;
if (a >= n) return (a == n);
if (a > n/2) a = n-a;
for (d = 1; d <= a; d++) {
if (r >= ULONG_MAX/n) {  // Possible overflow 
unsigned long nr, dr;  // reduced numerator / denominator 
g = gcd_ui(n, d);  nr = n/g;  dr = d/g;
g = gcd_ui(r, dr);  r = r/g;  dr = dr/g;
if (r >= ULONG_MAX/nr) return 0;  // Unavoidable overflow
r *= nr;
r /= dr;
n--;
} else {
r *= n--;
r /= d;
}
}
return r;
}



// Function for extracting entries of the emission probability matrix
// when the OPGPs are known
double Qentry(int OPGP,double Kaa,double Kab, double Kbb,int elem){
switch(OPGP){
case 1:
if(elem == 1)
return Kbb;
else if ((elem == 2)|(elem == 3))  
return Kab;
else if (elem == 4)
return Kaa;
case 2:
if(elem == 3)
return Kbb;
else if ((elem == 1)|(elem == 4))
return Kab;
else if (elem == 2)
return Kaa;
case 3:
if(elem == 2) 
return Kbb;
else if ((elem == 1)|(elem == 4))
return Kab;
else if (elem == 3)
return Kaa;
case 4:
if(elem == 4) 
return Kbb;
else if ((elem == 2)|(elem == 3))
return Kab;
else if (elem == 1)
return Kaa;
case 5:
if ((elem == 1)|(elem == 2))
return Kab;
else if ((elem == 3)|(elem == 4))
return Kaa;
case 6:
if ((elem == 1)|(elem == 2))
return Kaa;
else if ((elem == 3)|(elem == 4))
return Kab;
case 7:
if ((elem == 1)|(elem == 2))
return Kbb;
else if ((elem == 3)|(elem == 4))
return Kab;
case 8:
if ((elem == 1)|(elem == 2))
return Kab;
else if ((elem == 3)|(elem == 4))
return Kbb;
case 9:
if ((elem == 1)|(elem == 3))
return Kab;
else if ((elem == 2)|(elem == 4))
return Kaa;
case 10:
if ((elem == 1)|(elem == 3))
return Kaa;
else if ((elem == 2)|(elem == 4))
return Kab;
case 11:
if ((elem == 1)|(elem == 3))
return Kbb;
else if ((elem == 2)|(elem == 4))
return Kab;
case 12:
if ((elem == 1)|(elem == 3))
return Kab;
else if ((elem == 2)|(elem == 4))
return Kbb;
case 13:
return Kaa;
case 14:
return Kab;
case 15:
return Kab;
case 16:
return Kbb;
} // end of Switch
return -1;
}


// Function for extracting entries of the emission probability matrix
// when the OPGP are considered the baseline (and so phase is unknown and the r.f's are sex-specific)
double Qentry_up(int config,double Kaa,double Kab, double Kbb,int elem){
switch(config){
case 1:
if(elem == 1)
return Kbb;
else if ((elem == 2)|(elem == 3))  
return Kab;
else if (elem == 4)
return Kaa;
case 2:
if ((elem == 1)|(elem == 2))
return Kab;
else if ((elem == 3)|(elem == 4))
return Kaa;
case 3:
if ((elem == 1)|(elem == 2))
return Kbb;
else if ((elem == 3)|(elem == 4))
return Kab;
case 4:
if ((elem == 1)|(elem == 3))
return Kab;
else if ((elem == 2)|(elem == 4))
return Kaa;
case 5:
if ((elem == 1)|(elem == 3))
return Kbb;
else if ((elem == 2)|(elem == 4))
return Kab;
} // end of Switch
return -1;
}

// Function for returning a specified enetry of the transition matrix for a given recombination fraction value
double Tmat(int s1, int s2, double rval){
int sSum = s1 + s2*4;
if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
return (1-rval)*(1-rval);
else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
return rval*rval;
else
return (1-rval)*rval;
}

// Function for returning a specified enetry of the transition matrix for a given recombination fraction value
// when the r.f.'s are sex-specific 
double Tmat_ss(int s1, int s2, double r_f, double r_m){
int sSum = s1 + s2*4;
if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
return (1-r_f)*(1-r_m);
else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
return r_f*r_m;
else if((sSum == 1)|(sSum == 4)|(sSum == 11)|(sSum == 14))
return (1-r_f)*r_m;
else 
return r_f*(1-r_m);
}


////////////////////
*/

// First derivative function for rf's
double der_rf(int s1, int s2, double rval){
  //double er;
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15)){
    return -2*rval*(1-rval)*(1-2*rval);
  }
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12)){
    return 2*rval*rval*(1-2*rval);
  }
  else{
    return rval*(1-2*rval)*(1-2*rval);
  }
}

// Second derivative function for rf's
double der_rf2(int s1, int s2, double rval){
  //double er;
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15)){
    return -2*rval*(1-2*rval)*(1-6*ravl+6*rval*rval);
  }
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12)){
    return 4*rval*rval*(1-2*rval)*(1-3*rval);
  }
  else{
    return rval*(1-2*rval)*(1-6*rval);
  }
}


// First derivative functions for epsilon
double partial_der_epsilon(int geno, double ep, int a, int d){
  switch(geno){
  case 1:
    return binomial(a,d-a) *powl(ep,d-a) * powl(1-ep,a) * ((d-a)-d*ep);
  case 2:
    return 0;
  case 3:
    return binomial(a,d-a) * powl(ep,a) * powl(1-ep,d-a) * (a-d*ep);
  }
  return -1;
}


double der_epsilon(int OPGP, double epsilon, int a, int b, int elem){
  if((a == 0) & (b == 0))
    return 0;
  int d = a + b;
  switch(OPGP){
  case 1:
    if(elem == 1)
      return partial_der_epsilon(3, epsilon, a, d);
    else if ((elem == 2)|(elem == 3))  
      return 0;
    else if (elem == 4)
      return partial_der_epsilon(1, epsilon, a, d);
  case 2:
    if(elem == 3)
      return partial_der_epsilon(3, epsilon, a, d);
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 2)
      return partial_der_epsilon(1, epsilon, a, d);
  case 3:
    if(elem == 2) 
      return partial_der_epsilon(3, epsilon, a, d);
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 3)
      return partial_der_epsilon(1, epsilon, a, d);
  case 4:
    if(elem == 4) 
      return partial_der_epsilon(3, epsilon, a, d);
    else if ((elem == 2)|(elem == 3))
      return 0;
    else if (elem == 1)
      return partial_der_epsilon(1, epsilon, a, d);
  case 5:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return partial_der_epsilon(1, epsilon, a, d);
  case 6:
    if ((elem == 1)|(elem == 2))
      return partial_der_epsilon(1, epsilon, a, d);
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 7:
    if ((elem == 1)|(elem == 2))
      return partial_der_epsilon(3, epsilon, a, d);
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 8:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return partial_der_epsilon(3, epsilon, a, d);
  case 9:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return partial_der_epsilon(1, epsilon, a, d);
  case 10:
    if ((elem == 1)|(elem == 3))
      return partial_der_epsilon(1, epsilon, a, d);
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 11:
    if ((elem == 1)|(elem == 3))
      return partial_der_epsilon(3, epsilon, a, d);
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 12:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return partial_der_epsilon(3, epsilon, a, d);
  case 13:
    return partial_der_epsilon(1, epsilon, a, d);
  case 14:
    return 0;
  case 15:
    return 0;
  case 16:
    return partial_der_epsilon(3, epsilon, a, d);
  } // end of Switch
  return -999;
}

// Second derivative functions for epsilon
double partial_der_epsilon2(int geno, double ep, int a, int d){
  switch(geno){
  case 1:
    return binomial(a,d-a) * powl(ep,d-a) * powl(1-ep,a) * (-d*ep*(1-ep) + powl((d-a)-d*ep,2));
  case 2:
    return 0;
  case 3:
    return binomial(a,d-a) * powl(ep,a) * powl(1-ep,d-a) * (-d*ep*(1-ep) + powl(a-d*ep,2));
  }
  return -1;
}


double der_epsilon2(int OPGP, double epsilon, int a, int b, int elem){
  if((a == 0) & (b == 0))
    return 0;
  int d = a + b;
  switch(OPGP){
  case 1:
    if(elem == 1)
      return partial_der_epsilon2(3, epsilon, a, d);
    else if ((elem == 2)|(elem == 3))  
      return 0;
    else if (elem == 4)
      return partial_der_epsilon2(1, epsilon, a, d);
  case 2:
    if(elem == 3)
      return partial_der_epsilon2(3, epsilon, a, d);
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 2)
      return partial_der_epsilon2(1, epsilon, a, d);
  case 3:
    if(elem == 2) 
      return partial_der_epsilon2(3, epsilon, a, d);
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 3)
      return partial_der_epsilon2(1, epsilon, a, d);
  case 4:
    if(elem == 4) 
      return partial_der_epsilon2(3, epsilon, a, d);
    else if ((elem == 2)|(elem == 3))
      return 0;
    else if (elem == 1)
      return partial_der_epsilon2(1, epsilon, a, d);
  case 5:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return partial_der_epsilon2(1, epsilon, a, d);
  case 6:
    if ((elem == 1)|(elem == 2))
      return partial_der_epsilon2(1, epsilon, a, d);
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 7:
    if ((elem == 1)|(elem == 2))
      return partial_der_epsilon2(3, epsilon, a, d);
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 8:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return partial_der_epsilon2(3, epsilon, a, d);
  case 9:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return partial_der_epsilon2(1, epsilon, a, d);
  case 10:
    if ((elem == 1)|(elem == 3))
      return partial_der_epsilon2(1, epsilon, a, d);
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 11:
    if ((elem == 1)|(elem == 3))
      return partial_der_epsilon2(3, epsilon, a, d);
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 12:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return partial_der_epsilon2(3, epsilon, a, d);
  case 13:
    return partial_der_epsilon2(1, epsilon, a, d);
  case 14:
    return 0;
  case 15:
    return 0;
  case 16:
    return partial_der_epsilon2(3, epsilon, a, d);
  } // end of Switch
  return -999;
}



// rf's are equal and error parameter
SEXP der_fs_scaled_err_c(SEXP r, SEXP epsilon, SEXP depth_Ref, SEXP depth_Alt, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, snp2, snp_der, nInd_c, nSnps_c, *pOPGP, *pdepth_Ref, *pdepth_Alt;
  double *pr, *pKaa, *pKab, *pKbb, epsilon_c, delta;
  double alphaTilde[4], alphaDot[4], sum, sum_der, sum_der2, sum_der3, w_new, w_prev;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = INTEGER(OPGP);
  pdepth_Ref = INTEGER(depth_Ref);
  pdepth_Alt = INTEGER(depth_Alt);
  pKaa = REAL(Kaa);
  pKab = REAL(Kab);
  pKbb = REAL(Kbb);
  pr = REAL(r);
  epsilon_c = REAL(epsilon)[0];
  // Define the output variable
  SEXP ll, score, hessian, out;
  PROTECT(ll = allocVector(REALSXP, 1));
  PROTECT(score = allocVector(REALSXP, nSnps_c));
  PROTECT(hessian = allocVector(REALSXP, nSnps_c*(nSnps_c+1)/2));
  PROTECT(out = allocVector(VECSXP, 3));
  double *pll, *pscore, *phessian;
  pll = REAL(ll);
  pscore = REAL(score);
  phessian = REAL(hessian);
  //SEXP pout = PROTECT(allocVector(VECSXP, 3));
  double llval = 0, phi[4][nSnps_c], phi_prev[4][nSnps_c], score_c[nSnps_c];
  for(snp = 0; snp < nSnps_c; snp++){
    score_c[snp] = 0;
  }
  double psi[4][nSnps_c][nSnps_c], psi_prev[4][nSnps_c][nSnps_c], hessian_c[nSnps_c][nSnps_c];
  for(snp = 0; snp < nSnps_c; snp++){
    for(snp2 = 0; snp2 < nSnps_c; snp2++){
      hessian_c[snp][snp2] = 0;
    }
  }
  
  // Now compute the likelihood and score function
  for(ind = 0; ind < nInd_c; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
      alphaTilde[s1] = alphaDot[s1];
      // Compute the derivative for epsilon
      phi_prev[s1][nSnps_c-1] = 0.25 * der_epsilon(pOPGP[0], epsilon_c, pdepth_Ref[ind], pdepth_Alt[ind], s1+1);
      psi_prev[s1][nSnps_c-1][nSnps_c-1] = 0.25 * der_epsilon2(pOPGP[0], epsilon_c, pdepth_Ref[ind], pdepth_Alt[ind], s1+1);
    }
    
    // add contribution to likelihood
    w_prev = sum;
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        delta = Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1);
        alphaDot[s2] = sum * delta/w_prev;
        // add contribution to new weight
        w_new = w_new + alphaDot[s2];
        // Compute the derivatives
        // rf's
        sum_der = 0;
        sum_der2 = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum_der = sum_der + der_rf(s1, s2, pr[snp-1]) * alphaTilde[s1];
          sum_der2 = sum_der2 + der_rf2(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        // first derivate k=j-1
        phi[s2][snp-1] = sum_der * delta * 1/w_prev;
        // second derivate k=l=j-1
        psi[s2][snp-1][snp-1] = sum_der2 * delta * 1/w_prev;
        // second derivate k=1,...,j-2,l=j-1
        for(snp2 = 0; snp2 < snp-1; snp2++){
          sum_der = 0;
          for(s1 = 0; s1 < 4; s1++){
            sum_der += phi_prev[s2][snp2] * der_rf(s1, s2, pr[snp-1]);
          }
          psi[s2][snp-1][snp2] = sum_der * delta * 1/w_prev;
        }
        for(snp_der = 0; snp_der < snp-1; snp_der++){
          sum_der = 0;
          // rf first derivative
          for(s1 = 0; s1 < 4; s1++){
            sum_der += phi_prev[s1][snp_der] * Tmat(s1, s2, pr[snp-1]);
          }
          phi[s2][snp_der] = sum_der * delta * 1/w_prev;
          // rf and epsilon second derivative
          sum_der = 0;
          for(s1 = 0; s1 < 4; s1++){
            sum_der += (psi_prev[s1][nSnps_c-1][snp_der] * delta +   
              phi_prev[s1][snp_der] * der_epsilon(pOPGP[snp], epsilon_c, pdepth_Ref[ind + nInd_c*snp], pdepth_Alt[ind + nInd_c*snp], s2+1)) *
              Tmat(s1, s2, pr[snp-1];  
          }
          psi[s2][nSnps_c-1][snp_der] = sum_der/w_prev;
          // rf second derivative, k,l < j-2, k\neq l
          for(snp2 = 0; snp2 < snp_der+1; snp2++){
            sum_der = 0;
            for(s1 = 0; s1 < 4; s1++){
              sum_der += psi_prev[s1][snp_der][snp2] * Tmat(s1, s2, pr[snp-1]);
            }
            psi[s2][snp_der][snp2] = sum_der * delta * 1/w_prev;   
          }
        }
        // sequencing error parameter
        sum_der = 0;
        sun_der2 = 0;
        for(s1 = 0; s1 < 4; s1++){  
          sum_der += ((phi_prev[s1][nSnps_c-1] * delta + alphaTilde[s1] * 
            der_epsilon(pOPGP[snp], epsilon_c, pdepth_Ref[ind + nInd_c*snp], pdepth_Alt[ind + nInd_c*snp], s2+1)) * 
            Tmat(s1, s2, pr[snp-1]));
          sum_der2 += (psi_prev[s1][nSnps_c-1][nSnps_c-1] * delta + 
            2 * phi_prev[s1][nSnps_c-1] * der_epsilon(pOPGP[snp], epsilon_c, pdepth_Ref[ind + nInd_c*snp], pdepth_Alt[ind + nInd_c*snp], s2+1) +
            alphaTilde[s1] * der_epsilon2(pOPGP[snp], epsilon_c, pdepth_Ref[ind + nInd_c*snp], pdepth_Alt[ind + nInd_c*snp], s2+1)) *
            Tmat(s1, s2, pr[snp-1];
        }
        phi[s2][nSnps_c-1] = sum_der/w_prev;
        psi[s2][nSnps_c-1][nSnps_c-1] = sum_der2/w_prev;
      }
      // Add contribution to the likelihood
      llval = llval + log(w_prev);
      
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2];
        // update the derivative vectors
        for(snp_der = 0; snp_der < nSnps_c; snp_der++){
          phi_prev[s2][snp_der] = phi[s2][snp_der];
          for(snp2 = 0; snp2 < snp_der; snp2++){
            psi_prev[s2][snp_der][snp2] = psi[s2][snp_der][snp2];
          }
        }
      }
      w_prev = w_new;
    }
    llval = llval + log(w_prev);
    //Rprintf("Likelihood value: %f\n", llval);
    // add contributions to the score vector
    for(snp_der = 0; snp_der < nSnps_c; snp_der++){
      sum_der = 0;
      for(s2 = 0; s2 < 4; s2++)
        sum_der += phi[s2][snp_der]/w_prev;
      sum_der += score_c[snp_der];
      score_c[snp_der] = sum_der;
    }
    // add contributions to the hessian matrix
    for(snp_der = 0; snp_der < nSnps_c; snp_der++){
      for(snp2 = 0; snp2 < snp_der; snp2++){
        sum_der = 0;
        sum_der2 = 0;
        sum_der3 = 0;
        for(s2 = 0; s2 < 4; s2++){
          sum_der += -psi[s2][snp_der][snp2];
          sum_der2 += phi[s2][snp_der];
          sum_der3 += phi[s2][snp2];
        }
        sum_der = -sum_der/w_prev + sum_der2*sum_der3/(w_prev*w_prev);
        sum_der += hessian_c[snp_der];
        hessian_c[snp_der] = sum_der;
      }
    }  
  }
  // Compute the score for each parameter
  for(snp_der=0; snp_der < nSnps_c; snp_der++){
    pscore[snp_der] = score_c[snp_der];
  }
  // Compute the hessian
  for(snp_der=0; snp_der < nSnps_c; snp_der++){
    for(snp2=0; snp2 < snp_der; snp2++){
      phessian[snp_der + snp2*nSnps_c] = hessian_c[snp_der][snp2];
    }
  }
  
  // Clean up and return likelihood value
  //Rprintf("Likelihood value: %f\n", llval);
  SET_VECTOR_ELT(out, 0, ll);
  SET_VECTOR_ELT(out, 1, score);
  SET_VECTOR_ELT(out, 2, hessian);
  UNPROTECT(4);
  return out;
}



#include "weno.h"

double weno3_upwind(double vec0, double vec1, double vec2)
{


  // Initialize the working arrays
  double alpha0=0, alpha1=0;
  double beta0=0, beta1=0;
  double stencils0=0, stencils1=0;
  double w0=0, w1=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno3 upwind reconstruction of vec[]
  beta0 += vec1 * vec1;
  stencils0 += 1.5 * vec1;
  beta0 += -2 * vec0 * vec1;
  beta0 += vec0 * vec0;
  stencils0 += -0.5 * vec0;
  alpha0 = (1./3.) / (epsilon + beta0 * beta0);
  alphaSum += alpha0;
  beta1 += vec2 * vec2;
  stencils1 += 0.5 * vec2;
  beta1 += -2 * vec1 * vec2;
  beta1 += vec1 * vec1;
  stencils1 += 0.5 * vec1;
  alpha1 = (2./3.) / (epsilon + beta1 * beta1);
  alphaSum += alpha1;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;

  return (w0*stencils0 + w1*stencils1);
}

// This taken from the Shu review
double weno5_upwind(double vec0, double vec1, double vec2, double vec3, double vec4)
{


  // Initialize the working arrays
  double alpha0=0, alpha1=0, alpha2=0;
  double beta0=0, beta1=0, beta2=0;
  double stencils0=0, stencils1=0, stencils2=0;
  double w0=0, w1=0, w2=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno5 upwind reconstruction of vec[]
  beta0 = 13.0*pow(vec2-2.0*vec3+vec4, 2)/12.0
  beta0 += pow(vec2-4.0*vec3+vec4, 2)/4.0
  beta1 = 13.0*pow(vec1-2.0*vec2+vec3, 2)/12.0
  beta1 += pow(v1-v3, 2)/4.0
  beta2 = 13.0*pow(vec0-2.0*vec1+vec2, 2)/12.0
  beta2 += pow(vec0-4.0*vec1+3.0*vec2, 2)/4.0
  alpha0 = 0.3 / (epsilon + beta0 * beta0);
  alpha1 = 0.6 / (epsilon + beta1 * beta1);
  alpha2 = 0.1 / (epsilon + beta2 * beta2);
  alphaSum = alpha0 + alpha1 + alpha2;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;
  w2 = alpha2 / alphaSum;
  // q^r = sum_{j=0}^{k-1} q_{i-r+j}
  stencils0 += 1.0/3.0*vec2;
  stencils0 += 5.0/6.0*vec3;
  stencils0 += -1.0/6.0*vec4;
  stencils1 += -1.0/6.0*vec1;
  stencils1 += 5.0/6.0*vec2;
  stencils1 += 1.0/3.0*vec3;
  stencils2 += 1.0/3.0*vec0;
  stencils2 += -7.0/6.0*vec1;
  stencils2 += 11.0/6.0*vec2;

  return (w0*stencils0 + w1*stencils1 + w2*stencils2);
}

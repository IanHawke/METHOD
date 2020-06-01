//! Special relativistic magnetohydrodynamics model
/*!
    This script contains the function definitions for the srmhd model. The form
  of the quations has been taken from Giacomazzo et al's Spritz code paper,
  with modifications.
    For detailed documentation about the methods contained herein, see srmhd.h
  and model.h.
*/

#include "srmhd_vector_potential.h"
#include "weno.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>

SRMHD_Vector_Potential::SRMHD_Vector_Potential() : Model()
{
  this->Ncons = 9;
  this->Nprims = 9;
  this->Naux = 20;
}

SRMHD_Vector_Potential::SRMHD_Vector_Potential(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 9;
  this->Nprims = (this->data)->Nprims = 9;
  this->Naux = (this->data)->Naux = 20;

  // Solutions for C2P all cells
  solution = (double *) malloc(sizeof(double)*2*data->Nx*data->Ny*data->Nz);

  smartGuesses = 0;

  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("Sx");
  this->data->consLabels.push_back("Sy");  this->data->consLabels.push_back("Sz");
  this->data->consLabels.push_back("tau"); this->data->consLabels.push_back("Ax");
  this->data->consLabels.push_back("Ay");  this->data->consLabels.push_back("Az");
  this->data->consLabels.push_back("Phi");

  this->data->primsLabels.push_back("rho"); this->data->primsLabels.push_back("vx");
  this->data->primsLabels.push_back("vy");  this->data->primsLabels.push_back("vz");
  this->data->primsLabels.push_back("p");   this->data->primsLabels.push_back("Ax");
  this->data->primsLabels.push_back("Ay");  this->data->primsLabels.push_back("Az");
  this->data->primsLabels.push_back("Phi");

  this->data->auxLabels.push_back("h");   this->data->auxLabels.push_back("W");
  this->data->auxLabels.push_back("e");   this->data->auxLabels.push_back("c");
  this->data->auxLabels.push_back("b0");  this->data->auxLabels.push_back("bx");
  this->data->auxLabels.push_back("by");  this->data->auxLabels.push_back("bz");
  this->data->auxLabels.push_back("bsq"); this->data->auxLabels.push_back("vsq");
  this->data->auxLabels.push_back("BS");  this->data->auxLabels.push_back("Bsq");
  this->data->auxLabels.push_back("Ssq"); this->data->auxLabels.push_back("Bx");
  this->data->auxLabels.push_back("By");  this->data->auxLabels.push_back("Bz");
  this->data->auxLabels.push_back("Ex");  this->data->auxLabels.push_back("Ey");
  this->data->auxLabels.push_back("Ez");  this->data->auxLabels.push_back("Phi_staggered");
}

SRMHD_Vector_Potential::~SRMHD_Vector_Potential()
{
  free(solution);
}


//! Generates the net numerical flux given the current state
/*!
    We are using the flux vector splitting method described in Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws`. For the form of the fluxes see Relativistic Magneto..., Anton '10
  with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
  interfaces, John Muddle.
    Note: We are assuming that all primitive and auxiliary variables are up-to-date
  at the time of this function execution.
*/
void SRMHD_Vector_Potential::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  // Generate flux vector
  // Fx: flux in x-direction
  if (dir == 0) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(1, i, j, k)];

          // Sx
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                 (prims[ID(4, i, j, k)] + aux[ID(8, i, j, k)] * 0.5) -
                                 aux[ID(5, i, j, k)] * aux[ID(13, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // Sy
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(1, i, j, k)] -
                                 aux[ID(6, i, j, k)] * aux[ID(13, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // Sz
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(1, i, j, k)] -
                                 aux[ID(7, i, j, k)] * aux[ID(13, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)] +
                                 aux[ID(8, i, j, k)] * 0.5) * prims[ID(1, i, j, k)] -
                                 aux[ID(4, i, j, k)] * aux[ID(13, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // Ax
          f[ID(5, i, j, k)] = cons[ID(8, i, j, k)];

          // Ay
          f[ID(6, i, j, k)] = 0.0;

          // Az
          f[ID(7, i, j, k)] = 0.0;

          // Phi
          f[ID(8, i, j, k)] = cons[ID(5, i, j, k)];

        }
      } // End k loop
    } // End j loop
  } // End i loop

  // Fy: flux in y-direction
  else if (dir==1) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(2, i, j, k)];

          // Sx
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(2, i, j, k)] -
                              aux[ID(5, i, j, k)] * aux[ID(14, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sy
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                              prims[ID(4, i, j, k)] + aux[ID(8, i, j, k)] * 0.5 -
                              aux[ID(6, i, j, k)] * aux[ID(14, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sz
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(2, i, j, k)] -
                              aux[ID(7, i, j, k)] * aux[ID(14, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)] +
                              aux[ID(8, i, j, k)] * 0.5) * prims[ID(2, i, j, k)] -
                              aux[ID(4, i, j, k)] * aux[ID(14, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Ax
          f[ID(5, i, j, k)] = 0.0;

          // Ay
          f[ID(6, i, j, k)] = cons[ID(8, i, j, k)];

          // Az
          f[ID(7, i, j, k)] = 0.0;

          // Phi
          f[ID(8, i, j, k)] = cons[ID(6, i, j, k)];

        }
      } // End k loop
    } // End j loop
  } // End i loop

  // Fz: flux in z-direction
  else {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(3, i, j, k)];

          // Sx
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(3, i, j, k)] -
                              aux[ID(5, i, j, k)] * aux[ID(15, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sy
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(3, i, j, k)] -
                              aux[ID(6, i, j, k)] * aux[ID(15, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sz
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(3, i, j, k)] +
                              prims[ID(4, i, j, k)] + aux[ID(8, i, j, k)] * 0.5 -
                              aux[ID(7, i, j, k)] * aux[ID(15, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)] +
                              aux[ID(8, i, j, k)] * 0.5) * prims[ID(3, i, j, k)] -
                              aux[ID(4, i, j, k)] * aux[ID(15, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Ax
          f[ID(5, i, j, k)] = 0.0;
          // Ay
          f[ID(6, i, j, k)] = 0.0;
          // Az
          f[ID(7, i, j, k)] = cons[ID(8, i, j, k)];
          // Phi
          f[ID(8, i, j, k)] = cons[ID(7, i, j, k)];

        }
      } // End k loop
    } // End j loop
  } // End i loop
}



//! Single cell source required for divergence cleaning
/*!
    See Anton 2010, `Relativistic Magnetohydrodynamcis: Renormalized Eignevectors
  and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD_Vector_Potential::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  double xi = 1.5 / this->data->dt;
  for (int var(0); var < Ncons; var++) {
    if (var == 8) {
      // phi
      source[var] = -xi * cons[8];
    }
    else {
      source[var] = 0;
    }
  }
}

//! Source required for generalized Lorenz gauge
/*!
    See Giacomazzo's Spritz paper
*/
void SRMHD_Vector_Potential::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  // Order of weno scheme
  int order(3);

  // Wave speed
  double alpha(1);
  // Number of variables, which is 3 for electric fields
  int nvars(3);

  // First do the generalized Lorenz gauge damping
  double xi = 1.5 / d->dt;
  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int k(0); k < this->data->Nz; k++) {
        for (int var(0); var < Ncons; var++) {
          if (var == 8) {
            // phi
            source[this->data->id(var, i, j, k)] = -xi * cons[this->data->id(8, i, j, k)];
          }
          else {
            source[this->data->id(var, i, j, k)] = 0;
          }
        }
      }
    }
  }

  // Up and downwind fluxes, which give electric fields
  double *fplus, *fminus, *frecon;
  fplus = (double *) malloc(sizeof(double) * nvars * d->Nx * d->Ny * d->Nz);
  fminus = (double *) malloc(sizeof(double) * nvars * d->Nx * d->Ny * d->Nz);
  frecon = (double *) malloc(sizeof(double) * nvars * d->Nx * d->Ny * d->Nz);

  // Lax-Friedrichs approximation of flux
  double Ex(0), Ey(0), Ez(0);
  // We are going to do a double reconstruction
  // E is the flux term for B, but components and signs need care
  // dt Bx + dx ( 0 ) + dy ( Ez) + dz (-Ey)
  // dt By + dx (-Ez) + dy ( 0 ) + dz ( Ex)
  // dt Bz + dx ( Ey) + dy (-Ex) + dz ( 0 )
  // So Ex(+) is Ex + \alpha (By - Bz)
  // So Ey(+) is Ey + \alpha (Bz - Bx)
  // So Ez(+) is Ez + \alpha (Bx - By)
  // Get the fluxes first
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double Bx = aux[ID(13, i, j, k)];
        double By = aux[ID(14, i, j, k)];
        double Bz = aux[ID(15, i, j, k)];
        double vx = prims[ID(1, i, j, k)];
        double vy = prims[ID(2, i, j, k)];
        double vz = prims[ID(3, i, j, k)];
        Ex = By * vz - Bz * vy;
        Ey =-Bx * vz + Bz * vx;
        Ez = Bx * vy - By * vx;
        fplus [ID(0, i, j, k)] = 0.5 * (Ex + alpha * (By - Bz));
        fminus[ID(0, i, j, k)] = 0.5 * (Ex - alpha * (By - Bz));
        fplus [ID(1, i, j, k)] = 0.5 * (Ey + alpha * (Bz - Bx));
        fminus[ID(1, i, j, k)] = 0.5 * (Ey - alpha * (Bz - Bx));
        fplus [ID(2, i, j, k)] = 0.5 * (Ez + alpha * (Bx - By));
        fminus[ID(2, i, j, k)] = 0.5 * (Ez - alpha * (Bx - By));
      }
    }
  }
  for (int var(0); var < nvars; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          if (i >= order && i < d->Nx-order) {
            frecon[ID(var, i, j, k)] = weno5_upwind(
                                        fplus[ID(var, i, j-order  , k)],
                                        fplus[ID(var, i, j-order+1, k)],
                                        fplus[ID(var, i, j-order+2, k)],
                                        fplus[ID(var, i, j-order+3, k)],
                                        fplus[ID(var, i, j-order+4, k)]) +
                                       weno5_upwind(
                                        fminus[ID(var, i, j+order-1, k)],
                                        fminus[ID(var, i, j+order-2, k)],
                                        fminus[ID(var, i, j+order-3, k)],
                                        fminus[ID(var, i, j+order-4, k)],
                                        fminus[ID(var, i, j+order-5, k)]);
          }
          else {
            frecon[ID(var, i, j, k)] = 0.0;
          }
        }
      }
    }
  }
  // frecon now contains E at faces. We want it at corners.
  for (int var(0); var < nvars; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          if (i >= order && i < d->Nx-order) {
            double E_corner         = weno5_upwind(
                                        frecon[ID(var, i-order, j  , k)],
                                        frecon[ID(var, i-order+1, j, k)],
                                        frecon[ID(var, i-order+2, j, k)],
                                        frecon[ID(var, i-order+3, j, k)],
                                        frecon[ID(var, i-order+4, j, k)]) +
                                       weno5_upwind(
                                        frecon[ID(var, i+order-1, j, k)],
                                        frecon[ID(var, i+order-2, j, k)],
                                        frecon[ID(var, i+order-3, j, k)],
                                        frecon[ID(var, i+order-4, j, k)],
                                        frecon[ID(var, i+order-5, j, k)]);
            // This flux gives us the corner electric flux,
            // which is the source term for A
            source[ID(5+var, i, j, k)] += E_corner;
          }
        }
      }
    }
  }


}

//! Residual function to minimize in the format required by cminpack
/*!
    I know, its a horrible layout, alas we're stuck with it.

    Parameters
    ----------
    p : pointer to struct
      Struct contains additional arguments that are required (if any)
    n : int
      Size of system
    x : pointer to double
      The array containing the initial guess
    fvec : pointer to double
      The array containing the solution
    iflag : int
      Error flag
*/
int SRMHDresidual(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Retrieve additional arguments
  Args * args = (Args*) p;

  // Values must make sense
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  double Bsq(args->Bsq);
  double Ssq(args->Ssq);
  double BS(args->BS);
  double W(1 / sqrt(1 - x[0]));
  double rho(args->D / W);
  double h(x[1] / (rho * W * W));
  double pr((h - 1) * rho * (args->g - 1) / args->g);
  if (pr < 0 || rho < 0 || h < 0 || W < 1) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  // Values should be OK
  fvec[0] = (x[1] + Bsq) * (x[1] + Bsq) * x[0] - (2 * x[1] + Bsq) * BS * BS / (x[1] * x[1]) - Ssq;
  fvec[1] = x[1] + Bsq - pr - Bsq / (2 * W * W) - BS * BS / (2 * x[1] * x[1]) - args->D - args->tau;

  return 0;
}

void SRMHD_Vector_Potential::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(2);                     // Size of system
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.49011612e-8;   // Tolerance of rootfinder
  const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array


  // Update possible values
  // Need to have done A*->B* first
  // Ax, Ay, Az
  prims[5] = cons[5];
  prims[6] = cons[6];
  prims[7] = cons[7];

  // BS
  aux[10] = aux[13] * cons[1] + aux[14] * cons[2] + aux[15] * cons[3];
  // Bsq
  aux[11] = aux[13] * aux[13] + aux[14] * aux[14] + aux[15] * aux[15];
  // Ssq
  aux[12] = cons[1] * cons[1] + cons[2] * cons[2] + cons[3] * cons[3];


  // Set additional args for rootfind
  args.D = cons[0];
  args.g = d->gamma;
  args.BS = aux[10];
  args.Bsq = aux[11];
  args.Ssq = aux[12];
  args.tau = cons[4];

  sol[0] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
  sol[1] = prims[0] * aux[0] /
           (1 - sol[0]);

  // Solve residual = 0
  info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                    tol, wa, lwa);
  // If root find fails, add failed cell to the list
  if (info!=1) {
    //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
    throw std::runtime_error("C2P could not converge.\n");
  }
  // W
  aux[1] = 1 / sqrt(1 - sol[0]);
  // rho
  prims[0] = cons[0] / aux[1];
  // h
  aux[0] = sol[1] / (prims[0] * aux[1] * aux[1]);
  // p
  prims[4] = (aux[0] - 1) * prims[0] *
                             (d->gamma - 1) / d->gamma;
  // e
  aux[2] = prims[4] / (prims[0] * (d->gamma - 1));
  // vx, vy, vz
  prims[1] = (cons[5] * aux[10] + cons[1] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  prims[2] = (cons[6] * aux[10] + cons[2] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  prims[3] = (cons[7] * aux[10] + cons[3] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  // vsq
  aux[9] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
  // c
  aux[3] = sqrt(aux[2] * d->gamma * (d->gamma - 1) /   aux[0]);
  // b0
  aux[4] = aux[1] * (aux[13] * prims[1] + aux[14] * prims[2] + aux[15] * prims[3]);
  // bx, by, bz
  aux[5] = aux[13] / aux[1] + aux[4] * prims[1];
  aux[6] = aux[14] / aux[1] + aux[4] * prims[2];
  aux[7] = aux[15] / aux[1] + aux[4] * prims[3];
  // bsq
  aux[8] = (aux[13] * aux[13] + aux[14] * aux[14] + aux[15] * aux[15] +
                           aux[4] * aux[4]) / (aux[1] * aux[1]);


}


//! Solve for the primitive and auxiliary variables
/*!
    Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
  an N=2 rootfind using cminpack library.

  Initial inputs will be the current values of the conserved vector and the
  old values for the prims and aux vectors.
  Output is the current values of cons, prims and aux.
*/
void SRMHD_Vector_Potential::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(2);                     // Size of system
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.49011612e-7;   // Tolerance of rootfinder
  const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array
  std::vector<Failed> fails;          // Vector of failed structs. Stores location of failed cons2prims cells.

  // Do the A -> B first. Note the staggered grid effects at the top boundaries.
  for (int i(0); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz-1; k++) {
        // 1d (x derivs)
        aux[ID(13, i, j, k)] = 0;
        aux[ID(14, i, j, k)] = -(cons[ID(7, i, j, k)] - cons[ID(7, i+1, j, k)]) / d->dx;
        aux[ID(15, i, j, k)] =  (cons[ID(6, i, j, k)] - cons[ID(6, i+1, j, k)]) / d->dx;
        // 2d (y derivs)
        if (d->Ny > 1) {
          aux[ID(13, i, j, k)] += (cons[ID(7, i, j, k)] - cons[ID(7, i+1, j, k)]) / d->dy;
          aux[ID(15, i, j, k)] -= (cons[ID(5, i, j, k)] - cons[ID(5, i+1, j, k)]) / d->dy;
          // 3d (z derivs)
          if (d->Nz > 1) {
            aux[ID(13, i, j, k)] -= (cons[ID(6, i, j, k)] - cons[ID(6, i+1, j, k)]) / d->dz;
            aux[ID(14, i, j, k)] += (cons[ID(5, i, j, k)] - cons[ID(5, i+1, j, k)]) / d->dz;
          }
        }
      }
    }
  }

  // Loop through domain solving and setting the prim and aux vars
  for (int i(0); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz-1; k++) {
        // Update possible values
        // Ax, Ay, Az
        prims[ID(5, i, j, k)] = cons[ID(5, i, j, k)];
        prims[ID(6, i, j, k)] = cons[ID(6, i, j, k)];
        prims[ID(7, i, j, k)] = cons[ID(7, i, j, k)];

        // BS
        aux[ID(10, i, j, k)] = aux[ID(13, i, j, k)] * cons[ID(1, i, j, k)] +
                                  aux[ID(14, i, j, k)] * cons[ID(2, i, j, k)] +
                                  aux[ID(15, i, j, k)] * cons[ID(3, i, j, k)];
        // Bsq
        aux[ID(11, i, j, k)] = aux[ID(13, i ,j, k)] * aux[ID(13, i, j, k)] +
                                  aux[ID(14, i, j, k)] * aux[ID(14, i, j, k)] +
                                  aux[ID(15, i, j, k)] * aux[ID(15, i, j, k)];
        // Ssq
        aux[ID(12, i, j, k)] = cons[ID(1, i ,j, k)] * cons[ID(1, i, j, k)] +
                                  cons[ID(2, i, j, k)] * cons[ID(2, i, j, k)] +
                                  cons[ID(3, i, j, k)] * cons[ID(3, i, j, k)];


        // Set additional args for rootfind
        args.D = cons[ID(0, i, j, k)];
        args.g = d->gamma;
        args.BS = aux[ID(10, i, j, k)];
        args.Bsq = aux[ID(11, i, j, k)];
        args.Ssq = aux[ID(12, i, j, k)];
        args.tau = cons[ID(4, i, j, k)];

        sol[0] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                 prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                 prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        sol[1] = prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] /
                 (1 - sol[0]);

        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                          tol, wa, lwa);
        // If root find fails, add failed cell to the list
        if (info!=1) {
          Failed fail = {i, j, k};
          fails.push_back(fail);
        }
        else {
          // Now have the correct values for vsq and rho*h*Wsq
          solution[ID(0, i, j, k)] = sol[0];
          solution[ID(1, i, j, k)] = sol[1];
        }

      } // End k-loop
    } // End j-loop
  } // End i-loop





  // ################################## Smart guessing ########################### //
  // Are there any failures?
  if (fails.size() > 0) {
    int x, y, z;
    // Loop through any failed cells and try again, using the mean of successfull
    // surrounding cells solutions as an initial estimate
    for (Failed fail : fails) {
      x = fail.x;
      y = fail.y;
      z = fail.z;
      // Vector to contain successful neighbours
      std::vector<Failed> neighbours;
      if (x > 0) neighbours.push_back(Failed {x-1, y, z});
      if (y > 0) neighbours.push_back(Failed {x, y-1, z});
      if (z > 0) neighbours.push_back(Failed {x, y, z-1});
      if (x < d->Nx - 1) neighbours.push_back(Failed {x+1, y, z});
      if (y < d->Ny - 1) neighbours.push_back(Failed {x, y+1, z});
      if (z < d->Nz - 1) neighbours.push_back(Failed {x, y, z+1});

      sol[0] = 0;
      sol[1] = 0;
      for (Failed neighbour : neighbours) {
        sol[0] += solution[ID(0, neighbour.x, neighbour.y, neighbour.z)];
        sol[1] += solution[ID(1, neighbour.x, neighbour.y, neighbour.z)];
      }
      sol[0] /= neighbours.size();
      sol[1] /= neighbours.size();
      // Solve residual = 0
      info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                        tol, wa, lwa);
      if (info != 1) {
        printf("Smart guessing did not work, exiting\n");
        printf("(%d, %d, %d) failed\n", fail.x, fail.y, fail.z);
        // std::exit(1);
      }
      // else {
      //   smartGuesses++;
        // printf("Smart guessing worked!\n");
        solution[ID(0, x, y, z)] = sol[0];
        solution[ID(1, x, y, z)] = sol[1];
      // }
    }
  }


  for (int i(0); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz-1; k++) {
        // W
        aux[ID(1, i, j, k)] = 1 / sqrt(1 - solution[ID(0, i, j, k)]);
        // rho
        prims[ID(0, i, j, k)] = cons[ID(0, i, j, k)] / aux[ID(1, i, j, k)];
        // h
        aux[ID(0, i, j, k)] = solution[ID(1, i, j, k)] / (prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                 aux[ID(1, i, j, k)]);
        // p
        prims[ID(4, i, j, k)] = (aux[ID(0, i, j, k)] - 1) * prims[ID(0, i, j, k)] *
                                   (d->gamma - 1) / d->gamma;
        // e
        aux[ID(2, i, j, k)] = prims[ID(4, i, j, k)] / (prims[ID(0, i, j, k)] *
                                 (d->gamma - 1));
        // vx, vy, vz
        prims[ID(1, i, j, k)] = (aux[ID(13, i, j, k)] * aux[ID(10, i, j, k)] +
                                   cons[ID(1, i, j, k)] * solution[ID(1, i, j, k)]) / (solution[ID(1, i, j, k)] *
                                   (aux[ID(11, i, j, k)] + solution[ID(1, i, j, k)]));
        prims[ID(2, i, j, k)] = (aux[ID(14, i, j, k)] * aux[ID(10, i, j, k)] +
                                   cons[ID(2, i, j, k)] * solution[ID(1, i, j, k)]) / (solution[ID(1, i, j, k)] *
                                   (aux[ID(11, i, j, k)] + solution[ID(1, i, j, k)]));
        prims[ID(3, i, j, k)] = (aux[ID(15, i, j, k)] * aux[ID(10, i, j, k)] +
                                   cons[ID(3, i, j, k)] * solution[ID(1, i, j, k)]) / (solution[ID(1, i, j, k)] *
                                   (aux[ID(11, i, j, k)] + solution[ID(1, i, j, k)]));
        aux[ID(9, i, j, k)] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                 prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                                 prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        // c
        aux[ID(3, i, j, k)] = sqrt(aux[ID(2, i, j, k)] * d->gamma * (d->gamma -1) /
                              aux[ID(0, i, j, k)]);
        // b0
        aux[ID(4, i, j, k)] = aux[ID(1, i, j, k)] * (aux[ID(13, i, j, k)] * prims[ID(1, i, j, k)] +
                                                           aux[ID(14, i, j, k)] * prims[ID(2, i, j, k)] +
                                                           aux[ID(15, i, j, k)] * prims[ID(3, i, j, k)]);
        // bx, by, bz
        aux[ID(5, i, j, k)] = aux[ID(13, i, j, k)] / aux[ID(1, i, j, k)] +
                                 aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)];
        aux[ID(6, i, j, k)] = aux[ID(14, i, j, k)] / aux[ID(1, i, j, k)] +
                                 aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)];
        aux[ID(7, i, j, k)] = aux[ID(15, i, j, k)] / aux[ID(1, i, j, k)] +
                                 aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)];
        // bsq
        aux[ID(8, i, j, k)] = (aux[ID(13, i, j, k)] * aux[ID(13, i, j, k)] +
                                 aux[ID(14, i, j, k)] * aux[ID(14, i, j, k)] +
                                 aux[ID(15, i, j, k)] * aux[ID(15, i, j, k)] +
                                 aux[ID(4, i, j, k)] * aux[ID(4, i, j, k)]) /
                                 (aux[ID(1, i, j, k)] * aux[ID(1, i, j, k)]);
      } // End k-loop
    } // End j-loop
  } // End i-loop

}





//! Go from cons or prim containing A to aux containing B.
/*!
    See the Spritz paper
*/
void SRMHD_Vector_Potential::AToB(double *c_or_p, double *aux)
{

    // Syntax
    Data * d = this->data;

    // Note the staggered grid effects at the top boundaries.
    for (int i(0); i < d->Nx-1; i++) {
      for (int j(0); j < d->Ny-1; j++) {
        for (int k(0); k < d->Nz-1; k++) {
          // 1d (x derivs)
          aux[ID(13, i, j, k)] = 0;
          aux[ID(14, i, j, k)] = -(c_or_p[ID(7, i, j, k)] - c_or_p[ID(7, i+1, j, k)]) / d->dx;
          aux[ID(15, i, j, k)] =  (c_or_p[ID(6, i, j, k)] - c_or_p[ID(6, i+1, j, k)]) / d->dx;
          // 2d (y derivs)
          if (d->Ny > 1) {
            aux[ID(13, i, j, k)] += (c_or_p[ID(7, i, j, k)] - c_or_p[ID(7, i+1, j, k)]) / d->dy;
            aux[ID(15, i, j, k)] -= (c_or_p[ID(5, i, j, k)] - c_or_p[ID(5, i+1, j, k)]) / d->dy;
            // 3d (z derivs)
            if (d->Nz > 1) {
              aux[ID(13, i, j, k)] -= (c_or_p[ID(6, i, j, k)] - c_or_p[ID(6, i+1, j, k)]) / d->dz;
              aux[ID(14, i, j, k)] += (c_or_p[ID(5, i, j, k)] - c_or_p[ID(5, i+1, j, k)]) / d->dz;
            }
          }
        }
      }
    }

}


//! Generate to the conserved and auxiliary variables
/*!
    Relations have been taken from Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD_Vector_Potential::primsToAll(double *cons, double *prims, double *aux)
{


  // Syntax
  Data * d = this->data;

  // Do the A -> B first. Note the staggered grid effects at the top boundaries.
  this->AToB(prims, aux);

  for (int i(0); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz-1; k++) {
        // Ax, Ay, Az
        cons[ID(5, i, j, k)] = prims[ID(5, i, j, k)];
        cons[ID(6, i, j, k)] = prims[ID(6, i, j, k)];
        cons[ID(7, i, j, k)] = prims[ID(7, i, j, k)];

        // Bsq
        aux[ID(11, i, j, k)] = aux[ID(13, i, j, k)] * aux[ID(13, i, j, k)] +
                                     aux[ID(14, i, j, k)] * aux[ID(14, i, j, k)] +
                                     aux[ID(15, i, j, k)] * aux[ID(15, i, j, k)];

        // phi
        cons[ID(8, i, j, k)] = 0;

        // vsq
        aux[ID(9, i, j, k)] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                    prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                                    prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        // W
        aux[ID(1, i, j, k)] = 1.0 / sqrt(1 - aux[ID(9, i, j, k)]);

        // b0
        aux[ID(4, i, j, k)] = aux[ID(1, i, j, k)] * (
                                    prims[ID(1, i, j, k)] * aux[ID(13, i, j, k)] +
                                    prims[ID(2, i, j, k)] * aux[ID(14, i, j, k)] +
                                    prims[ID(3, i, j, k)] * aux[ID(15, i, j, k)]);

        // bx, by, bz
        aux[ID(5, i, j, k)] = aux[ID(13, i, j, k)] / aux[ID(1, i, j, k)] +
                                    aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)];
        aux[ID(6, i, j, k)] = aux[ID(14, i, j, k)] / aux[ID(1, i, j, k)] +
                                    aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)];
        aux[ID(7, i, j, k)] = aux[ID(15, i, j, k)] / aux[ID(1, i, j, k)] +
                                    aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)];

        // bsq
        aux[ID(8, i, j, k)] = (aux[ID(13, i, j, k)] * aux[ID(13, i, j, k)] +
                                    aux[ID(14, i, j, k)] * aux[ID(14, i, j, k)] +
                                    aux[ID(15, i, j, k)] * aux[ID(15, i, j, k)] +
                                    aux[ID(4, i, j, k)] * aux[ID(4, i, j, k)]) /
                                    (aux[ID(1, i, j, k)] * aux[ID(1, i, j, k)]);

        // h
        aux[ID(0, i, j, k)] = 1 + prims[ID(4, i, j, k)] / prims[ID(0, i, j, k)] *
                                 (d->gamma / (d->gamma - 1));

        // e
        aux[ID(2, i, j, k)] = prims[ID(4, i, j, k)] / (prims[ID(0, i, j, k)] * (d->gamma - 1));

        // c
        aux[ID(3, i, j, k)] = sqrt(aux[ID(2, i, j, k)] * d->gamma * (d->gamma - 1) / aux[ID(0, i, j, k)]);

        // D
        cons[ID(0, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)];

        // Sx, Sy, Sz
        cons[ID(1, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] * prims[ID(1, i, j, k)] -
                                     aux[ID(4, i, j, k)] * aux[ID(5, i, j, k)];
        cons[ID(2, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] * prims[ID(2, i, j, k)] -
                                     aux[ID(4, i, j, k)] * aux[ID(6, i, j, k)];
        cons[ID(3, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] * prims[ID(3, i, j, k)] -
                                     aux[ID(4, i, j, k)] * aux[ID(7, i, j, k)];

        // Ssq
        aux[ID(12, i, j, k)] = cons[ID(1, i, j, k)] * cons[ID(1, i, j, k)] +
                                     cons[ID(2, i, j, k)] * cons[ID(2, i, j, k)] +
                                     cons[ID(3, i, j, k)] * cons[ID(3, i, j, k)];

        // BS
        aux[ID(10, i, j, k)] = aux[ID(13, i, j, k)] * cons[ID(1, i, j, k)] +
                                     aux[ID(14, i, j, k)] * cons[ID(2, i, j, k)] +
                                     aux[ID(15, i, j, k)] * cons[ID(3, i, j, k)];

        // tau
        cons[ID(4, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] - (prims[ID(4, i, j, k)] +
                                     aux[ID(8, i, j, k)] / 2.0) - aux[ID(4, i, j, k)] *
                                     aux[ID(4, i, j, k)] - cons[ID(0, i, j, k)];

      } // End k-loop
    } // End j-loop
  } // End i-loop


}

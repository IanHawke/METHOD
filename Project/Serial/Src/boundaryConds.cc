#include "boundaryConds.h"
#include <stdio.h>

void Outflow::apply(double * cons, double * prims, double * aux)
{

  // Syntax
  Data * d(this->data);

  // Cons
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->Ng, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            prims[ID(var, i, j, k)] = prims[ID(var, d->Ng, j, k)];
            // Right
            prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->nx + d->Ng - 1, j, k)];
          }
        }
      }
    }
  }
  if (aux) {
    // Aux
    for (int var(0); var < d->Naux; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            aux[ID(var, i, j, k)] = aux[ID(var, d->Ng, j, k)];
            // Right
            aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->nx + d->Ng - 1, j, k)];
          }
        }
      }
    }
  }
  if (d->Ny > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i, d->Ng, k)];
            // Back
            cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->ny + d->Ng - 1, k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, i, j, k)] = prims[ID(var, i, d->Ng, k)];
              // Back
              prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, i, j, k)] = aux[ID(var, i, d->Ng, k)];
              // Back
              aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
  }

  if (d->Nz > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            // Bottom
            cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->Ng)];
            // Top
            cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->nz + d->Ng - 1)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->Ng)];
              // Top
              prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->Ng)];
              // Top
              aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
  }
}



void OutflowRotatedBW::apply(double * cons, double * prims, double * aux)
{

  // Syntax
  Data * d(this->data);

  // Cons
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->Ng, j+d->Ng, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            prims[ID(var, i, j, k)] = prims[ID(var, d->Ng, j, k)];
            // Right
            prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->nx + d->Ng - 1, j, k)];
          }
        }
      }
    }
  }
  if (aux) {
    // Aux
    for (int var(0); var < d->Naux; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            aux[ID(var, i, j, k)] = aux[ID(var, d->Ng, j, k)];
            // Right
            aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->nx + d->Ng - 1, j, k)];
          }
        }
      }
    }
  }
  if (d->Ny > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i, d->Ng, k)];
            // Back
            cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->ny + d->Ng - 1, k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, i, j, k)] = prims[ID(var, i, d->Ng, k)];
              // Back
              prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, i, j, k)] = aux[ID(var, i, d->Ng, k)];
              // Back
              aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
  }


  // Ensure there are no artifacts from the boundaries
  // Cons

  // Bottom left
  {
    // Bottom quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-d->Ng; i++) {
        for (int j(0); j < d->Ng+1; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons[ID(var, i, j, k)] = cons[ID(var, i + d->Ng, j + d->Ng, k)];
          }
        }
      }
    }
    // Left quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Ng+1; i++) {
        for (int j(0); j < d->Ny-d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i + d->Ng, j + d->Ng, k)];
          }
        }
      }
    }
  }
  // Top right
  {
    // Top quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-d->Ng; i++) {
        for (int j(0); j < d->Ng+1; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = cons[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
          }
        }
      }
    }
    // Right quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Ng+1; i++) {
        for (int j(0); j < d->Ny-d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = cons[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
          }
        }
      }
    }
  }



  // Prims
  if (prims) {
    // Bottom left
    {
      // Bottom quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              prims[ID(var, i, j, k)] = prims[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
      // Left quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              prims[ID(var, i, j, k)] = prims[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
    }
    // Top right
    {
      // Top quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              prims[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = prims[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
      // Right quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = prims[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
    }
  }

  // Aux
  if (aux) {
    // Bottom left
    {
      // Bottom quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              aux[ID(var, i, j, k)] = aux[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
      // Left quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              aux[ID(var, i, j, k)] = aux[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
    }
    // Top right
    {
      // Top quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              aux[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = aux[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
      // Right quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = aux[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
    }
  }

  // Nothing fancy for Z-direction
  if (d->Nz > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            // Bottom
            cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->Ng)];
            // Top
            cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->nz + d->Ng - 1)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->Ng)];
              // Top
              prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->Ng)];
              // Top
              aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
  }
}


void Periodic::apply(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Cons
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            prims[ID(var, i, j, k)] = prims[ID(var, d->nx + i, j, k)];
            // Right
            prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->Ng + i, j, k)];
          }
        }
      }
    }
  }
  // Aux
  if (aux) {
    for (int var(0); var < d->Naux; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            aux[ID(var, i, j, k)] = aux[ID(var, d->nx + i, j, k)];
            // Right
            aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->Ng + i, j, k)];
          }
        }
      }
    }
  }

  if (d->Ny > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i, d->ny + j, k)];
            // Back
            cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->Ng + j, k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, i, j, k)] = prims[ID(var, i, d->ny + j, k)];
              // Back
              prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->Ng + j, k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, i, j, k)] = aux[ID(var, i, d->ny + j, k)];
              // Back
              aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->Ng + j, k)];
            }
          }
        }
      }
    }
  }

  if (d->Nz > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            // Bottom
            cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->nz + k)];
            // Top
            cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->Ng + k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->nz + k)];
              // Top
              prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->Ng + k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->nz + k)];
              // Top
              aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->Ng + k)];
            }
          }
        }
      }
    }
  }
}

// This applies periodic boundary conditions to the RHS directly.
// It assumes a staggered grid.
// This means the staggered variables have one more "point" to consider
// Cell number   0   1   2         4   5   6
//             | x | x | x | ... | x | x | x |
// Edge number 0   1   2   3     4   5   6   7
void PeriodicRHSStaggered::apply(double * rhs, double * prims, double * aux)
{
  // HACK: generalise
  const int n_staggered_vars(4), n_centered_vars(5);
  int staggered_vars[n_staggered_vars] = {5, 6, 7, 8};
  int centered_vars[n_centered_vars] = {0, 1, 2, 3, 4};
  // Syntax
  Data * d(this->data);

  // RHS
  for (int ivar(0); ivar < n_centered_vars; ivar++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          int var = centered_vars[ivar];
          // Left
          rhs[ID(var, i, j, k)] = rhs[ID(var, d->nx + i, j, k)];
          // Right
          rhs[ID(var, d->nx + d->Ng + i, j, k)] = rhs[ID(var, d->Ng + i, j, k)];
        }
      }
    }
  }
  for (int ivar(0); ivar < n_staggered_vars; ivar++) {
    for (int i(0); i < d->Ng - 1; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          int var = staggered_vars[ivar];
          // Left
          rhs[ID(var, i, j, k)] = rhs[ID(var, d->nx + i + 1, j, k)];
          // Right
          rhs[ID(var, d->nx + d->Ng + i + 1, j, k)] = rhs[ID(var, d->Ng + i, j, k)];
        }
      }
    }
  }

  if (d->Ny > 1) {
    // RHS
    for (int ivar(0); ivar < n_centered_vars; ivar++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            int var = centered_vars[ivar];
            // Front
            rhs[ID(var, i, j, k)] = rhs[ID(var, i, d->ny + j, k)];
            // Back
            rhs[ID(var, i, d->ny + d->Ng + j, k)] = rhs[ID(var, i, d->Ng + j, k)];
          }
        }
      }
    }
    for (int ivar(0); ivar < n_staggered_vars; ivar++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng - 1; j++) {
          for (int k(0); k < d->Nz; k++) {
            int var = staggered_vars[ivar];
            // Front
            rhs[ID(var, i, j, k)] = rhs[ID(var, i, d->ny + j + 1, k)];
            // Back
            rhs[ID(var, i, d->ny + d->Ng + j + 1, k)] = rhs[ID(var, i, d->Ng + j, k)];
          }
        }
      }
    }
  }

  if (d->Nz > 1) {
    // RHS
    for (int ivar(0); ivar < n_centered_vars; ivar++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            int var = centered_vars[ivar];
            // Bottom
            rhs[ID(var, i, j, k)] = rhs[ID(var, i, j, d->nz + k)];
            // Top
            rhs[ID(var, i, j, d->nz + d->Ng + k)] = rhs[ID(var, i, j, d->Ng + k)];
          }
        }
      }
    }
    for (int ivar(0); ivar < n_staggered_vars; ivar++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng - 1; k++) {
            int var = staggered_vars[ivar];
            // Bottom
            rhs[ID(var, i, j, k)] = rhs[ID(var, i, j, d->nz + k + 1)];
            // Top
            rhs[ID(var, i, j, d->nz + d->Ng + k + 1)] = rhs[ID(var, i, j, d->Ng + k)];
          }
        }
      }
    }
  }
}


void Flow::apply(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);
  // Cons
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            prims[ID(var, i, j, k)] = prims[ID(var, d->nx + i, j, k)];
            // Right
            prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->Ng + i, j, k)];
          }
        }
      }
    }
  }
  // Aux
  if (aux) {
    for (int var(0); var < d->Naux; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            aux[ID(var, i, j, k)] = aux[ID(var, d->nx + i, j, k)];
            // Right
            aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->Ng + i, j, k)];
          }
        }
      }
    }
  }

  if (d->Ny > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i, d->Ng, k)];
            // Back
            cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->ny + d->Ng - 1, k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, i, j, k)] = prims[ID(var, i, d->Ng, k)];
              // Back
              prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, i, j, k)] = aux[ID(var, i, d->Ng, k)];
              // Back
              aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
  }

  if (d->Nz > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            // Bottom
            cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->Ng)];
            // Top
            cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->nz + d->Ng - 1)];
          }
        }
      }
    }
    // Prims
    if (prims) {
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->Ng)];
              // Top
              prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->Ng)];
              // Top
              aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
  }
}

//
//
// void ConductingChannel::apply(double * cons, double * prims, double * aux)
// {
//   // Syntax
//   Data * d(this->data);
//   // Cons
//   for (int var(0); var < d->Ncons; var++) {
//     for (int i(0); i < d->Ng; i++) {
//       for (int j(0); j < d->Ny; j++) {
//         for (int k(0); k < d->Nz; k++) {
//           // Left
//           cons[ID(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
//           // Right
//           cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
//         }
//       }
//     }
//   }
//   // Prims
//   if (prims) {
//     for (int var(0); var < d->Nprims; var++) {
//       for (int i(0); i < d->Ng; i++) {
//         for (int j(0); j < d->Ny; j++) {
//           for (int k(0); k < d->Nz; k++) {
//             // Left
//             prims[ID(var, i, j, k)] = prims[ID(var, d->nx + i, j, k)];
//             // Right
//             prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->Ng + i, j, k)];
//           }
//         }
//       }
//     }
//   }
//   // Aux
//   if (aux) {
//     for (int var(0); var < d->Naux; var++) {
//       for (int i(0); i < d->Ng; i++) {
//         for (int j(0); j < d->Ny; j++) {
//           for (int k(0); k < d->Nz; k++) {
//             // Left
//             aux[ID(var, i, j, k)] = aux[ID(var, d->nx + i, j, k)];
//             // Right
//             aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->Ng + i, j, k)];
//           }
//         }
//       }
//     }
//   }
//
//   if (d->Ny > 1) {
//     // Cons
//     for (int var(0); var < 8; var++) {
//       for (int i(0); i < d->Nx; i++) {
//         for (int j(0); j < d->Ng; j++) {
//           for (int k(0); k < d->Nz; k++) {
//             // Front
//             cons[ID(var, i, j, k)] = cons[ID(var, i, d->Ng, k)];
//             // Back
//             cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->ny + d->Ng - 1, k)];
//           }
//         }
//       }
//     }
//     // Prims
//     if (prims) {
//       for (int var(0); var < d->Nprims; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ng; j++) {
//             for (int k(0); k < d->Nz; k++) {
//               // Front
//               prims[ID(var, i, j, k)] = prims[ID(var, i, d->Ng, k)];
//               // Back
//               prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->ny + d->Ng - 1, k)];
//             }
//           }
//         }
//       }
//     }
//     // Aux
//     if (aux) {
//       for (int var(0); var < d->Naux; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ng; j++) {
//             for (int k(0); k < d->Nz; k++) {
//               // Front
//               aux[ID(var, i, j, k)] = aux[ID(var, i, d->Ng, k)];
//               // Back
//               aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->ny + d->Ng - 1, k)];
//             }
//           }
//         }
//       }
//     }
//   }
//
//   if (d->Nz > 1) {
//     // Cons
//     for (int var(0); var < d->Ncons; var++) {
//       for (int i(0); i < d->Nx; i++) {
//         for (int j(0); j < d->Ny; j++) {
//           for (int k(0); k < d->Ng; k++) {
//             // Bottom
//             cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->Ng)];
//             // Top
//             cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->nz + d->Ng - 1)];
//           }
//         }
//       }
//     }
//     // Prims
//     if (prims) {
//       for (int var(0); var < d->Nprims; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ny; j++) {
//             for (int k(0); k < d->Ng; k++) {
//               // Bottom
//               prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->Ng)];
//               // Top
//               prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->nz + d->Ng - 1)];
//             }
//           }
//         }
//       }
//     }
//     // Aux
//     if (aux) {
//       for (int var(0); var < d->Naux; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ny; j++) {
//             for (int k(0); k < d->Ng; k++) {
//               // Bottom
//               aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->Ng)];
//               // Top
//               aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->nz + d->Ng - 1)];
//             }
//           }
//         }
//       }
//     }
//   }
// }

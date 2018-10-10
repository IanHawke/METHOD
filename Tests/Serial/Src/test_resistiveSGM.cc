#include "gtest/gtest.h"
#include "srmhd.h"
#include "simData.h"
#include "simulation.h"
#include "fluxVectorSplitting.h"
#include "resistiveSGM.h"

#define ID(variable, idx, jdx, kdx)  ((variable)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))

// ID for the SGM matrices
// Mx, My, and Mz matrix
#define IDM(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d.Nx)*(d.Ny)*(d.Nz) + (mdx)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))
// dfxdw, dfydw, dfzdw
#define IDFW(ldx, mdx, idx, jdx, kdx)  ((ldx)*(d.Nprims)*(d.Nx)*(d.Ny)*(d.Nz) + (mdx)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))
// dwdsb
#define IDWS(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d.Nx)*(d.Ny)*(d.Nz) + (mdx)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))

// Function required for the memory allocation check
int deref(double * in)
{
    return * in;
}



/******************************************************************************
      See PyMETH for description of tests. These are identical.
******************************************************************************/


namespace
{

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  TEST(RSGM, Reset)
  /*
    Ensure that the reset function sets K and the source vector to zero.
  */
  {
    Data d(20, 20, 20, 0.01, 2.01, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
    SRMHD model(&d);
    Simulation sim(&d);
    FVS fluxMethod(&d, &model);
    ResistiveSGM subgridModel(&d, &fluxMethod);

    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {

          // Set K as non-zero
          subgridModel.K[ID(0, i, j, k)] = 4.0;
          EXPECT_NEAR(subgridModel.K[ID(0, i, j, k)], 4.0, 1e-15);
          // Set source as nonzero
          for (int var(0); var<d.Ncons; var++) {
            d.source[ID(var, i, j, k)] = 4.0;
            EXPECT_NEAR(d.source[ID(var, i, j, k)], 4.0, 1e-15);
          }
        }
      }
    }

    subgridModel.reset(d.source);

    // Check has zero'd everything
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Set K as non-zero
          EXPECT_NEAR(subgridModel.K[ID(0, i, j, k)], 0.0, 1e-15);
          // Set source as nonzero
          for (int var(0); var<d.Ncons; var++) {
            EXPECT_NEAR(d.source[ID(var, i, j, k)], 0.0, 1e-15);
          }
        }
      }
    }
  }


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  TEST(RSGM, DataAssignment1D)
  /*
    Checks that, for 1-dimensional simulations, the variables are set correctly
  */
  {
    Data d(100, 0, 0, 0.01, 2.01, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
    SRMHD model(&d);
    Simulation sim(&d);
    FVS fluxMethod(&d, &model);
    ResistiveSGM subgridModel(&d, &fluxMethod);

    // Set mid point where x=1
    int mid(d.Nx/2-1);
    EXPECT_NEAR(d.x[mid], 1.0, 1e-15);

    // Set primitive variables to known values
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          d.prims[ID(0, i, j, k)] = 0.1;
          d.prims[ID(1, i, j, k)] = 0.0;
          d.prims[ID(2, i, j, k)] = 0.41 * d.x[i];
          d.prims[ID(3, i, j, k)] = 0.51;
          d.prims[ID(4, i, j, k)] = 0.66;
          d.prims[ID(5, i, j, k)] = 0.0;
          d.prims[ID(6, i, j, k)] = 0.22 * d.x[i];
          d.prims[ID(7, i, j, k)] = 0.33;
        }
      }
    }

    // Check element 54 is unchanged by d.x
    EXPECT_NEAR(d.prims[ID(2, mid, 0, 0)], 0.41, 1e-15);
    EXPECT_NEAR(d.prims[ID(6, mid, 0, 0)], 0.22, 1e-15);

    // Set global variables and direction-free matrices (testing so set factor=false)
    subgridModel.set_vars(NULL, d.prims, NULL);
    subgridModel.set_dwdsb(NULL, d.prims, NULL);

    // Set Ex and q manually
    // E = - v cross B    ------>    Ex = -(vy Bz - vvz By) x
    // q = partial_a E^a ------>    q = partial_x Ex
    double Ex( -1*(0.41*0.33 - 0.51*0.22) );
    double q(Ex);

    // Check values of E, q and alpha
    {
      EXPECT_NEAR(subgridModel.E[ID(0, mid, 0, 0)], Ex, 1e-15);
      EXPECT_NEAR(subgridModel.q[ID(0, mid, 0, 0)], q, 1e-15);
      EXPECT_NEAR(subgridModel.alpha[ID(0, mid, 0, 0)], 0.0000001382527749, 1e-14);
    }

    // Check that dwdsb has been set correctly
    {
      // First, check values of A
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 0, mid, 0, 0)], 57.750012326390994*0.0000001382527749, 1e-12);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 1, mid, 0, 0)], 41250.008804565*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 2, mid, 0, 0)], -27500.005869709996*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 0, mid, 0, 0)], -41250.008804565*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 1, mid, 0, 0)], 60.54511232639099*0.0000001382527749, 1e-12);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 2, mid, 0, 0)], 4.19265*0.0000001382527749, 1e-13);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 0, mid, 0, 0)], 27500.005869709996*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 1, mid, 0, 0)], 4.19265*0.0000001382527749, 1e-13);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 2, mid, 0, 0)], 64.038987326391*0.0000001382527749, 1e-12);
      // Second, check values of B
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 0, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 1, mid, 0, 0)], -63114.76360705499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 2, mid, 0, 0)], 52202.88593900499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 0, mid, 0, 0)], 63750.01360705499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 1, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 2, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 0, mid, 0, 0)], -51250.01093900499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 1, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 2, mid, 0, 0)], 0.0, 1e-15);
      // Third, check values of C
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 0, mid, 0, 0)], -125000.02668049998*0.0000001382527749, 1e-8);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 1, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 2, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 0, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 1, mid, 0, 0)], -131050.02668049998*0.0000001382527749, 1e-8);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 2, mid, 0, 0)], -9075.0*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 0, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 1, mid, 0, 0)], -9075.0*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 2, mid, 0, 0)], -138612.52668049998*0.0000001382527749, 1e-8);
      // Finally, check values of D
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 0, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 1, mid, 0, 0)], -1167.1752187800998*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 2, mid, 0, 0)], -1488.2627721411*0.0000001382527749, 1e-10);
      // Just in case, check the rest is zero
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 0, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 1, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 2, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 0, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 1, mid, 0, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 2, mid, 0, 0)], 0.0, 1e-15);
    }
  }


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  TEST(RSGM, DataAssignment2D)
  /*
    Checks that, for 1-dimensional simulations, the variables are set correctly
  */
  {
    Data d(10, 10, 0, 0.1, 2.1, 0.1, 2.1, 0, 1, 0.4, 0.1, 4, 2, 50);
    SRMHD model(&d);
    Simulation sim(&d);
    FVS fluxMethod(&d, &model);
    ResistiveSGM subgridModel(&d, &fluxMethod);

    // Set mid point where x=1
    int mid(d.Ny/2-1);
    EXPECT_NEAR(d.y[mid], 1.0, 1e-15);

    // Set primitive variables to known values
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          d.prims[ID(0, i, j, k)] = 0.1;
          d.prims[ID(1, i, j, k)] = 0.51;
          d.prims[ID(2, i, j, k)] = 0.0;
          d.prims[ID(3, i, j, k)] = 0.41 * d.y[j];
          d.prims[ID(4, i, j, k)] = 0.66;
          d.prims[ID(5, i, j, k)] = 0.33;
          d.prims[ID(6, i, j, k)] = 0.0;
          d.prims[ID(7, i, j, k)] = 0.22 * d.y[j];
        }
      }
    }

    // Check element 54 is unchanged by d.x
    EXPECT_NEAR(d.prims[ID(3, 0, mid, 0)], 0.41, 1e-15);
    EXPECT_NEAR(d.prims[ID(7, 0, mid, 0)], 0.22, 1e-15);

    // Set global variables and direction-free matrices (testing so set factor=false)
    subgridModel.set_vars(NULL, d.prims, NULL);
    subgridModel.set_dwdsb(NULL, d.prims, NULL);

    // Set Ey and q manually
    // E = - v cross B    ------>    Ey = -(vz Bx - vx Bz) y
    // q = partial_a E^a ------>    q = partial_y Ey
    double Ey( -1*(0.41*0.33 - 0.51*0.22) );
    double q(Ey);

    // Check values of E, q and alpha
    {
      EXPECT_NEAR(subgridModel.E[ID(1, mid, mid, 0)], Ey, 1e-15);
      EXPECT_NEAR(subgridModel.q[ID(0, mid, mid, 0)], q, 1e-15);
      EXPECT_NEAR(subgridModel.alpha[ID(0, mid, mid, 0)], 0.0000001382527749, 1e-14);
    }

    // Check that dwdsb has been set correctly
    {
      // First, check values of A
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 0, mid, mid, 0)], 64.03898732639098*0.0000001382527749, 1e-12);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 1, mid, mid, 0)], 27500.005869709996*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 2, mid, mid, 0)], 4.1926499999999995*0.0000001382527749, 1e-13);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 0, mid, mid, 0)], -27500.005869709996*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 1, mid, mid, 0)], 57.75001232639099*0.0000001382527749, 1e-12);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 2, mid, mid, 0)], 41250.008804565*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 0, mid, mid, 0)], 4.192649999999999*0.0000001382527749, 1e-13);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 1, mid, mid, 0)], -41250.008804565*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 2, mid, mid, 0)], 60.545112326390985*0.0000001382527749, 1e-12);
      // Second, check values of B
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 0, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 1, mid, mid, 0)],-51250.01093900499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 2, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 0, mid, mid, 0)], 52202.88593900499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 1, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 2, mid, mid, 0)], -63114.76360705499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 0, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 1, mid, mid, 0)], 63750.01360705499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 2, mid, mid, 0)], 0.0, 1e-15);
      // Third, check values of C
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 0, mid, mid, 0)], -138612.52668049998*0.0000001382527749, 1e-8);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 1, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 2, mid, mid, 0)], -9075.0*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 0, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 1, mid, mid, 0)], -125000.02668049998*0.0000001382527749, 1e-8);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 2, mid, mid, 0)], 0.0, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 0, mid, mid, 0)], -9075.0*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 1, mid, mid, 0)], 0.0, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 2, mid, mid, 0)], -131050.02668049998*0.0000001382527749, 1e-8);
      // Finally, check values of D
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 0, mid, mid, 0)], -1488.2627721411*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 1, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 2, mid, mid, 0)], -1167.1752187800998*0.0000001382527749, 1e-10);
      // Just in case, check the rest is zero
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 0, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 1, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 2, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 0, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 1, mid, mid, 0)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 2, mid, mid, 0)], 0.0, 1e-15);
    }
  }


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  TEST(RSGM, DataAssignment3D)
  /*
    Checks that, for 1-dimensional simulations, the variables are set correctly
  */
  {
    Data d(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50);
    SRMHD model(&d);
    Simulation sim(&d);
    FVS fluxMethod(&d, &model);
    ResistiveSGM subgridModel(&d, &fluxMethod);

    // Set mid point where x=1
    int mid(d.Nz/2-1);
    EXPECT_NEAR(d.z[mid], 1.0, 1e-15);

    // Set primitive variables to known values
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          d.prims[ID(0, i, j, k)] = 0.1;
          d.prims[ID(1, i, j, k)] = 0.41 * d.z[k];
          d.prims[ID(2, i, j, k)] = 0.51;
          d.prims[ID(3, i, j, k)] = 0.0;
          d.prims[ID(4, i, j, k)] = 0.66;
          d.prims[ID(5, i, j, k)] = 0.22 * d.z[k];
          d.prims[ID(6, i, j, k)] = 0.33;
          d.prims[ID(7, i, j, k)] = 0.0;
        }
      }
    }

    // Check element 54 is unchanged by d.x
    EXPECT_NEAR(d.prims[ID(1, 0, 0, mid)], 0.41, 1e-15);
    EXPECT_NEAR(d.prims[ID(5, 0, 0, mid)], 0.22, 1e-15);

    // Set global variables and direction-free matrices (testing so set factor=false)
    subgridModel.set_vars(NULL, d.prims, NULL);
    subgridModel.set_dwdsb(NULL, d.prims, NULL);

    // Set Ey and q manually
    // E = - v cross B    ------>    Ey = -(vz Bx - vx Bz) y
    // q = partial_a E^a ------>    q = partial_y Ey
    double Ez( -1*(0.41*0.33 - 0.51*0.22) );
    double q(Ez);

    // Check values of E, q and alpha
    {
      EXPECT_NEAR(subgridModel.E[ID(2, mid, mid, mid)], Ez, 1e-15);
      EXPECT_NEAR(subgridModel.q[ID(0, mid, mid, mid)], q, 1e-15);
      EXPECT_NEAR(subgridModel.alpha[ID(0, mid, mid, mid)], 0.0000001382527749, 1e-14);
    }

    // Check that dwdsb has been set correctly
    {
      // First, check values of A
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 0, mid, mid, mid)], 60.545112326390985*0.0000001382527749, 1e-12);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 1, mid, mid, mid)], 4.192649999999999*0.0000001382527749, 1e-13);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 2, mid, mid, mid)], -41250.008804565*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 0, mid, mid, mid)], 4.1926499999999995*0.0000001382527749, 1e-13);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 1, mid, mid, mid)], 64.03898732639098*0.0000001382527749, 1e-12);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 2, mid, mid, mid)], 27500.005869709996*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 0, mid, mid, mid)], 41250.008804565*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 1, mid, mid, mid)], -27500.005869709996*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 2, mid, mid, mid)], 57.75001232639099*0.0000001382527749, 1e-12);
      // Second, check values of B
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 0, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 1, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 2, mid, mid, mid)], 63750.01360705499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 0, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 1, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 2, mid, mid, mid)],-51250.01093900499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 0, mid, mid, mid)], -63114.76360705499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 1, mid, mid, mid)], 52202.88593900499*0.0000001382527749, 1e-9);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 2, mid, mid, mid)], 0.0, 1e-15);
      // Third, check values of C
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 0, mid, mid, mid)], -131050.02668049998*0.0000001382527749, 1e-8);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 1, mid, mid, mid)], -9075.0*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 2, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 0, mid, mid, mid)], -9075.0*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 1, mid, mid, mid)], -138612.52668049998*0.0000001382527749, 1e-8);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 2, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 0, mid, mid, mid)], 0.0, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 1, mid, mid, mid)], 0.0, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 2, mid, mid, mid)], -125000.02668049998*0.0000001382527749, 1e-8);
      // Finally, check values of D
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 0, mid, mid, mid)], -1167.1752187800998*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 1, mid, mid, mid)], -1488.2627721411*0.0000001382527749, 1e-10);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 2, mid, mid, mid)], 0.0, 1e-15);
      // Just in case, check the rest is zero
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 0, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 1, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 2, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 0, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 1, mid, mid, mid)], 0.0, 1e-15);
      EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 2, mid, mid, mid)], 0.0, 1e-15);
    }
  }

  TEST(RSGM, Directionality)
  {

  }

  TEST(RSGM, RotationallyInvariant)
  {
    
  }


  // SLOWWWWWW
  // TEST(RSGM, MemoryAllocation)
  // /*
  //   Ensure constructor allocates memory for all the necessary arrays
  // */
  // {
  //   Data d(7, 5, 1, 0, 1, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
  //   SRMHD model(&d);
  //   Simulation sim(&d);
  //   FVS fluxMethod(&d, &model);
  //   ResistiveSGM subgridModel(&d, &fluxMethod);
  //
  //   for (int i(0); i<d.Nx; i++) {
  //     for (int j(0); j<d.Ny; j++) {
  //       for (int k(0); k<d.Nz; k++) {
  //         printf("(%d, %d, %d)\n", i, j, k);
  //         // dfxdw, dfydw, dfzdw
  //         for (int l(0); l<d.Ncons; l++) {
  //           for (int m(0); m<d.Nprims; m++) {
  //             ASSERT_EXIT((deref(&subgridModel.dfxdw[IDFW(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //             ASSERT_EXIT((deref(&subgridModel.dfydw[IDFW(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //             ASSERT_EXIT((deref(&subgridModel.dfzdw[IDFW(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           }
  //         }
  //
  //         // dwdsb
  //         for (int l(0); l<d.Nprims; l++) {
  //           for (int m(0); m<3; m++) {
  //             ASSERT_EXIT((deref(&subgridModel.dwdsb[IDWS(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           }
  //         }
  //
  //         // E
  //         for (int l(0); l<3; l++) {
  //             ASSERT_EXIT((deref(&subgridModel.E[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //         }
  //
  //         // q
  //         ASSERT_EXIT((deref(&subgridModel.q[ID(0, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //
  //         // K
  //         for (int l(0); l<3; l++) {
  //             ASSERT_EXIT((deref(&subgridModel.K[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //         }
  //
  //         // Mx, My, Mz
  //         for (int l(0); l<d.Ncons; l++) {
  //           for (int m(0); m<3; m++) {
  //             ASSERT_EXIT((deref(&subgridModel.Mx[IDM(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //             ASSERT_EXIT((deref(&subgridModel.My[IDM(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //             ASSERT_EXIT((deref(&subgridModel.Mz[IDM(l, m, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           }
  //         }
  //
  //         // Stiff flux x, y, z
  //         for (int l(0); l<3; l++) {
  //           ASSERT_EXIT((deref(&subgridModel.fx[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           ASSERT_EXIT((deref(&subgridModel.fy[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           ASSERT_EXIT((deref(&subgridModel.fz[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //         }
  //
  //         // Diffusion vector
  //         for (int l(0); l<3; l++) {
  //           ASSERT_EXIT((deref(&subgridModel.diffuX[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           ASSERT_EXIT((deref(&subgridModel.diffuY[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //           ASSERT_EXIT((deref(&subgridModel.diffuZ[ID(l, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //         }
  //
  //         // alpha
  //         ASSERT_EXIT((deref(&subgridModel.alpha[ID(0, i, j, k)]), exit(0)), ::testing::ExitedWithCode(0),".*");
  //       }
  //     }
  //   }
  //
  //
  // }
}

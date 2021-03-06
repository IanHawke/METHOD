#include "gtest/gtest.h"
#include "saveData.h"
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include <cstdlib>
#include <stdexcept>

Data data(100, 10, 2, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8);

TEST(Simulation, dataInitialisation)
{

  EXPECT_THROW( Simulation sim(&data), std::runtime_error);

  SRMHD model(&data);
  Simulation sim(&data);

  // Check standard data
  EXPECT_EQ(sim.data->nx, 100);
  EXPECT_EQ(sim.data->ny, 10);
  EXPECT_EQ(sim.data->nz, 2);
  EXPECT_EQ(sim.data->Ng, 4);
  EXPECT_EQ(sim.data->Nx, 108);
  EXPECT_EQ(sim.data->Ny, 18);
  EXPECT_EQ(sim.data->Nz, 10);
  EXPECT_EQ(sim.data->xmin, 0.0);
  EXPECT_EQ(sim.data->xmax, 1.0);
  EXPECT_EQ(sim.data->ymin, -0.5);
  EXPECT_EQ(sim.data->ymax, 0.5);
  EXPECT_EQ(sim.data->zmin, -0.1);
  EXPECT_EQ(sim.data->zmax, 0.1);
  EXPECT_EQ(sim.data->endTime, 0.8);
  EXPECT_EQ(sim.data->cfl, 0.5);
  EXPECT_EQ(sim.data->gamma, 5.0/3);
  EXPECT_EQ(sim.data->sigma, 1000.0);
  EXPECT_EQ(sim.data->memSet, 1);
  EXPECT_EQ(sim.data->Ncons, 9);
  EXPECT_EQ(sim.data->Nprims, 8);
  EXPECT_EQ(sim.data->Naux, 13);
  EXPECT_EQ(sim.data->alphaX, 1.0);
  EXPECT_EQ(sim.data->alphaY, 1.0);
  EXPECT_EQ(sim.data->alphaZ, 1.0);
  EXPECT_EQ(sim.data->t, 0.0);
  EXPECT_NEAR(sim.data->dt, 0.002886751346, 1e-11);
  EXPECT_EQ(sim.data->dx, 0.01);
  EXPECT_EQ(sim.data->dy, 0.1);
  EXPECT_EQ(sim.data->dz, 0.1);
  EXPECT_EQ(sim.data->iters, 0);

  // Check domain
  /* Note: Ng = 4 */
  EXPECT_EQ(sim.data->x[4], 0.005);
  EXPECT_EQ(sim.data->x[1], -0.025);
  EXPECT_EQ(sim.data->y[4], -0.45);
  EXPECT_EQ(sim.data->y[1], -0.75);
  EXPECT_EQ(sim.data->x[99+4], 0.995);


}

//! Check that the fields dont change if the system if homogenous
TEST(Simulation, equilibriumSimulation)
{
  Data data(30, 30, 10, 0, 1, 0, 1, 0, 1, 0.1);
  SRMHD model(&data);
  FVS fluxMethod(&data, &model);
  Simulation sim(&data);
  OTVortexSingleFluid init(&data);
  Periodic bcs(&data);
  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  SaveData save(&data);

  for (int i(0); i < data.Nx; i++) {
    for (int j(0); j < data.Ny; j++) {
      for (int k(0); k < data.Nz; k++) {
        data.prims[data.id(0, i, j, k)] = 0.125;
        data.prims[data.id(1, i, j, k)] = 0.1;
        data.prims[data.id(2, i, j, k)] = 0.1;
        data.prims[data.id(3, i, j, k)] = 0.1;
        data.prims[data.id(4, i, j, k)] = 0.1;
        data.prims[data.id(5, i, j, k)] = 0.125;
        data.prims[data.id(6, i, j, k)] = 0.125;
        data.prims[data.id(7, i, j, k)] = 0.125;
      }
    }
  }

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  sim.evolve();


  for (int i(0); i < data.Nx; i++) {
    for (int j(0); j < data.Ny; j++) {
      for (int k(0); k < data.Nz; k++) {
        EXPECT_NEAR(data.prims[data.id(0, i, j, k)], 0.125, 1e-15);
        EXPECT_NEAR(data.prims[data.id(1, i, j, k)], 0.1, 1e-15);
        EXPECT_NEAR(data.prims[data.id(2, i, j, k)], 0.1, 1e-15);
        EXPECT_NEAR(data.prims[data.id(3, i, j, k)], 0.1, 1e-15);
        EXPECT_NEAR(data.prims[data.id(4, i, j, k)], 0.1, 1e-15);
        EXPECT_NEAR(data.prims[data.id(5, i, j, k)], 0.125, 1e-15);
        EXPECT_NEAR(data.prims[data.id(6, i, j, k)], 0.125, 1e-15);
        EXPECT_NEAR(data.prims[data.id(7, i, j, k)], 0.125, 1e-15);
      }
    }
  }
}

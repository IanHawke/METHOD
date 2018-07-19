#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "backwardsRK.h"
#include "SSP2.h"
#include "SSP3.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "saveData.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))


using namespace std;

int main(int argc, char *argv[]) {


  const double MU(1000);
  // Set up domain
  int nx(120);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.4);
  double cfl(0.2);
  int Ng(4);
  double gamma(4.0/3.0);
  double sigma(10000);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(60);
  bool output(false);
  int safety(500);
  int tpb(1);
  int bpg(19200);

  char * ptr(0);
  double tmp(0);
  //! Overwrite any variables that have been passed in as main() arguments
  for (int i(0); i < argc; i++) {
    if (strcmp(argv[i], "nx") == 0) {
      nx = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "ny") == 0) {
      ny = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "nz") == 0) {
      nz = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "safety") == 0) {
      safety = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "output") == 0) {
      output = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "frameSkip") == 0) {
      frameSkip = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "gammanum") == 0) {
      tmp = (double)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "gammaden") == 0 && tmp!=0) {
      gamma = tmp/(double)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "endTime") == 0 && tmp!=0) {
      endTime = (double)strtol(argv[i+1], &ptr, 10);
    }
  }


  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip, tpb, bpg);

  // Choose particulars of simulation
  TwoFluidEMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  BrioWuTwoFluid init(&data);

  Outflow bcs(&data);

  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();
  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  // printf("Parallel speedup of %6.3fx\n",3.45/timeTaken);

  return 0;

}

#include "simData.h"
#include <stdexcept>
#include <cstdio>

Data::Data(int nx, int ny, int nz,
           double xmin, double xmax,
           double ymin, double ymax,
           double zmin, double zmax,
           double endTime, double cfl, int Ng,
           double gamma, double sigma,
           double cp,
           double mu1, double mu2, int frameSkip)
           :
           nx(nx), ny(ny), nz(nz),
           xmin(xmin), xmax(xmax),
           ymin(ymin), ymax(ymax),
           zmin(zmin), zmax(zmax),
           endTime(endTime), cfl(cfl), Ng(Ng),
           gamma(gamma), sigma(sigma),
           memSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(cp),
           mu1(mu1), mu2(mu2), frameSkip(frameSkip)
{

  this->Nx = nx + 2 * Ng;
  this->Ny = ny + 2 * Ng;
  this->Nz = nz + 2 * Ng;
  dims = 3;

  // Catch 2D case
  if (nz == 0) {
    this->Nz = 1;
    zmin = -1e20;
    zmax = 1e20;
    dims = 2;
  }
  // Catch 1D case
  if (ny == 0) {
    this->Nz = this->Ny = 1;
    zmin = ymin = -1e20;
    zmax = ymax = 1e20;
    dims = 1;
  }
  // Ensure there is some Resistivity
  if (this->sigma < 0.0) {
    throw std::invalid_argument("Conductivity must be non-negative, sigma >= 0.\n");
  }
  // Ensure charges are correct way round
  if (this->mu1 > 0.0 or this->mu2 < 0.0) {
    throw std::invalid_argument("Species 1 must have negative charge, mu1 < 0, and species 2 must have positive charge, mu2 > 0.\n");
  }

  // Determine the specs of the GPU(s) and thus set details in simData
  cudaGetDeviceCount(&GPUcount);
  cudaGetDeviceProperties(&prop, 0);
  printf("totGlobMem = %zu\n", prop.totalGlobalMem);
  printf("Shared mem per multiprocessor = %zu\n", prop.sharedMemPerMultiprocessor);
  printf("GPU name: %s\n", prop.name);
  printf("Shared mem per block = %zu\n", prop.sharedMemPerBlock);
  printf("Max threads per multiprocessor = %i\n", prop.maxThreadsPerMultiProcessor);
  printf("Number of multiprocessors = %i\n", prop.multiProcessorCount);
  printf("Global L1 cahche supported = %i\n", prop.globalL1CacheSupported);
  printf("Local L1 cahche supported = %i\n", prop.localL1CacheSupported);
  printf("Shared mem per multiprocessor = %zu\n", prop.sharedMemPerMultiprocessor);
  printf("L2 cache size = %i\n", prop.l2CacheSize);
  printf("Total global memory = %ld\n", prop.totalGlobalMem);
  printf("Execute kernels concurrently? (1/0) = %d\n", prop.concurrentKernels);
  printf("Compute Capability (major) %d\n", prop.major);

}
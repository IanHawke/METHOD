#include "fluxVectorSplitting.h"
#include "cudaErrorCheck.h"
#include <iostream>
#include <cassert>
// Macro for getting array index
#define ID(variable, idx, jdx, kdx)  ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDZ(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDY(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx) + (kdx)*(d->Ny))
#define IDX(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx) + (jdx)*(d->Nz)*(d->Nx) + (kdx)*(d->Nx))

/*
  For loading data into contiguous arrays we need to change the indexing macro.
    In its original form, the data is contiguous in the z-direction, with y and
  then x as the next fastest moving index---so IDZ(var, i, j, k) is identical to
  the normal indexer ID.
    To get data contiguous in y-direction, in a loop assign the elements
  IDY(var, i, j, k) = ID(var, i, j, k), now data is contiguous in y with z then x
  as next fastest moving index.
    Similarly, IDX(var, i, j, k) = ID(var, i, j, k) with arrange data contiguously
  in the x-direction with z and y the next fastest moving index.

  To transform back, apply the same indexing trick.

  Example:
    Flux reconstruction in y direction...

      f = fluxVector(dir = y);

    Rearrange so data is contiguous in y direction...

      for var in Ncons, i in Nx, j in Ny, k in Nz:
        fcontig[IDY(var, i, j, k)] = f[ID(var, i, j, k)]

    Reconstruct flux...

      frecon = fluxRecon(dir = y);

    Copy back data into original form...

      for var in Ncons, i in Nx, j in Ny, k in Nz:
        fnet[ID(var, i, j, k)] = frecon[IDY(var, i, j, k)]
*/





__global__
static void fluxRecon(double * cons, double * f, int stream, int width, double delta, int dir, long unsigned int Ntot)
{
  // Order of weno scheme
  const int order(2);

  // Up and downwind fluxes
  extern __shared__  double ftmp [];
  double * fplus = ftmp;
  double * fminus = ftmp + blockDim.x;
  double * frec = ftmp + 2 * blockDim.x;
  const int lID(threadIdx.x + blockIdx.x * (blockDim.x - 2*order)); // In this stream
  const int gID(lID + stream * (width - 2*order));                  // GlobalID

  // Load data into shared memory whilst applying Lax-Friedrichs approximation of flux
  if (lID < width) {
    double tempf = f[lID];
    double tempc = cons[lID];             //    USE THIS IN FUTURE TO SAVE MEMORY ACCESS
    fplus[threadIdx.x] = 0.5 * (tempf + tempc);
    fminus[threadIdx.x] = 0.5 * (tempf - tempc);
  }

  __syncthreads();

  if (threadIdx.x >= order && threadIdx.x <= blockDim.x-order+1 && lID < width) {
    frec[threadIdx.x] = weno3_upwind(fplus[threadIdx.x-order],
                                     fplus[threadIdx.x-order+1],
                                     fplus[threadIdx.x-order+2]) +
                        weno3_upwind(fminus[threadIdx.x+order-1],
                                     fminus[threadIdx.x+order-2],
                                     fminus[threadIdx.x+order-3]);
  }

  //! Now we are going to use the device array 'f' as the differenced, reconstructed flux vector, i.e. fnet in serial code
  __syncthreads();

  if (threadIdx.x > order && threadIdx.x <= blockDim.x-order && lID < width && gID < Ntot) {
    f[lID] = frec[threadIdx.x+1] / delta - frec[threadIdx.x] / delta;
  }
}



FVS::FVS(Data * data, Model * model, Bcs * bcs) : FluxMethod(data, model, bcs)
{
  // Syntax
  Data * d(this->data);

  // Order of weno scheme
  int order(2);

  // Total number of cells to send
  long unsigned int Ntot(d->Ncons * d->Nx * d->Ny * d->Nz);

  // Define thread set up
  TpB = 512;
  BpG = 40;
  // Resulting size of stream...
  Cwidth = BpG * (TpB - 2*order);
  originalWidth = width = Cwidth + 2*order;

  // ...means we need this many streams
  Nstreams = (Ntot / Cwidth) + 1;

  // Corresponding size of memcpys
  inMemsize = sizeof(double) * width;
  outMemsize = sizeof(double) * Cwidth;

  // Size of dynamically allocd __shared__ memory in device
  sharedMemUsagePerBlock = TpB * 3 * sizeof(double);

  assert(sharedMemUsagePerBlock <= d->prop.sharedMemPerBlock);

  printf("BPG = %d, TPB = %d\n", BpG, TpB);
  printf("Width = %d, Cwidth = %d, Ntot = %lu\n", width, Cwidth, Ntot);
  printf("Shared mem usage = %lu\n", sharedMemUsagePerBlock);

  // Allocate device arrays for each stream
  cons_d = new double*[Nstreams];
  flux_d = new double*[Nstreams];

  for (int i(0); i < Nstreams; i++) {
    gpuErrchk( cudaMalloc((void **)&cons_d[i], inMemsize) );
    gpuErrchk( cudaMalloc((void **)&flux_d[i], inMemsize) );
  }
  gpuErrchk( cudaHostAlloc((void **)&cons_h, Ntot * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&flux_h, Ntot * sizeof(double), cudaHostAllocPortable) );

  // Create streams
  stream = new cudaStream_t[Nstreams];
  printf("Created %d streams\n\n\n", Nstreams);
  for (int i(0); i<Nstreams; i++) {
    gpuErrchk( cudaStreamCreate(&stream[i]) );
  }

  if (d->Nz > 1) {
    gpuErrchk( cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&fz, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable) );
  }
  else if (d->Ny > 1) {
    gpuErrchk( cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable) );
  }
  else {
    gpuErrchk( cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable) );
  }
}

FVS::~FVS()
{
  // Syntax
  Data * d(this->data);
  for (int i(0); i < Nstreams; i++) {
    gpuErrchk( cudaFree(flux_d[i]) );
    gpuErrchk( cudaFree(cons_d[i]) );
  }
  gpuErrchk( cudaFreeHost(flux_h) );
  gpuErrchk( cudaFreeHost(cons_h) );
  delete [] flux_d;
  delete [] cons_d;
  delete [] stream;

  if (d->Nz > 1) {
    gpuErrchk( cudaFreeHost(fx) );
    gpuErrchk( cudaFreeHost(fy) );
    gpuErrchk( cudaFreeHost(fz) );
  }
  else if (d->Ny > 1) {
    gpuErrchk( cudaFreeHost(fx) );
    gpuErrchk( cudaFreeHost(fy) );
  }
  else {
    gpuErrchk( cudaFreeHost(fx) );
  }
}


void FVS::fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir)
{
  // Syntax
  Data * d(this->data);



  // Order of weno scheme
  int order(2);
  double delta;
  // Total number of data points for each vector
  int Ntot(d->Ncons * d->Nx * d->Ny * d->Nz);
  // Get flux vector
  this->model->fluxVector(cons, prims, aux, f, dir);

  // int count(0);
  // for (int var(0); var<d->Ncons; var++) {
  //   for (int i(0); i<d->Nx; i++) {
  //     for (int j(0); j<d->Ny; j++) {
  //         cons[ID(var, i, j, 0)] = f[ID(var, i, j, 0)] = count++;
  //     }
  //   }
  // }

  // Data must be loaded into device contiguously, so will have to rearrange
  if (dir==0) {
    delta = d->dx;
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            flux_h[IDX(var, i, j, k)] = f   [ID(var, i, j, k)];
            cons_h[IDX(var, i, j, k)] = cons[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  else if (dir==1) {
    delta = d->dy;
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            flux_h[IDY(var, i, j, k)] = f   [ID(var, i, j, k)];
            cons_h[IDY(var, i, j, k)] = cons[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  else {
    delta = d->dz;
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            flux_h[IDZ(var, i, j, k)] = f   [ID(var, i, j, k)];
            cons_h[IDZ(var, i, j, k)] = cons[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  // Data is now contiguous, send to GPU
  int lb, rb; // Left and right boundary of data sent to device
  // Set/Reset width and memsize...
  Cwidth = BpG * (TpB - 2*order);
  width = Cwidth + 2*order;

  // Corresponding size of memcpys
  inMemsize = sizeof(double) * width;
  outMemsize = sizeof(double) * Cwidth;

  // Call parallel reconstruction
  for (int i(0); i<Nstreams; i++) {
    // printf("Running stream %d\n", i);
    // First determine where in the contiguous array the left boundary of this stream corresponds to
    lb = i*(width - 2 * order);
    rb = lb + width;
    if (i == Nstreams-1) {
      rb = Ntot;
      // Final stream so only do remaining cells
      width = rb - lb;
      Cwidth = width - 2*order;
      inMemsize = sizeof(double) * width;
      outMemsize = sizeof(double) * Cwidth;
    }
    // Copy stream's data to device
    gpuErrchk( cudaMemcpyAsync(cons_d[i], cons_h + lb, inMemsize, cudaMemcpyHostToDevice, stream[i]) );
    gpuErrchk( cudaMemcpyAsync(flux_d[i], flux_h + lb, inMemsize, cudaMemcpyHostToDevice, stream[i]) );

    fluxRecon<<<BpG, TpB, sharedMemUsagePerBlock, stream[i]>>>(cons_d[i], flux_d[i], i, originalWidth, delta, dir, Ntot);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaStreamSynchronize(stream[i]) );
    gpuErrchk( cudaMemcpyAsync(flux_h+lb+order, flux_d[i]+order, outMemsize, cudaMemcpyDeviceToHost, stream[i]) );
    gpuErrchk( cudaStreamSynchronize(stream[i]) );
  }

  // Data must be loaded back into original order on the host
  if (dir==0) {
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            frecon[ID(var, i, j, k)] = flux_h[IDX(var, i, j, k)];
          }
        }
      }
    }
  }
  else if (dir==1) {
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            frecon[ID(var, i, j, k)] = flux_h[IDY(var, i, j, k)];
          }
        }
      }
    }
  }
  else {
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            frecon[ID(var, i, j, k)] = flux_h[IDZ(var, i, j, k)];
          }
        }
      }
    }
  }

}

void FVS::F(double * cons, double * prims, double * aux, double * f, double * fnet)
{
  // Syntax
  Data * d(this->data);


  // 3D domain, loop over all cells determining the net flux
  if (d->Ny > 1 && d->Nz > 1) {

    // Determine flux vectors
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    this->fluxReconstruction(cons, prims, aux, f, fy, 1);
    this->fluxReconstruction(cons, prims, aux, f, fz, 2);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
        for (int j(0); j < d->Ny-1; j++) {
          for (int k(0); k < d->Nz-1; k++) {
            fnet[ID(var, i, j, k)] = fx[ID(var, i, j, k)] + fy[ID(var, i, j, k)] + fz[ID(var, i, j, k)];
          }
        }
      }
    }
  }

  // 2D domain, loop over x- and y-directions determining the net flux
  else if (d->Ny > 1) {
    // printf("\nCalling x-direction reconstruction\n");
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    // printf("\nCalling y-direction reconstruction\n");
    this->fluxReconstruction(cons, prims, aux, f, fy, 1);
    // printf("Y-direction exited\n");
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
        for (int j(0); j < d->Ny-1; j++) {
          fnet[ID(var, i, j, 0)] = fx[ID(var, i, j, 0)] + fy[ID(var, i, j, 0)];
        }
      }
    }
    // printf("\nfnet = %19.16f: xdir = %19.16f, ydir = %19.16f\n", fnet[ID(6, 38, 5, 0)], fx[ID(6, 38, 5, 0)], fy[ID(6, 38, 5, 0)]);
  }

  // Otherwise, domain is 1D only loop over x direction
  else {
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
          fnet[ID(var, i, 0, 0)] = fx[ID(var, i, 0, 0)];
      }
    }
  }
  // printf("\nfnet = %19.16f\n", fnet[ID(6, 38, 5, 0)]);
  // printf("culprit: x = %d, y = %d\n", IDX(6, 38, 5, 0), IDY(6, 38, 5, 0));
  // exit(1);
}

#include "SSP3.h"
#include "cudaErrorCheck.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "twoFluidEMHD.h"
#include "deviceArguments.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <omp.h>


// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDCons(var, idx, jdx, kdx) ( (var) + (idx)*(d->Ncons)*(d->Nz)*(d->Ny) + (jdx)*(d->Ncons)*(d->Nz) + (kdx)*(d->Ncons)  )
#define IDPrims(var, idx, jdx, kdx) ( (var) + (idx)*(d->Nprims)*(d->Nz)*(d->Ny) + (jdx)*(d->Nprims)*(d->Nz) + (kdx)*(d->Nprims)  )
#define IDAux(var, idx, jdx, kdx) ( (var) + (idx)*(d->Naux)*(d->Nz)*(d->Ny) + (jdx)*(d->Naux)*(d->Nz) + (kdx)*(d->Naux)  )

// Device function for stage one of IMEX rootfind
__global__
static void stageOne(double * sol, double * cons, double * prims, double * aux, double * source,
              double * wa, double * fvec, const double dt, const double gam, const double tol, const int stream,
              const int origWidth, const int streamWidth, const int Ncons, const int Nprims, const int Naux, const int lwa,
              const double gamma, const double sigma, const double mu1, const double mu2, const double cp,
              ModelType modType_t);

// Device function for stage two of IMEX rootfind
__global__
static void stageTwo(double * sol, double * cons, double * prims, double * aux, double * source,
              double * cons1, double * source1, double * flux1, double * wa, double * fvec,
              const double dt, const double gam, const double tol, const int stream,
              const int origWidth, const int streamWidth, const int Ncons, const int Nprims, const int Naux, const int lwa,
              const double gamma, const double sigma, const double mu1, const double mu2, const double cp,
              ModelType modType_t);
// Device function for stage three of IMEX rootfind
__global__
static void stageThree(double * sol, double * cons, double * prims, double * aux, double * source,
              double * cons1, double * source1, double * flux1, double * flux2, double * wa, double * fvec,
              const double dt, const double gam, const double tol, const int stream,
              const int origWidth, const int streamWidth, const int Ncons, const int Nprims, const int Naux, const int lwa,
              const double gamma, const double sigma, const double mu1, const double mu2, const double cp,
              ModelType modType_t);



//! Device residual functions for stage one of IMEX SSP3
__device__
int IMEX3Residual1Parallel(void *p, int n, const double *x, double *fvec, int iflag);
//! Device residual functions for stage two of IMEX SSP3
__device__
int IMEX3Residual2Parallel(void *p, int n, const double *x, double *fvec, int iflag);
//! Device residual functions for stage three of IMEX SSP3
__device__
int IMEX3Residual3Parallel(void *p, int n, const double *x, double *fvec, int iflag);

//! BackwardsRK parameterized constructor
SSP3::SSP3(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              SSP2(data, model, bc, fluxMethod)

{
  Data * d(this->data);
  this->args = IMEX3Arguments(data);

  // Need work arrays
  cudaHostAlloc((void **)&x, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fvec, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&wa, sizeof(double) * lwa,
                cudaHostAllocPortable);
  // Interstage results
  cudaHostAlloc((void **)&U3, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U3guess, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&source3, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux3, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempprims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempaux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);
}

SSP3::~SSP3()
{

  // Clean up your mess
  cudaFreeHost(x);
  cudaFreeHost(fvec);
  cudaFreeHost(wa);
  cudaFreeHost(U3);
  cudaFreeHost(U3guess);
  cudaFreeHost(source3);
  cudaFreeHost(flux3);
  cudaFreeHost(tempprims);
  cudaFreeHost(tempaux);
}

//! Single step functions
void SSP3::step(double * cons, double * prims, double * aux, double dt)
{
  int i, j, k, var;
  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;

  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons ; var++) U1[ID(var, i, j, k)]  = cons[ID(var, i, j, k)];
        #if MATCH_SERIAL
          for (var = 0; var < d->Nprims ; var++) holdPrims[ID(var, i, j, k)]  = prims[ID(var, i, j, k)]; // holdVars ---> see constructor
          for (var = 0; var < d->Naux ; var++) holdAux[ID(var, i, j, k)]  = aux[ID(var, i, j, k)]; // holdVars ---> see constructor
        #endif
      }
    }
  }
  #if MATCH_SERIAL
    callStageOne(U1, holdPrims, holdAux, source1, dt);
    this->fluxMethod->F(U1, holdPrims, holdAux, d->f, flux1);
  #else
    callStageOne(U1, prims, aux, source1, dt);
    this->fluxMethod->F(U1, prims, aux, d->f, flux1);
  #endif

  this->bcs->apply(U1);
  this->bcs->apply(flux1);

  //########################### STAGE TWO #############################//
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons ; var++) U2[ID(var, i, j, k)]  = cons[ID(var, i, j, k)];
      }
    }
  }

  #if MATCH_SERIAL
    callStageTwo(U2, holdPrims, holdAux, source2, U1, source1, flux1, dt);
    this->model->getPrimitiveVars(U2, holdPrims, holdAux);
    this->fluxMethod->F(U2, holdPrims, holdAux, d->f, flux2);
  #else
    callStageTwo(U2, prims, aux, source2, U1, source1, flux1, dt);
    this->fluxMethod->F(U2, prims, aux, d->f, flux2);
  #endif
  this->bcs->apply(flux2);


  //########################### STAGE THREE ##############################//
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons ; var++) U3[ID(var, i, j, k)]  = cons[ID(var, i, j, k)];
        #if MATCH_SERIAL
          for (var = 0; var < d->Nprims ; var++) holdPrims[ID(var, i, j, k)]  = prims[ID(var, i, j, k)]; // holdVars ---> see constructor
          for (var = 0; var < d->Naux ; var++) holdAux[ID(var, i, j, k)]  = aux[ID(var, i, j, k)]; // holdVars ---> see constructor
        #endif
      }
    }
  }
  #if MATCH_SERIAL
    callStageThree(U3, holdPrims, holdAux, source3, source1, flux1, U2, flux2, dt);
    this->bcs->apply(U3, holdPrims, holdAux);
    this->model->getPrimitiveVars(U3, holdPrims, holdAux);
    this->model->sourceTerm(U3, holdPrims, holdAux, source1);
    this->fluxMethod->F(U3, holdPrims, holdAux, d->f, flux3);
  #else
    callStageThree(U3, prims, aux, source3, source1, flux1, U2, flux2, dt);
    this->fluxMethod->F(U3, prims, aux, d->f, flux3);
  #endif
  this->bcs->apply(flux3);


  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(is); i < ie; i++) {
      for (int j(js); j < je; j++) {
        for (int k(ks); k < ke; k++) {
          cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - dt *
                    (flux1[ID(var, i, j, k)] + flux2[ID(var, i, j, k)] + 4*flux3[ID(var, i, j, k)]) / 6.0 +
                     dt * (source1[ID(var, i, j, k)] + source2[ID(var, i, j, k)] + 4*source3[ID(var, i, j, k)]) / 6.0;
        }
      }
    }
  }
  this->model->getPrimitiveVars(cons, prims, aux);
  this->bcs->apply(cons, prims, aux);
}

void SSP3::callStageOne(double * cons, double * prims, double * aux, double * source, double dt)
{
  int i, j, k, var;

  Data * d(this->data);
  //########################### STAGE ONE #############################//
  // First need to copy data to the device
  // A single cell requires all cons, prims and aux for the step. Rearrange so
  // we can copy data in contiguous way
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons; var++)  args.cons_h [IDCons(var, i, j, k) ] = cons[ID(var, i, j, k)];
        for (var = 0; var < d->Nprims; var++) args.prims_h[IDPrims(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (var = 0; var < d->Naux; var++)   args.aux_h[IDAux(var, i, j, k)    ] = aux[ID(var, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < d->Nstreams; i++) {

    // Which cell is at the left bound?
    int lcell(i * d->tpb * d->bpg);
    // Which cell is at the right bound?
    int rcell(lcell + d->tpb * d->bpg);
    if (rcell > d->Ncells) rcell = d->Ncells; // Dont overshoot
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));
    // Send stream's data
    gpuErrchk( cudaMemcpyAsync(args.cons_d[i], args.cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_d[i], args.prims_h + lcell*d->Nprims, inMemsize*d->Nprims, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_d[i], args.aux_h + lcell*d->Naux, inMemsize*d->Naux, cudaMemcpyHostToDevice, args.stream[i]) );

    int sharedMem((d->Ncons + d->Ncons + d->Ncons + d->Nprims + d->Naux) * sizeof(double) * d->tpb);
    // Call kernel and operate on data
    stageOne <<< d->bpg, d->tpb, sharedMem, args.stream[i] >>>
            (args.sol_d[i], args.cons_d[i], args.prims_d[i], args.aux_d[i],
            args.source_d[i], args.wa_d[i], args.fvec_d[i], dt, args.gam, tol, i, d->tpb * d->bpg,
            width, d->Ncons, d->Nprims, d->Naux, lwa,
            d->gamma, d->sigma, d->mu1, d->mu2, d->cp,
            model->modType_t);
    cudaStreamSynchronize(args.stream[i]);
    gpuErrchk( cudaPeekAtLastError() );

    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(args.sol_h + lcell*d->Ncons, args.sol_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_h + lcell*d->Nprims, args.prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_h + lcell*d->Naux, args.aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.source_h + lcell*d->Ncons, args.source_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );


  // Rearrange data back into arrays
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons; var++)  cons[ID(var, i, j, k)]   = args.sol_h[IDCons(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  source[ID(var, i, j, k)] = args.source_h[IDCons(var, i, j, k)];
        for (var = 0; var < d->Nprims; var++) prims[ID(var, i, j, k) ] = args.prims_h[IDPrims(var, i, j, k)];
        for (var = 0; var < d->Naux; var++)   aux[ID(var, i, j, k)]    = args.aux_h[IDAux(var, i, j, k)];
      }
    }
  }
}

__global__
void stageOne(double * sol, double * cons, double * prims, double * aux, double * source,
              double * wa, double * fvec, const double dt, const double gam, const double tol, const int stream,
              const int origWidth, const int streamWidth, const int Ncons, const int Nprims, const int Naux, const int lwa,
              const double gamma, const double sigma, double mu1, const double mu2, const double cp,
              ModelType modType_t)
{
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  const int gID(lID + stream * origWidth);        //!< global index (in domain)
  extern __shared__ double sharedMem[];
  double * CONS = &sharedMem[tID * (Ncons + Ncons + Ncons + Nprims + Naux)];
  double * PRIMS = &CONS[Ncons];
  double * AUX = &PRIMS[Nprims];
  double * SOURCE = &AUX[Naux];
  double * SOL = &SOURCE[Ncons];
  int info;

  if (lID < streamWidth)
  {
    // Load data into shared memory
    for (int i(0); i < Ncons; i++) CONS[i] = cons[i + lID * Ncons];
    for (int i(0); i < Nprims; i++) PRIMS[i] = prims[i + lID * Nprims];
    for (int i(0); i < Naux; i++) AUX[i] = aux[i + lID * Naux];
    for (int i(0); i < Ncons; i++) SOURCE[i] = source[i + lID * Ncons];

    Model_D * model_d;
    // Store pointers to devuce arrays in the structure
    // to be passed into the residual function
    TimeIntAndModelArgs * args = new TimeIntAndModelArgs(dt, gamma, sigma, mu1, mu2, cp, gam, SOL,
                                                         CONS,  PRIMS, AUX, SOURCE);
    args->gID = gID;
    // Need to instantiate the correct device model
    switch (modType_t)
    {
      case ModelType::SRMHD:
        model_d = new SRMHD_D(args);
        break;
      case ModelType::SRRMHD:
        model_d = new SRRMHD_D(args);
        break;
      case ModelType::TFEMHD:
        model_d = new TFEMHD_D(args);
        break;
    }


    // First load initial guess (current value of cons)
    for (int i(0); i < Ncons; i++) SOL[i] = cons[i + lID * Ncons];

    // Rootfind
    if ((info = __cminpack_func__(hybrd1)(IMEX3Residual1Parallel, model_d, Ncons, SOL, &fvec[lID * Ncons], tol, &wa[lwa * lID], lwa)) != 1)
    {
      printf("IMEX failed stage 1 for gID %d: info %d\n", gID, info);
    }

    // #if MATCH_SERIAL==0
      for (int i(0); i < Nprims; i++) prims[i + lID * Nprims] = PRIMS[i];
      for (int i(0); i < Naux; i++) aux[i + lID * Naux] =  AUX[i];
      for (int i(0); i < Ncons; i++) source[i + lID * Ncons] = SOURCE[i];
    // #endif
      for (int i(0); i < Ncons; i++) sol[i + lID * Ncons] = SOL[i];

    // Clean up
    delete args;
    delete model_d;
  }
}


  //! Residual function to minimize for stage one of IMEX SSP2
  /*!
    Root of this function gives the values for U^(1).

    Parameters
    ----------
    p : pointer to BackwardsRK2 object
      The integrator object contains the argument object with the constar, primstar
      etc. arrays and the model object required for the single cell source term
      method.
    n : int
      Size of system
    x : pointer to double
      The array containing the guess
    fvec : pointer to double
      The array containing the residual as a result of the guess x
    iflag : int
      Error flag
  */
  __device__
  int IMEX3Residual1Parallel(void *p, int n, const double *x, double *fvec, int iflag)
  {
    // Cast void pointer
    Model_D * mod = (Model_D *)p;

    // First determine the prim and aux vars due to guess x
    mod->getPrimitiveVarsSingleCell((double *)x, mod->args->prims, mod->args->aux);
    // Determine the source contribution due to the guess x
    mod->sourceTermSingleCell((double *)x, mod->args->prims, mod->args->aux, mod->args->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - mod->args->cons[i] - mod->args->dt * mod->args->gam * mod->args->source[i];
      if (mod->args->source[i] != mod->args->source[i] || x[i] != x[i] || fvec[i] != fvec[i])
      {
        for (int j(0); j<n; j++) fvec[j] = 1e6;
        return 0;
      }
    }

    return 0;
  }

void SSP3::callStageTwo(double * cons, double * prims, double * aux, double * source, double * cons1, double * source1, double * flux1, double dt)
{
  int i, j, k, var;

  Data * d(this->data);
  //########################### STAGE TWO #############################//
  // First need to copy data to the device
  // A single cell requires all cons, prims and aux for the step. Rearrange so
  // we can copy data in contiguous way
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        #if MATCH_SERIAL
          for (var = 0; var < d->Ncons; var++)  args.cons_h [IDCons(var, i, j, k)  ] = cons[ID(var, i, j, k)];
          for (var = 0; var < d->Nprims; var++) args.prims_h[IDPrims(var, i, j, k) ] = prims[ID(var, i, j, k)];
        #endif
        for (var = 0; var < d->Naux; var++)   args.aux_h[IDAux(var, i, j, k)     ] = aux[ID(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  args.cons1_h[IDCons(var, i, j, k)  ] = cons1[ID(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  args.source1_h[IDCons(var, i, j, k)] = source1[ID(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  args.flux1_h[IDCons(var, i, j, k)]   = flux1[ID(var, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < d->Nstreams; i++) {

    // Which cell is at the left bound?
    int lcell(i * d->tpb * d->bpg);
    // Which cell is at the right bound?
    int rcell(lcell + d->tpb * d->bpg);
    if (rcell > d->Ncells) rcell = d->Ncells; // Dont overshoot
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));

    // Send stream's data
    #if MATCH_SERIAL
      gpuErrchk( cudaMemcpyAsync(args.cons_d[i], args.cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.prims_d[i], args.prims_h + lcell*d->Nprims, inMemsize*d->Nprims, cudaMemcpyHostToDevice, args.stream[i]) );
    #endif
    gpuErrchk( cudaMemcpyAsync(args.aux_d[i], args.aux_h + lcell*d->Naux, inMemsize*d->Naux, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.cons1_d[i], args.cons1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.source1_d[i], args.source1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.flux1_d[i], args.flux1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );

    int sharedMem((d->Ncons + d->Ncons + d->Ncons + d->Nprims + d->Naux) * sizeof(double) * d->tpb);

    // Call kernel and operate on data
    stageTwo <<< d->bpg, d->tpb, sharedMem, args.stream[i] >>>
            (args.sol_d[i], args.cons_d[i], args.prims_d[i], args.aux_d[i], args.source_d[i], args.cons1_d[i],
            args.source1_d[i], args.flux1_d[i], args.wa_d[i], args.fvec_d[i], dt, args.gam, tol, i, d->tpb * d->bpg,
            width, d->Ncons, d->Nprims, d->Naux, lwa,
            d->gamma, d->sigma, d->mu1, d->mu2, d->cp,
            model->modType_t);

    cudaStreamSynchronize(args.stream[i]);
    gpuErrchk( cudaPeekAtLastError() );
    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(args.sol_h + lcell*d->Ncons, args.sol_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_h + lcell*d->Nprims, args.prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_h + lcell*d->Naux, args.aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.source_h + lcell*d->Ncons, args.source_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );

  // Rearrange data back into arrays
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons; var++)  cons[ID(var, i, j, k)]    = args.sol_h[IDCons(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  source[ID(var, i, j, k)]  = args.source_h[IDCons(var, i, j, k)];
        for (var = 0; var < d->Nprims; var++) prims[ID(var, i, j, k)  ] = args.prims_h[IDPrims(var, i, j, k)];
        for (var = 0; var < d->Naux; var++)   aux[ID(var, i, j, k)]     = args.aux_h[IDAux(var, i, j, k)];
      }
    }
  }
}

__global__
void stageTwo(double * sol, double * cons, double * prims, double * aux, double * source,
              double * cons1, double * source1, double * flux1, double * wa, double * fvec,
              const double dt, const double gam, const double tol, const int stream,
              const int origWidth, const int streamWidth, const int Ncons, const int Nprims, const int Naux, const int lwa,
              const double gamma, const double sigma, const double mu1, const double mu2, const double cp,
              ModelType modType_t)
{
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  const int gID(lID + stream * origWidth);        //!< global index (in domain)
  int info;
  extern __shared__ double sharedMem[];
  double * CONS = &sharedMem[tID * (Ncons + Ncons + Ncons +  Nprims + Naux)];
  double * PRIMS = &CONS[Ncons];
  double * AUX = &PRIMS[Nprims];
  double * SOURCE = &AUX[Naux];
  double * SOL = &SOURCE[Ncons];

  if (lID < streamWidth)
  {

    for (int i(0); i < Ncons; i++) CONS[i] = cons[i + lID * Ncons];
    for (int i(0); i < Nprims; i++) PRIMS[i] = prims[i + lID * Nprims];
    for (int i(0); i < Naux; i++) AUX[i] = aux[i + lID * Naux];
    for (int i(0); i < Ncons; i++) SOURCE[i] = source[i + lID * Ncons];

    Model_D * model_d;

    // Store pointers to devuce arrays in the structure
    // to be passed into the residual function
    TimeIntAndModelArgs * args = new TimeIntAndModelArgs(dt, gamma, sigma, mu1, mu2, cp, gam, SOL,
                                                         CONS, PRIMS,
                                                         AUX, SOURCE,
                                                         &cons1[lID * Ncons], &source1[lID * Ncons],
                                                         &flux1[lID * Ncons]);
    args->gID = gID;
    // Need to instantiate the correct device model
    switch (modType_t)
    {
      case ModelType::SRMHD:
        model_d = new SRMHD_D(args);
        break;
      case ModelType::SRRMHD:
        model_d = new SRRMHD_D(args);
        break;
      case ModelType::TFEMHD:
        model_d = new TFEMHD_D(args);
        break;
    }

    // First load initial guess (current value of cons)
    for (int i(0); i < Ncons; i++) SOL[i] = CONS[i];

    if ((info = __cminpack_func__(hybrd1)(IMEX3Residual2Parallel, model_d, Ncons, SOL, &fvec[lID * Ncons], tol,  &wa[lwa * lID], lwa)) != 1)
    {
      printf("IMEX failed stage 2 for gID %d: info %d\n", gID, info);
    }

    // #if MATCH_SERIAL==0
      for (int i(0); i < Nprims; i++) prims[ i + lID * Nprims] = PRIMS[i];
      for (int i(0); i < Naux; i++  ) aux[   i + lID * Naux]   = AUX[i];
      for (int i(0); i < Ncons; i++ ) source[i + lID * Ncons]  = SOURCE[i];
    // #endif
    for (int i(0); i < Ncons; i++ ) sol[   i + lID * Ncons ] = SOL[i];

    // Clean up
    delete args;
    delete model_d;
  }
}



__device__
int IMEX3Residual2Parallel(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  Model_D * mod = (Model_D *)p;

  // First determine the prim and aux vars due to guess x
  mod->getPrimitiveVarsSingleCell((double *)x, mod->args->prims, mod->args->aux);
  // Determine the source contribution due to the guess x
  mod->sourceTermSingleCell((double *)x, mod->args->prims, mod->args->aux, mod->args->source);

  // Set residual
  for (int i(0); i < n; i++)
  {
    fvec[i] = x[i] - mod->args->cons[i] + mod->args->dt * (mod->args->flux1[i] - (1 - 2*mod->args->gam) * mod->args->source1[i] - mod->args->gam * mod->args->source[i]);
    if (mod->args->source[i] != mod->args->source[i] || x[i] != x[i] || fvec[i] != fvec[i])
    {
      for (int j(0); j<n; j++) fvec[j] = 1e6;
      return 0;
    }
  }

return 0;
}



void SSP3::callStageThree(double * cons, double * prims, double * aux, double * source, double * source1, double * flux1, double * cons2, double * flux2, double dt)
{
  int i, j, k, var;

  Data * d(this->data);
  //########################### STAGE THREE #############################//
  // First need to copy data to the device
  // A single cell requires all cons, prims and aux for the step. Rearrange so
  // we can copy data in contiguous way
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        #if MATCH_SERIAL
          for (var = 0; var < d->Ncons; var++)  args.cons_h [IDCons(var, i, j, k)  ] = cons[ID(var, i, j, k)];
          for (var = 0; var < d->Nprims; var++) args.prims_h[IDPrims(var, i, j, k) ] = prims[ID(var, i, j, k)];
          for (var = 0; var < d->Ncons; var++)  args.source1_h[IDCons(var, i, j, k)] = source1[ID(var, i, j, k)];
          for (var = 0; var < d->Ncons; var++)  args.flux1_h[IDCons(var, i, j, k)]   = flux1[ID(var, i, j, k)];
        #endif
        for (var = 0; var < d->Naux; var++)   args.aux_h[IDAux(var, i, j, k)     ] = aux[ID(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  args.cons1_h[IDCons(var, i, j, k)  ] = cons2[ID(var, i, j, k)]; // NOTE: Storing cons2 in host cons1 array to save on memory
        for (var = 0; var < d->Ncons; var++)  args.flux2_h[IDCons(var, i, j, k)]   = flux2[ID(var, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < d->Nstreams; i++) {

    // Which cell is at the left bound?
    int lcell(i * d->tpb * d->bpg);
    // Which cell is at the right bound?
    int rcell(lcell + d->tpb * d->bpg);
    if (rcell > d->Ncells) rcell = d->Ncells; // Dont overshoot
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));

    // Send stream's data
    #if MATCH_SERIAL
      gpuErrchk( cudaMemcpyAsync(args.cons_d[i], args.cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.prims_d[i], args.prims_h + lcell*d->Nprims, inMemsize*d->Nprims, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.source1_d[i], args.source1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.flux1_d[i], args.flux1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    #endif
    gpuErrchk( cudaMemcpyAsync(args.aux_d[i], args.aux_h + lcell*d->Naux, inMemsize*d->Naux, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.cons1_d[i], args.cons1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.flux2_d[i], args.flux2_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );

    int sharedMem((d->Ncons + d->Ncons + d->Ncons + d->Nprims + d->Naux) * sizeof(double) * d->tpb);

    // Call kernel and operate on data
    stageThree <<< d->bpg, d->tpb, sharedMem, args.stream[i] >>>
            (args.sol_d[i], args.cons_d[i], args.prims_d[i], args.aux_d[i], args.source_d[i], args.cons1_d[i],
            args.source1_d[i], args.flux1_d[i], args.flux2_d[i], args.wa_d[i], args.fvec_d[i], dt, args.gam, tol, i, d->tpb * d->bpg,
            width, d->Ncons, d->Nprims, d->Naux, lwa,
            d->gamma, d->sigma, d->mu1, d->mu2, d->cp,
            model->modType_t);

    cudaStreamSynchronize(args.stream[i]);
    gpuErrchk( cudaPeekAtLastError() );
    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(args.sol_h + lcell*d->Ncons, args.sol_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_h + lcell*d->Nprims, args.prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_h + lcell*d->Naux, args.aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.source_h + lcell*d->Ncons, args.source_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );

  // Rearrange data back into arrays
  #pragma omp parallel for
  for (i = 0; i < d->Nx; i++) {
    for (j = 0; j < d->Ny; j++) {
      for (k = 0; k < d->Nz; k++) {
        for (var = 0; var < d->Ncons; var++)  cons[ID(var, i, j, k)]    = args.sol_h[IDCons(var, i, j, k)];
        for (var = 0; var < d->Ncons; var++)  source[ID(var, i, j, k)]  = args.source_h[IDCons(var, i, j, k)];
        for (var = 0; var < d->Nprims; var++) prims[ID(var, i, j, k)  ] = args.prims_h[IDPrims(var, i, j, k)];
        for (var = 0; var < d->Naux; var++)   aux[ID(var, i, j, k)]     = args.aux_h[IDAux(var, i, j, k)];
      }
    }
  }
}

__global__
void stageThree(double * sol, double * cons, double * prims, double * aux, double * source,
              double * cons1, double * source1, double * flux1, double * flux2, double * wa, double * fvec,
              const double dt, const double gam, const double tol, const int stream,
              const int origWidth, const int streamWidth, const int Ncons, const int Nprims, const int Naux, const int lwa,
              const double gamma, const double sigma, const double mu1, const double mu2, const double cp,
              ModelType modType_t)
{
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  const int gID(lID + stream * origWidth);        //!< global index (in domain)
  int info;
  extern __shared__ double sharedMem[];
  double * CONS = &sharedMem[tID * (Ncons + Ncons + Ncons +  Nprims + Naux)];
  double * PRIMS = &CONS[Ncons];
  double * AUX = &PRIMS[Nprims];
  double * SOURCE = &AUX[Naux];
  double * SOL = &SOURCE[Ncons];

  if (lID < streamWidth)
  {

    for (int i(0); i < Ncons; i++) CONS[i] = cons[i + lID * Ncons];
    for (int i(0); i < Nprims; i++) PRIMS[i] = prims[i + lID * Nprims];
    for (int i(0); i < Naux; i++) AUX[i] = aux[i + lID * Naux];
    for (int i(0); i < Ncons; i++) SOURCE[i] = source[i + lID * Ncons];
    Model_D * model_d;

    // Store pointers to devuce arrays in the structure
    // to be passed into the residual function
    TimeIntAndModelArgs * args = new TimeIntAndModelArgs(dt, gamma, sigma, mu1, mu2, cp, gam, SOL,
                                                         CONS, PRIMS,
                                                         AUX, SOURCE,
                                                         &cons1[lID * Ncons], &source1[lID * Ncons],
                                                         &flux1[lID * Ncons], &flux2[lID * Ncons]);
    args->gID = gID;
    // Need to instantiate the correct device model
    switch (modType_t)
    {
      case ModelType::SRMHD:
        model_d = new SRMHD_D(args);
        break;
      case ModelType::SRRMHD:
        model_d = new SRRMHD_D(args);
        break;
      case ModelType::TFEMHD:
        model_d = new TFEMHD_D(args);
        break;
    }


    // First load initial guess (current value of cons)
    for (int i(0); i < Ncons; i++) SOL[i] = args->cons1[i] + dt*(args->flux1[i]+ args->flux2[i])/4.0;

    if ((info = __cminpack_func__(hybrd1)(IMEX3Residual3Parallel, model_d, Ncons, SOL, &fvec[lID * Ncons], tol,  &wa[lwa * lID], lwa)) != 1)
    {
      printf("IMEX failed stage 3 for gID %d: info %d\n", gID, info);
    }

    // #if MATCH_SERIAL==0
      for (int i(0); i < Nprims; i++) prims[ i + lID * Nprims] = PRIMS[i];
      for (int i(0); i < Naux; i++  ) aux[   i + lID * Naux]   = AUX[i];
      for (int i(0); i < Ncons; i++ ) source[i + lID * Ncons]  = SOURCE[i];
    // #endif
    for (int i(0); i < Ncons; i++ ) sol[   i + lID * Ncons ] = SOL[i];

    // Clean up
    delete args;
    delete model_d;
  }
}


__device__
int IMEX3Residual3Parallel(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  Model_D * mod = (Model_D *)p;

  // First determine the prim and aux vars due to guess x
  mod->getPrimitiveVarsSingleCell((double *)x, mod->args->prims, mod->args->aux);
  // Determine the source contribution due to the guess x
  mod->sourceTermSingleCell((double *)x, mod->args->prims, mod->args->aux, mod->args->source);

  // Set residual
  for (int i(0); i < n; i++)
  {
    fvec[i] = x[i] - mod->args->cons[i] + mod->args->dt * (mod->args->flux1[i] + mod->args->flux2[i]) / 4.0 - mod->args->dt * (0.2071067811865475 * mod->args->source1[i] + mod->args->gam * mod->args->source[i]);

    if (mod->args->source[i] != mod->args->source[i] || x[i] != x[i] || fvec[i] != fvec[i])
    {
      for (int j(0); j<n; j++) fvec[j] = 1e6;
      return 0;
    }
  }

return 0;
}

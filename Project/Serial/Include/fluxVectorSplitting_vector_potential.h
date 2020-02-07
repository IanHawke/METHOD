#ifndef FLUXVECTORSPLITTING_VECTOR_POTENTIAL_H
#define FLUXVECTORSPLITTING_VECTOR_POTENTIAL_H

#include "flux.h"
#include "weno.h"

//! <b> Flux vector splitting method </b>
/*!
  @par
    This calculates the RHS for the MHD+Vector Potential. Far too specialized,
    but there's too many tricks. See FVS files and the Spritz paper.
  @par
*/
class FVS_vector_potential : public FluxMethod
{
  public:

    //! Constructor
    /*!
        Calls the base class constructor to store pointers to Data and Model
      classes.

      @param[in] *data pointer to Data class
      @param[in] *model pointer to Model class
    */
    FVS_vector_potential(Data * data, Model * model) : FluxMethod(data, model) { }

    //! Flux reconstruction
    /*!
        Reconstructs the fluxes at the center of the cells to the faces upwind
      and downwind and computes the difference, giving an approximation of the
      net flux (in the specified direction) at the cell faces. Method uses a
      second order WENO reconstruction.

      @param[in] *cons pointer to conserved vector
      @param[in] *prims pointer to primitive vector
      @param[in] *aux pointer to auxiliary vector
      @param[in] *f pointer to a flux work array to store the initial flux vector
      @param[out] *frecon pointer to the array containing the reconstructed values
      of the fluxes at the cell faces
      @param[in] dir the direction in which to determine flux reconstruction with
      (0, 1, 2) = (x, y, z)
      @param[in] vars size of the vector to reconstruct (saves time when using subgrid models).
      Default values is -1, which autos to Ncons.
      @note This is an approximation of the net flux at the cell faces, not
      the cell centers, and so the final approximation of the flux through a cell
      required differencing the values obtained here at either end of the cell.
      This is performed in F
      @sa F
    */
    void fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir, int vars=-1);

    //! Numerical flux function
    /*!
        For a given state described by cons prims and aux arrays, determines an
      approximation of the net flux of the conserved quantities through all cells
      by taking difference of the reconstructed values at the cell faces.

      @param[in] *cons pointer to conserved vector
      @param[in] *prims pointer to primitive vector
      @param[in] *aux pointer to auxiliary vector
      @param[in] *f pointer to a flux work array to store the initial flux vector
      @param[out] *fnet pointer to the array containing the net flux through every cell
      @sa fluxReconstruction
    */
    void F(double * cons, double * prims, double * aux, double * f, double * fnet);

};

#endif

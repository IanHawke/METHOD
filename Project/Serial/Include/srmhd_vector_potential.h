#ifndef SRMHD_VECTOR_POTENTIAL_H
#define SRMHD_VECTOR_POTENTIAL_H

#include "srmhd.h"
#include "model.h"


/*
This is the human readable description of this models variables.

  SRMHD has nine conserved variables:
    D, Sx, Sy, Sz, tau, Ax, Ay, Az, phi
  Nine primitive variables:
    rho, vx, vy, vz, p, Ax, Ay, Az, phi
  Twenty auxiliary variables:
    h, W, e, c, b0, bx, by, bz, bsq, vsq, BS, Bsq, Ssq, Bx, By, Bz, Ex, Ey, Ez, Phi_staggered

This form uses a staggered grid. A* are defined at cell edges; everything else
at cell centres.
*/



//! <b>Special Relativistic MagnetHydroDynamics </b>
/*!
  @par
    The single fluid, special relativistic, ideal limit of the MHD equations.
  Ideal fluid, so resistivity does not play a part, hence no
  electric field evolution. <br>

  @note
  Model has nine conserved variables: <br>
   \f$\ \ \ D\f$, \f$S_x\f$, \f$S_y\f$, \f$S_z\f$, \f$\tau\f$, \f$B_x\f$, \f$B_y\f$, \f$B_z\f$, \f$\phi\f$ <br>
  Eight primitive variables: <br>
    \f$\ \ \ \rho\f$, \f$v_x\f$, \f$v_y\f$, \f$v_z\f$, \f$p\f$, \f$B_x\f$, \f$B_y\f$, \f$B_z\f$ <br>
  Thirteen auxiliary variables:<br>
    \f$\ \ \ h\f$, \f$W\f$, \f$e\f$, \f$c\f$, \f$b_0\f$, \f$b_x\f$, \f$b_y\f$, \f$b_z\f$, \f$b^2\f$, \f$v^2\f$, \f$B \cdot S\f$, \f$B^2\f$, \f$S^2\f$ <br>

  @par
    The SRMHD fluid equations are derived from the consideration of the conservation
  of the rest-mass and stress-energy tensor, and the conservation of the Maxwell
  dual tensor for a perfect magneto-fluid. That is, starting from
  \f{align}{
    \partial_\mu N^\mu &= 0 \\
    \partial_mu T^{\mu \nu}&= 0 \\
    \partial_mu {^*}F^{\mu \nu} &= 0
  \f}
  with \f$N^\mu = \rho u^mu\f$ the rest-mass density, \f$T^{\mu \nu} = \rho h^* u^\mu u^\nu + \eta^{\mu \nu} p^* - b^mu b^nu\f$
  as the stress-energy tensor for a perfect magneto-fluid, and the Maxwell tensor
  given by \f$ {^*}F^{\mu \nu} = u^\mu b^\nu - u^\nu b^\mu\f$.
  @par
    In addition, \f$ u^\mu, h^*=1+e+p/\rho+b^2/\rho, \rho, \eta, b^\mu \text{ and } p*=p+b^2/2\f$
  are the fluid four-velocity, specific enthalpy including the magnetic contribution,
  mass-energy density, mostly positive flat space-time metric, four-vector magnetic
  field in the fluid rest-frame and the total pressure including the magnetic pressure.
  @par
    Following this through, we arrive at the conservative form of the special relativistic,
  ideal limit, single fluid equations of motion:
  \f{align}{
    \partial_t
    \begin{pmatrix}
      D \\ S^j \\ \tau \\ B^k \\ \phi
    \end{pmatrix}
    + \partial_i
    \begin{pmatrix}
      Dv^i \\ S^j v^i + p^* \delta^{ij} - b^j B^i / W \\ \tau v^i + p^* v^i - b^0 B^i/ W \\ v^i B^k - v^k B^i + \delta^{ij}\phi \\ B^i
    \end{pmatrix}
    =
    \begin{pmatrix}
      0 \\ 0 \\ 0 \\ 0 \\ -\phi / c_p^2
    \end{pmatrix}.
  \f}
  @par
    Here, we sum over the \f$i\f$ coordinate directions, and note the following
  relations:
  \f{align}{
    D &= \rho W \\
    S^j &= \rho h^* W^2 v^j - b^0 b^j \\
    \tau &= \rho h^* W^2 - p^* - (b^0)^2 - D \\
    u^\mu &= W(c, v^i) \\
    W &= 1 / \sqrt(1 - v_i v^i) \\
    b^0 &= W B_i v^i \\
    b^i &= B^i / W + b^0 v^i \\
    b^2 &= B_i B^i / W^2 + (B_i v^i)^2 \\
    c_p &= const.
  \f}
  @par
    We have also included the additional scalar field \f$\phi\f$ such that any
  errors in the evolution of the magnetic fields that break the divergence constraint
  set by Maxwell's equations, namely \f$\nabla \cdot B = 0\f$, are driven to zero,
  on a timescale set by the constant parameter \f$c_p\f$.
  See Dedner et al. 2002.


  @sa Model
*/
class SRMHD_Vector_Potential : public Model
{

  public:

    int smartGuesses;     //!< Number of smart guess required

    double * solution;    //!< Pointer to array to hold solution of C2P for every cell. Size is 2*Nx*Ny*Nz


    SRMHD_Vector_Potential();     //!< Default constructor

    //! Parameterized constructor
    /*!
      @par
        Stores a pointer to the Data class for reference in its methods

      @param *data Pointer to Data class containing global simulation data
    */
    SRMHD_Vector_Potential(Data * data);

    ~SRMHD_Vector_Potential();     //!< Destructor


    //! Single cell source term contribution
    /*!
      @par
        Models that can posess a stiff source term and hence (semi-)implicit time
      integrators will require a source contribution (and cons2prims method) that
      applies to a single cell.
      @par
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}\f$
      @param *source pointer to source vector work array. Size is \f$N_{cons}\f$
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::sourceTermSingleCell
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
      @par
        Non-zero flux for cons[8], phi, as a result of implementing divergence
      cleaning. For details see Muddle.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param *source pointer to source vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @sa Model::sourceTerm
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell cons2prims conversion
    /*!
      @par
        For the same reason as outlined in sourceTermSingleCell, some models will
      require a single celled primitive conversion method.
      @par
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}\f$
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);

    //! Spectral decomposition
    /*!
      @par
        Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
      Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
      an N=2 rootfind using cminpack library.
      @par
        Initial inputs will be the current values of the conserved vector and the
      <b> old </b> values for the prims and aux vectors.
      Output will be the current values of cons, prims and aux.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @sa Model::getPrimitiveVars
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
      @par
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Anton 2010.
      Riemann Solver`

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
      @par
        Method determines the value for the flux vector.
      @par
        For the form of the fluxes see Anton 2010 with the inclusion of divergence
      cleaning from Dedner et al. 2002.
      interfaces, John Muddle.

      @note We are assuming that all primitive and auxiliary variables are up-to-date
      at the time of this function execution.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param *f pointer to flux vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param dir direction in which to generate flux vector. \f$(x, y, z) = (0, 1, 2)\f$
      @sa Model::fluxVector
    */
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);

    //! Finalise the simulation variables
    /*!
      @par
        Mostly, this probably wont be needed, but if there is any final steps to finish
      off a timestep, this can be done here.
    */
    void finalise(double *cons, double *prims, double *aux) { };

    //! A (vector potential) to B (magnetic field)
    /*!
      @par
        Uses corner differencing
    */
    void AToB(double *c_or_p, double *aux);

};


#endif

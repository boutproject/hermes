/*
    Copyright B.Dudson, J.Leddy, University of York, September 2016
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes.

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.

*/

class Hermes;

#ifndef __HERMES_H__
#define __HERMES_H__

#include <bout/physicsmodel.hxx>

#include <bout/constants.hxx>
#include <bout/invert/laplace3d.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/invert/laplacexz.hxx>
#include <invert_laplace.hxx>

#include "radiation.hxx"

class Hermes : public PhysicsModel {
public:
  virtual ~Hermes() {}

protected:
  int init(bool restarting);
  int rhs(BoutReal t);

  int precon(BoutReal t, BoutReal gamma, BoutReal delta);

private:
  // Equilibrium current
  Field2D Jpar0;

  // Evolving variables
  Field3D Ne, Pe; // Electron density and pressure
  Field3D VePsi;  // Combination of Ve and psi
  Field3D Vort;   // Vorticity
  Field3D NVi;    // Parallel momentum

  FieldGroup EvolvingVars;

  // Auxilliary variables
  Field3D Te;           // Electron temperature
  BoutReal Ti;          // Ion temperature
  Field3D Ve, Vi, Jpar; // Electron and ion parallel velocities
  Field3D psi;          // Electromagnetic potential (-A_||)
  Field3D phi;          // Electrostatic potential

  // Limited variables
  Field3D Telim;

  // Collisional damping
  Field3D nu, kappa_epar, Dn;

  BoutReal flux_limit_alpha;  // Flux limiter. < 0 disables
  BoutReal kappa_limit_alpha; // Heat flux limiter from SOLPS
  BoutReal eta_limit_alpha;   // Momentum flux limiter from SOLPS

  // Neutral gas model
  int neutral_model; // Neutrals 0 = None
  bool neutral_friction;
  BoutReal gamma_ratio;        // Ratio of specific heats
  BoutReal neutral_viscosity;  // Neutral gas viscosity
  BoutReal neutral_bulk;       // Neutral gas bulk viscosity
  BoutReal neutral_conduction; // Neutral gas thermal conduction
  BoutReal frecycle;           // Recycling fraction
  BoutReal Eionize;            // Ionisation energy loss
  bool excitation;             // Include electron impact excitation?

  Field3D Nn;  // Neutral gas density (evolving)
  Field3D Pn;  // Neutral gas pressure (evolving)
  Field3D NVn; // Neutral gas momentum (evolving)
  Field3D Tn;  // Neutral gas temperature

  Field2D Nn2D;  // Neutral gas density (evolving)
  Field2D Pn2D;  // Neutral gas pressure (evolving)
  Vector2D Vn2D; // Neutral gas velocity
  Field2D Tn2D;
  Field2D DivV2D; // Divergence of gas velocity

  // Transformation to cylindrical coordinates

  // Grad x = Txr * Grad R + Txz * Grad Z
  // Grad y = Tyr * Grad R + Tyz * Grad Z
  Field2D Txr, Txz;
  Field2D Tyr, Tyz;

  // Grad R = Urx * Grad x + Ury * Grad y
  // Grad Z = Uzx * Grad x + Uzy * Grad y
  Field2D Urx, Ury;
  Field2D Uzx, Uzy;

  // Atomic rates and transfer channels

  Field3D Wi;                    // Power transfer from electrons to ions
  UpdatedRadiatedPower hydrogen; // Atomic rates (H.Willett)
  Field3D Dnn;                   // Neutral gas diffusion rate
  Field3D S;                     // Plasma particle sink, neutral source
  Field3D Q;                     // Power transfer from plasma to neutrals
  Field3D Qe, Qi;                // Power transfer from electrons and ions
  Field3D F;                     // Ion-neutral friction
  Field3D Rp;                    // Radiation from the plasma
  Field3D Rn;                    // Radiation from the neutrals
  BoutReal Lmax;                 // Maximum neutral mean free path [m]
  Field3D lambda_int, lambda;    // for mean free path calculation
  Field3D Nn0;                   // scaling for neutral density approximation

  void neutral_rates(const Field3D &Ne, const Field3D &Te,
                     const Field3D &Vi, // Plasma quantities
                     const Field3D &Nn, const Field3D &Tn,
                     const Field3D &Vnpar, // Neutral gas
                     Field3D &S, Field3D &F, Field3D &Q,
                     Field3D &R, // Transfer rates
                     Field3D &Riz, Field3D &Rrc, Field3D &Rcx); // Rates

  // Impurity radiation
  BoutReal carbon_fraction;
  Field3D Rzrad;             // Radiated power
  RadiatedPower *carbon_rad; // Carbon cooling curve

  // Switches
  bool evolve_plasma; // Should plasma be evolved?

  bool electromagnetic; // Include magnetic potential psi
  bool FiniteElMass;    // Finite Electron Mass

  bool j_diamag; // Diamagnetic current: Vort <-> Pe
  bool j_par;    // Parallel current:    Vort <-> Psi
  bool parallel_flow;
  bool pe_par;             // Parallel pressure gradient: Pe <-> Psi
  bool resistivity;        // Resistivity: Psi -> Pe
  bool thermal_force;      // Force due to temperature gradients
  bool electron_viscosity; // Electron parallel viscosity
  bool electron_neutral;   // Include electron-neutral collisions
  bool poloidal_flows; // Include y derivatives in diamagnetic and ExB drifts
  bool thermal_flux;   // Include parallel and perpendicular energy flux from Te
                       // gradients
  bool thermal_conduction; // Braginskii electron heat conduction

  bool classical_diffusion; // Collisional diffusion, including viscosity

  // Anomalous perpendicular diffusion coefficients
  BoutReal anomalous_D;   // Density diffusion
  BoutReal anomalous_chi; // Electron thermal diffusion
  BoutReal anomalous_nu;  // Momentum diffusion (kinematic viscosity)

  bool ion_velocity; // Include Vi terms

  bool rotation;          ///< Include toroidal rotation?
  BoutReal rotation_rate; ///< Rotation rate (if rotation=true). Input in Hz
  Vector2D Omega_vec;     ///< Rotation vector (in vertical direction)
  Vector2D bxGradR;       ///< b cross Grad(R) for centrifugal drift

  bool phi3d; // Use a 3D solver for phi

  bool staggered; // Use staggered differencing along B

  bool boussinesq; // Use a fixed density (Nnorm) in the vorticity equation

  bool sinks;          // Sink terms for running 2D drift-plane simulations
  bool sheath_closure; // Sheath closure sink on vorticity (if sinks = true)
  bool drift_wave;     // Drift-wave closure (if sinks=true)

  bool radial_buffers; // Radial buffer regions
  int radial_inner_width; // Number of points in the inner radial buffer
  int radial_outer_width; // Number of points in the outer radial buffer
  BoutReal radial_buffer_D; // Diffusion in buffer region

  Field2D sink_invlpar; // Parallel inverse connection length (1/L_{||}) for
                        // sink terms
  Field2D alpha_dw;

  // Sheath heat transmission factor
  int sheath_model;       // Sets boundary condition model
  BoutReal sheath_gamma;  // Heat transmission
  BoutReal neutral_gamma; // Heat transmission for neutrals
  BoutReal neutral_vwall; // Scale velocity at the wall
  bool sheath_yup, sheath_ydown;

  Field2D wall_flux;  // Particle flux to wall (diagnostic)
  Field2D wall_power; // Power flux to wall (diagnostic)

  // Outflowing boundaries for neutrals
  bool outflow_ydown; // Allow neutral outflows?

  // Fix density in SOL
  bool sol_fix_profiles;
  FieldGenerator *sol_ne, *sol_te; // Generating functions

  // Output switches for additional information
  bool verbose;    // Outputs additional fields, mainly for debugging
  bool output_ddt; // Output time derivatives

  // Parallel diffusion

  BoutReal vort_mupar;
  BoutReal ne_pardiff; // Parallel diffusion of density
  BoutReal pe_kappa;   // Constant heat flux coefficient (normalised)

  // Collision times for electrons and ions
  Field3D tau_e, tau_i;

  // Perpendicular diffusion
  BoutReal Dne, Dvort, Dve, Dvi, Dte;

  BoutReal ion_neutral; // Ion-neutral collision rate

  BoutReal numdiff, hyper, hyperpar; ///< Numerical dissipation
  BoutReal neut_numdiff;             // Numdiff for neutral gas
  BoutReal ExBdiff;
  bool ExBpar;      // Include parallel terms in ExBdiff
  BoutReal ADpar;   // Added Dissipation method in the parallel direction
  bool ADpar_phine; // Include 1/Ne factor in phi ADpar term
  bool ADpar_bndry; // Boundaries included in ADpar?
  int low_pass_z;   // Fourier filter in Z
  BoutReal z_hyper_viscos, x_hyper_viscos, y_hyper_viscos; // 4th-order derivatives
  bool low_n_diffuse; // Diffusion in parallel direction at low density
  bool low_n_diffuse_perp; // Diffusion in perpendicular direction at low density
  BoutReal ne_hyper_z, pe_hyper_z; // Hyper-diffusion

  bool vepsi_dissipation; // Dissipation term in VePsi equation
  
  BoutReal resistivity_multiply; ///< Factor in front of nu

  // Sources and profiles

  bool ramp_mesh; // Use Ne,Pe in the grid file for starting ramp target
  BoutReal ramp_timescale;    // Length of time for the initial ramp
  Field2D NeTarget, TeTarget, PeTarget; // For adaptive sources

  bool adapt_source;           // Use a PI controller to feedback profiles
  bool adapt_fix_form;         // Fix the form of the source
  bool core_sources;           // Sources only in the core
  bool energy_source;          // Add the same amount of energy to each particle
  bool temperature_feedback;   // Feedback power source on temperature rather than pressure
  BoutReal source_p, source_i; // Proportional-Integral controller
  Field2D Sn, Spe;             // Sources in density and Pe
  BoutReal total_Sn, total_Spe; // Sum over all cells
  Field3D NeSource, PeSource;   // These are the actual source added

  bool source_vary_g11; // Multiply source by g11
  Field2D g11norm;

  BoutReal density_error_lasttime, density_error_last, density_error_integral;
  BoutReal pe_error_lasttime, pe_error_last, pe_error_integral;

  bool density_inflow; // Does incoming density have momentum?

  bool pe_bndry_flux;   // Allow flux of pe through radial boundaries
  bool ne_bndry_flux;   // Allow flux of ne through radial boundaries
  bool vort_bndry_flux; // Allow flux of vorticity through radial boundaries

  // Stress tensor components
  BoutReal tau_e0, tau_i0;

  // Diffusion parameters

  // Normalisation parameters
  BoutReal Tnorm, Nnorm, Bnorm;
  BoutReal AA, Cs0, rho_s0, Omega_ci, R0;
  BoutReal mi_me, beta_e;
  Field2D hthe;

  // Group of fields for communication
  FieldGroup evars; // Evolving variables

  // Curvature, Grad-B drift
  Vector3D Curlb_B; // Curl(b/B)

  // Numerical dissipation operators
  const Field3D D(const Field3D &f, BoutReal d); // Diffusion operator
  const Field3D H(const Field3D &f, BoutReal h); // Hyper-diffusion operator
  const Field3D H(const Field3D &f, BoutReal h,
                  const Field3D &mask); // Hyper-diffusion operator

  // Perturbed parallel gradient operators
  const Field3D Grad_parP(const Field3D &f);
  const Field3D Grad_parP_CtoL(const Field3D &f);
  const Field3D Grad_parP_LtoC(const Field3D &f);
  const Field3D Div_parP_LtoC(const Field3D &f, const Field3D &v);

  // Utility functions
  const Field2D CumSumY2D(const Field2D &f, bool reverse);
  const Field3D CumSumY3D(const Field3D &f, bool reverse);
  const BoutReal bcast_lasty(const BoutReal f);

  // Electromagnetic solver for finite electron mass case
  bool split_n0_psi; // Split the n=0 component of Apar (psi)?
  // Laplacian *aparSolver;
  LaplaceXZ *aparSolver;
  LaplaceXY *aparXY; // Solves n=0 component
  Field2D psi2D;     // Axisymmetric Psi

  // Solvers for the electrostatic potential

  bool split_n0;        // Split solve into n=0 and n~=0?
  LaplaceXY *laplacexy; // Laplacian solver in X-Y (n=0)
  Field2D phi2D;        // Axisymmetric phi

  bool newXZsolver;
  Laplacian *phiSolver; // Old Laplacian in X-Z
  LaplaceXZ *newSolver; // New Laplacian in X-Z
};

/// Fundamental constants

const BoutReal e0 = 8.854e-12;      // Permittivity of free space
const BoutReal mu0 = 4.e-7 * PI;    // Permeability of free space
const BoutReal qe = 1.602e-19;      // Electron charge
const BoutReal Me = 9.109e-31;      // Electron mass
const BoutReal Mp = 1.67262158e-27; // Proton mass

#endif // __HERMES_H__

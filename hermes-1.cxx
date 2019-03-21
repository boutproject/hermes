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
#include "hermes-1.hxx"

#include <derivs.hxx>
#include <field_factory.hxx>
#include <initialprofiles.hxx>

#include <invert_parderiv.hxx>

#include "div_ops.hxx"
#include "loadmetric.hxx"

#include <bout/constants.hxx>

#include <bout/assert.hxx>

BoutReal floor(const BoutReal &var, const BoutReal &f) {
  if (var < f)
    return f;
  return var;
}

const Field3D ceil(const Field3D &var, BoutReal f) {
  Field3D result;
  result.allocate();
  
  for (int jx=0;jx<mesh->ngx;jx++)
    for (int jy=0;jy<mesh->ngy;jy++)
      for (int jz=0;jz<mesh->ngz;jz++) {
        if (var(jx, jy, jz) < f) {
          result(jx, jy, jz) = var(jx,jy,jz);
        } else {
          result(jx, jy, jz) = f;
        }
      }
  return result;
}

int Hermes::init(bool restarting) {
  Options *opt = Options::getRoot();

  // Switches in model section
  Options *optsc = opt->getSection("Hermes");

  OPTION(optsc, evolve_plasma, true);

  OPTION(optsc, electromagnetic, true);
  OPTION(optsc, FiniteElMass, true);

  OPTION(optsc, j_diamag, true);
  OPTION(optsc, j_par, true);
  OPTION(optsc, parallel_flow, true);
  OPTION(optsc, pe_par, true);
  OPTION(optsc, resistivity, true);
  OPTION(optsc, thermal_flux, true);
  OPTION(optsc, thermal_force, true);
  OPTION(optsc, electron_viscosity, true);
  OPTION(optsc, electron_neutral, true);

  OPTION(optsc, poloidal_flows, false);
  OPTION(optsc, ion_velocity, true);

  OPTION(optsc, rotation, false);
  OPTION(optsc, rotation_rate, 0.0); // Input in Hz

  OPTION(optsc, thermal_conduction, true);

  OPTION(optsc, neutral_friction, false);
  OPTION(optsc, frecycle, 0.9);
  OPTION(optsc, Eionize, 30); // Energy loss per ionisation [eV]
  OPTION(optsc, Lmax, 1.0);

  OPTION(optsc, phi3d, false);

  OPTION(optsc, pe_kappa, -1.0);

  OPTION(optsc, ne_bndry_flux, true);
  OPTION(optsc, pe_bndry_flux, true);
  OPTION(optsc, vort_bndry_flux, false);

  OPTION(optsc, ramp_mesh, true);
  OPTION(optsc, ramp_timescale, 1e4);

  OPTION(optsc, energy_source, false);

  OPTION(optsc, ion_neutral, 0.0);

  OPTION(optsc, staggered, false);

  OPTION(optsc, boussinesq, false);

  OPTION(optsc, sinks, false);
  OPTION(optsc, sheath_closure, true);
  OPTION(optsc, drift_wave, false);

  OPTION(optsc, classical_diffusion, false);
  OPTION(optsc, anomalous_D, -1);
  OPTION(optsc, anomalous_chi, -1);
  OPTION(optsc, anomalous_nu, -1);

  OPTION(optsc, flux_limit_alpha, -1);
  OPTION(optsc, kappa_limit_alpha, -1);
  OPTION(optsc, eta_limit_alpha, -1);

  // Numerical dissipation terms
  OPTION(optsc, numdiff, -1.0);
  OPTION(optsc, neut_numdiff, numdiff);
  OPTION(optsc, hyper, -1);
  OPTION(optsc, hyperpar, -1);
  OPTION(optsc, ExBdiff, -1);
  OPTION(optsc, ExBpar, false);
  OPTION(optsc, ADpar, -1);
  OPTION(optsc, ADpar_phine, false);
  OPTION(optsc, ADpar_bndry, false);
  OPTION(optsc, low_pass_z, -1);
  OPTION(optsc, x_hyper_viscos, -1.0);
  OPTION(optsc, y_hyper_viscos, -1.0);
  OPTION(optsc, z_hyper_viscos, -1.0);
  
  OPTION(optsc, ne_hyper_z, -1.0);
  OPTION(optsc, pe_hyper_z, -1.0);

  OPTION(optsc, low_n_diffuse, true);
  OPTION(optsc, low_n_diffuse_perp, false);

  OPTION(optsc, vepsi_dissipation, false);
  
  OPTION(optsc, resistivity_multiply, 1.0);

  // Sheath boundary conditions
  OPTION(optsc, sheath_model, 0);
  OPTION(optsc, sheath_gamma, 6.5);
  OPTION(optsc, neutral_gamma, 5. / 4);
  OPTION(optsc, neutral_vwall, 1. / 3); // 1/3rd Franck-Condon energy at wall
  OPTION(optsc, sheath_yup, true);      // Apply sheath at yup?
  OPTION(optsc, sheath_ydown, true);    // Apply sheath at ydown?

  OPTION(optsc, outflow_ydown, false); // Allow outflowing neutrals?

  // Fix profiles in SOL
  OPTION(optsc, sol_fix_profiles, false);
  if (sol_fix_profiles) {
    sol_ne = FieldFactory::get()->parse("sol_ne", optsc);
    sol_te = FieldFactory::get()->parse("sol_te", optsc);
  }

  OPTION(optsc, radial_buffers, false);
  OPTION(optsc, radial_inner_width, 4);
  OPTION(optsc, radial_outer_width, 4);
  OPTION(optsc, radial_buffer_D, 1.0);

  // Output additional information
  OPTION(optsc, verbose, false);    // Save additional fields
  OPTION(optsc, output_ddt, false); // Save time derivatives

  // Normalisation
  OPTION(optsc, Tnorm, 100);  // Reference temperature [eV]
  OPTION(optsc, Nnorm, 1e19); // Reference density [m^-3]
  OPTION(optsc, Bnorm, 1.0);  // Reference magnetic field [T]

  OPTION(optsc, AA, 2.0); // Ion mass (2 = Deuterium)

  output.write("Normalisation Te=%e, Ne=%e, B=%e\n", Tnorm, Nnorm, Bnorm);
  SAVE_ONCE4(Tnorm, Nnorm, Bnorm, AA); // Save

  Cs0 = sqrt(qe * Tnorm / (AA * Mp)); // Reference sound speed [m/s]
  Omega_ci = qe * Bnorm / (AA * Mp);  // Ion cyclotron frequency [1/s]
  rho_s0 = Cs0 / Omega_ci;

  mi_me = AA * Mp / Me;
  
  // beta_e = Electron pressure / magnetic pressure
  beta_e = qe * Tnorm * Nnorm / (SQ(Bnorm) / (2.*mu0));

  output.write("\tmi_me=%e, beta_e=%e\n", mi_me, beta_e);
  SAVE_ONCE2(mi_me, beta_e);

  output.write("\t Cs=%e, rho_s=%e, Omega_ci=%e\n", Cs0, rho_s0, Omega_ci);
  SAVE_ONCE3(Cs0, rho_s0, Omega_ci);

  // Collision times
  BoutReal Coulomb = 6.6 - 0.5 * log(Nnorm * 1e-20) + 1.5 * log(Tnorm);
  tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3. / 2));
  tau_i0 = sqrt(AA) / (4.80e-8 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3. / 2));

  output.write("\ttau_e0=%e, tau_i0=%e\n", tau_e0, tau_i0);

  // Get a typical major radius (for flux limiter)
  Field2D Rxy;
  if (mesh->get(Rxy, "Rxy")) {
    R0 = 1.0; // Set to 1 meter
  } else {
    R0 = max(Rxy, true); // Maximum over all processors
  }
  output.write("\t R0 = %e m\n", R0);
  R0 /= rho_s0; // Normalise R0

  if (anomalous_D > 0.0) {
    // Normalise
    anomalous_D /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous D_perp = %e\n", anomalous_D);
  }
  if (anomalous_chi > 0.0) {
    // Normalise
    anomalous_chi /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous chi_perp = %e\n", anomalous_chi);
  }
  if (anomalous_nu > 0.0) {
    // Normalise
    anomalous_nu /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous nu_perp = %e\n", anomalous_nu);
  }

  if (ramp_mesh) {
    Jpar0 = 0.0;
  } else {
    // Read equilibrium current density
    // GRID_LOAD(Jpar0);
    // Jpar0 /= qe*Nnorm*Cs0;
    Jpar0 = 0.0;
  }

  string source;
  FieldFactory fact(mesh);

  if (sinks) {
    optsc->get("sink_invlpar", source, "0.05"); // 20 m
    sink_invlpar = fact.create2D(source);
    sink_invlpar *= rho_s0; // Normalise
    SAVE_ONCE(sink_invlpar);

    if (drift_wave) {
      alpha_dw = fact.create2D("Hermes:alpha_dw");
      SAVE_ONCE(alpha_dw);
    }
  }

  // Get switches from each variable section
  Options *optne = opt->getSection("Ne");
  optne->get("D", Dne, -1.0);
  optne->get("source", source, "0.0");
  Sn = fact.create2D(source);
  Sn /= Omega_ci;

  // Inflowing density carries momentum
  OPTION(optne, density_inflow, false);

  Options *optvort = opt->getSection("Vort");
  optvort->get("D", Dvort, -1.0);
  optvort->get("mupar", vort_mupar, -1);

  Options *optve = opt->getSection("VePsi");
  optve->get("D", Dve, -1.0);

  Options *optpe = opt->getSection("Pe");
  optpe->get("D", Dte, -1.0);
  optpe->get("source", source, "0.0");
  Spe = fact.create2D(source);
  Spe /= Omega_ci;

  OPTION(optsc, core_sources, false);
  if (core_sources) {
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      if (!mesh->periodicY(x)) {
        // Not periodic, so not in core
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          Sn(x, y) = 0.0;
          Spe(x, y) = 0.0;
        }
      }
    }
  }

  // Multiply sources by g11
  OPTION(optsc, source_vary_g11, false);
  if (source_vary_g11) {
    // Average metric tensor component
    g11norm = mesh->g11 / averageY(mesh->g11);
  }

  // Mid-plane power flux q_||
  string midplane_power;
  OPTION(optpe, midplane_power, "0.0");
  // Midplane power specified in Watts per m^2
  Field2D qfact;
  GRID_LOAD(qfact); // Factor to multiply to get volume source
  Field2D qmid = fact.create2D(midplane_power) * qfact;
  // Normalise from W/m^3
  qmid /= qe * Tnorm * Nnorm * Omega_ci;
  Spe += (2. / 3) * qmid;

  Options *optnvi = opt->getSection("NVi");
  optnvi->get("D", Dvi, -1.0);

  // Add variables to solver
  solver->add(Ne, "Ne");
  solver->add(Pe, "Pe");
  EvolvingVars.add(Ne, Pe);

  if (output_ddt) {
    SAVE_REPEAT2(ddt(Ne), ddt(Pe));
  }

  if (j_par || j_diamag) {
    // Have a source of vorticity
    solver->add(Vort, "Vort");
    EvolvingVars.add(Vort);
    if (output_ddt) {
      SAVE_REPEAT(ddt(Vort));
    }

    // Save the potential
    SAVE_REPEAT(phi);
  } else {
    Vort = 0.0;
  }

  if (electromagnetic || FiniteElMass) {
    solver->add(VePsi, "VePsi");
    EvolvingVars.add(VePsi);
    if (output_ddt) {
      SAVE_REPEAT(ddt(VePsi));
    }
  } else {
    // If both electrostatic and zero electron mass,
    // then Ohm's law has no time-derivative terms,
    // but is calculated from other evolving quantities
    VePsi = 0.0;
  }

  if (ion_velocity) {
    solver->add(NVi, "NVi");
    EvolvingVars.add(NVi);
    if (output_ddt) {
      SAVE_REPEAT(ddt(NVi));
    }
  } else {
    NVi = 0.0;
  }

  // No Ti evolution, but Ti can be used in parallel viscosity, sink
  output.write("  Fixed (constant) ion temperature\n");
  optsc->get("Ti", Ti, Tnorm); // Ion temperature [eV]
  Ti /= Tnorm;                 // Normalise

  OPTION(optsc, adapt_source, false);

  if (adapt_source) {
    // Adaptive sources to match profiles

    // PI controller, including an integrated difference term
    OPTION(optsc, source_p, 1e-2);
    OPTION(optsc, source_i, 1e-6);

    // Save the actual source added
    SAVE_REPEAT2(NeSource, PeSource);

    OPTION(optsc, adapt_fix_form, false);
    if (adapt_fix_form) {
      // Fix the form of the source, changing amplitude
      // Sum over the grid, for later normalisation
      total_Sn = total_Spe = 0.0;
      for (int i = 0; i < mesh->ngx; i++) {
        for (int j = 0; j < mesh->ngy; j++) {
          total_Sn += Sn(i, j);
          total_Spe += Spe(i, j);
        }
      }
      BoutReal local_values[2], values[2];
      local_values[0] = total_Sn;
      local_values[1] = total_Spe;
      MPI_Allreduce(local_values, values, 2, MPI_DOUBLE, MPI_SUM,
                    BoutComm::get());
      total_Sn = values[0];
      total_Spe = values[1];

      density_error_lasttime = -1.0; // signal no value
      solver->addToRestart(density_error_integral, "density_error_integral");

      pe_error_lasttime = -1.0;
      solver->addToRestart(pe_error_integral, "pe_error_integral");

      if (!restarting) {
        density_error_integral = -1. / source_i;
        pe_error_integral = -1. / source_i;
      }

      OPTION(optsc, temperature_feedback, false); // Feedback on Te rather than Pe
      
    } else {
      // Evolving the source profiles

      Field2D Snsave = copy(Sn);
      Field2D Spesave = copy(Spe);
      SOLVE_FOR2(Sn, Spe);
      Sn = Snsave;
      Spe = Spesave;
    }
  } else {
    SAVE_ONCE2(Sn, Spe);
  }

  /////////////////////////////////////////////////////////
  // Load metric tensor from the mesh, passing length and B
  // field normalisations
  TRACE("Loading metric tensor");

  bool loadmetric;
  OPTION(optsc, loadmetric, true);
  if (loadmetric) {
    // Load Rxy, Bpxy etc. to create orthogonal metric
    LoadMetric(rho_s0, Bnorm);
  } else {
    // To use non-orthogonal metric
    // Normalise
    mesh->dx /= rho_s0 * rho_s0 * Bnorm;
    mesh->Bxy /= Bnorm;
    // Metric is in grid file - just need to normalise
    mesh->g11 /= (Bnorm * Bnorm * rho_s0 * rho_s0);
    mesh->g22 *= (rho_s0 * rho_s0);
    mesh->g33 *= (rho_s0 * rho_s0);
    mesh->g12 /= Bnorm;
    mesh->g13 /= Bnorm;
    mesh->g23 *= (rho_s0 * rho_s0);

    mesh->J *= Bnorm / rho_s0;

    mesh->g_11 *= (Bnorm * Bnorm * rho_s0 * rho_s0);
    mesh->g_22 /= (rho_s0 * rho_s0);
    mesh->g_33 /= (rho_s0 * rho_s0);
    mesh->g_12 *= Bnorm;
    mesh->g_13 *= Bnorm;
    mesh->g_23 /= (rho_s0 * rho_s0);

    mesh->geometry(); // Calculate other metrics
  }

  /////////////////////////////////////////////////////////
  // Toroidal rotation
  // Note that this must be done after coordinates are a

  if (rotation) {
    // Normalise rotation rate from Hz to radians per cyclotron time

    rotation_rate *= TWOPI / Omega_ci;
    output.write("Normalised rotation rate = %e\n", rotation_rate);

    // Calculate vectors needed, mainly Grad(Z) and Grad(R)
    Field2D Rxy;
    if (mesh->get(Rxy, "Rxy")) {
      throw BoutException("Rotation requires an Rxy variable");
    }
    Field2D Zxy;
    if (mesh->get(Zxy, "Zxy")) {
      throw BoutException("Rotation requires an Rxy variable");
    }

    // Normalise Rxy and Zxy
    Rxy /= rho_s0;
    Zxy /= rho_s0;

    // Rotation rate vector for Coriolis force
    Omega_vec = rotation_rate * Grad(Zxy);

    // Magnetic field unit vector
    Vector2D b0vec;
    b0vec.covariant = false;
    b0vec.x = b0vec.z = 0;
    b0vec.y = 1. / mesh->g_22;

    // b cross Grad(R) for centrifugal force
    bxGradR = b0vec ^ Grad(Rxy);

    // Fill guard cells by communicating and applying boundary condition
    mesh->communicate(Omega_vec, bxGradR);
    Omega_vec.x.applyBoundary("neumann");
    Omega_vec.y.applyBoundary("neumann");
    Omega_vec.z.applyBoundary("neumann");

    bxGradR.x.applyBoundary("neumann");
    bxGradR.y.applyBoundary("neumann");
    bxGradR.z.applyBoundary("neumann");

    SAVE_ONCE2(Omega_vec, bxGradR);
  }

  /////////////////////////////////////////////////////////
  // Neutral models

  TRACE("Initialising neutral models");
  OPTION(optsc, neutral_model, 0);
  switch (neutral_model) {
  case 1: {
    // 2D (X-Z) diffusive model
    // Neutral gas dynamics
    solver->add(Nn, "Nn");
    solver->add(Pn, "Pn");
    EvolvingVars.add(Nn, Pn);

    Dnn = 0.0; // Neutral gas diffusion

    S = 0;
    F = 0;
    Q = 0;
    Rp = 0;
    Rn = 0;

    SAVE_REPEAT(Dnn);
    break;
  }
  case 2: {

    Nn = 0.0;         // initialise to 0
    Tn = 3.5 / Tnorm; // 3.5eV  for Franck-Condon energy

    S = 0;
    F = 0;
    Q = 0;
    Rp = 0;
    SAVE_REPEAT4(Nn0, lambda, lambda_int, Nn);
    break;
  }
  case 3: {
    /*! 2D (X-Y) full velocity model
     *
     * Vn2D is covariant (the default), so has components
     * V_x, V_y, V_z
     */
    Options *options = Options::getRoot()->getSection("neutral");
    options->get("gamma_ratio", gamma_ratio, 5. / 3);
    options->get("viscosity", neutral_viscosity, 1e-2);
    options->get("bulk", neutral_bulk, 1e-2);
    options->get("conduction", neutral_conduction, 1e-2);

    // Evolve 2D density, pressure, and velocity
    solver->add(Nn2D, "Nn");
    solver->add(Pn2D, "Pn");
    solver->add(Vn2D, "Vn");

    DivV2D.setBoundary("Pn"); // Same boundary condition as Pn
    SAVE_REPEAT(DivV2D);

    // Load necessary metrics for non-orth calculation
    Field2D etaxy, cosbeta;
    if (mesh->get(etaxy, "etaxy")) {
      etaxy = 0.0;
    }
    cosbeta = sqrt(1. - SQ(etaxy));

    // Calculate transformation to Cartesian coordinates
    Field2D Rxy, Zxy, hthe, Bpxy;

    if (mesh->get(Rxy, "Rxy")) {
      throw BoutException("Fluid neutrals model requires Rxy");
    }
    if (mesh->get(Zxy, "Zxy")) {
      throw BoutException("Fluid neutrals model requires Zxy");
    }
    if (mesh->get(hthe, "hthe")) {
      throw BoutException("Fluid neutrals model requires hthe");
    }
    if (mesh->get(Bpxy, "Bpxy")) {
      throw BoutException("Fluid neutrals model requires Bpxy");
    }

    // Normalise
    Rxy /= rho_s0;
    Zxy /= rho_s0;
    hthe /= rho_s0;
    Bpxy /= Bnorm;

    // Axisymmetric neutrals simplifies things considerably...

    Urx.allocate();
    Ury.allocate();
    Uzx.allocate();
    Uzy.allocate();

    Txr.allocate();
    Txz.allocate();
    Tyr.allocate();
    Tyz.allocate();

    for (int i = 0; i < mesh->ngx; i++)
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        // Central differencing of coordinates
        BoutReal dRdtheta, dZdtheta;
        if (j == mesh->ystart) {
          dRdtheta = (Rxy(i, j + 1) - Rxy(i, j)) / (mesh->dy(i, j));
          dZdtheta = (Zxy(i, j + 1) - Zxy(i, j)) / (mesh->dy(i, j));
        } else if (j == mesh->yend) {
          dRdtheta = (Rxy(i, j) - Rxy(i, j - 1)) / (mesh->dy(i, j));
          dZdtheta = (Zxy(i, j) - Zxy(i, j - 1)) / (mesh->dy(i, j));
        } else {
          dRdtheta = (Rxy(i, j + 1) - Rxy(i, j - 1)) / (2. * mesh->dy(i, j));
          dZdtheta = (Zxy(i, j + 1) - Zxy(i, j - 1)) / (2. * mesh->dy(i, j));
        }

        // Match to hthe, 1/|Grad y|
        BoutReal h = sqrt(SQ(dRdtheta) + SQ(dZdtheta));
        BoutReal grady = 1.0 / hthe(i, j);
        dRdtheta = dRdtheta / grady / h;
        dZdtheta = dZdtheta / grady / h;

        BoutReal dRdpsi, dZdpsi;
        if (i == 0) {
          // One-sided differences
          dRdpsi = (Rxy(i + 1, j) - Rxy(i, j)) / (mesh->dx(i, j));
          dZdpsi = (Zxy(i + 1, j) - Zxy(i, j)) / (mesh->dx(i, j));
          // dRdpsi = (-0.5*Rxy(i+2,j) + 2.*Rxy(i+1,j) -
          // 1.5*Rxy(i,j))/(mesh->dx(i,j));
          // dZdpsi = (-0.5*Zxy(i+2,j) + 2.*Zxy(i+1,j) -
          // 1.5*Zxy(i,j))/(mesh->dx(i,j));
        } else if (i == (mesh->ngx - 1)) {
          // One-sided differences
          dRdpsi = (Rxy(i, j) - Rxy(i - 1, j)) / (mesh->dx(i, j));
          dZdpsi = (Zxy(i, j) - Zxy(i - 1, j)) / (mesh->dx(i, j));
        } else {
          dRdpsi = (Rxy(i + 1, j) - Rxy(i - 1, j)) / (2. * mesh->dx(i, j));
          dZdpsi = (Zxy(i + 1, j) - Zxy(i - 1, j)) / (2. * mesh->dx(i, j));
        }

        // Match to Bp, |Grad psi|. NOTE: this only works if
        // X and Y are orthogonal.
        BoutReal dldpsi =
            sqrt(SQ(dRdpsi) + SQ(dZdpsi)) * cosbeta(i, j); // ~ 1/(R*Bp)
        dRdpsi /= dldpsi * Bpxy(i, j) * Rxy(i, j);
        dZdpsi /= dldpsi * Bpxy(i, j) * Rxy(i, j);

        Urx(i, j) = dRdpsi;
        Ury(i, j) = dRdtheta;
        Uzx(i, j) = dZdpsi;
        Uzy(i, j) = dZdtheta;

        // Poloidal (R,Z) transformation Jacobian
        BoutReal J = dRdpsi * dZdtheta - dZdpsi * dRdtheta;

        Txr(i, j) = dZdtheta / J;
        Txz(i, j) = -dRdtheta / J;
        Tyr(i, j) = -dZdpsi / J;
        Tyz(i, j) = dRdpsi / J;
      }

    Urx.applyBoundary("neumann");
    Ury.applyBoundary("neumann");
    Uzx.applyBoundary("neumann");
    Uzy.applyBoundary("neumann");

    Txr.applyBoundary("neumann");
    Txz.applyBoundary("neumann");
    Tyr.applyBoundary("neumann");
    Tyz.applyBoundary("neumann");

    SAVE_ONCE4(Urx, Ury, Uzx, Uzy);
    SAVE_ONCE4(Txr, Txz, Tyr, Tyz);

    // Atomic processes
    S = 0;
    F = 0;
    Q = 0;
    Rp = 0;

    break;
  }
  case 4: {
    // 3D model, diffusive in X-Z and fluid in Y
    solver->add(Nn, "Nn");
    solver->add(Pn, "Pn");
    solver->add(NVn, "NVn");
    EvolvingVars.add(Nn, Pn, NVn);

    if (output_ddt) {
      SAVE_REPEAT3(ddt(Nn), ddt(Pn), ddt(NVn));
    }

    Dnn = 0.0; // Neutral gas diffusion

    S = 0;
    F = 0;
    Q = 0;
    Rp = 0;

    SAVE_REPEAT(Dnn);
    break;
  }
  }
  // Shared save-repeats
  if (neutral_model > 0) {
    SAVE_REPEAT4(S, F, Q, Rp);
    SAVE_ONCE(neutral_model);
  }

  // Load hthe - only needed for neutral_model=2
  GRID_LOAD(hthe);

  /////////////////////////////////////////////////////////
  // Radiation due to hydrogen excitation

  OPTION(optsc, excitation, false); // Include electron impact excitation?

  /////////////////////////////////////////////////////////
  // Impurities
  TRACE("Impurities");

  OPTION(optsc, carbon_fraction, -1.);
  if (carbon_fraction > 0.0) {
    SAVE_REPEAT(Rzrad);
    SAVE_ONCE(carbon_fraction);
    carbon_rad = new HutchinsonCarbonRadiation();
  }

  /////////////////////////////////////////////////////////
  // Read profiles from the mesh
  TRACE("Reading profiles");

  Field2D NeMesh, TeMesh;
  if (mesh->get(NeMesh, "Ne0")) {
    // No Ne0. Try Ni0
    if (mesh->get(NeMesh, "Ni0")) {
      output << "WARNING: Neither Ne0 nor Ni0 found in mesh input\n";
    }
  }
  NeMesh *= 1e20; // Convert to m^-3

  NeMesh /= Nnorm; // Normalise

  if (mesh->get(TeMesh, "Te0")) {
    // No Te0
    output << "WARNING: Te0 not found in mesh\n";
    // Try to read Ti0
    if (mesh->get(TeMesh, "Ti0")) {
      // No Ti0 either
      output << "WARNING: No Te0 or Ti0. Setting TeMesh to 0.0\n";
      TeMesh = 0.0;
    }
  }

  TeMesh /= Tnorm; // Normalise

  NeTarget = NeMesh;
  TeTarget = TeMesh;
  PeTarget = NeMesh * TeMesh;

  if (!restarting && !ramp_mesh) {
    bool startprofiles;
    OPTION(optsc, startprofiles, true);
    if (startprofiles) {
      Ne += NeMesh; // Add profiles in the mesh file

      Pe += NeMesh * TeMesh;
    }

    // Check for negatives
    if (min(Ne, true) < 0.0) {
      throw BoutException("Starting density is negative");
    }
    if (max(Ne, true) < 1e-5) {
      throw BoutException("Starting density is too small");
    }

    if (min(Pe, true) < 0.0) {
      throw BoutException("Starting pressure is negative");
    }
    if (max(Pe, true) < 1e-5) {
      throw BoutException("Starting pressure is too small");
    }

    mesh->communicate(Ne, Pe);
  }

  /////////////////////////////////////////////////////////
  // Read curvature components
  TRACE("Reading curvature");

  Curlb_B.covariant = false; // Contravariant
  mesh->get(Curlb_B, "bxcv");
  if (mesh->ShiftXderivs) {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }

  Curlb_B.x /= Bnorm;
  Curlb_B.y *= rho_s0 * rho_s0;
  Curlb_B.z *= rho_s0 * rho_s0;

  Curlb_B *= 2. / mesh->Bxy;

  if (j_par) {
    SAVE_REPEAT(Ve);

    if (electromagnetic)
      SAVE_REPEAT(psi);
  }

  OPTION(optsc, split_n0, false); // Split into n=0 and n~=0
  OPTION(optsc, split_n0_psi, split_n0);
  // Phi solver
  if (phi3d) {
#ifdef PHISOLVER
    phiSolver3D = Laplace3D::create();
#endif
  } else {
    if (split_n0) {
      // Create an XY solver for n=0 component
      laplacexy = new LaplaceXY(mesh);
      // Set coefficients for Boussinesq solve
      laplacexy->setCoefs(1. / SQ(mesh->Bxy), 0.0);
      phi2D = 0.0; // Starting guess
    }

    // Create an XZ solver
    OPTION(optsc, newXZsolver, false);
    if (newXZsolver) {
      // Test new LaplaceXZ solver
      newSolver = LaplaceXZ::create(mesh);
      // Set coefficients for Boussinesq solve
      newSolver->setCoefs(1. / SQ(mesh->Bxy), 0.0);
    } else {
      // Use older Laplacian solver
      phiSolver = Laplacian::create(opt->getSection("phiSolver"));
      // Set coefficients for Boussinesq solve
      phiSolver->setCoefC(1. / SQ(mesh->Bxy));
    }
    phi = 0.0;
  }

  // Apar (Psi) solver
  // aparSolver = Laplacian::create(opt->getSection("aparSolver"));
  aparSolver = LaplaceXZ::create(mesh, opt->getSection("aparSolver"));
  if (split_n0_psi) {
    // Use another XY solver for n=0 psi component
    aparXY = new LaplaceXY(mesh);
    psi2D = 0.0;
  }

  Ve.setBoundary("Ve");
  nu.setBoundary("nu");
  phi.setBoundary("phi"); // For y boundaries
  Jpar.setBoundary("Jpar");

  nu = 0.0;
  kappa_epar = 0.0;
  Dn = 0.0;

  SAVE_REPEAT(Telim);
  if (verbose) {
    SAVE_REPEAT(Jpar);
    SAVE_REPEAT(kappa_epar);
    if (resistivity)
      SAVE_REPEAT(nu);
    // SAVE_REPEAT2(wall_flux, wall_power);
  }

  psi = phi = 0.0;

  // Preconditioner
  setPrecon((preconfunc)&Hermes::precon);

  return 0;
}

int Hermes::rhs(BoutReal time) {
  // printf("TIME = %e\r", time);

  if (!evolve_plasma) {
    Ne = 0.0;
    Pe = 0.0;
    Vort = 0.0;
    VePsi = 0.0;
    NVi = 0.0;
    sheath_model = 0;
  }

  // Communicate evolving variables
  mesh->communicate(EvolvingVars);

  Field3D Nelim = floor(Ne, 1e-5);

  Te = Pe / Nelim;
  Vi = NVi / Nelim;

  Telim = floor(Te, 0.1 / Tnorm);

  Field3D Pelim = Telim * Nelim;

  Field3D logPelim = log(Pelim);
  logPelim.applyBoundary("neumann");

  // Are there any currents? If not, then there is no source
  // for vorticity, phi = 0 and jpar = 0
  bool currents = j_par | j_diamag;

  // Local sound speed. Used for parallel advection operator
  Field3D sound_speed = sqrt(Telim * (5. / 3));

  //////////////////////////////////////////////////////////////
  // Calculate electrostatic potential phi
  //
  //

  TRACE("Electrostatic potential");

  if (!currents) {
    // Disabling electric fields
    // phi = 0.0; // Already set in initialisation
  } else {
    // Solve phi from Vorticity
    if (phi3d) {
#ifdef PHISOLVER
      phiSolver3D->setCoefC(Ne / SQ(mesh->Bxy));
      // phi.setBoundaryTo(3.*Te);
      if (mesh->lastX()) {
        for (int i = mesh->xend + 1; i < mesh->ngx; i++)
          for (int j = mesh->ystart; j <= mesh->yend; j++)
            for (int k = 0; k < mesh->ngz - 1; k++) {
              phi(i, j, k) = 3. * Te(i, j, k);
            }
      }
      phi = phiSolver3D->solve(Vort, phi);
#endif
    } else {
      // Phi flags should be set in BOUT.inp
      // phiSolver->setInnerBoundaryFlags(INVERT_DC_GRAD);
      // phiSolver->setOuterBoundaryFlags(INVERT_SET);

      // Sheath multiplier Te -> phi (2.84522 for Deuterium)
      BoutReal sheathmult = log(0.5 * sqrt(mi_me / PI));

      if (boussinesq) {

        if (split_n0) {
          ////////////////////////////////////////////
          // Boussinesq, split
          // Split into axisymmetric and non-axisymmetric components
          Field2D Vort2D = Vort.DC(); // n=0 component

          // Set the boundary to 2.8*Te
          // phi2D.setBoundaryTo(sheathmult * floor(Te.DC(), 0.0));
          phi2D.setBoundaryTo(sheathmult * Telim.DC());

          phi2D = laplacexy->solve(Vort2D, phi2D);
          
          // Solve non-axisymmetric part using X-Z solver
          if (newXZsolver) {
            newSolver->setCoefs(1. / SQ(mesh->Bxy), 0.0);
            phi = newSolver->solve(Vort - Vort2D, phi);
          } else {
            phiSolver->setCoefC(1. / SQ(mesh->Bxy));
            // phi = phiSolver->solve((Vort-Vort2D)*SQ(mesh->Bxy), phi);
            phi = phiSolver->solve((Vort - Vort2D) * SQ(mesh->Bxy),
                                   sheathmult * (Telim - Telim.DC()));
          }
          phi += phi2D; // Add axisymmetric part
        } else {
          ////////////////////////////////////////////
          // Boussinesq, non-split
          // Solve all components using X-Z solver

          if (newXZsolver) {
            // Use the new LaplaceXY solver
            // newSolver->setCoefs(1./SQ(mesh->Bxy), 0.0); // Set when
            // initialised
            phi = newSolver->solve(Vort, phi);
          } else {
            // Use older Laplacian solver
            // phiSolver->setCoefC(1./SQ(mesh->Bxy)); // Set when initialised
            phi = phiSolver->solve(Vort * SQ(mesh->Bxy), sheathmult * Telim);

            /*
            Field2D phidc = phi.DC();
            // Get gradient at outer boundary
            BoutReal gradient[mesh->ngy];
            for (int y=0;y<mesh->ngy;y++) {
              gradient[y] = phidc(mesh->xend+1,y) - phidc(mesh->xend,y);
            }
            // Broadcast to other processors sharing this X boundary
            MPI_Comm xcomm = mesh->getXcomm();
            int nxpe;
            MPI_Comm_size(xcomm, &nxpe);
            MPI_Bcast( gradient, mesh->ngy, MPI_DOUBLE, nxpe-1, xcomm);
            // Now subtract the gradient from the mean electric field
            // so that the gradient is zero at the outer edge
            // This imposes zero value and zero gradient at the outer
            // boundary

            for (int i=0;i<mesh->ngx;i++) {
              BoutReal xdist = mesh->GlobalNx - mesh->XGLOBAL(i) - 2.5;
              //output.write("%d: %e\n", i, xdist);
              for (int j=0;j<mesh->ngy;j++) {
                BoutReal change = xdist * gradient[j];
                for (int k=0;k<mesh->ngz-1;k++) {
                  phi(i,j,k) += change;
                }
              }
            }
            */
          }
        }
      } else {
        ////////////////////////////////////////////
        // Non-Boussinesq
        //
        phiSolver->setCoefC(Nelim / SQ(mesh->Bxy));
        phi = phiSolver->solve(Vort * SQ(mesh->Bxy) / Nelim, sheathmult * Telim);
      }
    }
    phi.applyBoundary(time);
    mesh->communicate(phi);
  }

  //////////////////////////////////////////////////////////////
  // Collisions and stress tensor
  TRACE("Collisions");

  tau_e = (Cs0 / rho_s0) * tau_e0 * (Telim ^ 1.5) /
          Nelim; // Normalised electron collision time

  tau_i = (Cs0 / rho_s0) * tau_i0 * pow(Ti, 1.5) /
          Nelim; // Normalised ion-ion collision time

  {
    // Smooth tau_e by averaging over nearby grid cells in X and Z
    Field3D tt = copy(tau_e);
    for (int i = 1; i < mesh->ngx - 1; i++) {
      for (int j = 0; j < mesh->ngy; j++) {
        for (int k = 0; k < mesh->ngz - 1; k++) {
          int kp = (k + 1) % (mesh->ngz - 1);
          int km = (k - 2 + mesh->ngz) % (mesh->ngz - 1);

          tau_e(i, j, k) =
              0.2 * (tt(i, j, k) + tt(i + 1, j, k) + tt(i - 1, j, k) +
                     tt(i, j, kp) + tt(i, j, km));
        }
      }
    }
  }

  // Collisional damping (normalised)
  if (resistivity || (!electromagnetic && !FiniteElMass)) {
    // Need to calculate nu if electrostatic and zero electron mass
    nu = resistivity_multiply / (1.96 * tau_e * mi_me);

    if (electron_neutral && (neutral_model != 0)) {
      /*
       * Include electron-neutral collisions. These can dominate
       * the resistivity at low temperatures (~1eV)
       *
       * This assumes a fixed cross-section, independent of energy
       *
       */
      BoutReal a0 = PI * SQ(5.29e-11); // Cross-section [m^2]

      // Electron thermal speed
      Field3D vth_e = sqrt(Telim);

      // Electron-neutral collision rate
      Field3D nu_ne = vth_e * Nnorm * Nn * a0 * rho_s0;

      // Add collision rate to the electron-ion rate
      nu += nu_ne;
    }
  }

  if (thermal_conduction || sinks) {
    // Braginskii expression for parallel conduction
    // kappa ~ n * v_th^2 * tau
    kappa_epar = 3.2 * mi_me * Telim * Nelim * tau_e;

    if (flux_limit_alpha > 0) {
      /* Apply a flux limiter
       *
       * Free streaming heat conduction calculated as:
       *     kappa_fs = alpha * n V_th R0
       *
       * then the effective heat flux is
       * kappa = kappa_par * kappa_fs / (kappa_par + kappa_fs)
       */

      Field3D kappa_fs = flux_limit_alpha * Nelim * sqrt(mi_me * Telim) * R0;
      kappa_epar = kappa_epar * kappa_fs / (kappa_epar + kappa_fs);
    }

    if (kappa_limit_alpha > 0.0) {
      /*
       * A better flux limiter, as used in SOLPS. This doesn't require
       * an input length scale R0
       *
       * Calculate the heat flux from Spitzer-Harm and flux limit
       *
       * Typical value of alpha ~ 0.2 for electrons
       *
       * R.Schneider et al. Contrib. Plasma Phys. 46, No. 1-2, 3 – 191 (2006)
       * DOI 10.1002/ctpp.200610001
       */

      // Spitzer-Harm heat flux
      Field3D q_SH = kappa_epar * Grad_par(Te);
      Field3D q_fl = kappa_limit_alpha * Nelim * Telim * sqrt(mi_me * Telim);

      kappa_epar = kappa_epar / (1. + abs(q_SH / q_fl));

      // Values of kappa on cell boundaries are needed for fluxes
      mesh->communicate(kappa_epar);
      kappa_epar.applyBoundary("dirichlet");
    }
  }

  nu.applyBoundary(time);

  //////////////////////////////////////////////////////////////
  // Calculate perturbed magnetic field psi
  TRACE("Calculating psi");

  if (!currents) {
    // No magnetic fields or currents
    psi = 0.0;
    Jpar = 0.0;
    Ve = 0.0; // Ve will be set after the sheath boundaries below
  } else {
    // Calculate electomagnetic potential psi from VePsi
    // VePsi = Ve - Vi + 0.5 * mi_me * beta_e * psi
    // where the first term comes from finite electron mass, and the second
    // from the parallel component of the electric field
    // Note that psi is -A_|| so Jpar = Delp2(psi)
    if (electromagnetic) {
      if (FiniteElMass) {
        // Solve Helmholtz equation for psi

        // aparSolver->setCoefA(-Ne*0.5*mi_me*beta_e);
        Field2D NDC = Ne.DC();
        aparSolver->setCoefs(1.0, -NDC * 0.5 * mi_me * beta_e);
        // aparSolver->setCoefs(1.0, -Ne*0.5*mi_me*beta_e);
        if (split_n0_psi) {
          // Solve n=0 component separately

          aparXY->setCoefs(1.0, -NDC * 0.5 * mi_me * beta_e);

          Field2D JDC = -NDC * VePsi.DC();
          aparXY->solve(JDC, psi2D);

          psi = aparSolver->solve(-NDC * VePsi - JDC, psi - psi2D) + psi2D;
        } else {
          psi = aparSolver->solve(-NDC * VePsi, psi);
          // psi = aparSolver->solve(-Ne*VePsi, psi);
        }

        Ve = VePsi - 0.5 * mi_me * beta_e * psi + Vi;

        Ve.applyBoundary(time);
        mesh->communicate(Ve, psi);

        Jpar = Ne * (Vi - Ve);
        Jpar.applyBoundary(time);

      } else {
        // Zero electron mass
        // No Ve term in VePsi, only electromagnetic term

        psi = VePsi / (0.5 * mi_me * beta_e);

        // Ve = (NVi - Delp2(psi)) / Nelim;
        Jpar = Div_Perp_Lap_FV(1.0, psi, true);
        // Jpar = Div_Perp_Lap_XYZ(1.0, psi, true);
        mesh->communicate(Jpar);

        Jpar.applyBoundary(time);
        Ve = (NVi - Jpar) / Nelim;
      }

      // psi -= psi.DC(); // Remove toroidal average, only keep fluctuations
    } else {
      // Electrostatic
      psi = 0.0;
      if (FiniteElMass) {
        // No psi contribution to VePsi
        Ve = VePsi + Vi;
      } else {
        // Zero electron mass and electrostatic.
        // Special case where Ohm's law has no time-derivatives
        mesh->communicate(phi);

        Ve = Vi + (Grad_parP_CtoL(phi) - Grad_parP_CtoL(Pe) / Ne) / nu;

        if (thermal_force) {
          Ve -= 0.71 * Grad_parP_CtoL(Te) / nu;
        }
      }

      Ve.applyBoundary(time);
      // Communicate auxilliary variables
      mesh->communicate(Ve);

      Jpar = NVi - Ne * Ve;
    }
    // Ve -= Jpar0 / Ne; // Equilibrium current density

    // Limit Ve to less than c
    // Ve = ceil(Ve, 3e8/Cs0);
  }

  //////////////////////////////////////////////////////////////
  // Sheath boundary conditions on Y up and Y down
  TRACE("Sheath boundaries");

  if (sheath_ydown) {
    switch (sheath_model) {
    case 0: { // Normal Bohm sheath
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = floor(Ne(r.ind, mesh->ystart, jz), 0.0);

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = -sqrt(tesheath); // Sound speed outwards

          if (Vi(r.ind, mesh->ystart, jz) < visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->ystart, jz);
          }

          // Sheath current
          // Note that phi/Te >= 0.0 since for phi < 0
          // vesheath is the electron saturation current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->ystart, jz), 0.0);
          BoutReal vesheath =
              -sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);

          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
            // Neumann conditions
            Ne(r.ind, jy, jz) = nesheath;
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->ystart, jz);
            Pelim(r.ind, jy, jz) = Pelim(r.ind, mesh->ystart, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->ystart, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->ystart, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->ystart, jz);
          }
        }
      }
      break;
    }
    case 1: {
      /*
        Loizu boundary conditions

        Temperature

        Grad_par(Te) = 0

        Density equation

        Grad_par(n) = -(n/Cs) Grad_par(Vi)
        -> n_p - n_m = - (n_p + n_m)/(2Cs) (Vi_p - Vi_m)

        Pressure

        Grad_par(Pe) = Te Grad_par(n)

        Potential

        Grad_par(phi) = -Cs Grad_par(Vi)
        -> phi_p - phi_m = -Cs (Vi_p - Vi_m)
       */

      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);

          // Zero gradient Te
          Te(r.ind, mesh->ystart - 1, jz) = Te(r.ind, mesh->ystart, jz);
          BoutReal Cs = sqrt(tesheath); // Sound speed

          // Ion velocity goes to the sound speed
          // Dirichlet boundary condition
          Vi(r.ind, mesh->ystart - 1, jz) =
              -2. * Cs - Vi(r.ind, mesh->ystart, jz);

          BoutReal g = 0.0;
          if (tesheath > 0.1 / Tnorm) {
            // Only divide by Cs if the temperature is greater than 0.1eV
            // to avoid divide-by-zero errors
            g = (Vi(r.ind, mesh->ystart - 1, jz) -
                 Vi(r.ind, mesh->ystart, jz)) /
                (-2. * Cs);
          }

          // Mixed boundary condition on n
          Ne(r.ind, mesh->ystart - 1, jz) =
              Ne(r.ind, mesh->ystart, jz) * (1 - g) / (1 + g);

          // Make sure nesheath doesn't go negative
          Ne(r.ind, mesh->ystart - 1, jz) = floor(
              Ne(r.ind, mesh->ystart - 1, jz), -Ne(r.ind, mesh->ystart, jz));

          // Density at the sheath
          BoutReal nesheath = 0.5 * (Ne(r.ind, mesh->ystart, jz) +
                                     Ne(r.ind, mesh->ystart - 1, jz));

          // Momentum
          NVi(r.ind, mesh->ystart - 1, jz) = NVi(r.ind, mesh->ystart, jz);
          if (NVi(r.ind, mesh->ystart - 1, jz) > 0.0) {
            // Limit flux to be <= 0
            NVi(r.ind, mesh->ystart - 1, jz) = -NVi(r.ind, mesh->ystart, jz);
          }
          //           NVi(r.ind,mesh->ystart+1,jz) =
          //           Ne(r.ind,mesh->ystart+1,jz) *
          //           Vi(r.ind,mesh->ystart+1,jz);

          // Pressure
          Pe(r.ind, mesh->ystart - 1, jz) =
              Pe(r.ind, mesh->ystart, jz) +
              tesheath * (Ne(r.ind, mesh->ystart - 1, jz) -
                          Ne(r.ind, mesh->ystart, jz));

          // Potential
          phi(r.ind, mesh->ystart - 1, jz) =
              phi(r.ind, mesh->ystart, jz) -
              Cs * (Vi(r.ind, mesh->ystart - 1, jz) -
                    Vi(r.ind, mesh->ystart, jz));
          BoutReal phisheath = 0.5 * (phi(r.ind, mesh->ystart, jz) +
                                      phi(r.ind, mesh->ystart - 1, jz));

          // Sheath current
          BoutReal phi_te = floor(phisheath / tesheath, 0.0);
          BoutReal vesheath =
              -Cs * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);

          BoutReal jsheath = nesheath * (-Cs - vesheath);

          // Electron velocity
          Ve(r.ind, mesh->ystart - 1, jz) =
              2. * vesheath - Ve(r.ind, mesh->ystart, jz);

          // Parallel velocity
          Jpar(r.ind, mesh->ystart - 1, jz) =
              2. * jsheath - Jpar(r.ind, mesh->ystart, jz);

          if (currents && !finite(Jpar(r.ind, mesh->ystart - 1, jz))) {
            output.write("JPAR: %d, %d: %e, %e, %e, %e\n", r.ind, jz, jsheath,
                         vesheath, Cs, nesheath);
            output.write(" -> %e, %e, %e\n", Ne(r.ind, mesh->ystart, jz),
                         Ne(r.ind, mesh->ystart - 1, jz), g);
            exit(1);
          }
          // Electron heat conductivity
          kappa_epar(r.ind, mesh->ystart - 1, jz) =
              kappa_epar(r.ind, mesh->ystart, jz);

          // Constant gradient on other cells
          for (int jy = mesh->ystart - 2; jy >= 0; jy--) {
            Vi(r.ind, jy, jz) =
                2. * Vi(r.ind, jy + 1, jz) - Vi(r.ind, jy + 2, jz);
            Ve(r.ind, jy, jz) =
                2. * Ve(r.ind, jy + 1, jz) - Ve(r.ind, jy + 2, jz);
            NVi(r.ind, jy, jz) =
                2. * NVi(r.ind, jy + 1, jz) - NVi(r.ind, jy + 2, jz);

            Ne(r.ind, jy, jz) =
                2. * Ne(r.ind, jy + 1, jz) - Ne(r.ind, jy + 2, jz);
            Te(r.ind, jy, jz) =
                2. * Te(r.ind, jy + 1, jz) - Te(r.ind, jy + 2, jz);
            Pe(r.ind, jy, jz) =
                2. * Pe(r.ind, jy + 1, jz) - Pe(r.ind, jy + 2, jz);

            phi(r.ind, jy, jz) =
                2. * phi(r.ind, jy + 1, jz) - phi(r.ind, jy + 2, jz);
            Vort(r.ind, jy, jz) =
                2. * Vort(r.ind, jy + 1, jz) - Vort(r.ind, jy + 2, jz);
            Jpar(r.ind, jy, jz) =
                2. * Jpar(r.ind, jy + 1, jz) - Jpar(r.ind, jy + 2, jz);
          }
        }
      }
      break;
    }
    case 2: { // Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->ystart, jz) -
                                     Ne(r.ind, mesh->ystart + 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = -sqrt(tesheath); // Sound speed outwards

          if (Vi(r.ind, mesh->ystart, jz) < visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->ystart, jz);
          }

          // Sheath current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->ystart, jz), 0.0);
          BoutReal vesheath =
              -sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);
          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
            // Neumann conditions
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);

            // Dirichlet conditions
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->ystart, jz);
            Pe(r.ind, jy, jz) =
                2. * nesheath * tesheath - Pe(r.ind, mesh->ystart, jz);

            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->ystart, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->ystart, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->ystart, jz);
          }
        }
      }
      break;
    }
    case 3: { // Insulating Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz - 1; jz++) {
          // Free density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->ystart, jz) -
                                     Ne(r.ind, mesh->ystart + 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);

          // Potential for a floating surface
          // BoutReal phisheath = log(0.5*sqrt(mi_me/PI)) * tesheath;

          // Zero-gradient potential
          // NOTE: This should probably not be zero-gradient,
          // since insulator could have electric field across it
          BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = -sqrt(tesheath); // Sound speed outwards

          if (Vi(r.ind, mesh->ystart, jz) < visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->ystart, jz);
          }

          // Sheath current set to zero, as insulating boundary
          BoutReal vesheath = visheath;

          // Apply boundary condition half-way between cells
          for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->ystart, jz);
            phi(r.ind, jy, jz) = 2. * phisheath - phi(r.ind, mesh->ystart, jz);
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->ystart, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->ystart, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->ystart, jz);
            Jpar(r.ind, jy, jz) = -Jpar(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->ystart, jz);
          }
        }
      }
      break;
    }
    }
  } else {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->ngz; jz++) {
        for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
          // Zero-gradient Te
          Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
          Telim(r.ind, jy, jz) = Telim(r.ind, mesh->ystart, jz);

          Pe(r.ind, jy, jz) = Ne(r.ind, jy, jz) * Te(r.ind, jy, jz);
          Pelim(r.ind, jy, jz) = Nelim(r.ind, jy, jz) * Telim(r.ind, jy, jz);

          // Grad_par(phi) = Grad_par(Pe)/Ne
          // phi(r.ind, jy, jz) = phi(r.ind, jy+1, jz) + 2.*( Pe(r.ind, jy, jz)
          // - Pe(r.ind, jy+1, jz) ) / (Ne(r.ind, jy+1, jz) + Ne(r.ind, jy, jz)
          // );

          /*
          phi(r.ind, jy, jz) = phi(r.ind, jy+1, jz) + phi(r.ind, jy+2, jz) -
          phi(r.ind, jy+3, jz);
          Ne(r.ind, jy, jz) = Ne(r.ind, jy+1, jz) + Ne(r.ind, jy+2, jz) -
          Ne(r.ind, jy+3, jz);
          Pe(r.ind, jy, jz) = Pe(r.ind, jy+1, jz) + Pe(r.ind, jy+2, jz) -
          Pe(r.ind, jy+3, jz);
          */
          // if(phi(r.ind, jy, jz) < phi(r.ind, mesh->ystart, jz)) {
          //  phi(r.ind, jy, jz) = phi(r.ind, mesh->ystart, jz);
          //}
        }
      }
    }
  }

  if (sheath_yup) {
    switch (sheath_model) {
    case 0: { // Normal Bohm sheath
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = floor(Ne(r.ind, mesh->yend, jz), 0.0);

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->yend, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = sqrt(tesheath); // Sound speed outwards

          if (Vi(r.ind, mesh->yend, jz) > visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->yend, jz);
          }

          // Sheath current
          // Note that phi/Te >= 0.0 since for phi < 0
          // vesheath is the electron saturation current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->yend, jz), 0.0);
          BoutReal vesheath =
              sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);
          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->yend + 1; jy < mesh->ngy; jy++) {
            // Neumann conditions
            Ne(r.ind, jy, jz) = nesheath;
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->yend, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->yend, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->yend, jz);
            Pelim(r.ind, jy, jz) = Pelim(r.ind, mesh->yend, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->yend, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->yend, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->yend, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->yend, jz);
          }
        }
      }
      break;
    }
    case 1: {
      /*
        Loizu boundary conditions

        Temperature

        Grad_par(Te) = 0

        Density equation

        Grad_par(n) = -(n/Cs) Grad_par(Vi)
        -> n_p - n_m = - (n_p + n_m)/(2Cs) (Vi_p - Vi_m)

        Pressure

        Grad_par(Pe) = Te Grad_par(n)

        Potential

        Grad_par(phi) = -Cs Grad_par(Vi)
        -> phi_p - phi_m = -Cs (Vi_p - Vi_m)
       */

      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);

          // Zero gradient Te
          Te(r.ind, mesh->yend + 1, jz) = Te(r.ind, mesh->yend, jz);
          BoutReal Cs = sqrt(tesheath); // Sound speed

          // Ion velocity goes to the sound speed
          // Dirichlet boundary condition
          Vi(r.ind, mesh->yend + 1, jz) = 2. * Cs - Vi(r.ind, mesh->yend, jz);

          BoutReal g = 0.0;
          if (tesheath > 0.1 / Tnorm) {
            // Only divide by Cs if the temperature is greater than 0.1eV
            // to avoid divide-by-zero errors
            g = (Vi(r.ind, mesh->yend + 1, jz) - Vi(r.ind, mesh->yend, jz)) /
                (2. * Cs);
          }

          // Mixed boundary condition on n
          Ne(r.ind, mesh->yend + 1, jz) =
              Ne(r.ind, mesh->yend, jz) * (1 - g) / (1 + g);

          // Make sure nesheath doesn't go negative
          Ne(r.ind, mesh->yend + 1, jz) =
              floor(Ne(r.ind, mesh->yend + 1, jz), -Ne(r.ind, mesh->yend, jz));
          // Density at the sheath
          BoutReal nesheath =
              0.5 * (Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend + 1, jz));

          // Momentum
          NVi(r.ind, mesh->yend + 1, jz) = NVi(r.ind, mesh->yend, jz);
          if (NVi(r.ind, mesh->yend + 1, jz) < 0.0) {
            // Limit flux to be >= 0
            NVi(r.ind, mesh->yend + 1, jz) = -NVi(r.ind, mesh->yend, jz);
          }
          //           NVi(r.ind,mesh->yend+1,jz) = Ne(r.ind,mesh->yend+1,jz) *
          //           Vi(r.ind,mesh->yend+1,jz);

          // Pressure
          Pe(r.ind, mesh->yend + 1, jz) =
              Pe(r.ind, mesh->yend, jz) +
              tesheath *
                  (Ne(r.ind, mesh->yend + 1, jz) - Ne(r.ind, mesh->yend, jz));

          // Potential
          phi(r.ind, mesh->yend + 1, jz) =
              phi(r.ind, mesh->yend, jz) -
              Cs * (Vi(r.ind, mesh->yend + 1, jz) - Vi(r.ind, mesh->yend, jz));
          BoutReal phisheath = 0.5 * (phi(r.ind, mesh->yend, jz) +
                                      phi(r.ind, mesh->yend + 1, jz));

          // Sheath current
          BoutReal phi_te = floor(phisheath / tesheath, 0.0);
          BoutReal vesheath =
              Cs * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);

          BoutReal jsheath = nesheath * (Cs - vesheath);

          // Electron velocity
          Ve(r.ind, mesh->yend + 1, jz) =
              2. * vesheath - Ve(r.ind, mesh->yend, jz);

          // Parallel velocity
          Jpar(r.ind, mesh->yend + 1, jz) =
              2. * jsheath - Jpar(r.ind, mesh->yend, jz);

          if (currents && !finite(Jpar(r.ind, mesh->yend + 1, jz))) {
            output.write("JPAR: %d, %d: %e, %e, %e, %e\n", r.ind, jz, jsheath,
                         vesheath, Cs, nesheath);
            output.write(" -> %e, %e, %e\n", Ne(r.ind, mesh->yend, jz),
                         Ne(r.ind, mesh->yend + 1, jz), g);
            exit(1);
          }
          // Electron heat conductivity
          kappa_epar(r.ind, mesh->yend + 1, jz) =
              kappa_epar(r.ind, mesh->yend, jz);

          // Constant gradient on other cells
          for (int jy = mesh->yend + 2; jy < mesh->ngy; jy++) {
            Vi(r.ind, jy, jz) =
                2. * Vi(r.ind, jy - 1, jz) - Vi(r.ind, jy - 2, jz);
            Ve(r.ind, jy, jz) =
                2. * Ve(r.ind, jy - 1, jz) - Ve(r.ind, jy - 2, jz);
            NVi(r.ind, jy, jz) =
                2. * NVi(r.ind, jy - 1, jz) - NVi(r.ind, jy - 2, jz);

            Ne(r.ind, jy, jz) =
                2. * Ne(r.ind, jy - 1, jz) - Ne(r.ind, jy - 2, jz);
            Te(r.ind, jy, jz) =
                2. * Te(r.ind, jy - 1, jz) - Te(r.ind, jy - 2, jz);
            Pe(r.ind, jy, jz) =
                2. * Pe(r.ind, jy - 1, jz) - Pe(r.ind, jy - 2, jz);

            phi(r.ind, jy, jz) =
                2. * phi(r.ind, jy - 1, jz) - phi(r.ind, jy - 2, jz);
            Vort(r.ind, jy, jz) =
                2. * Vort(r.ind, jy - 1, jz) - Vort(r.ind, jy - 2, jz);
            Jpar(r.ind, jy, jz) =
                2. * Jpar(r.ind, jy - 1, jz) - Jpar(r.ind, jy - 2, jz);
          }
        }
      }
      break;
    }
    case 2: { // Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->yend, jz) -
                                     Ne(r.ind, mesh->yend - 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->yend, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = sqrt(tesheath); // Sound speed outwards

          if (Vi(r.ind, mesh->yend, jz) > visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->yend, jz);
          }

          // Sheath current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->yend, jz), 0.0);
          BoutReal vesheath =
              sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);
          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->yend + 1; jy < mesh->ngy; jy++) {
            // Neumann conditions
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->yend, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = tesheath;

            // Dirichlet conditions
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->yend, jz);
            Pe(r.ind, jy, jz) =
                2. * nesheath * tesheath - Pe(r.ind, mesh->yend, jz);

            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->yend, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->yend, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->yend, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->yend, jz);
          }
        }
      }
      break;
    }
    case 3: { // Insulating Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->yend, jz) -
                                     Ne(r.ind, mesh->yend - 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);

          // Potential of floating surface
          // BoutReal phisheath = log(0.5*sqrt(mi_me/PI)) * tesheath;
          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->yend, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = sqrt(tesheath); // Sound speed outwards

          if (Vi(r.ind, mesh->yend, jz) > visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->yend, jz);
          }

          // Zero sheath current
          BoutReal vesheath = visheath;

          // Apply boundary condition half-way between cells
          for (int jy = mesh->yend + 1; jy < mesh->ngy; jy++) {
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->yend, jz);
            phi(r.ind, jy, jz) = 2. * phisheath - phi(r.ind, mesh->yend, jz);
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->yend, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->yend, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->yend, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->yend, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->yend, jz);
            Jpar(r.ind, jy, jz) = -Jpar(r.ind, mesh->yend, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->yend, jz);
          }
        }
      }
      break;
    }
    }
  }

  if (!currents) {
    // No currents, so reset Ve to be equal to Vi
    // VePsi also reset, so saved in restart file correctly
    Ve = Vi;
    VePsi = Ve;
  }

  //////////////////////////////////////////////////////////////
  // Fix profiles on lower Y in SOL region by applying
  // a Dirichlet boundary condition.
  // This is to remain consistent with no-flow boundary conditions
  // on velocity fields, and to avoid spurious fluxes of energy
  // through the boundaries.
  //
  // A generator is used (sol_ne, sol_te), and where sol_ne gives a negative
  // value, no boundary condition is applied. This is to allow
  // parts of the domain to be Dirichlet, and parts (e.g. PF) to be Neumann

  if (sol_fix_profiles) {
    TRACE("Fix profiles");
    for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone();
         idwn.next()) {

      BoutReal xnorm = mesh->GlobalX(idwn.ind);
      BoutReal ynorm =
          0.5 * (mesh->GlobalY(mesh->ystart) + mesh->GlobalY(mesh->ystart - 1));

      BoutReal neval = sol_ne->generate(xnorm, TWOPI * ynorm, 0.0, time);
      BoutReal teval = sol_te->generate(xnorm, TWOPI * ynorm, 0.0, time);

      if ((neval < 0.0) || (teval < 0.0))
        continue; // Skip, leave as previous boundary

      for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
        for (int jz = 0; jz < mesh->ngz - 1; jz++) {
          Ne(idwn.ind, jy, jz) = 2. * neval - Ne(idwn.ind, mesh->ystart, jz);
          Te(idwn.ind, jy, jz) = 2. * teval - Te(idwn.ind, mesh->ystart, jz);

          Pe(idwn.ind, jy, jz) = Ne(idwn.ind, jy, jz) * Te(idwn.ind, jy, jz);

          Telim(idwn.ind, jy, jz) = floor(Te(idwn.ind, jy, jz), 0.1 / Tnorm);

          // Zero gradient on Vi to allow flows through boundary
          Vi(idwn.ind, jy, jz) = Vi(idwn.ind, mesh->ystart, jz);

          NVi(idwn.ind, jy, jz) = Vi(idwn.ind, jy, jz) * Ne(idwn.ind, jy, jz);

          // At the boundary, Ve = Vi so no currents
          Ve(idwn.ind, jy, jz) =
              2. * Vi(idwn.ind, jy, jz) - Ve(idwn.ind, mesh->ystart, jz);
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////
  // Density
  // This is the electron density equation
  TRACE("density");

  ddt(Ne) = -Div_n_bxGrad_f_B_XPPM(Ne, phi, ne_bndry_flux, poloidal_flows,
                                   true); // ExB drift

  // Parallel flow
  if (parallel_flow) {
    if (currents) {
      // Parallel wave speed increased to electron sound speed
      ddt(Ne) -= Div_par_FV_FS(Ne, Ve, sqrt(mi_me) * sound_speed);
    } else {
      // Parallel wave speed is ion sound speed
      ddt(Ne) -= Div_par_FV_FS(Ne, Ve, sound_speed);
    }
  }

  if (j_diamag) {
    // Diamagnetic drift, formulated as a magnetic drift

    ddt(Ne) -= Div_f_v_XPPM(Ne, -Telim * Curlb_B,
                            ne_bndry_flux); // Grad-B, curvature drift
  }

  if (ramp_mesh && (time < ramp_timescale)) {
    ddt(Ne) += NeTarget / ramp_timescale;
  }

  if (classical_diffusion) {
    // Classical perpendicular diffusion

    Dn = (Telim + Ti) / (tau_e * mi_me * SQ(mesh->Bxy));
    ddt(Ne) += Div_Perp_Lap_FV(Dn, Ne, ne_bndry_flux);
    ddt(Ne) -= Div_Perp_Lap_FV(0.5 * Ne / (tau_e * mi_me * SQ(mesh->Bxy)), Te,
                               ne_bndry_flux);
  }
  if (anomalous_D > 0.0) {
    ddt(Ne) += Div_Perp_Lap_FV(anomalous_D, Ne.DC(), ne_bndry_flux);
    // ddt(Ne) += Div_Perp_Lap_XYZ(anomalous_D, Ne, ne_bndry_flux);
  }

  // Source
  if (adapt_source) {
    Field2D NeErr = averageY(Ne.DC() - NeTarget);

    if (adapt_fix_form) {
      // Fix the form of the source function

      // Average the error by integrating over the domain
      // weighted by the source function
      BoutReal local_error = 0.0;
      for (int i = 0; i < mesh->ngx; i++) {
        for (int j = 0; j < mesh->ngy; j++) {
          local_error += NeErr(i, j) * Sn(i, j);
        }
      }
      BoutReal error;
      MPI_Allreduce(&local_error, &error, 1, MPI_DOUBLE, MPI_SUM,
                    BoutComm::get());
      error /= total_Sn;

      // All processors now share the same error

      // PI controller, using crude integral of the error
      if (density_error_lasttime < 0.0) {
        // First time
        density_error_lasttime = time;
        density_error_last = error;
      }
      // Integrate using Trapezium rule
      if (time > density_error_lasttime) { // Since time can decrease
        density_error_integral += (time - density_error_lasttime) * 0.5 *
                                  (error + density_error_last);

        density_error_last = error;
        density_error_lasttime = time;
      }
      if (density_error_integral > 0.0) {
        density_error_integral = 0.0;
      }

      NeSource = Sn * (-source_p * error - source_i * density_error_integral);

    } else {
      if (core_sources) {
        // Sources only in core (periodic Y) domain
        // Try to keep NeTarget

        ddt(Sn) = 0.0;
        for (int x = mesh->xstart; x <= mesh->xend; x++) {
          if (!mesh->periodicY(x))
            continue; // Not periodic, so skip

          for (int y = mesh->ystart; y <= mesh->yend; y++) {
            Sn(x, y) -= source_p * NeErr(x, y);
            ddt(Sn)(x, y) = -source_i * NeErr(x, y);

            if (Sn(x, y) < 0.0) {
              Sn(x, y) = 0.0;
              if (ddt(Sn)(x, y) < 0.0)
                ddt(Sn)(x, y) = 0.0;
            }
          }
        }

        NeSource = Sn;
      } else {
        // core_sources = false
        NeSource = Sn * where(Sn, NeTarget, Ne);
        NeSource -= source_p * NeErr / NeTarget;

        ddt(Sn) = -source_i * NeErr;
      }
    }
  } else {
    // Add source. Ensure that sink will go to zero as Ne -> 0

    NeSource = Sn * where(Sn, 1.0, Ne);
  }

  if (source_vary_g11) {
    NeSource *= g11norm;
  }

  ddt(Ne) += NeSource;

  if (ExBdiff > 0.0) {

    if (ExBpar) {
      ddt(Ne) += ExBdiff *
                 Div_Perp_Lap_XYZ(SQ(mesh->dx) * mesh->g_11, Ne, ne_bndry_flux);
    } else {
      ddt(Ne) += ExBdiff * Div_Perp_Lap_FV_Index(1.0, Ne, ne_bndry_flux);
    }
  }

  if (ADpar > 0.0) {
    // ddt(Ne) += ADpar * AddedDissipation(1.0, Pe, Ne, false);

    ddt(Ne) += ADpar * AddedDissipation(Ne, Pe, Nelim, ADpar_bndry);

    if (ADpar_phine) {
      ddt(Ne) -= ADpar * AddedDissipation(Ne, phi, Nelim, ADpar_bndry);
    } else {
      ddt(Ne) -= ADpar * AddedDissipation(1.0, phi, Ne, ADpar_bndry);
    }
  }

  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions
    ddt(Ne) +=
        Div_par_diffusion(SQ(mesh->dy) * mesh->g_22 * 1e-4 / Nelim, Ne, false);
  }
  if (low_n_diffuse_perp) {
    ddt(Ne) += Div_Perp_Lap_FV_Index(1e-4 / Nelim, Ne, ne_bndry_flux);
  }

  if (ne_hyper_z > 0.) {
    ddt(Ne) -= ne_hyper_z * SQ(SQ(mesh->dz)) * D4DZ4(Ne);
  }

  ///////////////////////////////////////////////////////////
  // Vorticity
  // This is the current continuity equation

  TRACE("vorticity");

  ddt(Vort) = 0.0;

  if (currents) {
    // Only evolve vorticity if any diamagnetic or parallel currents
    // are included.

    if (j_par) {
      // Parallel current
      // ddt(Vort) += Div_par(Jpar);
      ddt(Vort) += 0.5 * (Div_par(Jpar) + Ne * Div_par(Vi - Ve) +
                          (Vi - Ve) * Grad_par(Ne));
    }

    if (j_diamag) {
      // Electron diamagnetic current
      // Note: This term is central differencing so that it balances
      // the corresponding compression term in the pressure equation

      ddt(Vort) += Div(Pe * Curlb_B);
    }

    if (rotation) {
      // Including solid-body rotation in toroidal direction
      // Adds a current because Coriolis and centrifugal drifts
      // depend on charge and mass

      ddt(Vort) += Div_f_v_XPPM(
          Ne,
          -SQ(rotation_rate) * bxGradR / mesh->Bxy // centrifugal force
              + 2. * Vi * Omega_vec / mesh->Bxy    // Coriolis force
          ,
          vort_bndry_flux);
    }

    // Advection of vorticity by ExB
    if (boussinesq) {
      // If the Boussinesq approximation is made then the
      // form of this becomes simple, similar to the terms in
      // the density and pressure equations
      ddt(Vort) -=
          Div_n_bxGrad_f_B_XPPM(Vort, phi, vort_bndry_flux, poloidal_flows);
    } else {
      // When the Boussinesq approximation is not made,
      // then the changing ion density introduces a number
      // of other terms.

      Field3D tmp = Ne * Delp2(phi) / mesh->Bxy;
      mesh->communicate(tmp);
      tmp.applyBoundary("neumann");

      ddt(Vort) -= 0.5 * Div_n_bxGrad_f_B_XPPM(Vort + tmp, phi, vort_bndry_flux,
                                               poloidal_flows);

      // Ion density
      Field3D dnidt = -Div_n_bxGrad_f_B_XPPM(Ne, phi); // ExB drift of ions
      // dnidt -= Div_parP_LtoC(Ne,Vi); // Parallel flow

      mesh->communicate(dnidt);
      dnidt.applyBoundary("neumann");
      ddt(Vort) +=
          Div_Perp_Lap_FV((0.5 * dnidt) / SQ(mesh->Bxy), phi, vort_bndry_flux);
    }

    if (classical_diffusion) {
      /// Field3D mu = 0.75*(1.+1.6*SQ(neoclassical_q))/tau_i;
      Field3D mu = 0.3 * Ti / (tau_i * SQ(mesh->Bxy));
      ddt(Vort) += Div_Perp_Lap_FV(mu, Vort, vort_bndry_flux);
    }

    if (anomalous_nu > 0.0) {
      // Perpendicular anomalous momentum diffusion
      ddt(Vort) += Div_Perp_Lap_FV(anomalous_nu, Vort.DC(), vort_bndry_flux);
    }

    // Sink of vorticity due to ion-neutral friction
    // ddt(Vort) += Sn*where(Sn, 0.0, Vort);
    if (ion_neutral > 0.0)
      ddt(Vort) -= ion_neutral * Vort;

    if (Dvort > 0.0)
      ddt(Vort) += Dvort * Div_Perp_Lap_FV(1.0, Vort, vort_bndry_flux);

    if (ExBdiff > 0.0) {

      if (ExBpar) {
        ddt(Vort) += ExBdiff * Div_Perp_Lap_XYZ(SQ(mesh->dx) * mesh->g_11, Vort,
                                                vort_bndry_flux);
        // ddt(Vort) += ExBdiff * Div_par_diffusion(SQ(mesh->dy)*mesh->g_22,
        // Vort);
      } else {
        ddt(Vort) +=
            ExBdiff * Div_Perp_Lap_FV_Index(1.0, Vort, vort_bndry_flux);
      }
    }

    if (ADpar > 0.0) {
      if (ADpar_phine) {
        ddt(Vort) -= ADpar * AddedDissipation(Ne, phi, Nelim, ADpar_bndry);
      } else {
        ddt(Vort) -= ADpar * AddedDissipation(1.0, phi, Ne, ADpar_bndry);
      }
    }
    if (z_hyper_viscos > 0) {
      // Form of hyper-viscosity to suppress zig-zags in Z
      ddt(Vort) -= z_hyper_viscos * SQ(SQ(mesh->dz)) * D4DZ4(Vort);
    }
    if (x_hyper_viscos > 0) {
      // Form of hyper-viscosity to suppress zig-zags in X
      ddt(Vort) -= x_hyper_viscos * D4DX4_FV_Index(Vort);
    }

    if (y_hyper_viscos > 0) {
      // Form of hyper-viscosity to suppress zig-zags in Y
      ddt(Vort) -= y_hyper_viscos * D4DY4_FV_Index(Vort, false);
    }
  }

  ///////////////////////////////////////////////////////////
  // Ohm's law
  // VePsi = Ve - Vi + 0.5*mi_me*beta_e*psi
  TRACE("Ohm's law");

  ddt(VePsi) = 0.0;

  Field3D NelimVe = Nelim;

  if (currents && (electromagnetic || FiniteElMass)) {
    // Evolve VePsi except for electrostatic and zero electron mass case

    if (resistivity) {
      ddt(VePsi) -= mi_me * nu * (Ve - Vi);
    }

    // Parallel electric field
    if (j_par) {
      ddt(VePsi) += mi_me * Grad_parP_CtoL(phi);
    }

    // Parallel electron pressure
    if (pe_par) {
      // ddt(VePsi) -= mi_me*(  Te*Grad_parP_CtoL(log(Ne)) +
      // 1.71*Grad_parP_CtoL(Te) );
      ddt(VePsi) -= mi_me * Grad_parP_CtoL(Pelim) / NelimVe;
    }

    if (thermal_force) {
      ddt(VePsi) -= mi_me * 0.71 * Grad_parP_CtoL(Te);
    }

    if (electron_viscosity) {
      // Electron parallel viscosity (Braginskii)
      Field3D ve_eta = 0.973 * mi_me * tau_e * Telim;

      if (flux_limit_alpha > 0) {
        // Limit to free streaming value
        Field3D ve_eta_fs = flux_limit_alpha * sqrt(mi_me * Telim) * R0;
        ve_eta = (ve_eta * ve_eta_fs) / (ve_eta + ve_eta_fs);
      }
      if (eta_limit_alpha > 0.) {
        // SOLPS-style flux limiter
        // Values of alpha ~ 0.5 typically
        Field3D q_cl = ve_eta * Grad_par(Ve);           // Collisional value
        Field3D q_fl = eta_limit_alpha * Pelim * mi_me; // Flux limit

        ve_eta = ve_eta / (1. + abs(q_cl / q_fl));

        mesh->communicate(ve_eta);
        ve_eta.applyBoundary("neumann");
      }
      ddt(VePsi) += Div_par_diffusion(ve_eta, Ve);
    }

    if (FiniteElMass) {
      // Finite Electron Mass. Small correction needed to conserve energy
      ddt(VePsi) -= Vi * Grad_par(Ve - Vi); // Parallel advection
      ddt(VePsi) -= bracket(phi, Ve - Vi);  // ExB advection
      // Should also have ion polarisation advection here
    }

    if (numdiff > 0.0) {
      ddt(VePsi) += sqrt(mi_me) * numdiff * Div_par_diffusion_index(Ve);
    }

    if (hyper > 0.0) {
      ddt(VePsi) -= hyper * mi_me * nu * Delp2(Jpar) / Nelim; 
    }

    if (hyperpar > 0.0) {
      ddt(VePsi) -= D4DY4_FV(SQ(SQ(mesh->dy)), Ve - Vi);
    }

    if (vepsi_dissipation) {
      // Adds dissipation term like in other equations
      // Maximum speed either electron sound speed or Alfven speed
      Field3D max_speed = Bnorm*mesh->Bxy / sqrt(SI::mu0* AA*SI::Mp*Nnorm*Nelim) / Cs0; // Alfven speed (normalised by Cs0)
      Field3D elec_sound = sqrt(mi_me)*sound_speed; // Electron sound speed
      for (int jx=0;jx<mesh->ngx;jx++)
        for (int jy=0;jy<mesh->ngy;jy++)
          for (int jz=0;jz<mesh->ngz;jz++) {
            if (elec_sound(jx, jy, jz) > max_speed(jx, jy, jz)) {
              max_speed(jx, jy, jz) = elec_sound(jx, jy, jz);
            }
      }
      
      ddt(VePsi) -= Div_par_FV_FS(Ve-Vi, 0.0, max_speed);
    }
  }

  ///////////////////////////////////////////////////////////
  // Ion velocity
  if (ion_velocity) {
    TRACE("Ion velocity");

    ddt(NVi) = -Div_n_bxGrad_f_B_XPPM(NVi, phi, ne_bndry_flux, poloidal_flows); // ExB drift
    
    ddt(NVi) -= Div_par_FV_FS(NVi, Vi, sound_speed, false); // Parallel flow

    // Ignoring polarisation drift for now
    if (pe_par) {
      ddt(NVi) -= Grad_parP_CtoL(Pe);
    }

    if (rotation) {
      // Including solid-body rotation in toroidal direction
      // Adds Coriolis and centrifugal drifts

      ddt(NVi) -= Div_f_v_XPPM(
          NVi,
          -SQ(rotation_rate) * bxGradR / mesh->Bxy // centrifugal force
              + 2. * Vi * Omega_vec / mesh->Bxy    // Coriolis force
          ,
          ne_bndry_flux);
    }

    // Ion-neutral friction
    if (ion_neutral > 0.0)
      ddt(NVi) -= ion_neutral * NVi;

    // Perpendicular viscosity
    if (Dvi > 0.0)
      ddt(NVi) += Dvi * Div_Perp_Lap_x3(Ne, Vi, false);

    if (numdiff > 0.0) {
      ddt(NVi) += numdiff * Div_par_diffusion_index(Vi);
      // ddt(NVi) += Div_par_diffusion(SQ(mesh->dy)*mesh->g_22*numdiff, Vi);
    }

    if (density_inflow) {
      // Particles arrive in cell at rate NeSource
      // This should come from a flow through the cell edge
      // with a flow velocity, and hence momentum

      ddt(NVi) += NeSource * (NeSource / Ne) * mesh->dy * sqrt(mesh->g_22);
    }

    if (classical_diffusion) {
      ddt(NVi) += Div_Perp_Lap_FV(Vi * Dn, Ne, ne_bndry_flux);
    }
    if (anomalous_D > 0.0) {
      ddt(NVi) +=
          Div_Perp_Lap_FV(Vi.DC() * anomalous_D, Ne.DC(), ne_bndry_flux);
      // ddt(NVi) += Div_Perp_Lap_XYZ(Vi*anomalous_D, Ne, ne_bndry_flux);
    }
    if (anomalous_nu > 0.0) {
      // Perpendicular anomalous momentum diffusion
      ddt(NVi) += Div_Perp_Lap_FV(anomalous_nu * Ne.DC(), Vi.DC(), ne_bndry_flux);
    }
    
    if (ExBdiff > 0.0) {

      if (ExBpar) {
        ddt(NVi) += ExBdiff * Div_Perp_Lap_XYZ(SQ(mesh->dx) * mesh->g_11, NVi,
                                               ne_bndry_flux);
      } else {
        ddt(NVi) += ExBdiff * Div_Perp_Lap_FV_Index(1.0, NVi, ne_bndry_flux);
      }
    }

    if (ADpar > 0.0) {
      ddt(NVi) += ADpar * AddedDissipation(Ne, Pe, NVi, ADpar_bndry);
    }

    if (hyperpar > 0.0) {
      ddt(NVi) -= D4DY4_FV(SQ(SQ(mesh->dy)), Vi) / mi_me;
    }

    if (low_n_diffuse) {
      // Diffusion which kicks in at very low density, in order to
      // help prevent negative density regions
      ddt(NVi) += Div_par_diffusion(
          Vi * SQ(mesh->dy) * mesh->g_22 * 1e-4 / Nelim, Ne, true);
    }
  }

  ///////////////////////////////////////////////////////////
  // Pressure equation
  TRACE("Electron pressure");
  ddt(Pe) = 0.0;

  // Divergence of heat flux due to ExB advection
  ddt(Pe) -=
      Div_n_bxGrad_f_B_XPPM(Pe, phi, pe_bndry_flux, poloidal_flows, true);

  if (parallel_flow) {
    if (currents) {
      // Like Ne term, parallel wave speed increased
      ddt(Pe) -= Div_par_FV_FS(Pe, Ve, sqrt(mi_me) * sound_speed);
    } else {
      ddt(Pe) -= Div_par_FV_FS(Pe, Ve, sound_speed);
    }
  }

  if (j_diamag) { // Diamagnetic flow
    // Magnetic drift (curvature) divergence.
    ddt(Pe) -= (5. / 3) * Div_f_v_XPPM(Pe, -Telim * Curlb_B, pe_bndry_flux);

    // This term energetically balances diamagnetic term
    // in the vorticity equation
    ddt(Pe) -= (2. / 3) * Pe * (Curlb_B * Grad(phi));
  }

  // Parallel heat conduction
  if (thermal_conduction) {
    ddt(Pe) += (2. / 3) * Div_par_diffusion(kappa_epar, Te);
  }

  if (thermal_flux) {
    // Parallel heat convection
    ddt(Pe) += (2. / 3) * 0.71 * Div_parP_LtoC(Te, Jpar);
  }

  if (currents && resistivity) {
    // Ohmic heating
    ddt(Pe) += nu * Jpar * (Jpar - Jpar0) / Nelim;
  }

  if (pe_hyper_z > 0.0) {
    ddt(Pe) -= pe_hyper_z * SQ(SQ(mesh->dz)) * D4DZ4(Pe);
  }

  ///////////////////////////////////
  // Heat transmission through sheath
  wall_power = 0.0; // Diagnostic output
  if (sheath_yup) {
    TRACE("sheath yup heat transmission");
    switch (sheath_model) {
    case 0:
    case 2:
    case 3: {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Temperature and density at the sheath entrance
          BoutReal tesheath = floor(
              0.5 * (Te(r.ind, mesh->yend, jz) + Te(r.ind, mesh->yend + 1, jz)),
              0.0);
          BoutReal nesheath = floor(
              0.5 * (Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend + 1, jz)),
              0.0);

          // Sound speed (normalised units)
          BoutReal Cs = sqrt(tesheath);

          // Heat flux
          BoutReal q = (sheath_gamma - 3) * tesheath * nesheath * Cs;

          // Multiply by cell area to get power
          BoutReal flux = q * (mesh->J(r.ind, mesh->yend) +
                               mesh->J(r.ind, mesh->yend + 1)) /
                          (sqrt(mesh->g_22(r.ind, mesh->yend)) +
                           sqrt(mesh->g_22(r.ind, mesh->yend + 1)));

          // Divide by volume of cell, and 2/3 to get pressure
          BoutReal power =
              flux / (mesh->dy(r.ind, mesh->yend) * mesh->J(r.ind, mesh->yend));
          ddt(Pe)(r.ind, mesh->yend, jz) -= (2. / 3) * power;
          wall_power(r.ind, mesh->yend) += power;
        }
      }
      break;
    }
    }
  }
  if (sheath_ydown) {
    TRACE("sheath ydown heat transmission");
    switch (sheath_model) {
    case 0:
    case 2:
    case 3: {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->ngz; jz++) {
          // Temperature and density at the sheath entrance
          BoutReal tesheath = floor(0.5 * (Te(r.ind, mesh->ystart, jz) +
                                           Te(r.ind, mesh->ystart - 1, jz)),
                                    0.0);
          BoutReal nesheath = floor(0.5 * (Ne(r.ind, mesh->ystart, jz) +
                                           Ne(r.ind, mesh->ystart - 1, jz)),
                                    0.0);

          // Sound speed (normalised units)
          BoutReal Cs = sqrt(tesheath);

          // Heat flux
          BoutReal q =
              (sheath_gamma - 3) * tesheath * nesheath * Cs; // NB: positive

          // Multiply by cell area to get power
          BoutReal flux = q * (mesh->J(r.ind, mesh->ystart) +
                               mesh->J(r.ind, mesh->ystart - 1)) /
                          (sqrt(mesh->g_22(r.ind, mesh->ystart)) +
                           sqrt(mesh->g_22(r.ind, mesh->ystart - 1)));

          // Divide by volume of cell, and 2/3 to get pressure
          BoutReal power = flux / (mesh->dy(r.ind, mesh->ystart) *
                                   mesh->J(r.ind, mesh->ystart));
          ddt(Pe)(r.ind, mesh->ystart, jz) -= (2. / 3) * power;
          wall_power(r.ind, mesh->ystart) += power;
        }
      }
      break;
    }
    }
  }
  // Transfer and source terms
  if (thermal_force) {
    ddt(Pe) -= (2. / 3) * 0.71 * Jpar * Grad_parP(Te);
  }

  if (pe_par) {
    // This term balances energetically the pressure term
    // in Ohm's law
    ddt(Pe) -= (2. / 3) * Pelim * Div_par(Ve);
  }
  if (ramp_mesh && (time < ramp_timescale)) {
    ddt(Pe) += PeTarget / ramp_timescale;
  }

  //////////////////////
  // Classical diffusion

  if (classical_diffusion) {
    // nu_rho2 = nu_ei * rho_e^2 in normalised units
    Field3D nu_rho2 = Telim / (tau_e * mi_me * SQ(mesh->Bxy));

    ddt(Pe) += (2. / 3) * (Div_Perp_Lap_FV(nu_rho2, Pe, pe_bndry_flux) +
                           (11. / 12) * Div_Perp_Lap_FV(nu_rho2 * Nelim, Te,
                                                        pe_bndry_flux));
  }

  //////////////////////
  // Anomalous diffusion

  if (anomalous_D > 0.0) {
    ddt(Pe) += Div_Perp_Lap_FV(anomalous_D * Te.DC(), Ne.DC(), ne_bndry_flux);
    // ddt(Pe) += Div_Perp_Lap_XYZ(anomalous_D*Te, Ne, ne_bndry_flux);
  }
  if (anomalous_chi > 0.0) {
    // ddt(Pe) += Div_Perp_Lap_XYZ(anomalous_chi*Ne, Te, pe_bndry_flux);
    ddt(Pe) += Div_Perp_Lap_FV(anomalous_chi * Ne.DC(), Te.DC(), pe_bndry_flux);
  }

  if (Dte > 0.0)
    ddt(Pe) += Dte * Div_Perp_Lap_x3(1.0, Te, pe_bndry_flux);

  //////////////////////
  // Sources

  if (adapt_source) {
    // Add source. Ensure that sink will go to zero as Pe -> 0
    Field2D PeErr;

    if (!temperature_feedback) {
      // Feedback on the pressure. This is the default
      PeErr = averageY(Pe.DC() - PeTarget);
    } else {
      // Feedback on temperature
      PeErr = averageY(Te.DC() - TeTarget);
    }

    if (adapt_fix_form) {
      // Fix the form of the source function

      // Average the error by integrating over the domain
      // weighted by the source function
      BoutReal local_error = 0.0;
      for (int i = 0; i < mesh->ngx; i++) {
        for (int j = 0; j < mesh->ngy; j++) {
          local_error += PeErr(i, j) * Spe(i, j);
        }
      }
      BoutReal error;
      MPI_Allreduce(&local_error, &error, 1, MPI_DOUBLE, MPI_SUM,
                    BoutComm::get());
      error /= total_Spe;

      // All processors now share the same error

      // PI controller, using crude integral of the error
      if (pe_error_lasttime < 0.0) {
        // First time
        pe_error_lasttime = time;
        pe_error_last = error;
      }
      // Integrate using Trapezium rule
      if (time > pe_error_lasttime) { // Since time can decrease
        pe_error_integral +=
            (time - pe_error_lasttime) * 0.5 * (error + pe_error_last);

        pe_error_lasttime = time;
        pe_error_last = error;
      }
      if (pe_error_integral > 0.0) {
        pe_error_integral = 0.0;
      }

      PeSource = Spe * (-source_p * error - source_i * pe_error_integral);

    } else {

      if (core_sources) {
        // Sources only in core

        ddt(Spe) = 0.0;
        for (int x = mesh->xstart; x <= mesh->xend; x++) {
          if (!mesh->periodicY(x))
            continue; // Not periodic, so skip

          for (int y = mesh->ystart; y <= mesh->yend; y++) {
            Spe(x, y) -= source_p * PeErr(x, y);
            ddt(Spe)(x, y) = -source_i * PeErr(x, y);

            if (Spe(x, y) < 0.0) {
              Spe(x, y) = 0.0;
              if (ddt(Spe)(x, y) < 0.0)
                ddt(Spe)(x, y) = 0.0;
            }
          }
        }

        if (energy_source) {
          // Add the same amount of energy to each particle
          PeSource = Spe * Nelim / Nelim.DC();
        } else {
          PeSource = Spe;
        }
      } else {

        Spe -= source_p * PeErr / PeTarget;
        ddt(Spe) = -source_i * PeErr;

        if (energy_source) {
          // Add the same amount of energy to each particle
          PeSource = Spe * Nelim / Nelim.DC();
        } else {
          PeSource = Spe * where(Spe, PeTarget, Pe);
        }
      }
    }
  } else {
    // Not adapting sources

    if (energy_source) {
      // Add the same amount of energy to each particle
      PeSource = Spe * Nelim / Nelim.DC();
    } else {
      // Add the same amount of energy per volume
      // If no particle source added, then this can lead to
      // a small number of particles with a lot of energy!
      PeSource = Spe * where(Spe, 1.0, Pe);
    }
  }

  if (source_vary_g11) {
    PeSource *= g11norm;
  }

  ddt(Pe) += PeSource;

  //////////////////////
  // Numerical dissipation

  if (ExBdiff > 0.0) {
    if (ExBpar) {
      ddt(Pe) += ExBdiff *
                 Div_Perp_Lap_XYZ(SQ(mesh->dx) * mesh->g_11, Pe, pe_bndry_flux);
    } else {
      ddt(Pe) += ExBdiff * Div_Perp_Lap_FV_Index(1.0, Pe, pe_bndry_flux);
    }
  }

  if (ADpar > 0.0) {
    // ddt(Pe) += ADpar * AddedDissipation(1.0, Pe, Pe, false);

    ddt(Pe) += ADpar * AddedDissipation(Ne, Pe, Pelim, ADpar_bndry);
    if (ADpar_phine) {
      ddt(Pe) -= ADpar * AddedDissipation(Ne, phi, Telim * Nelim, ADpar_bndry);
    } else {
      ddt(Pe) -= ADpar * AddedDissipation(1.0, phi, Pe, ADpar_bndry);
    }

    ddt(Pe) += ADpar * AddedDissipation(1.0, Te, 1.0, ADpar_bndry);
  }

  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions
    ddt(Pe) += Div_par_diffusion(Te * SQ(mesh->dy) * mesh->g_22 * 1e-4 / Nelim,
                                 Ne, false);
  }

  ///////////////////////////////////////////////////////////
  // Radial buffer regions for turbulence simulations

  if (radial_buffers) {
    /// Radial buffer regions

    // Calculate Z averages
    Field2D PeDC = Pe.DC();
    Field2D NeDC = Ne.DC();
    Field2D VortDC = Vort.DC();
    
    if ((mesh->XGLOBAL(mesh->xstart) - mesh->xstart) < radial_inner_width) {
      // This processor contains points inside the inner radial boundary
      
      int imax = mesh->xstart + radial_inner_width - 1
        - (mesh->XGLOBAL(mesh->xstart) - mesh->xstart);
      if (imax > mesh->xend) {
        imax = mesh->xend;
      }
      
      int imin =mesh->xstart;
      if (!mesh->firstX()) {
        --imin; // Calculate in guard cells, for radial fluxes
      }
      int ncz = mesh->ngz - 1;
      
      for (int i = imin; i <= imax; ++i) {
        // position inside the boundary (0 = on boundary, 0.5 = first cell)
        BoutReal pos = static_cast<BoutReal>(mesh->XGLOBAL(i) - mesh->xstart) + 0.5;
        
        // Diffusion coefficient which increases towards the boundary
        BoutReal D = radial_buffer_D * (1. - pos / radial_inner_width);
        
        for (int j = mesh->ystart; j <= mesh->yend; ++j) {
          BoutReal dx = mesh->dx(i,j);
          BoutReal dx_xp = mesh->dx(i+1,j);
          BoutReal J = mesh->J(i,j);
          BoutReal J_xp = mesh->J(i+1,j);

          // Calculate metric factors for radial fluxes
          BoutReal rad_flux_factor = 0.25*(J + J_xp)*(dx + dx_xp);
          BoutReal x_factor = rad_flux_factor / (J * dx);
          BoutReal xp_factor = rad_flux_factor / (J_xp * dx_xp);
          
          for (int k = 0; k < ncz; ++k) {
            // Relax towards constant value on flux surface
            ddt(Pe)(i,j,k) -= D*(Pe(i,j,k) - PeDC(i,j));
            ddt(Ne)(i,j,k) -= D*(Ne(i,j,k) - NeDC(i,j));
            ddt(Vort)(i,j,k) -= D*(Vort(i,j,k) - VortDC(i,j));

            // Radial fluxes
            BoutReal f = D * (Ne(i + 1, j, k) - Ne(i, j, k));
            ddt(Ne)(i, j, k) += f * x_factor;
            ddt(Ne)(i + 1, j, k) -= f * xp_factor;

            f = D * (Pe(i + 1, j, k) - Pe(i, j, k));
            ddt(Pe)(i, j, k) += f * x_factor;
            ddt(Pe)(i + 1, j, k) -= f * xp_factor;

            f = D * (Vort(i + 1, j, k) - Vort(i, j, k));
            ddt(Vort)(i, j, k) += f * x_factor;
            ddt(Vort)(i + 1, j, k) -= f * xp_factor;
          }
        }
      }
    }
    // Number of points in outer guard cells
    int nguard = mesh->ngx-mesh->xend-1;
    
    if (mesh->GlobalNx - nguard - mesh->XGLOBAL(mesh->xend) <= radial_outer_width) {
      
      // Outer boundary
      int imin = mesh->GlobalNx - nguard - radial_outer_width - mesh->XGLOBAL(0);
      if (imin < mesh->xstart) {
        imin = mesh->xstart;
      }
      int ncz = mesh->ngz - 1;
      for (int i = imin; i <= mesh->xend; ++i) {

        // position inside the boundary
        BoutReal pos = static_cast<BoutReal>(mesh->GlobalNx - nguard -  mesh->XGLOBAL(i)) - 0.5;
        
        // Diffusion coefficient which increases towards the boundary
        BoutReal D = radial_buffer_D * (1. - pos / radial_outer_width);
        
        for (int j = mesh->ystart; j <= mesh->yend; ++j) {
          BoutReal dx = mesh->dx(i,j);
          BoutReal dx_xp = mesh->dx(i+1,j);
          BoutReal J = mesh->J(i,j);
          BoutReal J_xp = mesh->J(i+1,j);

          // Calculate metric factors for radial fluxes
          BoutReal rad_flux_factor = 0.25*(J + J_xp)*(dx + dx_xp);
          BoutReal x_factor = rad_flux_factor / (J * dx);
          BoutReal xp_factor = rad_flux_factor / (J_xp * dx_xp);
          
          for (int k = 0; k < ncz; ++k) {
            ddt(Pe)(i,j,k) -= D*(Pe(i,j,k) - PeDC(i,j));
            ddt(Ne)(i,j,k) -= D*(Ne(i,j,k) - NeDC(i,j));
            ddt(Vort)(i,j,k) -= D*(Vort(i,j,k) - VortDC(i,j));
            //ddt(Vort)(i,j,k) -= D*Vort(i,j,k);
            
            BoutReal f = D * (Vort(i + 1, j, k) - Vort(i, j, k));
            ddt(Vort)(i, j, k) += f * x_factor;
            ddt(Vort)(i + 1, j, k) -= f * xp_factor;
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////
  // Neutral gas
  if (neutral_model > 0) {
    TRACE("Neutral gas model");
    // Include a neutral gas model

    // Set 3D variables Nn to neutral density,
    // Tn to neutral temperature
    // Nelim to electron density
    switch (neutral_model) {
    case 1: {
      Nn = floor(Nn, 1e-8);
      Tn = Pn / Nn;
      Tn = floor(Tn, 0.01 / Tnorm);
      break;
    }
    case 2: {
      // Lower limit for neutral density (mainly for first time through when Nn
      // = 0)
      Nn = floor(Nn, 1e-8);
      // Franck-Condon energy for Tn
      break;
    }
    case 3: {
      // Navier-Stokes for axisymmetric neutral gas profiles
      // Nn2D, Pn2D and Tn2D are unfloored
      Tn2D = Pn2D / Nn2D;

      // Nn and Tn are 3D floored fields, used for rate coefficients
      Nn = floor(Nn2D, 1e-8);
      Tn = Pn2D / Nn;
      Tn = floor(Tn, 0.01 / Tnorm);
      break;
    }
      // Case 4 handled later
    }

    Nelim = floor(Ne, 1e-19); // Smaller limit for rate coefficients

    // Calculate atomic rates

    switch (neutral_model) {
    case 1: {
      TRACE("neutral model 1");
      // Neutral gas dynamics

      // Set a floor of room temperature

      // Calculate atomic processes
      for (int i = 0; i < mesh->ngx; i++)
        for (int j = 0; j < mesh->ngy; j++)
          for (int k = 0; k < mesh->ngz; k++) {
            // Charge exchange frequency, normalised to ion cyclotron frequency
            BoutReal sigma_cx =
                Nelim(i, j, k) * Nnorm *
                hydrogen.chargeExchange(Telim(i, j, k) * Tnorm) / Omega_ci;

            // Ionisation frequency, normalised to ion cyclotron frequency
            BoutReal sigma_iz = Nelim(i, j, k) * Nnorm * Nn(i, j, k) *
                                hydrogen.ionisation(Telim(i, j, k) * Tnorm) /
                                Omega_ci;

            // Neutral thermal velocity
            BoutReal vth_n = sqrt(Tn(i, j, k)); // Normalised to Cs0

            // Neutral-neutral mean free path
            BoutReal a0 = PI * SQ(5.29e-11);
            BoutReal lambda_nn = 1. / (Nnorm * Nn(i, j, k) * a0); // meters
            if (lambda_nn > Lmax) {
              // Limit maximum mean free path
              lambda_nn = Lmax;
            }

            lambda_nn /= rho_s0; // Normalised length to rho_s0
            // Neutral-Neutral collision rate, normalised to ion cyclotron
            // frequency
            BoutReal sigma_nn = vth_n / lambda_nn;

            // Total neutral collision frequency, normalised to ion cyclotron
            // frequency
            BoutReal sigma = sigma_cx + sigma_nn + sigma_iz;

            // Neutral gas diffusion
            Dnn(i, j, k) = SQ(vth_n) / sigma;

            // Rates
            BoutReal R_rc = SQ(Nelim(i, j, k)) *
                            hydrogen.recombination(Nelim(i, j, k) * Nnorm,
                                                   Telim(i, j, k) * Tnorm) *
                            Nnorm / Omega_ci; // Time normalisation
            BoutReal R_iz = Nelim(i, j, k) * Nn(i, j, k) *
                            hydrogen.ionisation(Telim(i, j, k) * Tnorm) *
                            Nnorm / Omega_ci; // Time normalisation
            BoutReal R_cx = sigma_cx * Nn(i, j, k);

            // Plasma sink / neutral source
            S(i, j, k) = R_rc - R_iz;

            // Power transfer from plasma to neutrals

            Q(i, j, k) = R_cx * (3. / 2) * (Telim(i, j, k) - Tn(i, j, k));

            // Power transfer due to ionisation and recombination
            Q(i, j, k) +=
                (3. / 2) * (Telim(i, j, k) * R_rc - Tn(i, j, k) * R_iz);

            // Ion-neutral friction
            F(i, j, k) = R_cx    // Charge-Exchange
                         + R_rc; // Recombination

            // Radiated power from plasma
            // Factor of 1.09 so that recombination becomes an energy source at
            // 5.25eV
            Rp(i, j, k) = (1.09 * Telim(i, j, k) - 13.6 / Tnorm) * R_rc +
                          (Eionize / Tnorm) * R_iz; // Ionisation energy
          }

      // Neutral density
      ddt(Nn) = +S + Div_Perp_Lap_x3(Dnn, Nn, true);

      ddt(Ne) -= S;

      // Neutral pressure
      ddt(Pn) = (2. / 3) * Q + Div_Perp_Lap_x3(Dnn, Pn, true)
          //+ Div_Perp_Lap_x3(Tn*Dnn, Nn, true)  // Density diffusion
          //+ Div_Perp_Lap_x3(Nn*Dnn, Tn, true)  // Temperature diffusion
          ;
      ddt(Pe) -= (2. / 3) * (Q + Rp);

      if (neutral_friction) {
        // Vorticity
        if (boussinesq) {
          ddt(Vort) -= Div_Perp_Lap_FV(F / (Nelim * SQ(mesh->Bxy)), phi,
                                       vort_bndry_flux);
        } else {
          ddt(Vort) -= Div_Perp_Lap_FV(F / SQ(mesh->Bxy), phi, vort_bndry_flux);
        }
      }
      break;
    }
    case 2: {
      TRACE("neutral model 2");
      // Pseudo-recycling with Nn exponential approximation
      //
      // Calculate scaling of Nn exponential (such that n_lost = n_recy on each
      // field line)
      BoutReal nnexp, nlost, nnexp_tot;
      BoutReal vth_n, sigma_cx, sigma_iz, fluxout;
      BoutReal Nn0max = 10.0; // max(Nelim,true);
      lambda = 0.0;
      Nn0 = 0.0;
      // Approximate neutral density at t=0 to be exponential away from plate
      // with max density equal to ion max density
      if (time < 1e-10) {
        Field2D ll;
        ll = CumSumY2D(hthe * mesh->dy / Lmax, true);
        Nn = max(Nelim) * exp(-ll);
      }
      // calculate iz and cx mean free paths
      for (int i = mesh->xstart; i <= mesh->xend; i++) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          for (int k = 0; k < mesh->ngz; k++) {
            vth_n = sqrt(Tn(i, j, k)) * Cs0; // in m/s
            sigma_cx =
                Nelim(i, j, k) * Nnorm *
                hydrogen.chargeExchange(Telim(i, j, k) * Tnorm); // in seconds
            sigma_iz =
                Nelim(i, j, k) * Nnorm * Nn(i, j, k) *
                hydrogen.ionisation(Telim(i, j, k) * Tnorm); // in seconds
            // mean-free path for cx and iz is sqrt(l_cx*l_iz), and each are
            // vth/<sig-v>
            lambda(i, j, k) =
                vth_n / sqrt(sigma_cx * sigma_iz); // min( vth_n /
                                                   // sqrt(sigma_cx*sigma_iz),
                                                   // Lmax); // in meters
          }
        }
      }

      // Sum in y dy/lambda (in meters) for exponential (if lambda is constant,
      // results in exp(-y/lambda) )
      lambda_int = CumSumY3D(
          hthe * mesh->dy / lambda,
          true); // length in poloidal plane - appropriate for neutrals

      // int(S dV) = fr * N_lost (ie. volume integral of ionization density
      // source = recycling fraction * particle loss rate)
      for (int k = 0; k < mesh->ngz; k++) {
        for (int i = mesh->xstart; i <= mesh->xend; i++) {
          // plasma density flux out to divertor plate
          fluxout = 0.5 * (NVi(i, mesh->yend, k) +
                           NVi(i, mesh->yend + 1, k)); // [d^-2][t^-1]

          // number of particles lost per dt (flux out times perp volume) [t^-1]
          nlost = bcast_lasty(
              fluxout * 0.5 *
              (mesh->J(i, mesh->yend) * mesh->dx(i, mesh->yend) * mesh->dz /
                   sqrt(mesh->g_22(i, mesh->yend)) +
               mesh->J(i, mesh->yend + 1) * mesh->dx(i, mesh->yend + 1) *
                   mesh->dz / sqrt(mesh->g_22(i, mesh->yend + 1))));

          // Integrate ionization rate over volume to get volume loss rate
          // (simple integration using trap rule)
          nnexp = 0.0;
          for (int j = mesh->ystart; j <= mesh->yend; j++) {
            sigma_iz = hydrogen.ionisation(Telim(i, j, k) * Tnorm) * Nnorm /
                       Omega_ci; // ionization rate [d^3]/[t]
            BoutReal dV = mesh->J(i, j) * mesh->dx(i, j) * mesh->dy(i, j) *
                          mesh->dz; // volume element
            nnexp += Nelim(i, j, k) * sigma_iz * exp(-lambda_int(i, j, k)) *
                     dV; // full integral of density source [d^3]/[t]
          }
          MPI_Allreduce(&nnexp, &nnexp_tot, 1, MPI_DOUBLE, MPI_SUM,
                        mesh->getYcomm(i)); // add up all y-procs

          // neutral density factor
          for (int j = 0; j < mesh->ngy; j++) {
            Nn0(i, j, k) =
                frecycle * nlost / nnexp_tot; // ( [d^-3] ) Max neutral density
                                              // can't exceed max plasma density
          }
        }
      }
      Nn0 = -floor(-Nn0, -max(Ne, true));

      // Approximate neutral density as exponential
      Nn = Nn0 * exp(-lambda_int);

      // Calculate neutral density and source terms
      for (int i = mesh->xstart; i <= mesh->xend; i++)
        for (int j = mesh->ystart; j <= mesh->yend; j++)
          for (int k = 0; k < mesh->ngz; k++) {
            // Rates
            BoutReal R_rc = SQ(Nelim(i, j, k)) *
                            hydrogen.recombination(Nelim(i, j, k) * Nnorm,
                                                   Telim(i, j, k) * Tnorm) *
                            Nnorm / Omega_ci; // Time normalisation
            BoutReal R_iz = Nelim(i, j, k) * Nn(i, j, k) *
                            hydrogen.ionisation(Telim(i, j, k) * Tnorm) *
                            Nnorm / Omega_ci; // Time normalisation
            BoutReal R_cx = Nelim(i, j, k) * Nn(i, j, k) *
                            hydrogen.chargeExchange(Telim(i, j, k) * Tnorm) *
                            Nnorm / Omega_ci;

            // Ionisation plasma source, recombination sink
            S(i, j, k) = R_rc - R_iz;

            // Ionisation power transfer to plasma from neutrals
            Q(i, j, k) = -Tn(i, j, k) * R_iz;
            // Recombination plasma power sink
            Q(i, j, k) += Telim(i, j, k) * R_rc;
            // Power transfer from plasma to neutrals due to charge exchange
            Q(i, j, k) += R_cx * (Telim(i, j, k) - Tn(i, j, k));

            // Ion-neutral friction due to charge exchange
            F(i, j, k) = R_cx;
            // Friction due to recombination
            F(i, j, k) += R_rc;

            // Radiated power from plasma
            // Factor of 1.09 so that recombination becomes an energy source at
            // 5.25eV
            Rp(i, j, k) = (1.09 * Telim(i, j, k) - 13.6 / Tnorm) * R_rc +
                          (Eionize / Tnorm) * R_iz; // Ionisation energy
          }

      // Add sources to plasma fields
      ddt(Ne) -= S;
      ddt(Pe) -= Q + Rp;
      if (neutral_friction) {
        // Vorticity
        if (boussinesq) {
          ddt(Vort) -= Div_Perp_Lap_FV(F / (Nelim * SQ(mesh->Bxy)), phi,
                                       vort_bndry_flux);
        } else {
          ddt(Vort) -= Div_Perp_Lap_FV(F / SQ(mesh->Bxy), phi, vort_bndry_flux);
        }
      }
      break;
    } // End case 2
    case 3: {
      TRACE("Neutral model 3");
      //////////////////////////////////////////////////////
      // 2D (X-Y) full velocity model
      //
      // Evolves density Nn2D, velocity vector Vn2D and pressure Pn2D
      //
      mesh->communicate(Nn2D, Vn2D, Pn2D);

      if (outflow_ydown) {
        // Outflowing boundaries at ydown. If flow direction is
        // into domain then zero value is set. If flow is out of domain
        // then Neumann conditions are set

        for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone();
             idwn.next()) {

          if (Vn2D.y(idwn.ind, mesh->ystart) < 0.0) {
            // Flowing out of domain
            Vn2D.y(idwn.ind, mesh->ystart - 1) = Vn2D.y(idwn.ind, mesh->ystart);
          } else {
            // Flowing into domain
            Vn2D.y(idwn.ind, mesh->ystart - 1) =
                -Vn2D.y(idwn.ind, mesh->ystart);
          }
          // Neumann boundary condition on X and Z components
          Vn2D.x(idwn.ind, mesh->ystart - 1) = Vn2D.x(idwn.ind, mesh->ystart);
          Vn2D.z(idwn.ind, mesh->ystart - 1) = Vn2D.z(idwn.ind, mesh->ystart);

          // Neumann conditions on density and pressure
          Nn2D(idwn.ind, mesh->ystart - 1) = Nn2D(idwn.ind, mesh->ystart);
          Pn2D(idwn.ind, mesh->ystart - 1) = Pn2D(idwn.ind, mesh->ystart);
        }
      }

      // Density
      ddt(Nn2D) = -Div(Vn2D, Nn2D);

      Field2D Nn2D_floor = floor(Nn2D, 1e-2);
      // Velocity
      ddt(Vn2D) = -Grad(Pn2D) / Nn2D_floor;

      //////////////////////////////////////////////////////
      // Momentum advection

      // Convert to cylindrical coordinates for velocity
      // advection term. This is to avoid Christoffel symbol
      // terms in curvilinear geometry
      Field2D vr = Txr * Vn2D.x + Tyr * Vn2D.y; // Grad R component
      Field2D vz = Txz * Vn2D.x + Tyz * Vn2D.y; // Grad Z component

      // Advect as scalars (no Christoffel symbols needed)
      ddt(vr) = -V_dot_Grad(Vn2D, vr);
      ddt(vz) = -V_dot_Grad(Vn2D, vz);

      // Convert back to field-aligned coordinates
      ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
      ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

      //////////////////////////////////////////////////////
      // Viscosity
      // This includes dynamic ( neutral_viscosity)
      // and bulk/volume viscosity ( neutral_bulk )

      ddt(vr) = Laplace_FV(neutral_viscosity, vr);
      ddt(vz) = Laplace_FV(neutral_viscosity, vz);

      ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
      ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

      DivV2D = Div(Vn2D);
      DivV2D.applyBoundary(time);
      mesh->communicate(DivV2D);

      // ddt(Vn2D) += Grad( (neutral_viscosity/3. + neutral_bulk) * DivV2D ) /
      // Nn2D_floor;

      //////////////////////////////////////////////////////
      // Pressure
      ddt(Pn2D) =
          -Div(Vn2D, Pn2D) -
          (gamma_ratio - 1.) * Pn2D * DivV2D * floor(Nn2D, 0) / Nn2D_floor +
          Laplace_FV(neutral_conduction, Pn2D / Nn2D);

      ///////////////////////////////////////////////////////////////////
      // Boundary condition on fluxes

      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        // Calculate flux of ions into target from Ne and Vi boundary
        // This calculation is supposed to be consistent with the flow
        // of plasma from Div_par_FV(Ne, Ve)

        BoutReal flux_ion = 0.0;
        for (int jz = 0; jz < mesh->ngz - 1; jz++) {
          flux_ion -=
              frecycle * 0.5 *
              (Ne(r.ind, mesh->ystart, jz) + Ne(r.ind, mesh->ystart - 1, jz)) *
              0.5 * (Ve(r.ind, mesh->ystart, jz) +
                     Ve(r.ind, mesh->ystart - 1, jz)); // Flux through surface
                                                       // [m^-2 s^-1], should be
                                                       // positive since Ve <
                                                       // 0.0
        }
        flux_ion /= mesh->ngz - 1; // Average over Z

        // Flow of neutrals inwards
        BoutReal flow = flux_ion * (mesh->J(r.ind, mesh->ystart) +
                                    mesh->J(r.ind, mesh->ystart - 1)) /
                        (sqrt(mesh->g_22(r.ind, mesh->ystart)) +
                         sqrt(mesh->g_22(r.ind, mesh->ystart - 1)));

        // Rate of change of neutrals in final cell
        BoutReal dndt = flow / (mesh->J(r.ind, mesh->ystart) *
                                mesh->dy(r.ind, mesh->ystart));
        ddt(Nn2D)(r.ind, mesh->ystart) += dndt;
        ddt(Pn2D)(r.ind, mesh->ystart) +=
            dndt * (3.5 / Tnorm); // Franck-Condon energy

        // Loss of thermal energy to the target.
        // This depends on the reflection coefficient
        // and is controlled by the option neutral_gamma
        //         q = neutral_gamma * n * T * cs

        // Density at the target
        BoutReal Nnout =
            0.5 * (Nn2D(r.ind, mesh->ystart) + Nn2D(r.ind, mesh->ystart - 1));
        if (Nnout < 0.0)
          Nnout = 0.0;
        // Temperature at the target
        BoutReal Tnout =
            0.5 * (Tn2D(r.ind, mesh->ystart) + Tn2D(r.ind, mesh->ystart - 1));
        if (Tnout < 0.0)
          Tnout = 0.0;

        // gamma * n * T * cs
        BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
        // Multiply by cell area to get power
        BoutReal heatflux = q * (mesh->J(r.ind, mesh->ystart) +
                                 mesh->J(r.ind, mesh->ystart - 1)) /
                            (sqrt(mesh->g_22(r.ind, mesh->ystart)) +
                             sqrt(mesh->g_22(r.ind, mesh->ystart - 11)));

        // Divide by volume of cell, and multiply by 2/3 to get pressure
        ddt(Pn2D)(r.ind, mesh->ystart) -=
            (2. / 3) * heatflux /
            (mesh->dy(r.ind, mesh->ystart) * mesh->J(r.ind, mesh->ystart));
      }

      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        // Calculate flux of ions into target from Ne and Vi boundary
        // This calculation is supposed to be consistent with the flow
        // of plasma from Div_par_FV(Ne, Ve)

        BoutReal flux_ion = 0.0;
        for (int jz = 0; jz < mesh->ngz - 1; jz++) {
          flux_ion +=
              frecycle * 0.5 *
              (Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend + 1, jz)) *
              0.5 *
              (Ve(r.ind, mesh->yend, jz) +
               Ve(r.ind, mesh->yend + 1,
                  jz)); // Flux through surface [m^-2 s^-1], should be positive
        }
        flux_ion /= mesh->ngz - 1; // Average over Z

        // Flow of neutrals inwards
        BoutReal flow = flux_ion * (mesh->J(r.ind, mesh->yend) +
                                    mesh->J(r.ind, mesh->yend + 1)) /
                        (sqrt(mesh->g_22(r.ind, mesh->yend)) +
                         sqrt(mesh->g_22(r.ind, mesh->yend + 1)));

        // Rate of change of neutrals in final cell
        BoutReal dndt =
            flow / (mesh->J(r.ind, mesh->yend) * mesh->dy(r.ind, mesh->yend));
        ddt(Nn2D)(r.ind, mesh->yend) += dndt;
        ddt(Pn2D)(r.ind, mesh->yend) +=
            dndt * (3.5 / Tnorm); // Franck-Condon energy

        // Loss of thermal energy to the target.
        // This depends on the reflection coefficient
        // and is controlled by the option neutral_gamma
        //         q = neutral_gamma * n * T * cs

        // Density at the target
        BoutReal Nnout =
            0.5 * (Nn2D(r.ind, mesh->yend) + Nn2D(r.ind, mesh->yend + 1));
        if (Nnout < 0.0)
          Nnout = 0.0;
        // Temperature at the target
        BoutReal Tnout =
            0.5 * (Tn2D(r.ind, mesh->yend) + Tn2D(r.ind, mesh->yend + 1));
        if (Tnout < 0.0)
          Tnout = 0.0;

        // gamma * n * T * cs
        BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
        // Multiply by cell area to get power
        BoutReal heatflux =
            q * (mesh->J(r.ind, mesh->yend) + mesh->J(r.ind, mesh->yend + 1)) /
            (sqrt(mesh->g_22(r.ind, mesh->yend)) +
             sqrt(mesh->g_22(r.ind, mesh->yend + 1)));

        // Divide by volume of cell, and multiply by 2/3 to get pressure
        ddt(Pn2D)(r.ind, mesh->yend) -=
            (2. / 3) * heatflux /
            (mesh->dy(r.ind, mesh->yend) * mesh->J(r.ind, mesh->yend));
      }

      // Exchange of parallel momentum. This could be done
      // in a couple of ways, but here we use the fact that
      // Vn2D is covariant and b = e_y / (JB) to write:
      //
      // V_{||n} = b dot V_n = Vn2D.y / (JB)
      Field2D Vnpar = Vn2D.y / (mesh->J * mesh->Bxy);

      /////////////////////////////////////////////////////
      // Atomic processes

      Field3D Riz, Rrc, Rcx;
      neutral_rates(Ne, Telim, Vi, Nn, Tn, Vnpar, S, F, Q, Rp, Riz, Rrc, Rcx);

      Field3D Fvort = Rrc + Rcx; // Friction for vorticity

      // Loss of momentum in the X and Z directions
      ddt(Vn2D).x -= (Rcx.DC() + Riz.DC()) * Vn2D.x / Nn2D_floor;
      ddt(Vn2D).z -= (Rcx.DC() + Riz.DC()) * Vn2D.z / Nn2D_floor;

      // Particles
      ddt(Nn2D) += S.DC(); // Average over toroidal angle z
      ddt(Ne) -= S;

      // Energy
      ddt(Pn2D) += (2. / 3) * Q.DC();
      ddt(Pe) -= (2. / 3) * (Q + Rp);

      // Momentum. Note need to turn back into covariant form
      ddt(Vn2D).y += F.DC() * (mesh->J * mesh->Bxy) / Nn2D_floor;
      ddt(NVi) -= F;

      if (neutral_friction) {
        // Vorticity
        if (boussinesq) {
          ddt(Vort) -= Div_Perp_Lap_FV(Fvort / (Nelim * SQ(mesh->Bxy)), phi,
                                       vort_bndry_flux);
        } else {
          ddt(Vort) -=
              Div_Perp_Lap_FV(Fvort / SQ(mesh->Bxy), phi, vort_bndry_flux);
        }
      }

      // Density evolution
      for (int i = 0; i < mesh->ngx; i++)
        for (int j = 0; j < mesh->ngy; j++) {
          if ((Nn2D(i, j) < 1e-8) && (ddt(Nn2D)(i, j) < 0.0)) {
            ddt(Nn2D)(i, j) = 0.0;
          }
        }

      break;
    } // End case 3
    case 4: {
      //////////////////////////////////////////////////////
      // 3D model, diffusive in X-Z, fluid in Y
      //
      // Evolves neutral density Nn, pressure Pn, and
      // parallel momentum NVn
      //
      TRACE("Neutral model 4");

      Nn = floor(Nn, 1e-8);
      Field3D Nnlim =
          floor(Nn, 1e-5); // Used where division by neutral density is needed
      Tn = Pn / Nn;
      Tn = floor(Tn, 0.01 / Tnorm);

      Field3D Vn = NVn / Nnlim; // Neutral parallel velocity

      Field3D Pnlim = Nn * Tn;
      Pnlim.applyBoundary("neumann");

      /////////////////////////////////////////////////////
      // Boundary conditions
      TRACE("Neutral boundary conditions");

      if (sheath_ydown) {
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->ngz; jz++) {
            // Free boundary (constant gradient) density
            BoutReal nnwall = 0.5 * (3. * Nn(r.ind, mesh->ystart, jz) -
                                     Nn(r.ind, mesh->ystart + 1, jz));
            if (nnwall < 0.0)
              nnwall = 0.0;

            BoutReal tnwall = Tn(r.ind, mesh->ystart, jz);

            Nn(r.ind, mesh->ystart - 1, jz) =
                2 * nnwall - Nn(r.ind, mesh->ystart, jz);

            // Zero gradient temperature, heat flux added later
            Tn(r.ind, mesh->ystart - 1, jz) = tnwall;

            // Set pressure consistent at the boundary
            Pn(r.ind, mesh->ystart - 1, jz) =
                2. * nnwall * tnwall - Pn(r.ind, mesh->ystart, jz);
            // No flow into wall
            Vn(r.ind, mesh->ystart - 1, jz) = -Vn(r.ind, mesh->ystart, jz);
            NVn(r.ind, mesh->ystart - 1, jz) = -NVn(r.ind, mesh->ystart, jz);
          }
        }
      }

      if (sheath_yup) {
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->ngz; jz++) {
            // Free boundary (constant gradient) density
            BoutReal nnwall = 0.5 * (3. * Nn(r.ind, mesh->yend, jz) -
                                     Nn(r.ind, mesh->yend - 1, jz));
            if (nnwall < 0.0)
              nnwall = 0.0;

            BoutReal tnwall = Tn(r.ind, mesh->yend, jz);

            Nn(r.ind, mesh->yend + 1, jz) =
                2 * nnwall - Nn(r.ind, mesh->yend, jz);

            // Zero gradient temperature, heat flux added later
            Tn(r.ind, mesh->yend + 1, jz) = tnwall;

            // Set pressure consistent at the boundary
            Pn(r.ind, mesh->yend + 1, jz) =
                2. * nnwall * tnwall - Pn(r.ind, mesh->yend, jz);
            // No flow into wall
            Vn(r.ind, mesh->yend + 1, jz) = -Vn(r.ind, mesh->yend, jz);
            NVn(r.ind, mesh->yend + 1, jz) = -NVn(r.ind, mesh->yend, jz);
          }
        }
      }

      /////////////////////////////////////////////////////
      // Atomic processes
      TRACE("Atomic processes");

      Field3D Riz, Rrc, Rcx;
      neutral_rates(Ne, Telim, Vi, Nn, Tn, Vn, S, F, Q, Rp, Riz, Rrc, Rcx);

      // Neutral cross-field diffusion coefficient
      BoutReal neutral_lmax = 0.1 / rho_s0;
      Field3D Rnn = Nn * sqrt(Tn) / neutral_lmax; // Neutral-neutral collisions
      Dnn = Pnlim / (Riz + Rcx + Rnn);
      mesh->communicate(Dnn);
      Dnn.applyBoundary("dirichlet_o2");

      // Logarithms used to calculate perpendicular velocity
      // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
      //
      // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
      //               = Grad(logTn + logNn)
      // Field3D logNn = log(Nn);
      // Field3D logTn = log(Tn);

      Field3D logPnlim = log(Pnlim);
      logPnlim.applyBoundary("neumann");

      /////////////////////////////////////////////////////
      // Neutral density
      TRACE("Neutral density");

      ddt(Nn) = -Div_par_FV(Nn, Vn) // Advection
                + S                 // Source from recombining plasma
                + Div_Perp_Lap_XYZ(Dnn * Nn, logPnlim,
                                   false) // Perpendicular diffusion
          ;

      ddt(Ne) -= S; // Sink of plasma density

      /////////////////////////////////////////////////////
      // Neutral momentum
      TRACE("Neutral momentum");

      ddt(NVn) = -Div_par_FV(NVn, Vn) // Momentum flow
                 + F                  // Friction with plasma
                 - Grad_par(Pnlim)    // Pressure gradient
                 + Div_Perp_Lap_XYZ(Dnn * NVn, logPnlim,
                                    false) // Perpendicular diffusion

                 + Div_par_diffusion(Dnn * Nn, Vn, false) // Parallel viscosity
          ;

      if (neut_numdiff > 0.0) {
        ddt(NVn) += neut_numdiff * Div_par_diffusion_index(Vn);
        // ddt(NVn) += Div_par_diffusion(SQ(mesh->dy)*mesh->g_22*neut_numdiff,
        // Vn);
      }

      ddt(NVi) -= F;

      if (neutral_friction) {
        // Vorticity
        if (boussinesq) {
          ddt(Vort) -= Div_Perp_Lap_FV((Rcx + Rrc) / (Nelim * SQ(mesh->Bxy)),
                                       phi, vort_bndry_flux);
        } else {
          ddt(Vort) -= Div_Perp_Lap_FV((Rcx + Rrc) / SQ(mesh->Bxy), phi,
                                       vort_bndry_flux);
        }
      }

      /////////////////////////////////////////////////////
      // Neutral pressure
      TRACE("Neutral pressure");

      ddt(Pn) =
          -Div_par_FV(Pn, Vn)              // Advection
          - (2. / 3) * Pnlim * Div_par(Vn) // Compression
          + (2. / 3) * Q +
          Div_Perp_Lap_XYZ(Dnn * Pn, logPnlim, false) // Perpendicular diffusion
          + Div_Perp_Lap_XYZ(Dnn * Nn, Tn, false)     // Conduction
          + Div_par_diffusion(Dnn * Nn, Tn)           // Parallel conduction
          ;

      ddt(Pe) -= (2. / 3) * (Q + Rp);

      ///////////////////////////////////////////////////////////////////
      // Boundary condition on fluxes
      TRACE("Neutral boundary fluxes");
      wall_flux = 0.0;

      if (sheath_ydown) {
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          // Calculate flux of ions into target from Ne and Vi boundary
          // This calculation is supposed to be consistent with the flow
          // of plasma from Div_par_FV(Ne, Ve)

          for (int jz = 0; jz < mesh->ngz - 1; jz++) {
            BoutReal flux_ion =
                -0.5 * (Ne(r.ind, mesh->ystart, jz) +
                        Ne(r.ind, mesh->ystart - 1, jz)) *
                0.5 * (Ve(r.ind, mesh->ystart, jz) +
                       Ve(r.ind, mesh->ystart - 1, jz)); // Flux through surface
                                                         // [m^-2 s^-1], should
                                                         // be positive since Ve
                                                         // < 0.0

            // Flow of neutrals inwards
            BoutReal flow = frecycle * flux_ion *
                            (mesh->J(r.ind, mesh->ystart) +
                             mesh->J(r.ind, mesh->ystart - 1)) /
                            (sqrt(mesh->g_22(r.ind, mesh->ystart)) +
                             sqrt(mesh->g_22(r.ind, mesh->ystart - 1)));

            // Rate of change of neutrals in final cell
            BoutReal dndt = flow / (mesh->J(r.ind, mesh->ystart) *
                                    mesh->dy(r.ind, mesh->ystart));

            ddt(Nn)(r.ind, mesh->ystart, jz) += dndt;
            ddt(Pn)(r.ind, mesh->ystart, jz) +=
                dndt * (3.5 / Tnorm); // Franck-Condon energy
            ddt(NVn)(r.ind, mesh->ystart, jz) +=
                dndt * neutral_vwall * sqrt(3.5 / Tnorm);

            // Loss of thermal energy to the target.
            // This depends on the reflection coefficient
            // and is controlled by the option neutral_gamma
            //         q = neutral_gamma * n * T * cs

            // Density at the target
            BoutReal Nnout = 0.5 * (Nn(r.ind, mesh->ystart, jz) +
                                    Nn(r.ind, mesh->ystart - 1, jz));
            if (Nnout < 0.0)
              Nnout = 0.0;
            // Temperature at the target
            BoutReal Tnout = 0.5 * (Tn(r.ind, mesh->ystart, jz) +
                                    Tn(r.ind, mesh->ystart - 1, jz));
            if (Tnout < 0.0)
              Tnout = 0.0;

            // gamma * n * T * cs
            BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
            // Multiply by cell area to get power
            BoutReal heatflux = q * (mesh->J(r.ind, mesh->ystart) +
                                     mesh->J(r.ind, mesh->ystart - 1)) /
                                (sqrt(mesh->g_22(r.ind, mesh->ystart)) +
                                 sqrt(mesh->g_22(r.ind, mesh->ystart - 11)));

            // Divide by volume of cell, and multiply by 2/3 to get pressure
            ddt(Pn)(r.ind, mesh->ystart, jz) -=
                (2. / 3) * heatflux /
                (mesh->dy(r.ind, mesh->ystart) * mesh->J(r.ind, mesh->ystart));
          }
        }
      }

      if (sheath_yup) {
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          // Calculate flux of ions into target from Ne and Vi boundary
          // This calculation is supposed to be consistent with the flow
          // of plasma from Div_par_FV(Ne, Ve)

          for (int jz = 0; jz < mesh->ngz - 1; jz++) {
            BoutReal flux_ion =
                frecycle * 0.5 *
                (Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend + 1, jz)) *
                0.5 * (Ve(r.ind, mesh->yend, jz) +
                       Ve(r.ind, mesh->yend + 1, jz)); // Flux through surface
                                                       // [m^-2 s^-1], should be
                                                       // positive

            // Flow of neutrals inwards
            BoutReal flow = flux_ion * (mesh->J(r.ind, mesh->yend) +
                                        mesh->J(r.ind, mesh->yend + 1)) /
                            (sqrt(mesh->g_22(r.ind, mesh->yend)) +
                             sqrt(mesh->g_22(r.ind, mesh->yend + 1)));

            // Rate of change of neutrals in final cell
            BoutReal dndt = flow / (mesh->J(r.ind, mesh->yend) *
                                    mesh->dy(r.ind, mesh->yend));
            ddt(Nn)(r.ind, mesh->yend, jz) += dndt;
            ddt(Pn)(r.ind, mesh->yend, jz) +=
                dndt * (3.5 / Tnorm); // Franck-Condon energy
            ddt(NVn)(r.ind, mesh->yend, jz) -=
                dndt * neutral_vwall * sqrt(3.5 / Tnorm);

            // Loss of thermal energy to the target.
            // This depends on the reflection coefficient
            // and is controlled by the option neutral_gamma
            //         q = neutral_gamma * n * T * cs

            // Density at the target
            BoutReal Nnout = 0.5 * (Nn(r.ind, mesh->yend, jz) +
                                    Nn(r.ind, mesh->yend + 1, jz));
            if (Nnout < 0.0)
              Nnout = 0.0;
            // Temperature at the target
            BoutReal Tnout = 0.5 * (Tn(r.ind, mesh->yend, jz) +
                                    Tn(r.ind, mesh->yend + 1, jz));
            if (Tnout < 0.0)
              Tnout = 0.0;

            // gamma * n * T * cs
            BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
            // Multiply by cell area to get power
            BoutReal heatflux = q * (mesh->J(r.ind, mesh->yend) +
                                     mesh->J(r.ind, mesh->yend + 1)) /
                                (sqrt(mesh->g_22(r.ind, mesh->yend)) +
                                 sqrt(mesh->g_22(r.ind, mesh->yend + 1)));

            // Divide by volume of cell, and multiply by 2/3 to get pressure
            ddt(Pn)(r.ind, mesh->yend, jz) -=
                (2. / 3) * heatflux /
                (mesh->dy(r.ind, mesh->yend) * mesh->J(r.ind, mesh->yend));
          }
        }
      }
      break;
    } // End case 4
    } // End neutral_model switch
  }

  //////////////////////////////////////////////////////////////
  // Impurities

  if (carbon_fraction > 0.0) {
    TRACE("Carbon impurity radiation");
    Rzrad = carbon_rad->power(Te * Tnorm, Ne * Nnorm,
                              Ne * (Nnorm * carbon_fraction)); // J / m^3 / s
    Rzrad /= SI::qe * Tnorm * Nnorm * Omega_ci;                // Normalise

    ddt(Pe) -= (2. / 3) * Rzrad;
  }

  //////////////////////////////////////////////////////////////
  // Parallel closures for 2D simulations

  if (sinks) {
    // Sink terms for 2D simulations

    // Field3D nsink = 0.5*Ne*sqrt(Telim)*sink_invlpar;   // n C_s/ (2L)  //
    // Sound speed flow to targets
    Field3D nsink = 0.5 * sqrt(Ti) * Ne * sink_invlpar;
    nsink = floor(nsink, 0.0);

    ddt(Ne) -= nsink;

    Field3D conduct = (2. / 3) * kappa_epar * Te * SQ(sink_invlpar);
    conduct = floor(conduct, 0.0);
    ddt(Pe) -= conduct      // Heat conduction
               + Te * nsink // Advection
        ;

    if (sheath_closure) {
      ///////////////////////////
      // Sheath dissipation closure

      Field3D phi_te = floor(phi / Telim, 0.0);
      Field3D jsheath = Nelim * sqrt(Telim) *
                        (1 - (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te));

      ddt(Vort) += jsheath * sink_invlpar;
    } else {
      ///////////////////////////
      // Vorticity closure
      ddt(Vort) -= Div_Perp_Lap_FV(nsink / SQ(mesh->Bxy), phi, vort_bndry_flux);
      // ddt(Vort) -= (nsink / Ne)*Vort;
      // Note: If nsink = n * nu then this reduces to
      // ddt(Vort) -= nu * Vort
    }

    if (drift_wave) {
      // Include a drift-wave closure in the core region
      // as in Hasegawa-Wakatani and SOLT models

      Field2D Tedc = Telim.DC(); // Zonal average
      Field3D Q = phi - Telim * log(Nelim);

      // Drift-wave operator
      Field3D Adw = alpha_dw * (Tedc ^ 1.5) * (Q - Q.DC());

      ddt(Ne) += Adw;
      ddt(Vort) += Adw;
    }

    // Electron and ion parallel dynamics not evolved
  }

  if (low_pass_z >= 0) {
    // Low pass Z filtering, keeping up to and including low_pass_z
    ddt(Ne) = lowPass(ddt(Ne), low_pass_z);
    ddt(Pe) = lowPass(ddt(Pe), low_pass_z);

    if (currents) {
      ddt(Vort) = lowPass(ddt(Vort), low_pass_z);
      if (electromagnetic || FiniteElMass) {
        ddt(VePsi) = lowPass(ddt(VePsi), low_pass_z);
      }
    }
    if (ion_velocity) {
      ddt(NVi) = lowPass(ddt(NVi), low_pass_z);
    }
  }

  if (!evolve_plasma) {
    ddt(Ne) = 0.0;
    ddt(Pe) = 0.0;
    ddt(Vort) = 0.0;
    ddt(VePsi) = 0.0;
    ddt(NVi) = 0.0;
  }

  return 0;
}

/*!
 * Preconditioner. Solves the heat conduction
 *
 * @param[in] t  The simulation time
 * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J). Related
 * to timestep
 * @param[in] delta   Not used here
 */
int Hermes::precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  static InvertPar *inv = NULL;
  if (!inv) {
    // Initialise parallel inversion class
    inv = InvertPar::Create();
    inv->setCoefA(1.0);
  }
  if (thermal_conduction) {
    // Set the coefficient in front of Grad2_par2
    inv->setCoefB(-(2. / 3) * gamma * kappa_epar);
    Field3D dT = ddt(Pe);
    dT.applyBoundary("neumann");
    ddt(Pe) = inv->solve(dT);
  }

  if (neutral_model == 1) {
    // Neutral gas diffusion
    // Solve (1 - gamma*Dnn*Delp2)^{-1}
    static Laplacian *inv = NULL;
    if (!inv) {
      inv = Laplacian::create();
      // Zero value outer boundary

      inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);

      inv->setCoefA(1.0);
    }

    // output << min(Dnn) << ", " << max(Dnn) << endl;

    inv->setCoefD(-gamma * Dnn);

    ddt(Nn) = inv->solve(ddt(Nn));

    ddt(Pn) = inv->solve(ddt(Pn));
  }
  return 0;
}

const Field3D Hermes::D(const Field3D &f, BoutReal d) { // Diffusion operator
  if (d < 0.0)
    return 0.0;

  return Div_par_diffusion(d * SQ(mesh->dy * mesh->g_22), f)
      //+ d*SQ(mesh->dx)*D2DX2(f)
      //+ Div_Perp_Lap_FV(d*SQ(mesh->dx*mesh->g_11), f, false)
      //+ Div_Perp_Lap_x3(d*SQ(mesh->dx), f, false)
      ;
}

const Field3D Hermes::H(const Field3D &f, BoutReal h) { // Diffusion operator
  if (h < 0.0)
    return 0.0;

  return -D4DY4_FV(h * SQ(SQ(mesh->dy)), f)
      //- h*SQ(SQ(mesh->dx))*D4DX4(f)
      //- h*SQ(SQ(mesh->dz))*D4DZ4(f)
      ;
}

const Field3D Hermes::H(const Field3D &f, BoutReal h,
                        const Field3D &mask) { // Diffusion operator
  if (h < 0.0)
    return 0.0;

  return -D4DY4_FV(h * mask * SQ(SQ(mesh->dy)), f) -
         h * mask * SQ(SQ(mesh->dx)) * D4DX4(f) -
         h * mask * SQ(SQ(mesh->dz)) * D4DZ4(f);
}

const Field3D Hermes::Grad_parP(const Field3D &f) {
  return Grad_par(f); //+ 0.5*beta_e*bracket(psi, f, BRACKET_ARAKAWA);
}

const Field3D Hermes::Grad_parP_CtoL(const Field3D &f) {
  if (staggered) {
    return Grad_par_CtoL(f); //+ 0.5*beta_e*bracket(psi, f, BRACKET_ARAKAWA);
  }
  return Grad_par(f); //+ 0.5*beta_e*bracket(psi, f, BRACKET_ARAKAWA);
}

const Field3D Hermes::Grad_parP_LtoC(const Field3D &f) {
  if (staggered) {
    return Grad_par_LtoC(f); //+ 0.5*beta_e*bracket(psi, f, BRACKET_ARAKAWA);
  }
  return Grad_par(f); // + 0.5*beta_e*bracket(psi, f, BRACKET_ARAKAWA);
}

const Field3D Hermes::Div_parP_LtoC(const Field3D &f, const Field3D &v) {
  return Div_par_FV(
      f,
      v); //+ 0.5*beta_e*mesh->Bxy*bracket(psi, f/mesh->Bxy, BRACKET_ARAKAWA);
}

const Field2D Hermes::CumSumY2D(const Field2D &f, bool reverse) {
  // Cumulative sum in Y one xz-slice at a time
  //    -- reverse is option to sum from the end of Y
  Field2D result = 0.0;

  if (reverse) {
    for (int i = mesh->xstart; i <= mesh->xend; i++) {
      // All but last processor receive
      if (!mesh->lastY()) {
        mesh->wait(mesh->irecvYOutOutdest(&result(i, mesh->yend + 1), 1,
                                          mesh->ngx * i));
      }
      // Calculate sum (reversed)
      for (int j = mesh->yend; j >= mesh->ystart; j--) {
        result[i][j] = result[i][j + 1] + f[i][j];
      }
      // Send the value at yend to the next processor.
      mesh->sendYInOutdest(&result(i, mesh->ystart), 1, mesh->ngx * i);
    }
  } else {
    for (int i = mesh->xstart; i <= mesh->xend; i++) {
      // All but first processor receive
      if (!mesh->firstY()) {
        mesh->wait(mesh->irecvYInIndest(&result(i, mesh->ystart - 1), 1,
                                        mesh->ngx * i));
      }
      // Calculate sum
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        result[i][j] = result[i][j - 1] + f[i][j];
      }
      // Send the value at yend to the next processor.
      mesh->sendYOutIndest(&result(i, mesh->yend), 1, mesh->ngx * i);
    }
  }

  return result;
}

const Field3D Hermes::CumSumY3D(const Field3D &f, bool reverse) {
  // Cumulative sum in Y one xz-slice at a time
  //    -- reverse is option to sum from the end of Y
  Field3D result = 0.0;

  if (reverse) {
    for (int k = 0; k < mesh->ngz; k++) {
      for (int i = mesh->xstart; i <= mesh->xend; i++) {
        // All but last processor receive
        if (!mesh->lastY()) {
          mesh->wait(mesh->irecvYOutOutdest(&result(i, mesh->yend + 1, k), 1,
                                            mesh->ngx * i + mesh->ngz * k));
        }
        // Calculate sum (reversed)
        for (int j = mesh->yend; j >= mesh->ystart; j--) {
          result[i][j][k] = result[i][j + 1][k] + f[i][j][k];
        }
        // Send the value at yend to the next processor.
        mesh->sendYInOutdest(&result(i, mesh->ystart, k), 1,
                             mesh->ngx * i + mesh->ngz * k);
      }
    }
  } else {
    for (int k = 0; k < mesh->ngz; k++) {
      for (int i = mesh->xstart; i <= mesh->xend; i++) {
        // All but first processor receive
        if (!mesh->firstY()) {
          mesh->wait(mesh->irecvYInIndest(&result(i, mesh->ystart - 1, k), 1,
                                          mesh->ngx * i + mesh->ngz * k));
        }
        // Calculate sum
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          result[i][j][k] = result[i][j - 1][k] + f[i][j][k];
        }
        // Send the value at yend to the next processor.
        mesh->sendYOutIndest(&result(i, mesh->yend, k), 1,
                             mesh->ngx * i + mesh->ngz * k);
      }
    }
  }

  return result;
}

const BoutReal Hermes::bcast_lasty(const BoutReal f) {
  BoutReal myf;
  // All but last processor receive
  if (!mesh->lastY()) {
    mesh->wait(mesh->irecvYOutOutdest(&myf, 1, 0));
  } else {
    myf = f;
  }
  // Send the value at yend to the next processor.
  mesh->sendYInOutdest(&myf, 1, 0);

  return myf;
}

// Standard main() function
BOUTMAIN(Hermes);

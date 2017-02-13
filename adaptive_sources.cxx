#include "adaptive_sources.hxx"

#include <field_factory.hxx>

void AdaptiveSource::init(const string &prefix, const Field2D &start, Options *opt, Solver *solver, Mesh *mesh, bool restarting) {

  // The error integral needs to be preserved across restarts
  solver->addToRestart(error_integral, prefix+"_error_integral");

  // Save the starting source profile
  profile = start;

  // Options for the feedback controller
  OPTION(opt, controller_i, 1e-4);
  OPTION(opt, controller_p, 1e-2);
  OPTION(opt, source_core_only, false);
  
  if(controller_i < 1e-10) {
    throw BoutException("controller_i too small. Must be > 1e-10");
  }
  if(controller_p < 0.0) {
    throw BoutException("controller_p must be positive");
  }
  
  if(!restarting) {
    // Set this so that te starting source = profile
    error_integral = 0.0;
  }
  
  // Cache the cell volume, as this is used every time
  volume = mesh->J * mesh->dx * mesh->dy;
}

const Field2D AdaptiveSource::get(const Field2D &var, const Field2D &target, BoutReal time) {

  // Average error over the domain
  // Weighting by cell volume
  BoutReal error = 0.0, weight = 0.0;

  // Loop over flux surfaces
  for(int x=mesh->xstart;x<=mesh->xend;x++) {
    if(source_core_only && !mesh->periodicY(x)) // Not periodic, and only sources in core, so skip
      continue;
    
    for(int y=mesh->ystart; y<=mesh->yend; y++) {
      BoutReal vol = volume(x,y);
      error += (var(x,y) - target(x,y)) * vol;
      weight += vol;
    }
  }

  BoutReal sending[2] = {error, weight};
  BoutReal receiving[2];
  
  // Sum across all processors
  MPI_Allreduce(sending, receiving, 2, MPI_DOUBLE, MPI_SUM, BoutComm::get());
  // all processors now have the total error and weight
  // Calculate the average error over the domain
  error = receiving[0] / receiving[1]; 

  if(last_time < 0.0) {
    // First time
    last_time = time;
    last_error = error;
  }

  // Integrate using Trapezium rule
  if(time > last_time) { // Since time can decrease
    error_integral += (time - last_time)*
      0.5*(error + last_error);
  }
  last_error = error;
  last_time = time;
  
  BoutReal source_multiplier = 1.0 + controller_p * error + controller_i * error_integral;
  
  return source_multiplier * profile;
}


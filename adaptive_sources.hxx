
#include <options.hxx>
#include <bout/solver.hxx>
#include <field2d.hxx>

/// Provides an adaptive source
///
///
class AdaptiveSource {
public:
  AdaptiveSource() : last_time(-1.0), error_integral(0.0) {}
  
  /// Initialise the source
  /// adding variables to solver as needed
  ///
  /// @param[in] prefix
  void init(const string &prefix, const Field2D &start, Options *opt, Solver *solver, Mesh *mesh, bool restarting);
  
  /// Update the source if needed, and return at given time
  const Field2D get(const Field2D &var, const Field2D &target, BoutReal time);
private:

  Field2D profile; // Source profile
  
  Field2D volume; // Cell volume
  
  BoutReal last_time; // The previous time the source was updated
  BoutReal last_error; // The error at the last evaluation
  BoutReal error_integral; // Integral of error over time

  BoutReal controller_p; // Proportional to current error
  BoutReal controller_i; // Proportional to integral of error

  bool source_core_only;  // Source only in the core?
};

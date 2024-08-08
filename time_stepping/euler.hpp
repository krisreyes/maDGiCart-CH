#pragma once

#include "time_integrator.hpp"

/**
 * Explicit Euler time integrator.
 */
class Euler : public TimeIntegrator {
 public:
  
  Euler(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);

  void doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
      override;

 private:

  std::unique_ptr<SolutionState> stage_solution_;

};
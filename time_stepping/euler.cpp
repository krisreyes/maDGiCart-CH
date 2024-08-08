#include "euler.hpp"

#include "data_structures/exec_includes.hpp"


Euler::Euler(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) : TimeIntegrator(rhs, opts)
{
  stage_solution_ = rhs.createSolutionState();
}


void Euler::doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
{
  profile();

  rhs.evalRHSImpl(state, time, *stage_solution_);

  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto stage_sol = write_access(stage_solution_->getVec(ivec));
    auto residual = write_access(dstate_dt.getVec(ivec));

    maDGForAll(i, 0, residual.size(), {  //
      residual[i] = stage_sol[i];
    });
  }

  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    // auto residual = write_access(dstate_dt.getVec(ivec));
    auto sol = read_write_access(state.getVec(ivec));
    // auto res = read_access(stage_solution_->getVec(ivec));
    auto res = read_access(dstate_dt.getVec(ivec));

    maDGForAll(i, 0, sol.size(), {
      sol[i] += dt * res[i];
      // residual[i] = res[i];
    });
  }
}


static auto eulerinstance = FactoryRegistry<TimeIntegrator>::get().add(
    "euler",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) { return std::make_unique<Euler>(rhs, opts); });
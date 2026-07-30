#pragma once

#include "time_integrator.hpp"

#include <cstddef>
#include <complex>
#include <string>
#include <vector>

class FFTSemiImplicit : public TimeIntegrator {
 public:
  FFTSemiImplicit(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);
  ~FFTSemiImplicit() override;

  void solve(TimeIntegrableRHS& rhs, InitialConditions& initial_conditions) override;

  void doTimeStep(TimeIntegrableRHS&, SolutionState&, SolutionState&, double, double) override
  {
    Logger::get().FatalMessage("FFTSemiImplicit advances through solve(), not doTimeStep().");
  }

  int_t  getCurrentStep() const override { return iter_; }
  real_t getCurrentTime() const override { return time_; }
  real_t getTimeStepSize() const override { return dt_; }

 private:
  struct ModeSymbol {
    double k2 = 0.0;
    double k4 = 0.0;
    bool   keep_nonlinear = true;
  };

  int    dim_ = 0;
  int    nx_  = 0;
  int    ny_  = 0;
  int    nz_  = 0;
  size_t n_total_ = 0;

  double dx_ = 0.0;
  double dy_ = 0.0;
  double dz_ = 0.0;

  double mean_          = 0.0;
  double eps2_          = 0.0;
  double sigma_         = 0.0;
  double stabilization_ = 0.0;
  double max_abs_c_     = 0.0;
  double growth_factor_ = 1.0;
  double shrink_factor_ = 0.5;
  int    max_retries_   = 20;
  bool   dealias_       = true;
  bool   use_measure_   = false;
  std::string laplacian_type_;

  int_t  iter_ = 0;
  real_t time_ = 0.0;
  real_t dt_   = 0.0;

  std::vector<double> old_field_;
  std::vector<double> field_;
  std::vector<double> trial_field_;
  std::vector<double> nonlinear_;
  std::vector<ModeSymbol> symbols_;
  std::vector<std::complex<double>> spatial_;
  std::vector<std::complex<double>> spectral_;
  std::vector<std::complex<double>> nonlinear_hat_;
  std::vector<std::complex<double>> line_buffer_;

  SolutionState& mutableSolutionState()
  {
    return const_cast<SolutionState&>(TimeIntegrator::getCurrentSolutionState());
  }

  SolutionState& mutableResidualState()
  {
    return const_cast<SolutionState&>(TimeIntegrator::getCurrentResidual());
  }

  void initialize(TimeIntegrableRHS& rhs);
  void allocateFFT();
  void fft3D(std::vector<std::complex<double>>& data, bool inverse);
  void transformAxis(std::vector<std::complex<double>>& data, int axis, bool inverse);
  void buildSymbols();
  void copyStateToField(const TimeIntegrableRHS& rhs, const SolutionState& state, std::vector<double>& field) const;
  void copyFieldToState(const TimeIntegrableRHS& rhs, const std::vector<double>& field, SolutionState& state) const;

  bool takeStep(const std::vector<double>& input, double dt, std::vector<double>& output);
  void projectMean(std::vector<double>& field) const;
  bool fieldIsUsable(const std::vector<double>& field) const;
  double relativeUpdateNorm(const std::vector<double>& before, const std::vector<double>& after) const;
  double maxAbs(const std::vector<double>& field) const;
};

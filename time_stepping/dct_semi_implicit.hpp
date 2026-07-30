#pragma once

// -----------------------------------------------------------------------------------------------
// DctSemiImplicit
//
// Pseudospectral semi-implicit time integrator for the (local, source-augmented) Cahn-Hilliard
// equation:
//
//     du/dt = -eps2 * Laplacian^2(u) + Laplacian(u^3 - u) - sigma * (u - m)
//
// Boundary conditions:
//     x, y : periodic          -> complex FFT
//     z    : Neumann (du/dz=0) -> discrete cosine transform (DCT-II / DCT-III), matching a
//                                  cell-centered, ghost-symmetric FD grid.
//
// This is the "Path B" solver: same physics/time-stepping scheme as FFTSemiImplicit
// (fft_semi_implicit.cpp/.hpp), generalized to mixed periodic/Neumann boundaries via FFTW instead
// of the fully-periodic in-house radix FFT. FFTSemiImplicit is left untouched; this is a
// separate, independently-registered integrator.
//
// Registered names: "dct_semi_implicit" (long) and "dct_si" (short alias).
//
// Requires linking against FFTW3 (see time_stepping/CMakeLists.txt).
// -----------------------------------------------------------------------------------------------

#include "time_integrator.hpp"

#include <complex>
#include <cstddef>
#include <fftw3.h>
#include <string>
#include <vector>

class DctSemiImplicit : public TimeIntegrator {
 public:
  DctSemiImplicit(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);
  ~DctSemiImplicit() override;

  void solve(TimeIntegrableRHS& rhs, InitialConditions& initial_conditions) override;

  void doTimeStep(TimeIntegrableRHS&, SolutionState&, SolutionState&, double, double) override
  {
    Logger::get().FatalMessage("DctSemiImplicit advances through solve(), not doTimeStep().");
  }

  int_t  getCurrentStep() const override { return iter_; }
  real_t getCurrentTime() const override { return time_; }
  real_t getTimeStepSize() const override { return dt_; }

 private:
  struct ModeSymbol {
    double k2             = 0.0;
    double k4             = 0.0;
    bool   keep_nonlinear = true;
  };

  int    nx_ = 0, ny_ = 0, nz_ = 0;   // nz_ is the Neumann (cosine) axis
  size_t n_total_ = 0;

  double dx_ = 0.0, dy_ = 0.0, dz_ = 0.0;

  double mean_  = 0.0;   // target mass ratio m
  double eps2_  = 0.0;   // interfacial-width parameter
  double sigma_ = 0.0;   // mean-regressive coefficient

  double stabilization_ = 0.0;
  double max_abs_c_     = 0.0;
  double growth_factor_ = 1.0;
  double shrink_factor_ = 0.5;
  int    max_retries_   = 20;
  bool   dealias_       = true;
  bool   use_measure_   = false;
  std::string laplacian_type_;

  // Single forward-transform scale factor: nx_ * ny_ * (2 * nz_).
  // Used to normalize after the inverse transform and to scale the k=0 (mean) source term.
  double transform_scale_ = 1.0;

  int_t  iter_ = 0;
  real_t time_ = 0.0;
  real_t dt_   = 0.0;

  std::vector<double> old_field_;
  std::vector<double> field_;
  std::vector<double> trial_field_;
  std::vector<double> nonlinear_;
  std::vector<ModeSymbol> symbols_;

  std::vector<std::complex<double>> buffer_;         // in-place x/y FFT + DCT(z) workspace
  std::vector<std::complex<double>> linear_hat_;      // transform of u
  std::vector<std::complex<double>> nonlinear_hat_;   // transform of u^3

  std::vector<double> re_plane_, im_plane_;           // scratch for splitting/merging during DCT(z)

  fftw_plan plan_x_fwd_ = nullptr, plan_x_bwd_ = nullptr;
  fftw_plan plan_y_fwd_ = nullptr, plan_y_bwd_ = nullptr;
  fftw_plan plan_z_fwd_ = nullptr, plan_z_bwd_ = nullptr;

  SolutionState& mutableSolutionState()
  {
    return const_cast<SolutionState&>(TimeIntegrator::getCurrentSolutionState());
  }

  SolutionState& mutableResidualState() { return const_cast<SolutionState&>(TimeIntegrator::getCurrentResidual()); }

  void initialize(TimeIntegrableRHS& rhs);
  void allocateFFTW();
  void buildSymbols();

  // transform: FFT(x) -> FFT(y) -> DCT(z), and its inverse in reverse order
  void forwardTransform(std::vector<std::complex<double>>& data);
  void inverseTransform(std::vector<std::complex<double>>& data);
  void dctZ(std::vector<double>& real_plane, bool inverse);
  void splitComplex(const std::vector<std::complex<double>>& in, std::vector<double>& re, std::vector<double>& im)
      const;
  void mergeComplex(
      const std::vector<double>& re, const std::vector<double>& im, std::vector<std::complex<double>>& out) const;

  void copyStateToField(const TimeIntegrableRHS& rhs, const SolutionState& state, std::vector<double>& field) const;
  void copyFieldToState(const TimeIntegrableRHS& rhs, const std::vector<double>& field, SolutionState& state) const;

  bool   takeStep(const std::vector<double>& input, double dt, std::vector<double>& output);
  void   projectMean(std::vector<double>& field) const;
  bool   fieldIsUsable(const std::vector<double>& field) const;
  double relativeUpdateNorm(const std::vector<double>& before, const std::vector<double>& after) const;
  double maxAbs(const std::vector<double>& field) const;
};

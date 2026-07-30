#include "fft_semi_implicit.hpp"

#include "program_options/program_options.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

int
fftMode(int i, int n)
{
  return (i <= n / 2) ? i : i - n;
}

size_t
linearIndex(int i, int j, int k, int ny, int nz)
{
  return (static_cast<size_t>(i) * static_cast<size_t>(ny) + static_cast<size_t>(j)) * static_cast<size_t>(nz)
         + static_cast<size_t>(k);
}

double
axisSymbol(int mode, int n, double h, const std::string& laplacian_type)
{
  if (laplacian_type == "spectral") {
    const double length = static_cast<double>(n) * h;
    const double k     = 2.0 * pi * static_cast<double>(mode) / length;
    return k * k;
  }

  if (laplacian_type == "fd" || laplacian_type == "finite_difference") {
    const double s = std::sin(pi * static_cast<double>(mode) / static_cast<double>(n));
    return 4.0 * s * s / (h * h);
  }

  Logger::get().FatalMessage("Unknown fft_si_laplacian option: " + laplacian_type);
  return 0.0;
}

real_t
finalTimeTolerance(real_t time, real_t final_time)
{
  return 64.0 * std::numeric_limits<real_t>::epsilon()
         * std::max({real_t(1.0), std::fabs(time), std::fabs(final_time)});
}

int
smallestFactor(int n)
{
  if (n % 2 == 0) {
    return 2;
  }
  for (int f = 3; f * f <= n; f += 2) {
    if (n % f == 0) {
      return f;
    }
  }
  return n;
}

void
fft1DRecursive(std::complex<double>* data, int n, bool inverse)
{
  if (n <= 1) {
    return;
  }

  const int radix = smallestFactor(n);
  const int m     = n / radix;
  const double sign = inverse ? 1.0 : -1.0;

  if (radix == n) {
    std::vector<std::complex<double>> out(n);
    for (int k = 0; k < n; ++k) {
      std::complex<double> sum(0.0, 0.0);
      for (int j = 0; j < n; ++j) {
        const double angle = sign * 2.0 * pi * static_cast<double>(j * k) / static_cast<double>(n);
        sum += data[j] * std::complex<double>(std::cos(angle), std::sin(angle));
      }
      out[k] = sum;
    }
    std::copy(out.begin(), out.end(), data);
    return;
  }

  std::vector<std::complex<double>> work(n);
  for (int j = 0; j < radix; ++j) {
    for (int q = 0; q < m; ++q) {
      work[j * m + q] = data[j + radix * q];
    }
    fft1DRecursive(work.data() + j * m, m, inverse);
  }

  for (int k = 0; k < n; ++k) {
    const int km = k % m;
    std::complex<double> sum(0.0, 0.0);
    for (int j = 0; j < radix; ++j) {
      const double angle = sign * 2.0 * pi * static_cast<double>(j * k) / static_cast<double>(n);
      sum += work[j * m + km] * std::complex<double>(std::cos(angle), std::sin(angle));
    }
    data[k] = sum;
  }
}

}  // namespace


FFTSemiImplicit::FFTSemiImplicit(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts)
    : TimeIntegrator(rhs, opts)
{
  initialize(rhs);
}


FFTSemiImplicit::~FFTSemiImplicit()
{
}


void
FFTSemiImplicit::initialize(TimeIntegrableRHS& rhs)
{
  dim_ = Options::get().dimension();
  Logger::get().FatalAssert(dim_ == 2 || dim_ == 3, "FFTSemiImplicit only supports dimension 2 or 3.");

  if (Options::get().bc_x() != "periodic" || Options::get().bc_y() != "periodic"
      || (dim_ == 3 && Options::get().bc_z() != "periodic")) {
    Logger::get().FatalMessage("FFTSemiImplicit currently requires periodic boundaries in every active direction.");
  }

  Logger::get().FatalAssert(rhs.nEquations() == 1, "FFTSemiImplicit expects one scalar concentration equation.");

  nx_ = Options::get().domain_resolution_x();
  ny_ = Options::get().domain_resolution_y();
  nz_ = (dim_ == 3) ? Options::get().domain_resolution_z() : 1;

  Logger::get().FatalAssert(nx_ > 0 && ny_ > 0 && nz_ > 0, "FFT grid dimensions must be positive.");

  n_total_ = static_cast<size_t>(nx_) * static_cast<size_t>(ny_) * static_cast<size_t>(nz_);
  Logger::get().FatalAssert(
      static_cast<size_t>(rhs.dofsPerEquation()) == n_total_,
      "FFTSemiImplicit expects dofsPerEquation == nx * ny * nz. This is true for periodic Cartesian CH grids.");
  {
    auto idx = read_access(rhs.interiorIndices());
    Logger::get().FatalAssert(
        static_cast<size_t>(idx.size()) == n_total_,
        "FFTSemiImplicit expects every periodic grid point to be an interior index.");
  }

  const double lx = Options::get().domain_x_end() - Options::get().domain_x_begin();
  Logger::get().FatalAssert(lx > 0.0, "domain_x_end must be larger than domain_x_begin.");

  dx_ = lx / static_cast<double>(nx_);
  dy_ = dx_;
  dz_ = dx_;

  mean_          = Options::get().ch_m();
  eps2_          = Options::get().ch_eps2();
  sigma_         = Options::get().ch_sigma();
  stabilization_ = Options::get().fft_si_stabilization();
  max_abs_c_     = Options::get().fft_si_max_abs_c();
  growth_factor_ = Options::get().fft_si_growth_factor();
  shrink_factor_ = Options::get().fft_si_shrink_factor();
  max_retries_   = Options::get().fft_si_max_retries();
  dealias_       = Options::get().fft_si_dealias();
  use_measure_   = Options::get().fft_si_plan_measure();
  laplacian_type_ = Options::get().fft_si_laplacian();

  Logger::get().FatalAssert(eps2_ >= 0.0, "eps2 must be non-negative.");
  Logger::get().FatalAssert(sigma_ >= 0.0, "sigma must be non-negative.");
  Logger::get().FatalAssert(stabilization_ >= 0.0, "fft_si_stabilization must be non-negative.");
  Logger::get().FatalAssert(max_abs_c_ > 0.0, "fft_si_max_abs_c must be positive.");
  Logger::get().FatalAssert(growth_factor_ >= 1.0, "fft_si_growth_factor must be >= 1.");
  Logger::get().FatalAssert(shrink_factor_ > 0.0 && shrink_factor_ < 1.0, "fft_si_shrink_factor must be in (0, 1).");
  Logger::get().FatalAssert(max_retries_ >= 0, "fft_si_max_retries must be non-negative.");

  old_field_.resize(n_total_);
  field_.resize(n_total_);
  trial_field_.resize(n_total_);
  nonlinear_.resize(n_total_);
  symbols_.resize(n_total_);

  allocateFFT();
  buildSymbols();

  std::ostringstream ss;
  ss << "Using fft_semi_implicit integrator: grid=" << nx_ << "x" << ny_ << "x" << nz_
     << ", dx=" << std::scientific << dx_
     << ", eps2=" << eps2_
     << ", sigma=" << sigma_
     << ", stabilization=" << stabilization_
     << ", laplacian=" << laplacian_type_
     << ", dealias=" << (dealias_ ? "true" : "false");
  Logger::get().InfoMessage(ss.str());
}


void
FFTSemiImplicit::allocateFFT()
{
  spatial_.assign(n_total_, std::complex<double>(0.0, 0.0));
  spectral_.assign(n_total_, std::complex<double>(0.0, 0.0));
  nonlinear_hat_.assign(n_total_, std::complex<double>(0.0, 0.0));
  line_buffer_.resize(static_cast<size_t>(std::max({nx_, ny_, nz_})));

  if (use_measure_) {
    Logger::get().WarningMessage("fft_si_plan_measure ignored by built-in FFT backend.");
  }
}


void
FFTSemiImplicit::fft3D(std::vector<std::complex<double>>& data, bool inverse)
{
  transformAxis(data, 2, inverse);
  transformAxis(data, 1, inverse);
  transformAxis(data, 0, inverse);
}


void
FFTSemiImplicit::transformAxis(std::vector<std::complex<double>>& data, int axis, bool inverse)
{
  if (axis == 0) {
    for (int j = 0; j < ny_; ++j) {
      for (int k = 0; k < nz_; ++k) {
        for (int i = 0; i < nx_; ++i) {
          line_buffer_[i] = data[linearIndex(i, j, k, ny_, nz_)];
        }
        fft1DRecursive(line_buffer_.data(), nx_, inverse);
        for (int i = 0; i < nx_; ++i) {
          data[linearIndex(i, j, k, ny_, nz_)] = line_buffer_[i];
        }
      }
    }
    return;
  }

  if (axis == 1) {
    for (int i = 0; i < nx_; ++i) {
      for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
          line_buffer_[j] = data[linearIndex(i, j, k, ny_, nz_)];
        }
        fft1DRecursive(line_buffer_.data(), ny_, inverse);
        for (int j = 0; j < ny_; ++j) {
          data[linearIndex(i, j, k, ny_, nz_)] = line_buffer_[j];
        }
      }
    }
    return;
  }

  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      for (int k = 0; k < nz_; ++k) {
        line_buffer_[k] = data[linearIndex(i, j, k, ny_, nz_)];
      }
      fft1DRecursive(line_buffer_.data(), nz_, inverse);
      for (int k = 0; k < nz_; ++k) {
        data[linearIndex(i, j, k, ny_, nz_)] = line_buffer_[k];
      }
    }
  }
}


void
FFTSemiImplicit::buildSymbols()
{
  for (int i = 0; i < nx_; ++i) {
    const int    mi  = fftMode(i, nx_);
    const double sx2 = axisSymbol(mi, nx_, dx_, laplacian_type_);
    for (int j = 0; j < ny_; ++j) {
      const int    mj  = fftMode(j, ny_);
      const double sy2 = axisSymbol(mj, ny_, dy_, laplacian_type_);
      for (int k = 0; k < nz_; ++k) {
        const int    mk  = fftMode(k, nz_);
        const double sz2 = (dim_ == 3) ? axisSymbol(mk, nz_, dz_, laplacian_type_) : 0.0;
        const size_t id  = linearIndex(i, j, k, ny_, nz_);

        const double k2 = sx2 + sy2 + sz2;
        symbols_[id].k2 = k2;
        symbols_[id].k4 = k2 * k2;

        const bool keep_x = std::abs(mi) <= nx_ / 3;
        const bool keep_y = std::abs(mj) <= ny_ / 3;
        const bool keep_z = (dim_ == 2) || (std::abs(mk) <= nz_ / 3);
        symbols_[id].keep_nonlinear = keep_x && keep_y && keep_z;
      }
    }
  }
}


void
FFTSemiImplicit::copyStateToField(
    const TimeIntegrableRHS& rhs, const SolutionState& state, std::vector<double>& field) const
{
  auto idx = read_access(rhs.interiorIndices());
  auto sol = read_access(state.getVec(0));

  for (size_t p = 0; p < n_total_; ++p) {
    field[p] = static_cast<double>(sol[idx[p]]);
  }
}


void
FFTSemiImplicit::copyFieldToState(
    const TimeIntegrableRHS& rhs, const std::vector<double>& field, SolutionState& state) const
{
  auto idx = read_access(rhs.interiorIndices());
  auto sol = write_access(state.getVec(0));

  for (size_t p = 0; p < n_total_; ++p) {
    sol[idx[p]] = field[p];
  }
}


bool
FFTSemiImplicit::takeStep(const std::vector<double>& input, double dt, std::vector<double>& output)
{
  for (size_t p = 0; p < n_total_; ++p) {
    const double u = input[p];
    if (!std::isfinite(u)) {
      return false;
    }
    nonlinear_[p] = u * u * u;

    spatial_[p] = std::complex<double>(u, 0.0);
  }

  spectral_ = spatial_;
  fft3D(spectral_, false);

  for (size_t p = 0; p < n_total_; ++p) {
    nonlinear_hat_[p] = spectral_[p];
    spatial_[p] = std::complex<double>(nonlinear_[p], 0.0);
  }

  spectral_ = spatial_;
  fft3D(spectral_, false);

  for (size_t p = 0; p < n_total_; ++p) {
    if (dealias_ && !symbols_[p].keep_nonlinear) {
      spectral_[p] = std::complex<double>(0.0, 0.0);
    }

    const double k2 = symbols_[p].k2;
    const double k4 = symbols_[p].k4;
    const double stabilizer = stabilization_ * k2;
    const double denominator = 1.0 - dt * k2 + dt * eps2_ * k4 + dt * sigma_ + dt * stabilizer;

    std::complex<double> rhs = (1.0 + dt * stabilizer) * nonlinear_hat_[p] - dt * k2 * spectral_[p];

    if (p == 0) {
      rhs += dt * sigma_ * mean_ * static_cast<double>(n_total_);
    }

    spectral_[p] = rhs / denominator;
  }

  spatial_ = spectral_;
  fft3D(spatial_, true);

  const double inv_n = 1.0 / static_cast<double>(n_total_);
  for (size_t p = 0; p < n_total_; ++p) {
    output[p] = spatial_[p].real() * inv_n;
  }

  projectMean(output);
  return fieldIsUsable(output);
}


void
FFTSemiImplicit::projectMean(std::vector<double>& field) const
{
  double sum = 0.0;
  for (double v : field) {
    sum += v;
  }

  const double correction = sum / static_cast<double>(n_total_) - mean_;
  for (double& v : field) {
    v -= correction;
  }
}


bool
FFTSemiImplicit::fieldIsUsable(const std::vector<double>& field) const
{
  for (double v : field) {
    if (!std::isfinite(v) || std::fabs(v) > max_abs_c_) {
      return false;
    }
  }
  return true;
}


double
FFTSemiImplicit::relativeUpdateNorm(const std::vector<double>& before, const std::vector<double>& after) const
{
  double update2 = 0.0;
  double state2  = 0.0;

  for (size_t p = 0; p < n_total_; ++p) {
    const double du = after[p] - before[p];
    update2 += du * du;
    state2 += after[p] * after[p];
  }

  const double denom = std::max(1.0e-30, std::sqrt(state2 / static_cast<double>(n_total_)));
  return std::sqrt(update2 / static_cast<double>(n_total_)) / denom;
}


double
FFTSemiImplicit::maxAbs(const std::vector<double>& field) const
{
  double value = 0.0;
  for (double v : field) {
    value = std::max(value, std::fabs(v));
  }
  return value;
}


void
FFTSemiImplicit::solve(TimeIntegrableRHS& rhs, InitialConditions& initial_conditions)
{
  auto& state    = mutableSolutionState();
  auto& residual = mutableResidualState();

  initial_conditions.set(rhs, state);
  copyStateToField(rhs, state, field_);
  projectMean(field_);
  copyFieldToState(rhs, field_, state);

  const auto& time_opts = getTimeOptions();
  time_ = time_opts.t0_;
  dt_   = time_opts.dt_initial_;
  iter_ = time_opts.initial_step_;

  if (time_opts.use_cfl_time_step_) {
    Logger::get().FatalMessage("FFTSemiImplicit does not use CFL time stepping. Use time_step_size/adaptive limits.");
  }

  if (time_opts.use_adaptive_time_step_) {
    dt_ = std::min(std::max(dt_, time_opts.min_time_step_size_), time_opts.max_time_step_size_);
  }

  int consecutive_rejections = 0;

  while (true) {
    const real_t tol = finalTimeTolerance(time_, time_opts.tfinal_);
    const real_t remaining = time_opts.tfinal_ - time_;

    if (remaining <= tol) {
      time_ = time_opts.tfinal_;
      Logger::get().InfoMessage("fft_semi_implicit reached final time.");
      break;
    }

    if (time_opts.max_time_steps_ && iter_ >= time_opts.max_time_steps_) {
      Logger::get().InfoMessage("fft_semi_implicit reached max_time_steps.");
      break;
    }

    if (dt_ > remaining || remaining - dt_ <= tol) {
      dt_ = remaining;
    }

    old_field_ = field_;

    bool accepted = takeStep(old_field_, dt_, trial_field_);

    if (!accepted) {
      ++consecutive_rejections;
      if (consecutive_rejections > max_retries_) {
        Logger::get().FatalMessage("fft_semi_implicit exceeded fft_si_max_retries at time=" + std::to_string(time_));
      }

      dt_ *= shrink_factor_;
      if (dt_ < time_opts.min_time_step_size_) {
        Logger::get().FatalMessage("fft_semi_implicit dt fell below min_time_step_size: " + std::to_string(dt_));
      }

      Logger::get().WarningMessage(
          "fft_semi_implicit rejected step; retrying at dt=" + std::to_string(dt_));
      continue;
    }

    consecutive_rejections = 0;
    field_.swap(trial_field_);
    time_ += dt_;
    ++iter_;

    copyFieldToState(rhs, field_, state);
    rhs.evalRHSImpl(state, time_, residual);

    const double rel_update = relativeUpdateNorm(old_field_, field_);
    std::ostringstream ss;
    ss << "fft_semi_implicit accepted iter=" << iter_ << ", time=" << std::scientific << time_
       << ", dt=" << dt_ << ", rel_update=" << rel_update << ", max_abs_c=" << maxAbs(field_);
    Logger::get().TraceMessage(ss.str());

    notifyObservers(Event::TimeStepComplete);
    notifyObservers(Event::SolutionUpdate);

    if (time_opts.converged_rel_tol_ > 0.0 && rel_update < time_opts.converged_rel_tol_) {
      Logger::get().InfoMessage("fft_semi_implicit exiting due to converged_rel_tol on relative update.");
      break;
    }

    if (time_opts.use_adaptive_time_step_) {
      dt_ = std::min(time_opts.max_time_step_size_, dt_ * growth_factor_);
      dt_ = std::max(time_opts.min_time_step_size_, dt_);
    }
  }
}


static auto fft_semi_implicit_instance = FactoryRegistry<TimeIntegrator>::get().add(
    "fft_semi_implicit",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<FFTSemiImplicit>(rhs, opts);
    });

static auto fft_si_instance = FactoryRegistry<TimeIntegrator>::get().add(
    "fft_si",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<FFTSemiImplicit>(rhs, opts);
    });

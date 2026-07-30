#include "dct_semi_implicit.hpp"

#include "program_options/program_options.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

size_t
linearIndex(int i, int j, int k, int ny, int nz)
{
  return (static_cast<size_t>(i) * static_cast<size_t>(ny) + static_cast<size_t>(j)) * static_cast<size_t>(nz)
         + static_cast<size_t>(k);
}

// Periodic-axis mode unwrap: FFT bin index -> signed wavenumber index.
int
fftMode(int i, int n)
{
  return (i <= n / 2) ? i : i - n;
}

// Eigenvalue of -Laplacian_1D restricted to a single periodic axis, for FFT mode `mode`.
double
periodicAxisSymbol(int mode, int n, double h, const std::string& laplacian_type)
{
  if (laplacian_type == "spectral") {
    const double length = static_cast<double>(n) * h;
    const double k      = 2.0 * pi * static_cast<double>(mode) / length;
    return k * k;
  }

  if (laplacian_type == "fd" || laplacian_type == "finite_difference") {
    const double s = std::sin(pi * static_cast<double>(mode) / static_cast<double>(n));
    return 4.0 * s * s / (h * h);
  }

  Logger::get().FatalMessage("Unknown fft_si_laplacian option: " + laplacian_type);
  return 0.0;
}

// Eigenvalue of -Laplacian_1D restricted to a single Neumann axis (cell-centered, ghost-symmetric),
// for DCT-II mode `mode` (mode = 0 .. n-1, no wraparound needed -- cosine modes are already
// non-negative frequencies).
//
// NOTE: this is the standard DCT-II eigenvalue for the 3-point Neumann discrete Laplacian.
// Validate against a small scipy.fft.dctn/idctn reference before trusting production runs,
// per the plan's step 7.
double
neumannAxisSymbol(int mode, int n, double h, const std::string& laplacian_type)
{
  if (laplacian_type == "spectral") {
    const double length = static_cast<double>(n) * h;
    const double k       = pi * static_cast<double>(mode) / length;
    return k * k;
  }

  if (laplacian_type == "fd" || laplacian_type == "finite_difference") {
    const double s = std::sin(pi * static_cast<double>(mode) / (2.0 * static_cast<double>(n)));
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

}  // namespace


DctSemiImplicit::DctSemiImplicit(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts)
    : TimeIntegrator(rhs, opts)
{
  initialize(rhs);
}


DctSemiImplicit::~DctSemiImplicit()
{
  if (plan_x_fwd_) fftw_destroy_plan(plan_x_fwd_);
  if (plan_x_bwd_) fftw_destroy_plan(plan_x_bwd_);
  if (plan_y_fwd_) fftw_destroy_plan(plan_y_fwd_);
  if (plan_y_bwd_) fftw_destroy_plan(plan_y_bwd_);
  if (plan_z_fwd_) fftw_destroy_plan(plan_z_fwd_);
  if (plan_z_bwd_) fftw_destroy_plan(plan_z_bwd_);
}


void
DctSemiImplicit::initialize(TimeIntegrableRHS& rhs)
{
  Logger::get().FatalAssert(Options::get().dimension() == 3, "DctSemiImplicit currently requires dimension 3.");

  Logger::get().FatalAssert(
      Options::get().bc_x() == "periodic" && Options::get().bc_y() == "periodic",
      "DctSemiImplicit requires periodic boundaries in x and y (use fft_si for fully periodic problems).");
  Logger::get().FatalAssert(
      Options::get().bc_z() == "neumann",
      "DctSemiImplicit requires Neumann boundaries in z (use fft_si for fully periodic problems).");

  Logger::get().FatalAssert(rhs.nEquations() == 1, "DctSemiImplicit expects one scalar concentration equation.");

  nx_ = Options::get().domain_resolution_x();
  ny_ = Options::get().domain_resolution_y();
  nz_ = Options::get().domain_resolution_z();

  Logger::get().FatalAssert(nx_ > 0 && ny_ > 0 && nz_ > 0, "DCT grid dimensions must be positive.");

  n_total_ = static_cast<size_t>(nx_) * static_cast<size_t>(ny_) * static_cast<size_t>(nz_);
  Logger::get().FatalAssert(
      static_cast<size_t>(rhs.dofsPerEquation()) == n_total_,
      "DctSemiImplicit expects dofsPerEquation == nx * ny * nz.");
  {
    auto idx = read_access(rhs.interiorIndices());
    Logger::get().FatalAssert(
        static_cast<size_t>(idx.size()) == n_total_,
        "DctSemiImplicit expects every grid point to be an interior index.");
  }

  const double lx = Options::get().domain_x_end() - Options::get().domain_x_begin();
  Logger::get().FatalAssert(lx > 0.0, "domain_x_end must be larger than domain_x_begin.");

  // Same convention as FFTSemiImplicit: a single uniform grid spacing derived from the x extent,
  // applied to all three axes (dy = dz = dx). If your grid is not cubic-cell, replace this with
  // real domain_y_end/domain_z_end accessors if/when they exist in your Options class.
  dx_ = lx / static_cast<double>(nx_);
  dy_ = dx_;
  dz_ = dx_;

  mean_ = Options::get().ch_m();

  // Two supported parameterizations, resolved by precedence:
  //   1) (eps, gamma), if both are explicitly set positive: eps2 = eps^2, sigma = 1/(eps^2*gamma^2).
  //      Convenient for one-off/manual runs (see run_sim_fft_cpu.sh examples).
  //   2) (eps2, sigma) directly, otherwise. This is what build_ensemble.py's generated SLURM
  //      script actually passes -- it stores (m, epsilon, gamma[, sigma]) per simulation in
  //      params.csv, but its run_one_sim() template computes eps2 = epsilon^2 in bash and invokes
  //      the simulator with --eps2/--sigma, not --eps/--gamma. --eps/--gamma both default to 0.0,
  //      so that path is never accidentally taken by an ensemble run.
  const double eps_in   = Options::get().ch_eps();
  const double gamma_in = Options::get().ch_gamma();
  const bool   use_eps_gamma = (eps_in > 0.0) && (gamma_in > 0.0);

  if (use_eps_gamma) {
    eps2_  = eps_in * eps_in;
    sigma_ = 1.0 / (eps2_ * gamma_in * gamma_in);
  } else {
    eps2_  = Options::get().ch_eps2();
    sigma_ = Options::get().ch_sigma();
    Logger::get().FatalAssert(
        eps2_ > 0.0 && sigma_ > 0.0,
        "dct_semi_implicit needs either (--eps > 0 and --gamma > 0), or (--eps2 > 0 and --sigma > 0).");
  }

  stabilization_  = Options::get().fft_si_stabilization();
  max_abs_c_      = Options::get().fft_si_max_abs_c();
  growth_factor_  = Options::get().fft_si_growth_factor();
  shrink_factor_  = Options::get().fft_si_shrink_factor();
  max_retries_    = Options::get().fft_si_max_retries();
  dealias_        = Options::get().fft_si_dealias();
  use_measure_    = Options::get().fft_si_plan_measure();
  laplacian_type_ = Options::get().fft_si_laplacian();

  Logger::get().FatalAssert(stabilization_ >= 0.0, "fft_si_stabilization must be non-negative.");
  Logger::get().FatalAssert(max_abs_c_ > 0.0, "fft_si_max_abs_c must be positive.");
  Logger::get().FatalAssert(growth_factor_ >= 1.0, "fft_si_growth_factor must be >= 1.");
  Logger::get().FatalAssert(shrink_factor_ > 0.0 && shrink_factor_ < 1.0, "fft_si_shrink_factor must be in (0, 1).");
  Logger::get().FatalAssert(max_retries_ >= 0, "fft_si_max_retries must be non-negative.");

  transform_scale_ = static_cast<double>(nx_) * static_cast<double>(ny_) * 2.0 * static_cast<double>(nz_);

  old_field_.resize(n_total_);
  field_.resize(n_total_);
  trial_field_.resize(n_total_);
  nonlinear_.resize(n_total_);
  symbols_.resize(n_total_);

  allocateFFTW();
  buildSymbols();

  std::ostringstream ss;
  ss << "Using dct_semi_implicit integrator: grid=" << nx_ << "x" << ny_ << "x" << nz_
     << " (bc: periodic, periodic, neumann), dx=" << std::scientific << dx_ << ", m=" << mean_;
  if (use_eps_gamma) {
    ss << ", eps=" << eps_in << ", gamma=" << gamma_in << " -> eps2=" << eps2_ << ", sigma=" << sigma_;
  } else {
    ss << " (from --eps2/--sigma directly)" << ", eps2=" << eps2_ << ", sigma=" << sigma_;
  }
  ss
     << ", stabilization=" << stabilization_ << ", laplacian=" << laplacian_type_
     << ", dealias=" << (dealias_ ? "true" : "false")
     << ", plan_measure=" << (use_measure_ ? "true (FFTW_MEASURE)" : "false (FFTW_ESTIMATE)");
  Logger::get().InfoMessage(ss.str());
  if (use_measure_) {
    // Unlike fft_si (which warns and ignores this flag with its in-house backend), this FFTW
    // backend genuinely honors fft_si_plan_measure: plan creation will run real timed trials
    // over several candidate FFT strategies, which costs extra startup time per unique
    // (nx, ny, nz) combination. Worth knowing if you're sweeping many small ensemble runs.
    Logger::get().InfoMessage(
        "dct_semi_implicit: fft_si_plan_measure=true will spend extra startup time benchmarking "
        "FFTW strategies for this grid size.");
  }
}


void
DctSemiImplicit::allocateFFTW()
{
  buffer_.assign(n_total_, std::complex<double>(0.0, 0.0));
  linear_hat_.assign(n_total_, std::complex<double>(0.0, 0.0));
  nonlinear_hat_.assign(n_total_, std::complex<double>(0.0, 0.0));
  re_plane_.assign(n_total_, 0.0);
  im_plane_.assign(n_total_, 0.0);

  const unsigned flags = use_measure_ ? FFTW_MEASURE : FFTW_ESTIMATE;
  auto* data = reinterpret_cast<fftw_complex*>(buffer_.data());

  // --- x axis: transform dim stride = ny_*nz_ (x is the slowest index); howmany dims: j (stride
  //     nz_, count ny_) and k (stride 1, count nz_).
  {
    fftw_iodim dims[1]        = {{nx_, ny_ * nz_, ny_ * nz_}};
    fftw_iodim howmany_dims[2] = {{ny_, nz_, nz_}, {nz_, 1, 1}};
    plan_x_fwd_ = fftw_plan_guru_dft(1, dims, 2, howmany_dims, data, data, FFTW_FORWARD, flags);
    plan_x_bwd_ = fftw_plan_guru_dft(1, dims, 2, howmany_dims, data, data, FFTW_BACKWARD, flags);
  }

  // --- y axis: transform dim stride = nz_ (y is the middle index); howmany dims: i (stride
  //     ny_*nz_, count nx_) and k (stride 1, count nz_).
  {
    fftw_iodim dims[1]        = {{ny_, nz_, nz_}};
    fftw_iodim howmany_dims[2] = {{nx_, ny_ * nz_, ny_ * nz_}, {nz_, 1, 1}};
    plan_y_fwd_ = fftw_plan_guru_dft(1, dims, 2, howmany_dims, data, data, FFTW_FORWARD, flags);
    plan_y_bwd_ = fftw_plan_guru_dft(1, dims, 2, howmany_dims, data, data, FFTW_BACKWARD, flags);
  }

  // --- z axis: real-to-real DCT-II (forward, REDFT10) / DCT-III (backward, REDFT01), applied to
  //     a real plane of shape (nx_, ny_, nz_). Transform dim stride = 1, count nz_; howmany dims:
  //     i (stride ny_*nz_, count nx_) and j (stride nz_, count ny_).
  //
  //     NOTE: this plan is later executed (via dctZ()) on two *different* vectors -- re_plane_ and
  //     im_plane_ -- not just the one it was created with. FFTW's "new-array execute" only
  //     guarantees correct results across different arrays if they share the same alignment as
  //     the array used at plan-creation time (see FFTW docs on SIMD/alignment). std::vector<double>
  //     gives no such guarantee between two separate vectors, so we force FFTW_UNALIGNED here to
  //     avoid a correctness footgun that would otherwise depend on how the system's FFTW was built
  //     (e.g. AVX-enabled builds are far more likely to trip on this than plain SSE ones). The x/y
  //     plans below don't need this: they're always executed on the exact buffer_ pointer they
  //     were planned with, so their alignment is unconditionally correct.
  {
    auto* real_data = re_plane_.data();
    fftw_iodim dims[1]        = {{nz_, 1, 1}};
    fftw_iodim howmany_dims[2] = {{nx_, ny_ * nz_, ny_ * nz_}, {ny_, nz_, nz_}};
    fftw_r2r_kind fwd_kind = FFTW_REDFT10;
    fftw_r2r_kind bwd_kind = FFTW_REDFT01;
    const unsigned z_flags = flags | FFTW_UNALIGNED;
    plan_z_fwd_ = fftw_plan_guru_r2r(1, dims, 2, howmany_dims, real_data, real_data, &fwd_kind, z_flags);
    plan_z_bwd_ = fftw_plan_guru_r2r(1, dims, 2, howmany_dims, real_data, real_data, &bwd_kind, z_flags);
  }

  Logger::get().FatalAssert(
      plan_x_fwd_ && plan_x_bwd_ && plan_y_fwd_ && plan_y_bwd_ && plan_z_fwd_ && plan_z_bwd_,
      "DctSemiImplicit failed to create one or more FFTW plans.");
}


void
DctSemiImplicit::splitComplex(
    const std::vector<std::complex<double>>& in, std::vector<double>& re, std::vector<double>& im) const
{
  for (size_t p = 0; p < n_total_; ++p) {
    re[p] = in[p].real();
    im[p] = in[p].imag();
  }
}


void
DctSemiImplicit::mergeComplex(
    const std::vector<double>& re, const std::vector<double>& im, std::vector<std::complex<double>>& out) const
{
  for (size_t p = 0; p < n_total_; ++p) {
    out[p] = std::complex<double>(re[p], im[p]);
  }
}


void
DctSemiImplicit::dctZ(std::vector<double>& real_plane, bool inverse)
{
  fftw_plan plan = inverse ? plan_z_bwd_ : plan_z_fwd_;
  fftw_execute_r2r(plan, real_plane.data(), real_plane.data());
}


void
DctSemiImplicit::forwardTransform(std::vector<std::complex<double>>& data)
{
  auto* raw = reinterpret_cast<fftw_complex*>(data.data());
  fftw_execute_dft(plan_x_fwd_, raw, raw);
  fftw_execute_dft(plan_y_fwd_, raw, raw);

  splitComplex(data, re_plane_, im_plane_);
  dctZ(re_plane_, /*inverse=*/false);
  dctZ(im_plane_, /*inverse=*/false);
  mergeComplex(re_plane_, im_plane_, data);
}


void
DctSemiImplicit::inverseTransform(std::vector<std::complex<double>>& data)
{
  splitComplex(data, re_plane_, im_plane_);
  dctZ(re_plane_, /*inverse=*/true);
  dctZ(im_plane_, /*inverse=*/true);
  mergeComplex(re_plane_, im_plane_, data);

  auto* raw = reinterpret_cast<fftw_complex*>(data.data());
  fftw_execute_dft(plan_y_bwd_, raw, raw);
  fftw_execute_dft(plan_x_bwd_, raw, raw);
  // NOTE: caller is responsible for the final 1/transform_scale_ normalization.
}


void
DctSemiImplicit::buildSymbols()
{
  for (int i = 0; i < nx_; ++i) {
    const int    mi  = fftMode(i, nx_);
    const double sx2 = periodicAxisSymbol(mi, nx_, dx_, laplacian_type_);
    for (int j = 0; j < ny_; ++j) {
      const int    mj  = fftMode(j, ny_);
      const double sy2 = periodicAxisSymbol(mj, ny_, dy_, laplacian_type_);
      for (int k = 0; k < nz_; ++k) {
        // DCT-II modes are already non-negative frequency indices; no wraparound needed.
        const double sz2 = neumannAxisSymbol(k, nz_, dz_, laplacian_type_);
        const size_t id  = linearIndex(i, j, k, ny_, nz_);

        const double k2 = sx2 + sy2 + sz2;
        symbols_[id].k2 = k2;
        symbols_[id].k4 = k2 * k2;

        const bool keep_x = std::abs(mi) <= nx_ / 3;
        const bool keep_y = std::abs(mj) <= ny_ / 3;
        const bool keep_z = k <= (2 * nz_) / 3;
        symbols_[id].keep_nonlinear = keep_x && keep_y && keep_z;
      }
    }
  }
}


void
DctSemiImplicit::copyStateToField(
    const TimeIntegrableRHS& rhs, const SolutionState& state, std::vector<double>& field) const
{
  auto idx = read_access(rhs.interiorIndices());
  auto sol = read_access(state.getVec(0));

  for (size_t p = 0; p < n_total_; ++p) {
    field[p] = static_cast<double>(sol[idx[p]]);
  }
}


void
DctSemiImplicit::copyFieldToState(
    const TimeIntegrableRHS& rhs, const std::vector<double>& field, SolutionState& state) const
{
  auto idx = read_access(rhs.interiorIndices());
  auto sol = write_access(state.getVec(0));

  for (size_t p = 0; p < n_total_; ++p) {
    sol[idx[p]] = field[p];
  }
}


bool
DctSemiImplicit::takeStep(const std::vector<double>& input, double dt, std::vector<double>& output)
{
  for (size_t p = 0; p < n_total_; ++p) {
    const double u = input[p];
    if (!std::isfinite(u)) {
      return false;
    }
    nonlinear_[p] = u * u * u;
    buffer_[p]    = std::complex<double>(u, 0.0);
  }

  forwardTransform(buffer_);
  linear_hat_ = buffer_;

  for (size_t p = 0; p < n_total_; ++p) {
    buffer_[p] = std::complex<double>(nonlinear_[p], 0.0);
  }
  forwardTransform(buffer_);
  nonlinear_hat_ = buffer_;

  for (size_t p = 0; p < n_total_; ++p) {
    if (dealias_ && !symbols_[p].keep_nonlinear) {
      nonlinear_hat_[p] = std::complex<double>(0.0, 0.0);
    }

    const double k2          = symbols_[p].k2;
    const double k4          = symbols_[p].k4;
    const double stabilizer  = stabilization_ * k2;
    const double denominator = 1.0 - dt * k2 + dt * eps2_ * k4 + dt * sigma_ + dt * stabilizer;

    std::complex<double> rhs = (1.0 + dt * stabilizer) * linear_hat_[p] - dt * k2 * nonlinear_hat_[p];

    if (p == 0) {
      rhs += dt * sigma_ * mean_ * transform_scale_;
    }

    buffer_[p] = rhs / denominator;
  }

  inverseTransform(buffer_);

  const double inv_scale = 1.0 / transform_scale_;
  for (size_t p = 0; p < n_total_; ++p) {
    output[p] = buffer_[p].real() * inv_scale;
  }

  projectMean(output);
  return fieldIsUsable(output);
}


void
DctSemiImplicit::projectMean(std::vector<double>& field) const
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
DctSemiImplicit::fieldIsUsable(const std::vector<double>& field) const
{
  for (double v : field) {
    if (!std::isfinite(v) || std::fabs(v) > max_abs_c_) {
      return false;
    }
  }
  return true;
}


double
DctSemiImplicit::relativeUpdateNorm(const std::vector<double>& before, const std::vector<double>& after) const
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
DctSemiImplicit::maxAbs(const std::vector<double>& field) const
{
  double value = 0.0;
  for (double v : field) {
    value = std::max(value, std::fabs(v));
  }
  return value;
}


void
DctSemiImplicit::solve(TimeIntegrableRHS& rhs, InitialConditions& initial_conditions)
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
    Logger::get().FatalMessage("DctSemiImplicit does not use CFL time stepping. Use time_step_size/adaptive limits.");
  }

  if (time_opts.use_adaptive_time_step_) {
    dt_ = std::min(std::max(dt_, time_opts.min_time_step_size_), time_opts.max_time_step_size_);
  }

  int consecutive_rejections = 0;

  while (true) {
    const real_t tol       = finalTimeTolerance(time_, time_opts.tfinal_);
    const real_t remaining = time_opts.tfinal_ - time_;

    if (remaining <= tol) {
      time_ = time_opts.tfinal_;
      Logger::get().InfoMessage("dct_semi_implicit reached final time.");
      break;
    }

    if (time_opts.max_time_steps_ && iter_ >= time_opts.max_time_steps_) {
      Logger::get().InfoMessage("dct_semi_implicit reached max_time_steps.");
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
        Logger::get().FatalMessage("dct_semi_implicit exceeded fft_si_max_retries at time=" + std::to_string(time_));
      }

      dt_ *= shrink_factor_;
      if (dt_ < time_opts.min_time_step_size_) {
        Logger::get().FatalMessage("dct_semi_implicit dt fell below min_time_step_size: " + std::to_string(dt_));
      }

      Logger::get().WarningMessage("dct_semi_implicit rejected step; retrying at dt=" + std::to_string(dt_));
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
    ss << "dct_semi_implicit accepted iter=" << iter_ << ", time=" << std::scientific << time_ << ", dt=" << dt_
       << ", rel_update=" << rel_update << ", max_abs_c=" << maxAbs(field_);
    Logger::get().TraceMessage(ss.str());

    notifyObservers(Event::TimeStepComplete);
    notifyObservers(Event::SolutionUpdate);

    if (time_opts.converged_rel_tol_ > 0.0 && rel_update < time_opts.converged_rel_tol_) {
      Logger::get().InfoMessage("dct_semi_implicit exiting due to converged_rel_tol on relative update.");
      break;
    }

    if (time_opts.use_adaptive_time_step_) {
      dt_ = std::min(time_opts.max_time_step_size_, dt_ * growth_factor_);
      dt_ = std::max(time_opts.min_time_step_size_, dt_);
    }
  }
}


static auto dct_semi_implicit_instance = FactoryRegistry<TimeIntegrator>::get().add(
    "dct_semi_implicit",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<DctSemiImplicit>(rhs, opts);
    });

static auto dct_si_instance = FactoryRegistry<TimeIntegrator>::get().add(
    "dct_si",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<DctSemiImplicit>(rhs, opts);
    });

# FFT Semi-Implicit Solver On CCR

This branch keeps the full maDGiCart project intact and adds one new time
integrator:

```text
--time_integrator fft_semi_implicit
```

Short alias:

```text
--time_integrator fft_si
```

## What Each Directory Does

Use three separate things on CCR:

```text
madgicart-fft/             source tree checked out to this branch
madgicart-fft-gcc-build/   CMake build directory containing the maDGiCart binary
madg-gcc-amd64.sif         container image with compilers/libraries
```

The source tree controls the code. The build directory controls the executable.
The wrapper controls which source/build directories the container sees.

## Get The Branch On CCR

```bash
cd /projects/academic/kreyes3/BNL-M2DT/rlee/maDGiCart-github
git clone https://github.com/RachelLee57/maDGiCart-CH-rachel.git madgicart-fft
cd madgicart-fft
git checkout agent/fft-semi-implicit-integrator
```

If you upload the folder through OnDemand instead of cloning, make sure the
uploaded folder is named `madgicart-fft` or set `MADG_SRC` when running the
wrapper.

## Build

Build inside the same container environment used by the existing CPU workflow:

```bash
cd /projects/academic/kreyes3/BNL-M2DT/rlee/maDGiCart-github
mkdir -p madgicart-fft-gcc-build

apptainer exec --env LC_ALL=C \
  --bind "$PWD/madgicart-fft:/M2DT/maDGiCart/maDGiCart-CH" \
  --bind "$PWD/madgicart-fft-gcc-build:/M2DT/maDGiCart/maDGiCart-CH-gcc-build" \
  /projects/academic/kreyes3/BNL-M2DT/maDGiCart/madg-gcc-amd64.sif \
  bash -lc 'cd /M2DT/maDGiCart/maDGiCart-CH-gcc-build && cmake /M2DT/maDGiCart/maDGiCart-CH -DMADG_USE_SERIAL=On && make -j'
```

If CCR uses `singularity` instead of `apptainer`, replace `apptainer` with
`singularity`.

## Run

The branch includes a wrapper:

```bash
/projects/academic/kreyes3/BNL-M2DT/rlee/maDGiCart-github/madgicart-fft/scripts/run_sim_fft_cpu.sh --help | grep -i fft
```

That command should show `fft_semi_implicit` or `fft_si`. If it does not, the
source branch or build directory is not the FFT-enabled one.

Example smoke run:

```bash
/projects/academic/kreyes3/BNL-M2DT/rlee/maDGiCart-github/madgicart-fft/scripts/run_sim_fft_cpu.sh \
  --dimension 3 \
  --time_integrator fft_si \
  --m 0.0 \
  --eps2 2.42e-05 \
  --sigma 291.4488109293305 \
  --domain_resolution_x 32 \
  --domain_resolution_y 32 \
  --domain_resolution_z 16 \
  --domain_x_end 1.0 \
  --bc_x periodic \
  --bc_y periodic \
  --bc_z periodic \
  --time_step_size 1e-5 \
  --use_adaptive_time_step true \
  --min_time_step_size 1e-7 \
  --max_time_step_size 1e-3 \
  --max_time_steps 1000000 \
  --final_time 1e-2 \
  --converged_rel_tol 0.0 \
  --fft_si_laplacian fd \
  --fft_si_stabilization 2.0 \
  --fft_si_dealias true \
  --solution_output_file /scratch/$USER/fft_si_smoke/out
```

For notebook-exact spectral behavior, use:

```text
--fft_si_laplacian spectral
--fft_si_stabilization 0
```

For production sweeps, start with:

```text
--fft_si_laplacian fd
--fft_si_stabilization 2.0
--use_adaptive_time_step true
--time_step_size 1e-5
--max_time_step_size 1e-3
```

## Important Limitations

The FFT solver is periodic-only and expects one scalar Cahn-Hilliard equation.
It will stop with a fatal message for Neumann/non-periodic boundaries.

It uses a built-in CPU mixed-radix FFT backend, so no FFTW module is required.

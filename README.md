Repository for IMSRG Emulator Research.

# TODO:

- [ ] Add SINDy module
- [ ] Add more plotting routines
  - [ ] Come up with experiment for parametric error surface?
  - [ ] Add singular value distribution plot?
- [x] Implement rEPI module in `emulate.py` interactive module
- [ ] Extend parametric emulator to 2D parametric problem
- [ ] Push example of input argument piping
- [ ] Implement routine for storing emulated system steps on file

# Purpose

The purpose here is to streamline experimentation of data-driven emulator techniques for the IMSRG (and, in principle, any other dynamical system). The project rationale is this: data-driven emulator techniques begin in the same way, which is to collect a number of snapshots of the evolving system, and then perform some sort of optimization on an operator that drives the system. In this way, we want to design a general interface where the user can pick the emulation technique they want, pointing to the location of the data, and run the emulator interactively with relevant parameters. We also want to include some routines for analyzing the emulation efficacy.

# Docs

Refer to `davis9ja.github.io/imsrg_emu/`.

# How to run

Running an emulator is simple. Call the "executable" (as much as we can define an executable in Python), `emulate.py`, with a subcommand to determine the emulator technique to run, and a positional argument that points to the data path (depedent on emulator technique).

Call `python emulate.py -h` at anytime, and for any subcommand, for a description of all possible positional arguments and optional keyword arguments.

### Example: standard DMD

    python emulate.py standard path/to/csv/data/file  --nobs 20 --trunc 6 --t0 0.0 --t1 20.0 --dt 0.05

### Example: parametric DMD

    python emulate.py parametric path/to/data/list path/to/param/list <testParm> --emuType rKOI  --nobs 20 --trunc 6 --t0 0.0 --t1 20.0 --dt 0.05

# How to import to your own code

Export `imsrg_emu/` to your $PYTHONPATH

       import imsrg_emu as ie

       dmdrkoi = ie.dmd_rkoi.DMD_rKOI()
       dmdrkoi.fit(data_list, param_arr, nobs_t, r)
       dmdrkoi.interp_dmd(test_param)
       result = dmdrkoi.predict(s_range, ds)


       dmdstd = ie.dmd_std.DMD_STD()
       dmdstd.fit(data_mat, nobs)
       result = dmdstd.predict(s_range, ds)

*Note: right now, the parametric DMD is implemented for only a 1D trajectory in parametric space; e.g. varying pairing strength $`g`$ in the pairing model. Need further testing for exploring a parametric surface, manifold, etc.

**Important note: DMD results are very sensitive to decimal precision. Just keep as many digits as you can (e.g. double precision). Make sure you know how many digits are in your data files.

# Koopman operator theory

The Koopman operator linearizes the dynamical system characeterized by a first-order ODE,

<!---
```math
\frac{d}{dt}\mathbf{x}(t) = f(\mathbf{x}(t)).
```
-->

$$\frac{d}{dt}\mathbf{x}(t) = f(\mathbf{x}(t)).$$

The discrete-time Koopman operator propagates a measurement function $`g`$, which evaluates the system $`x`$ at a step $`x_k`$, foward in time to a measurement of the next step $`x_{k+1}`$. Methodologically, we define the Koopman operator $`\mathcal{K}`$ such that

```math
\mathcal{K}g(x_k) = g(x_{k+1}),
```

where the dynamical step from $`k`$ to $`k+1`$ is governed by a flow map $`\mathbf{F}`$ such that $`x_{k+1} = \mathbf{F}(x_k)`$. Thus, the discrete-time Koopman operator is data-driven in that it's defined in a *measurement* basis that measures the flow map $`\mathbf{F}`$. The eigenfunctions of the Koopman operator completely characterize this flow map, such that

```math
\mathcal{K}\phi(x_k) = \lambda \phi(x_{k+1}).
```

The eigenfunctions of $`\mathcal{K}`$ quantifies the flow map, which in turns allows us to reproduce the dynamics the system by simply evaluating the eigenfunctions for every relevant timestep.

Refer to Brunton et al. 2021 (arXiv:2102.12086v2) for more information.

# Dynamic Mode Decomposition (DMD)

The DMD is an algorithm for approximating the discrete-time Koopman operator. The DMD operator $`A`$ is defined as the best-fit operator which propagates a single snapshot of the evolving dynamical system forward in time, such that

```math
Ax_k = x_{k+1}.
```

In order to find this operator, we start by collecting *linear* measurements of the evolving dynamical system into a matrix of $`N`$ snapshot columns $`\mathbf{\chi}`$. The DMD operator is constructed in this finite measurement basis via the offset matrices $`\mathbf{X}`$ and $`\mathbf{X}'`$, which contain columns $`1`$ through $`N-1`$ in $`\mathbf{\chi}`$ and $`2`$ through $`N`$ in $`\mathbf{\chi}`$, respectively, so that

```math
A\mathbf{X} = \mathbf{X}'.
```
The very simple solution to $`A`$ is just

```math
A = \mathbf{X}'\mathbf{X}^\dagger,
```

where $`^\dagger`$ is the Moore-Penrose pseudoinverse.

Refer to Brunton et al. 2021 (arXiv:2102.12086v2) for more information.

# Reduced Koopman Operator Interpolation (rKOI)

rKOI algorithm as presented in Huhn et al. 2022:

1. Collect $`N`$ data matrices, $`\chi_1,\dots,\chi_N`$ for corresponding parametric realizations $`\mu_1,\dots,\mu_N`$.
2. Compute singular vector matrices $`U_1,\dots,U_N`$.       
3. Compute DMD operators $`A_1,\dots,A_N`$ and mode amplitudes $`b_1,\dots,b_N`$.
4. For test parametric realization $`\mu_\theta`$:
   1. Interpolate $`U_\theta`$ on $`U_1, \dots, U_N`$
   2. Interpolate $`A_\theta`$ on $`A_1, \dots, A_N`$   
   3. Interpolate $`b_\theta`$ on $`b_1, \dots, b_N`$
   4. Compute eigendecomp $`A_\theta W_\theta = W_\theta \Lambda_\theta`$
   5. Compute DMD modes $`\Phi_\theta = U_\theta W_\theta`$
   6. Compute DMD expansion for $`\mu_\theta`$ system

Refer to Huhn et al. 2022 (arXiv:2204.12006v1) for more information.

# Reduced Eigenpair Interpolation (rEPI)

1. Collect $`N`$ data matrices, $`\chi_1,\dots,\chi_N`$ for corresponding parametric realizations $`\mu_1,\dots,\mu_N`$.
2. Compute singular vector matrices $`U_1,\dots,U_N`$.
3. Compute DMD operator eigendecomps $`W_1,\dots,W_N`$ and $`\Lambda_1,\dots,\Lambda_N`$, and mode amplitudes $`b_1,\dots,b_N`$.
4. For test parametric realization $`\mu_\theta`$
   1. Interpolate $`U_\theta`$ on $`U_1, \dots, U_N`$
   2. Interpolate $`W_\theta`$ on $`W_1, \dots, W_N`$
   3. Interpolate $`\Lambda_\theta`$ on $`\Lambda_1, \dots, \Lambda_N`$   
   4. Interpolate $`b_\theta`$ on $`b_1, \dots, b_N`$
   5. Compute eigendecomp $`A_\theta W_\theta = W_\theta \Lambda_\theta`$
   6. Compute DMD modes $`\Phi_\theta = U_\theta W_\theta`$
   7. Compute DMD expansion for $`\mu_\theta`$ system

# Sparse Identification of Nonlinear Dynamics (SINDy)

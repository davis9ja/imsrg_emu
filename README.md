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

*Note: right now, the parametric DMD is implemented for only a 1D trajectory in parametric space; e.g. varying pairing strength $g$ in the pairing model. Need further testing for exploring a parametric surface, manifold, etc.

**Important note: DMD results are very sensitive to decimal precision. Just keep as many digits as you can (e.g. double precision). Make sure you know how many digits are in your data files.

*Another note: if you have trouble viewing equations, visit [github.com/davis9ja/imsrg_emu](github.com/davis9ja/imsrg_emu) for the same repo README rendered with Mathjax.

# Koopman operator theory

The Koopman operator linearizes the dynamical system characeterized by a first-order ODE,

$$\frac{d}{dt}\mathbf{x}(t) = f(\mathbf{x}(t)).$$

The discrete-time Koopman operator propagates a measurement function $g$, which evaluates the system $x$ at a step $x_k$, foward in time to a measurement of the next step $x_{k+1}$. Methodologically, we define the Koopman operator $\mathcal{K}$ such that

$$\mathcal{K}g(x_k) = g(x_{k+1}),$$

where the dynamical step from $k$ to $k+1$ is governed by a flow map $\mathbf{F}$ such that $x_{k+1} = \mathbf{F}(x_k)$. Thus, the discrete-time Koopman operator is data-driven in that it's defined in a *measurement* basis that measures the flow map $\mathbf{F}$. The eigenfunctions of the Koopman operator completely characterize this flow map, such that

$$\mathcal{K}\phi(x_k) = \lambda \phi(x_{k+1}).$$

The eigenfunctions of $\mathcal{K}$ quantifies the flow map, which in turns allows us to reproduce the dynamics the system by simply evaluating the eigenfunctions for every relevant timestep.

Refer to Brunton et al. 2021 (arXiv:2102.12086v2) for more information.

# Dynamic Mode Decomposition (DMD)

The DMD is an algorithm for approximating the discrete-time Koopman operator. The DMD operator $A$ is defined as the best-fit operator which propagates a single snapshot of the evolving dynamical system forward in time, such that

$$Ax_k = x_{k+1}.$$

In order to find this operator, we start by collecting *linear* measurements of the evolving dynamical system into a matrix of $N$ snapshot columns $\mathbf{\chi}$. The DMD operator is constructed in this finite measurement basis via the offset matrices $\mathbf{X}$ and $\mathbf{X}'$, which contain columns $1$ through $N-1$ in $\mathbf{\chi}$ and $2$ through $N$ in $\mathbf{\chi}$, respectively, so that

$$A\mathbf{X} = \mathbf{X}'.$$

The very simple solution to $A$ is just

$$A = \mathbf{X}'\mathbf{X}^\dagger,$$

where $^\dagger$ is the Moore-Penrose pseudoinverse.

Refer to Brunton et al. 2021 (arXiv:2102.12086v2) for more information.

# Enforcing consistencies in grid points for training parametric interpolation methods

DMD emulation as a parametric problem is a difficult one, and quite complicated for our case. However, we have found that parametric emulation can generate DMD operators, and evolved Hamiltonians, to an accuracy on the same order as a direct computation of the DMD expansion. We present two algorithms for solving the parametric DMD problem, with some caveats specific to our case. 

One caveat is that we enforce the same physical constaints on an emulated operator as we do the grid points (which are direct calculations of the DMD operator for varying parametric realizations). This is because constraints are enforced on objects which are *computed from* the DMD operator; therefore, they do not carry into the interpolation. These contraints include positive, real DMD eigenvalues, ordering the eigenvalues from largest to smallest, and forcing an eigenvalue ceiling of 1. The ordering is also imposed on the DMD eigenvectors. For the grid points ONLY, this ordering is imposed on the amplitudes. 

Another caveat is that we only impose ordering consistency on emulated values which are *computed from interpolated objects*. This includes the DMD eigenvalues and DMD eigenvectors. We DO NOT need to order the DMD mode amplitudes, because the DMD mode amplitudes are directly interpolated---in other words, the DMD mode amplitudes are not calculated from the emulated values.

The final caveat is that we impose a sign convention based on the first data point which appears in the interpolation training set. In principle, we could pick a random training point; the algorithm does not because it sequentially computes the DMD training information across the parameter range. The sign convention comes from the full DMD eigenvectors given by $\Phi$. We pick e.g. the first training point eigenvectors $\Phi_1$, and then calculate the sign overlap $O_i = \Phi_1 \cdot \Phi_i$ with all subsequent $\Phi_i$ where $i\neq 1$. The sign overlaps tell us where the $i$-th set of eigenvectors are negative (direction) with respect to $\Phi_1$. Multiplying $\Phi_i \cdot O_i$ ensures that the columns of $\Phi_i$ have the same sign as $\Phi_1$, which guarantees that the DMD mode amplitudes across all training and emulated points vary smoothly with the parameter.

All of these caveats contribute to ensuring that the spaces we interpolate in are *smooth with respect to the parameters*. Interpolation methods are most successful where the function that must be interpolated is smooth.

# Reduced Koopman Operator Interpolation (rKOI)

rKOI algorithm as presented in Huhn et al. 2022:

1. Collect $N$ data matrices, $\chi_1,\dots,\chi_N$ for corresponding parametric realizations $\mu_1,\dots,\mu_N$.
2. Compute singular vector matrices $U_1,\dots,U_N$.       
3. Compute DMD operators $A_1,\dots,A_N$ and mode amplitudes $b_1,\dots,b_N$.
4. For test parametric realization $\mu_\theta$:
   1. Interpolate $U_\theta$ on $U_1, \dots, U_N$
   2. Interpolate $A_\theta$ on $A_1, \dots, A_N$   
   3. Interpolate $b_\theta$ on $b_1, \dots, b_N$
   4. Compute eigendecomp $A_\theta W_\theta = W_\theta \Lambda_\theta$
   5. Compute DMD modes $\Phi_\theta = U_\theta W_\theta$
   6. Compute DMD expansion for $\mu_\theta$ system

Refer to Huhn et al. 2022 (arXiv:2204.12006v1) for more information.


# Reduced Eigenpair Interpolation (rEPI)

1. Collect $N$ data matrices, $\chi_1,\dots,\chi_N$ for corresponding parametric realizations $\mu_1,\dots,\mu_N$.
2. Compute singular vector matrices $U_1,\dots,U_N$.
3. Compute DMD operator eigendecomps $W_1,\dots,W_N$ and $\Lambda_1,\dots,\Lambda_N$, and mode amplitudes $b_1,\dots,b_N$.
4. For test parametric realization $\mu_\theta$
   1. Interpolate $U_\theta$ on $U_1, \dots, U_N$
   2. Interpolate $W_\theta$ on $W_1, \dots, W_N$
   3. Interpolate $\Lambda_\theta$ on $\Lambda_1, \dots, \Lambda_N$   
   4. Interpolate $b_\theta$ on $b_1, \dots, b_N$
   5. Compute eigendecomp $A_\theta W_\theta = W_\theta \Lambda_\theta$
   6. Compute DMD modes $\Phi_\theta = U_\theta W_\theta$
   7. Compute DMD expansion for $\mu_\theta$ system

# Sparse Identification of Nonlinear Dynamics (SINDy)

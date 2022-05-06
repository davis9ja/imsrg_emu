Repository for IMSRG Emulator Research.

# Purpose

The purpose here is to streamline experimentation of data-driven emulator techniques for the IMSRG (and, in principle, any other dynamical system). The project rationale is this: data-driven emulator techniques begin in the same way, which is to collect a number of snapshots of the evolving system, and then perform some sort of optimization on an operator that drives the system. In this way, we want to design a general interface where the user can pick the emulation technique they want, pointing to the location of the data, and run the emulator interactively with relevant parameters. We also want to include some routines for analyzing the emulation efficacy.

# How to run

Running an emulator is simple. Call the "executable" (as much as we can define an executable in Python), `emulator.py`, with a subcommand to determine the emulator technique to run, and a positional argument that points to the data path (depedent on emulator technique).

Call `python emulator.py -h` at anytime, and for any subcommand, for a description of all possible positional arguments and optional keyword arguments.

### Example: standard DMD

    python emulator.py standard path/to/csv/data/file  --nobs 20 --trunc 6 --t0 0.0 --t1 20.0 --dt 0.05

### Example: parametric DMD

    python emulator.py parametric /path/to/data/list --emuType rKOI  --nobs 20 --trunc 6 --t0 0.0 --t1 20.0 --dt 0.05

*Note: right now, the parametric DMD is implemented for only a 1D trajectory in parametric space; e.g. varying pairing strength $`g`$ in the pairing model. Need further testing for exploring a parametric surface, manifold, etc.

# Koopman operator theory

The Koopman operator linearizes the dynamical system characeterized by a first-order ODE,

```math
\frac{d}{dt}\mathbf{x}(t) = f(\mathbf{x}(t)).
```

The discrete-time Koopman operator propagates a measurement function $`g`$, which evaluates the system $`x`$ at a step $`x_k`$, foward in time to a measurement of the next step $`x_{k+1}`$. Methodologically, we define the Koopman operator $`\mathcal{K}`$ such that

```math
\mathcal{K}g(x_k) = g(x_{k+1}),
```

where the dynamical step from $`k`$ to $`k+1`$ is governed by a flow map $`\mathbf{F}`$ such that $`x_{k+1} = \mathbf{F}(x_k)`$. Thus, the discrete-time Koopman operator is data-driven in that it's defined in a *measurement* basis that measures the flow map $\mathbf{F}$. The eigenfunctions of the Koopman operator completely characterize this flow map, such that

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

In order to find this operator, we start by collecting *linear* measurements of the evolving dynamical system into a matrix of $`N`$ snapshot columns $`\mathbf{\Chi}`$. The DMD operator is constructed in this finite measurement basis via the offset matrices $`\mathbf{X}`$ and $`\mathbf{X}'`$, which contain columns $`1`$ through $`N-1`$ in $`\mathbf{\Chi}`$ and $`2`$ through $`N`$ in $`\mathbf{\Chi}`$, respectively, so that

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



# Reduced Eigenpair Interpolation (rEPI)

# Sparse Identification of Nonlinear Dynamics (SINDy)

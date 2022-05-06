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

The Koopman operator linearizes the dynamical system characeterized by the first-order ODE,

```math
\frac{d}{dt}\mathbf{x} = \mathbf{f}(\mathbf{x})
```

# Dynamic Mode Decomposition (DMD)

# Reduced Koopman Operator Interpolation (rKOI)

# Reduced Eigenpair Interpolation (rEPI)

# Sparse Identification of Nonlinear Dynamics (SINDy)

# CESDSOL
### Complex Equations — Simple Domain partial differential equations SOLver

CESDSOL is a C++ library providing capabilities to solve arbitrary partial differential equations with emphasis on highly nonlinear problems arising in field theory. It is battle-tested on numerous problems of finding stationary configurations for field-theoretic and condensed matter solitons, including gravitating ones, and black holes interacting with matter fields and of studying dynamics of soliton collisions.

## Features

- Stationary and explicit transient problems
- Arbitrary spatial dimension of problems
- Arbitrary number of equations
- Arbitrary order of derivatives, both spatial and temporal
- Arbitrarily-spaced cartesian grids
- Arbitrary boundary conditions
- Finite-difference discretization of arbitrary order of accuracy
- Field equations can be supported with discrete variable equations
- Variety of linear, nonlinear and ODE solvers
- Intel MKL as linear algebra backend
- Provided with Wolfram Language interface for generation of problem descriptor and reading of solver output

CESDSOL is developed in a very modular fashion allowing for extensions in different directions. Degrees of freedom include:

- New problem types
- New grid types
- New discretization types
- New math backends
- New solvers

## Building

Use `cmake . -DMainPath=%path to .cpp with your problem definition%` to build.

## Using

Refer to examples for guide.

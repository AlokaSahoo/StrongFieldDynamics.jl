# StrongFieldDynamics

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlokaSahoo.github.io/StrongFieldDynamics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlokaSahoo.github.io/StrongFieldDynamics.jl/dev/)
[![Build Status](https://github.com/AlokaSahoo/StrongFieldDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AlokaSahoo/StrongFieldDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AlokaSahoo/StrongFieldDynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlokaSahoo/StrongFieldDynamics.jl)

**StrongFieldDynamics.jl** is a Julia package designed for the theoretical study of atomic ionization processes in the presence of intense laser fields. The package implements the strong field approximation (SFA) using the partial wave expansion approach, enabling accurate and efficient simulations of laser-atom interactions.

## Overview

StrongFieldDynamics.jl provides a flexible framework to model and analyze the ionization dynamics of atomic systems subjected to ultrashort, high-intensity laser pulses. By leveraging the partial wave expansion, the package allows for the decomposition of the electron wavefunction into angular momentum components, facilitating detailed studies of angular distributions and multiphoton processes.

## Key Features

- **Partial Wave Expansion:** Efficiently expands the continuum electron wavefunction in terms of spherical harmonics, capturing angular effects in ionization.
- **Strong Field Approximation (SFA):** Implements SFA for modeling direct ionization and rescattering processes.
- **Atomic Potential Effects:** Includes several approaches to account for the influence of the atomic potential on the ionized electron, beyond the plain SFA.
- **Customizable Laser Pulses:** Supports a variety of laser pulse shapes, durations, and polarization schemes.
- **Observables Calculation:** Computes ionization rates, photoelectron spectra, and angular distributions.
- **Extensible Framework:** Modular design for easy extension to new atomic systems or laser configurations.

## Applications

- Study of multiphoton and tunneling ionization in atoms.
- Analysis of angular and energy-resolved photoelectron spectra.
- Investigation of strong-field phenomena such as above-threshold ionization (ATI) and rescattering.
- Benchmarking and comparison with experimental data.

## Future Scope

Planned extensions include:
- High Harmonic Generation (HHG) modeling.
- Inclusion of more complex atomic and molecular targets.
- Support for additional strong-field processes and advanced numerical techniques.

## Documentation

For detailed theory, usage examples, and API reference, see the [documentation](https://AlokaSahoo.github.io/StrongFieldDynamics.jl/stable/).
Also, detailed documentation is available in the docstrings of individual functions. Use `?function_name` in the Julia REPL for help.

## Installation

```julia
using Pkg
Pkg.add("StrongFieldDynamics")
```

## Getting Started

After installation, see the documentation for usage examples and further instructions.

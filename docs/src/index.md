```@meta
CurrentModule = StrongFieldDynamics
```

# StrongFieldDynamics

Documentation for [StrongFieldDynamics](https://github.com/AlokaSahoo/StrongFieldDynamics.jl).

**StrongFieldDynamics.jl** is a Julia package designed for the theoretical study of atomic ionization processes in the presence of intense laser fields. The package implements the strong field approximation (SFA) using the partial wave expansion approach, enabling accurate and efficient simulations of laser-atom interactions.

## Overview

The package provides a flexible framework to model and analyze the ionization dynamics of atomic systems subjected to ultrashort, high-intensity laser pulses. By leveraging the partial wave expansion, StrongFieldDynamics.jl allows for the decomposition of the electron wavefunction into angular momentum components, facilitating detailed studies of angular distributions and multiphoton processes.

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

## Documentation Structure

- [Theory](theory.md): Detailed background on the strong field approximation, partial wave expansion, and the physical models implemented.

## Getting Started

To get started, see the [installation instructions](https://github.com/AlokaSahoo/StrongFieldDynamics.jl#installation) and the usage examples in the documentation.

```@index
```

```@autodocs
Modules = [StrongFieldDynamics]
```

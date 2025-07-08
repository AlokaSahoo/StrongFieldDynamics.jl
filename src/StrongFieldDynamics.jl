"""
    StrongFieldDynamics

A Julia package for strong field atomic physics calculations using the Strong Field Approximation (SFA).

This package provides tools for:
- Computing electron wavefunctions in strong laser fields
- Calculating photoionization and high harmonic generation amplitudes
- Analyzing pulse shapes and field integrals
- Unit conversions for atomic physics calculations
- Visualization of strong field dynamics

# Main Components
- Amplitude calculations for various strong field processes
- Electron wavefunction computations in intense laser fields
- Pulse shape analysis and field integration utilities
- Atomic units and unit conversion functions
- Plotting utilities for visualization

# Example
```julia
using StrongFieldDynamics
```
"""
module StrongFieldDynamics

# Core type definitions and utility functions
include("types-fucntions.jl")

# Amplitude calculations for strong field processes
include("amplitudes.jl")

# Electron wavefunction computations in strong fields
include("electron-wavefunction.jl")

# Plotting and visualization utilities
include("plots.jl")

# Pulse shape analysis and field integration
include("pulse-shape-integral.jl")

# Unit conversions and atomic units
include("units.jl")

end
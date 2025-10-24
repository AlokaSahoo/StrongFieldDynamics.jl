using Base.Threads
using Distributed
using Printf
using ProgressMeter
using OffsetArrays
using SpecialFunctions
using SphericalHarmonics
using WignerSymbols

export ClebschGordan, Ylm
export PulseEnvelope, Sin2, Gaussian, Rectangular
export Pulse, AtomicElectron, ContinuumElectron, PartialWave, Settings
export ContinuumSolution, Bessel, Distorted
export Gauge, VelocityGauge, LengthGauge
export IonizationScheme, Hydrogenic, Atomic
export AngularDistribution, EnergyDistribution, MomentumDistribution
export compute_angular_distribution, compute_energy_distribution, compute_momentum_distribution
export Cartesian, spherical2cartesian



"""
    PulseEnvelope

Enumeration of different pulse envelope types for laser pulses.

# Variants
- `Sin2`: Sin-squared envelope, ideal for smooth turn-on/turn-off
- `Gaussian`: Gaussian envelope, commonly used in ultrafast optics
- `Rectangular`: Rectangular (flat-top) envelope, for constant field amplitude

# Usage
```julia
# Different envelope types
sin2_pulse = Pulse(I₀=1e14, λ=800, cycles=6, envelope=Sin2)
gauss_pulse = Pulse(I₀=1e14, λ=800, cycles=6, envelope=Gaussian)
rect_pulse = Pulse(I₀=1e14, λ=800, cycles=6, envelope=Rectangular)
```

# Physics
- **Sin2**: Envelope f(t) = sin²(ωt/2N) for smooth transitions
- **Gaussian**: Envelope f(t) = exp(-t²/2σ²) for natural pulse shapes
- **Rectangular**: Envelope f(t) = 1 for constant intensity during pulse
"""
@enum PulseEnvelope begin
    Sin2
    Gaussian
    Rectangular
end


"""
    Pulse

Represents a laser pulse with all necessary parameters for strong-field calculations.
All internal values are stored in atomic units.

# Fields
- `I₀::Float64`: Peak intensity (W/cm²) - stored in atomic units
- `A₀::Float64`: Peak vector potential amplitude (a.u.)
- `λ::Float64`: Wavelength (nm) - stored in atomic units
- `ω::Float64`: Angular frequency (a.u.)
- `cycles::Int64`: Number of optical cycles
- `Tp::Float64`: Pulse duration (a.u.)
- `Up::Float64`: Ponderomotive energy (a.u.)
- `f::Function`: Envelope function f(t)
- `cep::Float64`: Carrier-envelope phase (radians)
- `helicity::Int64`: Helicity (+1 or -1)
- `ϵ::Float64`: Ellipticity parameter (0=linear, 1=circular)
- `u::OffsetVector{Float64, Vector{Float64}}`: Polarization unit vector
- `Sv::Function`: Volkov-phase function
- `duration::Tuple{Float64, Float64}`: Time interval (t_start, t_end) for pulse duration

# Required Arguments
- `I₀`: Peak intensity (can include units, e.g., 1e14u"W/cm^2", or numeric value assuming W/cm²)
- `λ`: Wavelength (can include units, e.g., 800u"nm", or numeric value assuming nm)
- `cycles`: Number of optical cycles (positive integer)
- `envelope`: Envelope type (Sin2, Gaussian, Rectangular) or custom function f(t)
- `cep`: Carrier-envelope phase in radians (default: 0.0)
- `helicity`: Helicity (+1 or -1, default: +1)
- `ϵ`: Ellipticity parameter (0=linear, 1=circular, default: 0.0)
- `duration`: Time interval tuple (t_start, t_end) in a.u. (optional, auto-calculated if not provided)

# Optional Arguments
- `Sv`: Custom Volkov-phase function (default: t->0.0)

# Examples
```julia
# Basic pulse with automatic duration
pulse1 = Pulse(I₀=1e14, λ=800, cycles=6)

# Pulse with custom duration
pulse2 = Pulse(I₀=1e14, λ=800, cycles=6, duration=(-10.0, 50.0))

# Pulse with Gaussian envelope
pulse3 = Pulse(I₀=1e14, λ=800, cycles=6, envelope=Gaussian)
```
"""
struct Pulse
    I₀::Float64
    A₀::Float64
    λ::Float64
    ω::Float64
    np::Int64
    Tp::Float64
    Up::Float64
    f::Function
    cep::Float64
    helicity::Int64
    ϵ::Float64
    u::OffsetVector{Float64, Vector{Float64}}
    Sv::Function
    duration::Tuple{Float64, Float64}
    
    # Inner constructor for unit conversion - now only accepts keyword arguments
    function Pulse(; I₀, λ, cycles::Int64, envelope::Union{PulseEnvelope,Function}=Sin2, 
                   cep::Real=0.0, helicity::Int64=1, ϵ::Float64=0.0, Sv::Function=t->0.0, 
                   duration::Union{Tuple{Float64, Float64}, Nothing}=nothing)
        
        # Validate inputs
        if !(helicity == 1 || helicity == -1)
            throw(ArgumentError("helicity must be either +1 or -1"))
        end
        
        if cycles <= 0
            throw(ArgumentError("cycles must be positive"))
        end
        
        if !(0.0 <= ϵ <= 1.0)
            throw(ArgumentError("ellipticity ϵ must be between 0 and 1"))
        end
        
        # Convert inputs to atomic units
        # --- Conversion factors from CODATA2022 ---
        # 1 a.u. of intensity = 3.50944506e16 W/cm^2
        # 1 a.u. of length = BohrRadius (in meters) = 5.29177210544e-11 m
        # 1 a.u. of time = 2.4188843265857e-17 s
        
        # Speed of light in atomic units
        c_au = 1.0/0.0072973525643  # CODATA2022 value for speed of light in a.u.

        # Convert intensity to atomic units (from W/cm^2)
        I_au = I₀ / 3.50944506e16

        # Convert wavelength to atomic units (from nm)
        λ_au = (λ * 1e-9) / (5.29177210544e-11)

        # Calculate derived quantities using CODATA2022 constants
        ω = 2π * c_au / λ_au

        # Electric field amplitude in atomic units
        # E₀(a.u.) = sqrt(4 * π * I_au / c_au)
        # But in atomic units, E₀(a.u.) = sqrt(8π * I_au / c_au) is also used
        # For SFA, use E₀(a.u.) = sqrt(4 * I_au)
        E₀_au = sqrt(8π * I_au / c_au)

        # Vector potential amplitude A₀ = E₀/ω
        A₀ = E₀_au / ω  # sqrt(2* I_au) / ω  # E₀_au / ω

        # Pulse duration and ponderomotive energy
        Tp = 2π * cycles / ω
        Up = A₀^2 / 4

        # Initialize duration variable
        pulse_duration = duration

        # Define envelope functions with proper normalization
        f = if envelope isa PulseEnvelope
            if envelope == Sin2
                pulse_duration = pulse_duration === nothing ? (0.0, Tp) : pulse_duration
                t -> ( 0 ≤ t ≤ Tp) ? sin(ω * t / 2 / cycles)^2 : 0.0
            elseif envelope == Gaussian
                pulse_duration = pulse_duration === nothing ? (-3*Tp, 3*Tp) : pulse_duration
                σ = Tp / (2 * sqrt(2 * log(2)))
                t -> exp(-(t^2) / (2 * σ^2))
            elseif envelope == Rectangular
                pulse_duration = pulse_duration === nothing ? (0.0, Tp) : pulse_duration
                t -> ( 0 ≤ t ≤ Tp) ? 1.0 : 0.0
            else
                throw(ArgumentError("Unknown envelope type: $envelope. Use Sin2, Gaussian, Rectangular, or provide custom function"))
            end
        else
            # For custom envelope functions, use provided duration or default
            pulse_duration = pulse_duration === nothing ? (0.0, Tp) : pulse_duration
            envelope
        end

        # Calculate the polarization vector components in spherical basis
        # u-1, u0, u+1 components
        #------------------ Birger's part ----------------------#
        up = -1/sqrt(2*(1 + ϵ^2)) * (1 + helicity * ϵ)
        u0 = 0.0
        um = 1/sqrt(2*(1 + ϵ^2)) * (1 - helicity * ϵ)
        #------------------My corrctions  ----------------------#
        # up = -1/sqrt(2*(1 + ϵ^2)) * (1 - helicity * ϵ)
        # u0 = 0.0
        # um = 1/sqrt(2*(1 + ϵ^2)) * (1 + helicity * ϵ)
        #
        u = OffsetVector([um, u0, up], -1:1)

        new(I_au, A₀, λ_au, ω, cycles, Tp, Up, f, float(cep), helicity, ϵ, u, Sv, pulse_duration)
    end
end


"""
    ContinuumSolution

Enumeration of different approximations for the continuum electron wavefunction in strong-field ionization calculations.

# Variants
- `Bessel`
- `Distorted`

## `Bessel`
**Volkov states** - Free electron wavefunction in an oscillating electric field.

## `Distorted`
**Distorted waves** - Numerical solution including both Coulomb and laser field effects.

# Usage
```julia
# Create continuum electrons with different wavefunctions
bessel_electron = ContinuumElectron(1.0, sqrt(2.0), Bessel)       # SFA calculation
distorted_electron = ContinuumElectron(1.0, sqrt(2.0), Distorted) # Full calculation
```
"""
@enum ContinuumSolution begin
    Bessel
    Distorted
end

"""
    IonizationScheme

Enumeration of different approaches for modeling the initial bound state in strong-field ionization.

# Variants
- `Hydrogenic`: Hydrogen-like wavefunctions with effective nuclear charge
- `Atomic`: Full atomic structure calculation (Hartree-Fock or DFT)

# Usage
```julia
# Different ionization schemes
hydrogenic_calc = Settings(pulse, 1; ionization_scheme=Hydrogenic)
atomic_calc = Settings(pulse, 18; ionization_scheme=Atomic)
```
"""
@enum IonizationScheme begin
    Hydrogenic
    Atomic
end

"""
    Gauge

Enumeration of different gauge choices for the electromagnetic field in strong-field calculations.

# Variants
- `LengthGauge`: Length gauge (minimal coupling A·p)
- `VelocityGauge`: Velocity gauge (electric field E·r)

# Usage
```julia
# Different gauge choices
length_calc = Settings(pulse, 1; gauge=LengthGauge)
velocity_calc = Settings(pulse, 1; gauge=VelocityGauge)
```
"""
@enum Gauge begin
    LengthGauge
    VelocityGauge
end

"""
    ContinuumElectron

Defines a continuum (photo) electron for a particular energy.

# Fields
- `ε::Float64`
- `p::Float64`
- `solution::ContinuumSolution`
"""
struct ContinuumElectron
    ε::Float64
    p::Float64
    solution::ContinuumSolution
end


"""
    AtomicElectron

Represents a bound atomic electron with quantum numbers and radial wavefunction.

# Fields
- `Z::Int64`: Atomic number
- `n::Int64`: Principal quantum number
- `l::Int64`: Orbital angular momentum quantum number
- `j::Rational{Int64}`: Total angular momentum quantum number
- `ε::Float64`: Binding energy (positive value in a.u.)
- `r::Vector{Float64}`: Radial grid points (a.u.)
- `P::Vector{Float64}`: Radial wavefunction P(r) = r·R(r)
"""
struct AtomicElectron
    Z::Int64
    n::Int64
    l::Int64
    j::Rational{Int64}
    ε::Float64
    r::Vector{Float64}
    P::Vector{Float64}
end


"""
    PartialWave

Represents a partial wave component of the continuum electron wavefunction.

# Fields
- `ε::Float64`: Kinetic energy (a.u.)
- `l::Int64`: Orbital angular momentum quantum number
- `j::Rational{Int64}`: Total angular momentum quantum number
- `P::Vector{Float64}`: Radial wavefunction component
- `δ::Float64`: Phase shift (radians)
"""
struct PartialWave
    ε::Float64
    l::Int64
    j::Rational{Int64}
    P::Vector{Float64}
    δ::Float64
end


"""
    EnergyDistribution

Results from energy-resolved photoelectron distribution calculation at fixed angles.

# Fields
- `θ::Float64`: Polar angle (radians)
- `ϕ::Float64`: Azimuthal angle (radians)
- `energies::Vector{Float64}`: Energy grid (a.u.)
- `distribution::Vector{Float64}`: Differential ionization probability d²P/dΩdE
- `pulse::Pulse`: Laser pulse
"""
struct EnergyDistribution
    θ::Float64
    ϕ::Float64
    energies::Vector{Float64}
    distribution::Vector{Float64}
    pulse::Pulse
    
    function EnergyDistribution(θ::Float64, ϕ::Float64, energies::Vector{Float64}, distribution::Vector{Float64}, pulse::Pulse)
        length(energies) == length(distribution) || throw(ArgumentError("energies and distribution must have same length"))
        new(θ, ϕ, energies, distribution, pulse)
    end
end


"""
    AngularDistribution

Results from angle-resolved photoelectron distribution calculation at fixed energy.

# Fields
- `energy::Float64`: Fixed photoelectron energy (a.u.)
- `θ::Vector{Float64}`: Polar angle grid (radians)
- `ϕ::Vector{Float64}`: Azimuthal angle grid (radians)
- `distribution::Matrix{Float64}`: Angular distribution P(θ,ϕ) 
- `pulse::Pulse`: Laser pulse
"""
struct AngularDistribution
    energy::Float64
    θ::Vector{Float64}
    ϕ::Vector{Float64}
    distribution::Matrix{Float64}
    pulse::Pulse
    
    function AngularDistribution(energy::Float64, θ::Vector{Float64}, ϕ::Vector{Float64}, distribution::Matrix{Float64}, pulse::Pulse)
        size(distribution) == (length(θ), length(ϕ)) || 
            throw(ArgumentError("distribution size must match (length(θ), length(ϕ))"))
        new(energy, θ, ϕ, distribution, pulse)
    end
end


"""
    MomentumDistribution

Results from momentum-resolved photoelectron distribution calculation in spherical coordinates.

# Fields
- `p::Vector{Float64}`: Momentum magnitude grid (a.u.)
- `θ::Vector{Float64}`: Polar angle grid (radians, 0 to π)
- `φ::Vector{Float64}`: Azimuthal angle grid (radians, 0 to 2π)
- `distribution::Array{Float64,3}`: 3D momentum distribution P(p,θ,φ)
- `pulse::Pulse`: Laser pulse
"""
struct MomentumDistribution
    p::Vector{Float64}
    θ::Vector{Float64}
    φ::Vector{Float64}
    distribution::Array{Float64,3}
    pulse::Pulse
    
    function MomentumDistribution(p::Vector{Float64}, θ::Vector{Float64}, φ::Vector{Float64}, distribution::Array{Float64,3}, pulse::Pulse)
        size(distribution) == (length(p), length(θ), length(φ)) || 
            throw(ArgumentError("distribution size must match momentum grid dimensions"))
        new(p, θ, φ, distribution, pulse)
    end
end


"""
    compute_energy_distribution(Z::Int64, pulse::Pulse; settings=Settings(), θ::Real=pi/2, ϕ::Real=0.0,
                               energy_range::Tuple{Float64,Float64}=(0.0, 0.5), n_points::Int64=200) -> EnergyDistribution

Computes the energy-resolved photoelectron spectrum at specified angles using SFA.

# Arguments
- `Z::Int64`: Atomic number of the target atom
- `pulse::Pulse`: Laser pulse parameters
- `settings=Settings()`: Calculation settings (ionization scheme, continuum solution, gauge)
- `θ::Real=π/2`: Polar detection angle (radians)
- `ϕ::Real=0.0`: Azimuthal detection angle (radians)
- `energy_range::Tuple{Float64,Float64}=(0.0, 0.5)`: Energy range (a.u.)
- `n_points::Int64=200`: Number of energy points

# Returns
- `EnergyDistribution`: Structure containing energies, spectrum, and pulse parameters

# Example
```julia
pulse = Pulse(I₀=1e14, λ=800, cycles=8)
ed = compute_energy_distribution(18, pulse; energy_range=(0.0, 2.0), n_points=500)
```
"""
function compute_energy_distribution(Z::Int64, pulse::Pulse;
                                    settings=Settings(), 
                                    θ::Real=pi/2, ϕ::Real=0.0,
                                    energy_range::Tuple{Float64,Float64}=(0.0, 0.5),
                                    n_points::Int64=200)

    energies = range(energy_range[1], energy_range[2], length=n_points) |> collect
    distribution = zeros(Float64, n_points)
    
    a_electron = StrongFieldDynamics.compute_atomic_electron(Z, settings.ionization_scheme) ;

    p_electrons = [StrongFieldDynamics.ContinuumElectron(ep, sqrt(2*ep), settings.continuum_solution) for ep in energies]
    
    if Distributed.nworkers() > 1
        # Create a function for computing a single energy point
        compute_single_energy = function(i)
            return StrongFieldDynamics.probability(pulse, a_electron, p_electrons[i], float(θ), float(ϕ))
        end
        
        # Use pmap to distribute computation across workers
        distribution = @showprogress "Computing energy distribution..." pmap(compute_single_energy, 1:n_points)
    else
        # Use multithreading for single-node computation
        @showprogress Threads.@threads for i in eachindex(p_electrons)
            distribution[i] = StrongFieldDynamics.probability(pulse, a_electron, p_electrons[i], float(θ), float(ϕ))
        end
    end
    
    return EnergyDistribution(θ, ϕ, energies, distribution, pulse)
end


"""
    compute_angular_distribution(Z::Int64, pulse::Pulse; settings=Settings(), energy::Float64=1.0,
                                θ_range::Tuple{Float64,Float64}=(0.0, π), n_θ::Int=50,
                                ϕ_range::Tuple{Float64,Float64}=(0.0, 2π), n_ϕ::Int=200) -> AngularDistribution

Computes the angular distribution of photoelectrons at fixed energy using SFA.

# Arguments
- `Z::Int64`: Atomic number of the target atom
- `pulse::Pulse`: Laser pulse parameters
- `settings=Settings()`: Calculation settings (ionization scheme, continuum solution, gauge)
- `energy::Float64=1.0`: Fixed photoelectron energy (a.u.)
- `θ_range::Tuple{Float64,Float64}=(π/2, π/2)`: Polar angle range (radians)
- `n_θ::Int=1`: Number of polar angle points
- `ϕ_range::Tuple{Float64,Float64}=(0.0, 2π)`: Azimuthal angle range (radians)
- `n_ϕ::Int=200`: Number of azimuthal angle points

# Returns
- `AngularDistribution`: Structure containing angles, distribution, and pulse parameters

# Example
```julia
pulse = Pulse(I₀=5e13, λ=800, cycles=6, ϵ=0.5)
ad = compute_angular_distribution(36, pulse; energy=2.0, θ_range=(0.0, π), n_θ=100, n_ϕ=360)
```
"""
function compute_angular_distribution(Z::Int64, pulse::Pulse;
                                     settings=Settings(),
                                     energy::Float64=1.0,
                                     θ_range::Tuple{Float64,Float64}=(π/2, π/2),
                                     n_θ::Int=1,
                                     ϕ_range::Tuple{Float64,Float64}=(0.0, 2π),
                                     n_ϕ::Int=200)
    
    θs = range(θ_range[1], θ_range[2], length=n_θ) |> collect
    ϕs = range(ϕ_range[1], ϕ_range[2], length=n_ϕ) |> collect
    distribution = zeros(Float64, n_θ, n_ϕ)

    a_electron = StrongFieldDynamics.compute_atomic_electron(Z, settings.ionization_scheme) ;
    
    p_electron = StrongFieldDynamics.ContinuumElectron(energy, sqrt(2*energy), settings.continuum_solution)

    if Distributed.nworkers() > 1
        # Create a function for computing a single angular point
        compute_single_angle = function(i)
            θ_idx = ((i-1) ÷ n_ϕ) + 1
            ϕ_idx = ((i-1) % n_ϕ) + 1
            return StrongFieldDynamics.probability(pulse, a_electron, p_electron, float(θs[θ_idx]), float(ϕs[ϕ_idx]))
        end
        
        # Use pmap to distribute computation across workers
        results = @showprogress "Computing angular distribution..." pmap(compute_single_angle, 1:(n_θ*n_ϕ))
        
        # Reshape results into distribution array
        distribution = reshape(results, n_θ, n_ϕ)
    else
        # Use multithreading for single-node computation
        @showprogress Threads.@threads for i in 1:(n_θ*n_ϕ)
            θ_idx = ((i-1) ÷ n_ϕ) + 1
            ϕ_idx = ((i-1) % n_ϕ) + 1
            distribution[θ_idx, ϕ_idx] = StrongFieldDynamics.probability(pulse, a_electron, p_electron, float(θs[θ_idx]), float(ϕs[ϕ_idx]))
        end
    end

    return AngularDistribution(energy, θs, ϕs, distribution, pulse)
end


"""
    compute_momentum_distribution(Z::Int64, pulse::Pulse; settings=Settings(), energy_range::Tuple{Float64,Float64}=(0.0, 0.5),
                                 n_p::Int=50, n_theta::Int=1, n_phi::Int=50, coupled::Bool=true) -> MomentumDistribution

Computes the 3D momentum distribution of photoelectrons using SFA in spherical coordinates.

# Arguments
- `Z::Int64`: Atomic number of the target atom
- `pulse::Pulse`: Laser pulse parameters
- `settings=Settings()`: Calculation settings (ionization scheme, continuum solution, gauge)
- `energy_range::Tuple{Float64,Float64}=(0.0, 0.5)`: Energy range (a.u.)
- `n_p::Int=50`: Number of momentum magnitude points
- `n_theta::Int=1`: Number of polar angle points (typically 1 for θ=π/2)
- `n_phi::Int=50`: Number of azimuthal angle points

# Returns
- `MomentumDistribution`: Structure containing momentum grids, 3D distribution, and pulse parameters

# Example
```julia
pulse = Pulse(I₀=2e14, λ=800, cycles=4, helicity=-1)
md = compute_momentum_distribution(54, pulse; energy_range=(0.0, 5.0), n_p=100, n_phi=200)
```

# Notes
- For typical SFA calculations, `n_theta=1` with θ=π/2 is sufficient
- Higher `n_phi` values provide better angular resolution
- `coupled=true` includes interchannel coupling effects
"""
function compute_momentum_distribution(Z::Int64, pulse::Pulse;
                                      settings=Settings(), 
                                      energy_range::Tuple{Float64,Float64}=(0.0, 0.5), n_p::Int=50, 
                                      n_theta::Int=1, n_phi::Int=50, coupled::Bool=true)
    energies = range(energy_range[1], energy_range[2], length=n_p) |> collect
    θs = [float(pi/2)] |> collect
    # θs = range(0.0, π, length=n_theta) |> collect
    φs = range(0.0, 2π, length=n_phi) |> collect
    distribution = zeros(Float64, n_p, n_theta, n_phi)

    local_probability_func = coupled ? StrongFieldDynamics.probability : StrongFieldDynamics.probability_uncoupled
    # local_probability_func = StrongFieldDynamics.probability

    a_electron = StrongFieldDynamics.compute_atomic_electron(Z, settings.ionization_scheme) ;

    p_electrons = [ StrongFieldDynamics.ContinuumElectron(energy, sqrt(2*energy), settings.continuum_solution) for energy in energies ]

    if Distributed.nworkers() > 1
        # Create a function for computing a single momentum point
        compute_single_momentum = function(i)
            p_electron = p_electrons[i]
            row = zeros(Float64, n_theta, n_phi)
            for j in eachindex(θs)
                for k in eachindex(φs)
                    row[j, k] = local_probability_func(pulse, a_electron, p_electron, θs[j], φs[k])
                end
            end
            return row
        end
        
        # Use pmap to distribute computation
        results = @showprogress "Computing momentum distribution..." pmap(compute_single_momentum, 1:n_p)
        
        # Populate the distribution array from results
        for i in 1:n_p
            distribution[i, :, :] = results[i]
        end
    else    # Multithreading computation
        Threads.@threads :dynamic for i in eachindex(p_electrons)
            for j in eachindex(θs)
                for k in eachindex(φs)
                    distribution[i, j, k] = local_probability_func(pulse, a_electron, p_electrons[i], θs[j], φs[k])
                end
            end
        end
    end
    
    return MomentumDistribution(sqrt.(2*energies), θs, φs, distribution, pulse)
end


"""
    Settings

Configuration settings for strong-field ionization calculations.

# Fields
- `ionization_scheme::IonizationScheme`: Method for bound state calculation
- `continuum_solution::ContinuumSolution`: Continuum electron wavefunction approximation
- `gauge::Gauge`: Electromagnetic gauge choice

# Constructors
```julia
# Basic constructor with defaults
settings = Settings(pulse, 1)

# With specific options
settings = Settings(pulse, 18; 
                   ionization_scheme=Atomic,
                   continuum_solution=Distorted,
                   gauge=VelocityGauge)
```

# Validation
- Atomic number must be positive
"""
struct Settings
    ionization_scheme::IonizationScheme
    continuum_solution::ContinuumSolution
    gauge::Gauge

    function Settings(;ionization_scheme::IonizationScheme=Hydrogenic,
                     continuum_solution::ContinuumSolution=Bessel,
                     gauge::Gauge=VelocityGauge)
        
        # Warning and Error
        if gauge == LengthGauge
            throw(ArgumentError("LengthGauge computations are not implemented yet!!!"))
        end

        new(ionization_scheme, continuum_solution, gauge)
    end
end


"""
    Cartesian

A simple data structure to represent 3D Cartesian coordinates.

# Fields
- `x::Number`: x-coordinate
- `y::Number`: y-coordinate  
- `z::Number`: z-coordinate

# Example
```julia
# Create Cartesian coordinates
point = Cartesian(1.0, 2.0, 3.0)
println("Position: (", point.x, ", ", point.y, ", ", point.z, ")")

# Convert from spherical coordinates
r, θ, ϕ = 2.0, π/4, π/6
cartesian_point = spherical2cartesian(r, θ, ϕ)
```
"""
struct Cartesian
    x::Number
    y::Number
    z::Number
end

"""
    spherical2cartesian(r::Float64, θ::Float64, ϕ::Float64) -> Cartesian

Convert spherical coordinates to Cartesian coordinates.

# Arguments
- `r::Float64`: Radial distance (magnitude)
- `θ::Float64`: Polar angle (angle from z-axis, 0 to π radians)
- `ϕ::Float64`: Azimuthal angle (angle in xy-plane from x-axis, 0 to 2π radians)

# Returns
- `Cartesian`: A Cartesian coordinate struct with x, y, z components

# Mathematical Relations
The conversion follows standard spherical coordinate conventions:
- `x = r sin(θ) cos(ϕ)`
- `y = r sin(θ) sin(ϕ)`
- `z = r cos(θ)`

# Examples
```julia
# Convert point at unit sphere
r, θ, ϕ = 1.0, π/2, 0.0  # Point on positive x-axis
cart = spherical2cartesian(r, θ, ϕ)
# Returns: Cartesian(1.0, 0.0, 0.0)

# Convert general point
r, θ, ϕ = 2.0, π/4, π/3
cart = spherical2cartesian(r, θ, ϕ)
# Returns: Cartesian(0.707..., 1.224..., 1.414...)
"""
function spherical2cartesian(r::Float64, θ::Float64, ϕ::Float64)
    x = r * sin(θ) * cos(ϕ)
    y = r * sin(θ) * sin(ϕ)
    z = r * cos(θ)

    return Cartesian(x, y, z)
end

"""
    ClebschGordan(ja, ma, jb, mb, Jab, Mab)

Computes the Clebsch-Gordan coefficients
"""
function ClebschGordan(ja, ma, jb, mb, Jab, Mab)
    try  
        WignerSymbols.clebschgordan(Float64, ja, ma, jb, mb, Jab, Mab)
    catch er
        return 0.
    end
end


"""
    Ylm(l, m, θ, ϕ)

Computes the Spherical hamonics using the Associated Legender Polynomials as 
    Ylm(l, m, θ, ϕ) = Plm(l, m, cos(θ)) * exp(im*m*ϕ)
"""
function Ylm(l, m, θ, ϕ)
    m = convert(Int64, m)
    if l < 0 || abs(m) > l 
        return zero(ComplexF64) 
    else
        return SphericalHarmonics.sphericalharmonic(θ, ϕ, l, m)
    end

end
# function Ylm(l, m, θ, ϕ)
#     m = convert(Int64, m)
#     if abs(m) > l return zero(ComplexF64) end
#     if m ≥ 0
#         P = AssociatedLegendrePolynomials.legendre(LegendreSphereNorm(), l, m, cos(θ))
#         return P * exp(im * m * ϕ)
#     else
#         P = AssociatedLegendrePolynomials.legendre(LegendreSphereNorm(), l, -m, cos(θ))
#         return (-1)^m * conj(P * exp(im * (-m) * ϕ))
#     end
# end


# # Add display method for Pulse
function Base.show(io::IO, p::Pulse)
    # Convert atomic units back to conventional units for display
    λ_nm = p.λ * (5.29177210544e-11) * 1e9
    I_Wcm2 = p.I₀ * 3.51e16
    
    # Calculate additional derived parameters
    ν_Hz = 2.998e8 / (λ_nm * 1e-9)
    photon_eV = 1240.0 / λ_nm
    E₀_Vm = sqrt(2 * I_Wcm2 * 1e4 / (2.998e8 * 8.854e-12))
    E₀_au = E₀_Vm / 5.14e11
    
    # Calculate critical parameters for strong-field physics
    γ_Keldysh = sqrt(0.5 * 13.6) / (p.Up * 27.2114)
    
    # Determine polarization type
    pol_type = if p.ϵ ≈ 0
        "\033[36mLinear\033[0m"
    elseif abs(p.ϵ) ≈ 1
        p.helicity == 1 ? "\033[35mRight-circular\033[0m" : "\033[35mLeft-circular\033[0m"
    else
        handedness = p.helicity == 1 ? "Right" : "Left"
        "\033[33mElliptical $(handedness) (ϵ=$(round(p.ϵ, digits=2)))\033[0m"
    end
    
    # Strip ANSI codes for length calculation - define once at the top
    strip_ansi(s) = replace(s, r"\033\[[0-9;]*m" => "")
    
    # Format parameters with consistent width
    col_width = 50  # Width for each column
    
    # Left column: USER-DEFINED PARAMETERS (input parameters)
    cep_fraction = p.cep / π
    cep_display = if p.cep ≈ 0
        @sprintf("\033[93m%.2f\033[0m rad", p.cep)
    else
        @sprintf("\033[93m%.2f\033[0m rad (\033[93m%.2fπ\033[0m)", p.cep, cep_fraction)
    end
    
    helicity_str = p.helicity == 1 ? "\033[32m+1\033[0m" : "\033[31m-1\033[0m"
    
    left_params = [
        ("Peak Intensity ", @sprintf("\033[91m%.2e\033[0m W/cm² (\033[94m%.2e\033[0m a.u.)", I_Wcm2, p.I₀)),
        ("Wavelength (λ)", @sprintf("\033[91m%.1f\033[0m nm (\033[94m%.3f\033[0m a.u.)", λ_nm, p.λ)),
        ("Optical Cycles", @sprintf("\033[92m%d\033[0m", p.np)),
        ("CEP Phase", cep_display),
        ("Helicity (λ)", helicity_str),
        ("Ellipticity (ϵ)", @sprintf("\033[93m%.2f\033[0m", p.ϵ)),
        ("Polarization Type", pol_type)
    ]
    
    # Right column: DERIVED PARAMETERS (calculated from inputs)
    duration_fs = p.Tp * 0.024189
    period_fs = 2π / p.ω * 0.024189
    up_eV = p.Up * 27.2114
    excursion_nm = p.A₀/p.ω * 0.0529
    
    right_params = [
        ("Peak E-field (E₀)", @sprintf("\033[91m%.2e\033[0m V/m (\033[94m%.3f\033[0m a.u.)", E₀_Vm, E₀_au)),
        ("Vector Potential (A₀)", @sprintf("\033[94m%.4f\033[0m a.u.", p.A₀)),
        ("Frequency (ω)", @sprintf("\033[91m%.2e\033[0m Hz (\033[94m%.3f\033[0m a.u.)", ν_Hz, p.ω)),
        ("Photon Energy (ℏω)", @sprintf("\033[91m%.3f\033[0m eV", photon_eV)),
        ("Pulse Duration (Tp)", @sprintf("\033[91m%.1f\033[0m fs (\033[94m%.2f\033[0m a.u.)", duration_fs, p.Tp)),
        ("Optical Period (T₀)", @sprintf("\033[91m%.2f\033[0m fs (\033[94m%.2f\033[0m a.u.)", period_fs, 2π/p.ω)),
        ("Ponderomotive Energy (Up)", @sprintf("\033[91m%.2f\033[0m eV (\033[94m%.4f\033[0m a.u.)", up_eV, p.Up))
    ]
    
    # Additional derived parameters
    additional_params = [
        ("Keldysh Parameter (γ)", @sprintf("\033[93m%.3f\033[0m", γ_Keldysh)),
        ("Excursion Amplitude (α₀)", @sprintf("\033[91m%.2f\033[0m nm (\033[94m%.2f\033[0m a.u.)", excursion_nm, p.A₀/p.ω)),
        ("Spherical Components (u₋₁,₀,₊₁)", @sprintf("(\033[94m%.3f\033[0m, \033[94m%.3f\033[0m, \033[94m%.3f\033[0m)", p.u[-1], p.u[0], p.u[1]))
    ]
    
    # Helper function for formatting two-column rows
    format_row(left_label, left_val, right_label, right_val) = begin
        left_content = "  \033[1m$(left_label):\033[0m $(left_val)"
        right_content = "\033[1m$(right_label):\033[0m $(right_val)"
        
        left_clean = strip_ansi(left_content)
        right_clean = strip_ansi(right_content)
        
        left_pad = max(0, col_width - length(left_clean))
        right_pad = max(0, col_width - length(right_clean) - 2)
        
        left_spaces = " "^left_pad
        right_spaces = " "^right_pad
        
        "\033[96m│\033[0m$(left_content)$(left_spaces)\033[96m│\033[0m $(right_content)$(right_spaces) \033[96m│\033[0m"
    end
    
    # Calculate box width
    box_width = col_width * 2 + 3  # Two columns plus borders
    
    # Print header with column titles
    println(io, "\033[96m╭" * "─"^(col_width) * "┬" * "─"^(col_width) * "╮\033[0m")
    
    # Main title
    title = "\033[1m\033[96mLASER PULSE PARAMETERS\033[0m"
    title_clean = "LASER PULSE PARAMETERS"
    title_pad = (box_width - length(title_clean) - 2) ÷ 2
    title_extra = (box_width - length(title_clean) - 2) % 2
    println(io, "\033[96m│\033[0m" * " "^title_pad * title * " "^(title_pad + title_extra) * "\033[96m│\033[0m")
    
    # Column headers
    left_header = "\033[1m\033[93mUSER-DEFINED PARAMETERS\033[0m"
    right_header = "\033[1m\033[93mDERIVED PARAMETERS\033[0m"
    left_header_clean = "USER-DEFINED PARAMETERS"
    right_header_clean = "DERIVED PARAMETERS"
    
    left_header_pad = (col_width - length(left_header_clean) - 2) ÷ 2
    right_header_pad = (col_width - length(right_header_clean) - 2) ÷ 2
    
    println(io, "\033[96m├" * "─"^(col_width) * "┼" * "─"^(col_width) * "┤\033[0m")
    println(io, "\033[96m│\033[0m " * " "^left_header_pad * left_header * " "^(col_width - length(left_header_clean) - left_header_pad - 1) * 
                "\033[96m│\033[0m " * " "^right_header_pad * right_header * " "^(col_width - length(right_header_clean) - right_header_pad - 1) * "\033[96m│\033[0m")
    println(io, "\033[96m├" * "─"^(col_width) * "┼" * "─"^(col_width) * "┤\033[0m")
    
    # Print parameter rows
    max_rows = max(length(left_params), length(right_params))
    for i in 1:max_rows
        left_label, left_val = i <= length(left_params) ? left_params[i] : ("", "")
        right_label, right_val = i <= length(right_params) ? right_params[i] : ("", "")
        
        if left_label == "" && right_label == ""
            continue
        elseif left_label == ""
            # Only right column
            right_content = "\033[1m$(right_label):\033[0m $(right_val)"
            right_clean = strip_ansi(right_content)
            right_pad = max(0, col_width - length(right_clean) - 2)
            right_spaces = " "^right_pad
            println(io, "\033[96m│\033[0m" * " "^col_width * "\033[96m│\033[0m $(right_content)$(right_spaces) \033[96m│\033[0m")
        elseif right_label == ""
            # Only left column
            left_content = "  \033[1m$(left_label):\033[0m $(left_val)"
            left_clean = strip_ansi(left_content)
            left_pad = max(0, col_width - length(left_clean))
            left_spaces = " "^left_pad
            right_empty_spaces = " "^(col_width - 1)
            println(io, "\033[96m│\033[0m$(left_content)$(left_spaces)\033[96m│\033[0m$(right_empty_spaces) \033[96m│\033[0m")
        else
            println(io, format_row(left_label, left_val, right_label, right_val))
        end
    end
    
    # Add additional derived parameters section
    println(io, "\033[96m├" * "─"^(col_width) * "┴" * "─"^(col_width) * "┤\033[0m")
    println(io, "\033[96m│  \033[1m\033[93mADDITIONAL STRONG-FIELD PARAMETERS \033[0m" * " "^(box_width - 39) * "\033[96m│\033[0m")
    println(io, "\033[96m├" * "─"^(box_width-2) * "┤\033[0m")
    
    for (param_label, param_val) in additional_params
        content = "  \033[1m$(param_label):\033[0m $(param_val)"
        content_clean = strip_ansi(content)
        padding = box_width - length(content_clean) - 3
        padding_spaces = " "^padding
        println(io, "\033[96m│\033[0m$(content)$(padding_spaces) \033[96m│\033[0m")
    end
    
    println(io, "\033[96m╰" * "─"^(box_width-2) * "╯\033[0m")
end


# # Add display method for Settings
# function Base.show(io::IO, settings::Settings)
    
#     scheme_name = if settings.ionization_scheme == Hydrogenic
#         "Hydrogenic"
#     elseif settings.ionization_scheme == Atomic
#         "Atomic"
#     else
#         "Unknown Scheme"
#     end

#     gauge_name = if settings.gauge == LengthGauge
#         "Length Gauge"
#     elseif settings.gauge == VelocityGauge
#         "Velocity Gauge"
#     else
#         "Unknown Gauge"
#     end

#     # Format parameters
#     atom_name = if settings.atom <= 118
#         elements = ["He", "Li", "Ne", "Ar", "Kr", "Xe"]
#         if settings.atom <= length(elements)
#             elements[settings.atom]
#         else
#             "Element"
#         end
#     else
#         "Unknown"
#     end

#     println(io, "\033[96m╭────────────────────────────────────────────────────────────╮\033[0m")
#     println(io, "\033[96m│\033[0m                   \033[1m\033[93mCALCULATION SETTINGS\033[0m                   \033[96m│\033[0m")
#     println(io, "\033[96m├────────────────────────────────────────────────────────────┤\033[0m")
    
#     params = [
#         ("Atom", "\033[91m$(atom_name)\033[0m (Z = \033[91m$(settings.atom)\033[0m)"),
#         ("Ionization Scheme", "\033[92m$(scheme_name)\033[0m"),
#         ("Continuum Solution", "\033[92m$(settings.continuum_solution)\033[0m"),
#         ("Gauge Choice", "\033[92m$(gauge_name)\033[0m")
#     ]
    
#     for (label, value) in params
#         content = "  \033[1m$(label):\033[0m $(value)"
#         content_clean = replace(content, r"\033\[[0-9;]*m" => "")
#         padding = 62 - length(content_clean)
#         padding_spaces = " "^max(0, padding)
#         println(io, "\033[96m│\033[0m$(content)$(padding_spaces) \033[96m│\033[0m")
#     end
    
#     println(io, "\033[96m╰────────────────────────────────────────────────────────────╯\033[0m")
# end

using AssociatedLegendrePolynomials
using Distributed
using Printf
using ProgressMeter
using SpecialFunctions
using SphericalHarmonics
using WignerSymbols

export ClebschGordan, Ylm
export Pulse, AtomicElectron, ContinuumElectron, ContinuumSolution, PartialWave
export AngularDistribution, EnergyDistribution, MomentumDistribution
export compute_angular_distribution, compute_energy_distribution, compute_momentum_distribution


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
- `u::QuasiVector{Float64, Tuple{UnitRange{Int64}}}`: Polarization unit vector
- `Sv::Function`: Vector potential function A(t)

# Required Arguments
- `I₀`: Peak intensity (can include units, e.g., 1e14u"W/cm^2", or numeric value assuming W/cm²)
- `λ`: Wavelength (can include units, e.g., 800u"nm", or numeric value assuming nm)
- `cycles`: Number of optical cycles (positive integer)
- `envelope`: Envelope type (:sin2, :gauss, :flat) or custom function f(t)
- `cep`: Carrier-envelope phase in radians (default: 0.0)
- `helicity`: Helicity (+1 or -1, default: +1)
- `ϵ`: Ellipticity parameter (0=linear, 1=circular, default: 0.0)

# Optional Arguments
- `Sv`: Custom vector potential function A(t) (default: t->0.0)
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
    u::QuasiVector{Float64, Tuple{UnitRange{Int64}}}
    Sv::Function
    
    # Inner constructor for unit conversion - now only accepts keyword arguments
    function Pulse(; I₀, λ, cycles::Int64, envelope::Union{Symbol,Function}=:sin2, 
                   cep::Real=0.0, helicity::Int64=1, ϵ::Float64=0.0, Sv::Function=t->0.0)
        
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
        A₀ = E₀_au / ω

        # Pulse duration and ponderomotive energy
        Tp = 2π * cycles / ω
        Up = A₀^2 / 4

        # Define envelope functions with proper normalization
        f = if envelope isa Symbol
            if envelope == :sin2
                t -> ( 0 ≤ t ≤ Tp) ? sin(ω * t / 2 / cycles)^2 : 0.0
            elseif envelope == :gauss
                σ = Tp / (2 * sqrt(2 * log(2)))
                t -> exp(-(t^2) / (2 * σ^2))
            elseif envelope == :flat
                t -> (abs(t) <= Tp/2) ? 1.0 : 0.0
            else
                throw(ArgumentError("Unknown envelope type: $envelope. Use :sin2, :gauss, :flat, or provide custom function"))
            end
        else
            envelope
        end

        # Calculate the polarization vector components in spherical basis
        # u-1, u0, u+1 components
        up = -1/sqrt(2*(1 + ϵ^2)) * (1 + helicity * ϵ)
        u0 = 0.0
        um = 1/sqrt(2*(1 + ϵ^2)) * (1 - helicity * ϵ)
        u = QuasiVector([um, u0, up], -1:1)

        new(I_au, A₀, λ_au, ω, cycles, Tp, Up, f, float(cep), helicity, ϵ, u, Sv)
    end
end


"""
    ContinuumSolution

Bessel      -> ...
Coulomb     -> ...
Distorted   -> ...
"""
@enum ContinuumSolution begin
    Bessel
    Coulomb
    Distorted
end


"""
    ContinuumElectron

Defines a continuum (photo) electron for a particular energy.

`ε::Float64`
`p::Float64`
`solution::Symbol`
"""
struct ContinuumElectron
    ε::Float64
    p::Float64
    solution::Symbol
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

Results from energy-resolved photoelectron spectrum calculation at fixed angles.

# Fields
- `θ::Float64`: Polar angle (radians)
- `ϕ::Float64`: Azimuthal angle (radians)
- `energies::Vector{Float64}`: Energy grid (a.u.)
- `spectrum::Vector{Float64}`: Differential ionization probability d²P/dΩdE
- `pulse::Pulse`: Laser pulse
"""
struct EnergyDistribution
    θ::Float64
    ϕ::Float64
    energies::Vector{Float64}
    spectrum::Vector{Float64}
    pulse::Pulse
    
    function EnergyDistribution(θ::Float64, ϕ::Float64, energies::Vector{Float64}, spectrum::Vector{Float64}, pulse::Pulse)
        length(energies) == length(spectrum) || throw(ArgumentError("energies and spectrum must have same length"))
        new(θ, ϕ, energies, spectrum, pulse)
    end
end


"""
    AngularDistribution

Results from angle-resolved photoelectron distribution at fixed energy and theta.

# Fields
- `energy::Float64`: Fixed photoelectron energy (a.u.)
- `θ::Float64`: Fixed polar angle (radians)
- `ϕ::Vector{Float64}`: Azimuthal angle grid (radians)
- `distribution::Vector{Float64}`: Angular distribution P(ϕ) at fixed θ
- `pulse::Pulse`: Laser pulse
"""
struct AngularDistribution
    energy::Float64
    θ::Float64
    ϕ::Vector{Float64}
    distribution::Vector{Float64}
    pulse::Pulse
    
    function AngularDistribution(energy::Float64, θ::Float64, ϕ::Vector{Float64}, distribution::Vector{Float64}, pulse::Pulse)
        length(distribution) == length(ϕ) || 
            throw(ArgumentError("distribution length must match length(ϕ)"))
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
    compute_energy_distribution(a_electron::AtomicElectron, pulse::Pulse; θ::Real=pi/2, ϕ::Real=0.0,
                               energy_range::Tuple{Float64,Float64}=(1e-6, 0.5), n_points::Int64=200, coupled::Bool=true) -> EnergyDistribution

Computes the energy-resolved photoelectron spectrum at specified angles using SFA.

# Arguments
- `a_electron::AtomicElectron`: Initial bound state
- `pulse::Pulse`: Laser pulse parameters
- `θ::Real=π/2`: Polar detection angle (radians)
- `ϕ::Real=0.0`: Azimuthal detection angle (radians)
- `energy_range::Tuple{Float64,Float64}=(0.0, 0.5)`: Energy range (a.u.)
- `n_points::Int64=200`: Number of energy points
"""
function compute_energy_distribution(a_electron::AtomicElectron, pulse::Pulse; 
                                    θ::Real=pi/2, ϕ::Real=0.0,
                                    energy_range::Tuple{Float64,Float64}=(0.0, 0.5),
                                    n_points::Int64=200,
                                    coupled::Bool=true)

    energies = range(energy_range[1], energy_range[2], length=n_points) |> collect
    distribution = zeros(Float64, n_points)
    
    p_electrons = [StrongFieldDynamics.ContinuumElectron(ep, sqrt(2*ep), :Bessel) for ep in energies]

    if coupled
        @showprogress Threads.@threads for i in eachindex(p_electrons)
            distribution[i] = StrongFieldDynamics.probability(pulse, a_electron, p_electrons[i], float(θ), float(ϕ))
        end
    else
        @showprogress Threads.@threads for i in eachindex(p_electrons)
            distribution[i] = StrongFieldDynamics.probability_uncoupled(pulse, a_electron, p_electrons[i], float(θ), float(ϕ))
        end
    end
    
    return EnergyDistribution(θ, ϕ, energies, distribution, pulse)
end


"""
    compute_angular_distribution(electron::ContinuumElectron, bound::AtomicElectron,
                                pulse::Pulse; energy::Float64=1.0,
                                θ::Float64=0.0,
                                ϕ_range::Tuple{Float64,Float64}=(0.0, 2π),
                                n_ϕ::Int=100,
                                coupled::Bool=true) -> AngularDistribution

Computes the angular distribution of photoelectrons at fixed energy and theta using SFA.

# Arguments
- `a_electron::AtomicElectron`: Initial bound state
- `pulse::Pulse`: Laser pulse parameters
- `energy::Float64=1.0`: Fixed photoelectron energy (a.u.)
- `θ::Real=pi/2`: Fixed polar angle (radians)
- `ϕ_range::Tuple{Float64,Float64}=(0.0, 2π)`: Azimuthal angle range (radians)
- `n_ϕ::Int=200`: Number of azimuthal angle points
"""
function compute_angular_distribution(a_electron::AtomicElectron,
                                     pulse::Pulse; energy::Float64=1.0,
                                     θ::Real=pi/2,
                                     ϕ_range::Tuple{Float64,Float64}=(0.0, 2π),
                                     n_ϕ::Int=200,
                                     coupled::Bool=true)
    
    ϕs = range(ϕ_range[1], ϕ_range[2], length=n_ϕ) |> collect
    distribution = zeros(Float64, n_ϕ)
    
    p_electron = StrongFieldDynamics.ContinuumElectron(energy, sqrt(2*energy), :Bessel)

    if coupled
        @showprogress Threads.@threads for i in eachindex(ϕs)
            distribution[i] = StrongFieldDynamics.probability(pulse, a_electron, p_electron, float(θ), float(ϕs[i]))
        end
    else
        @showprogress Threads.@threads for i in eachindex(ϕs)
            distribution[i] = StrongFieldDynamics.probability_uncoupled(pulse, a_electron, p_electron, float(θ), float(ϕs[i]))
        end
    end
    return AngularDistribution(energy, θ, ϕs, distribution, pulse)
end


"""
    compute_momentum_distribution(electron::ContinuumElectron, bound::AtomicElectron,
                                 pulse::Pulse; p_max::Float64=2.0,
                                 n_p::Int=50, n_θ::Int=25, n_φ::Int=50) -> MomentumDistribution

Computes the 3D momentum distribution of photoelectrons using SFA in spherical coordinates.

# Arguments
- `a_electron::AtomicElectron`: Initial bound state
- `pulse::Pulse`: Laser pulse parameters
- `energy_range::Tuple{Float64,Float64}=(0.0, 0.5)`: Energy range (a.u.)
- `n_p::Int=50`: Number of momentum magnitude points
- `n_θ::Int=1`: Number of polar angle points
- `n_φ::Int=50`: Number of azimuthal angle points
"""
function compute_momentum_distribution(a_electron::AtomicElectron, pulse::Pulse; 
                                      energy_range::Tuple{Float64,Float64}=(0.0, 0.5), n_p::Int=50, 
                                      n_theta::Int=1, n_phi::Int=50, coupled::Bool=true)
    energies = range(energy_range[1], energy_range[2], length=n_p) |> collect
    # θ = range(0.0, π, length=n_theta) |> collect
    θs = [float(pi/2)] |> collect
    φs = range(0.0, 2π, length=n_phi) |> collect
    distribution = zeros(Float64, n_p, n_theta, n_phi)

    local_probability_func = coupled ? StrongFieldDynamics.probability : StrongFieldDynamics.probability_uncoupled

    p_electrons = [ StrongFieldDynamics.ContinuumElectron(energy, sqrt(2*energy), :Bessel) for energy in energies ]

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
        results = pmap(compute_single_momentum, 1:n_p)
        
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

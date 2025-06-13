using AssociatedLegendrePolynomials
using SpecialFunctions
using WignerSymbols

export ClebschGordan, Ylm
export Pulse, AtomicElectron, ContinuumElectron, ContinuumSolution, PartialWave
export AngularDistribution, EnergyDistribution, MomentumDistribution
export compute_angular_distribution, compute_energy_distribution, compute_momentum_distribution

struct Pulse
    I       ::Float64
    A₀      ::Float64               # 
    # A       ::Function              # A(t) = A₀ f(t) exp(-iωt)
    λ       ::Float64
    ω       ::Float64
    np      ::Int64                 # number of cycles
    Tp      ::Float64
    Up      ::Float64
    f       ::Function              # Carrier envelope 
    ϕ       ::Float64               # CEP phase
    pol     ::Int64                 # polarization direction
    ϵ       ::Float64               # elipticity
    u       ::QuasiVector{Float64, Tuple{UnitRange{Int64}}}       # Polarization unit vector
    Sv      ::Function
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

Defines a atomic (bound) electron with Z, n, l,

- `Z::Int64`            ... Atomic number of atom
- `n::Int64`            ... Principal quantum number
- `l::Int64`            ... Orbital angular momentum quantum number
- `j::Rational{Int64}`  ... Total angular momentum quantum number
- `ε::Float64`          ... Ionization potential (in Hartree atomic unit)
- `r::Vector{Float64}`  ... Radial grid
- `P::Vector{Float64}`  ... Large component
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

- `ε::Float64`
- `l::Int64`
- `j::Rational{Int64}`
- `P::Vector{Float64}`
- `δ::Float64`
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

Defines the output from the energy distribution calculation.

- `θ::Float64`                  ...
- `ϕ::Float64`                  ...
- `energies::Vector{Float64}`   ...

"""
struct EnergyDistribution
    θ::Float64
    ϕ::Float64
    energies::Vector{Float64}
end


"""
    AngularDistribution

Defines the output from the energy distribution calculation.

- `θ::Float64`                  ...
- `ϕ::Vector{Float64}`          ...
- `energy::Float64`             ...

"""
struct AngularDistribution
    θ::Float64
    ϕ::Vector{Float64}
    energy::Float64
end


"""
    MomentumDistribution

Defines the output from the momentum distribution calculation.

- `θ::Float64`                  ...
- `ϕ::Vector{Float64}`          ...
- `energies::Vector{Float64}`   ...

"""
struct MomentumDistribution
    θ::Float64
    ϕ::Vector{Float64}
    energies::Vector{Float64}
end


"""
    compute_energy_distribution()

Computes the energy distribution...
"""
function compute_energy_distribution(p_electron::ContinuumElectron, a_electron::AtomicElectron, pulse::Pulse, settings)
    
end


"""
    compute_momentum_distribution()

Computes the momentum distribution...
"""
function compute_momentum_distribution()
    println("Not implemented")
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
    if abs(m) > l return zero(ComplexF64) end
    if m ≥ 0
        P = AssociatedLegendrePolynomials.legendre(LegendreSphereNorm(), l, m, cos(θ))
        return P * exp(im * m * ϕ)
    else
        P = AssociatedLegendrePolynomials.legendre(LegendreSphereNorm(), l, -m, cos(θ))
        return (-1)^m * conj(P * exp(im * (-m) * ϕ))
    end
end

using AssociatedLegendrePolynomials
using SpecialFunctions
using WignerSymbols

export Pulse, AtomicElectron, ContinuumElectron, ContinuumSolution, ClebschGordan

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


struct ContinuumElectron
    ε::Float64
    p::Float64
    solution::Symbol
end


struct AtomicElectron
    ε::Float64
    n::Int64
    l::Int64
    j::Rational{Int64}
    r::Vector{Float64}
    P::Vector{Float64}
end


struct PartialWave
    ε::Float64
    l::Int64
    j::Rational{Int64}
    P::Vector{Float64}
    δ::Float64
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

using QuadGK

export sin2Sv


@doc raw"""
The pulse shape integral
"""
function F1_integral(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64 ; sign=1)

    integrand(t) = pulse.f(t) * exp(-im*(a_electron.ε + sign*pulse.ω)*t + im*pulse.Sv(t, θ, pulse, a_electron, p_electron))
    return pulse.A₀ * exp(-im * sign * pulse.ϕ) * quadgk(integrand, 0.0, pulse.Tp, rtol=1e-10)[1]
end


@doc raw"""
The pulse shape integral
"""
function F2_integral(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64)

    function A2(t::Float64)
        a = ( pulse.A₀ * pulse.f(t) * exp(-im * (pulse.ω * t + pulse.ϕ)) ) .* pulse.pol
        return (a' * a)
    end

    integrand(t) = A2(t) * exp( (-im * a_electron.ε * t) + im * pulse.Sv(t, θ, pulse, a_electron, p_electron) ) 

    return quadgk(integrand, 0.0, pulse.Tp, rtol=1e-10)[1]
end


"""
    sin2Sv(t::Float64, θ::Float64, pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

Defines a sin² envelope.
"""
function sin2Sv(t::Float64, θ::Float64, pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

    # Inonization potential
    εp = a_electron.ε
    # Laser pulse 
    ω  = pulse.ω ;   Up = pulse.Up ;    np = pulse.np;     ξ  = pulse.ϕ - pulse.pol * θ
    a  =  pulse.A₀ * p_electron.p * sin(θ) / (sqrt(2) * ω)

    ck = [-1/4, 1/2, -1/4]

    if t < 0 
        return εp * t
    elseif t ≤ pulse.Tp
        tau = t
    else
        tau = pulse.Tp
    end

    # Compute the additional terms evaluated at tau
    term1 = (3/8) * Up * tau
    term2 = - (Up * np) / (2 * ω) * sin((ω * tau) / np)
    term3 = (Up * np) / (16 * ω) * sin(2 * (ω * tau) / np)
    sum_term = sum((ck[k+2] / (1 + k/np)) * (sin(ω * (1 + k/np) * tau + ξ) - sin(ξ)) for k in -1:1)
    f_tau = term1 + term2 + term3 + a * sum_term

    return εp * t + f_tau
end


# function sin2Sv(t::Float64, θ::Float64, pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

#     # Inonization potential
#     εp = a_electron.ε

#     ω  = pulse.ω
#     Up = pulse.Up
#     np = pulse.np
#     ξ  = pulse.ϕ - θ

#     a  =  pulse.A₀ * p_electron.p * sin(θ) / (sqrt(2) * ω)

#     ck = [-1/4, 1/2, -1/4]

#     if t < 0 
#         return εp * t

#     elseif t ≤ pulse.Tp
#         sum_term = 0.0
#         for k in -1:1
#             coeff = ck[k + 2] / (1 + k / np)  # Adjust index: k=-1 -> 1, k=0 -> 2, k=1 -> 3
#             sum_term += coeff * (sin(ω * (1 + k / np) * t + ξ) - sin(ξ))
#         end
#         return εp * t + (3/8) * Up * t -
#                (Up * np / (2 * ω)) * sin(ω * t / np) +
#                (Up * np / (16 * ω)) * sin(2 * ω * t / np) +
#                a * sum_term

#     elseif t ≥ pulse.Tp
#         sum_term = 0.0
#         for k in -1:1
#             coeff = ck[k + 2] / (1 + k / np)  # Adjust index as before
#             sum_term += coeff * (sin(ω * (1 + k / np) * Tp + ξ) - sin(ξ))
#         end
#         return εp * t + (3/8) * Up * Tp -
#                (Up * np / (2 * ω)) * sin(ω * Tp / np) +
#                (Up * np / (16 * ω)) * sin(2 * ω * Tp / np) +
#                a * sum_term
#     end
# end





################################################################## Danish's Part #################################################################

# abstract type Sign end

# struct Plus <: Sign end
# struct Minus <: Sign end

# signval(::Plus) = 1
# signval(::Minus) = -1

# abstract type SaddleType end

# struct Type1 <: SaddleType end  # For F₁ (with ±ω term)
# struct Type2 <: SaddleType end  # For F₂ (no ω term)

# # Saddle Point Equation for two pulse shape integrals
# function saddle_point_eq(t, pulse::Pulse, momentum::Abstractmomentum, ionization_potential::AbstractIonizationPotential, s::Sign, ::Type1)
#     A = vector_potential(t, pulse)
#     E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
#     return E - (ionization_potential.value + signval(s) * pulse.omega)
# end

# function saddle_point_eq(t, pulse::Pulse, momentum::Abstractmomentum, ionization_potential::AbstractIonizationPotential, s::Sign, ::Type2)
#     A = vector_potential(t, pulse)
#     E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
#     return E - ionization_potential.value
# end

# # Corresponding phase for the pulse shape integrals
# function phase(momentum::Abstractmomentum, t, s::Sign, pulse::Pulse, ionization_potential::AbstractIonizationPotential, ::Type1)
#     S_V = volkov_phase(momentum, t, pulse, ionization_potential)
#     return -1im * (ionization_potential.value + signval(s) * pulse.omega) * t + 1im * S_V
# end

# function phase(momentum::Abstractmomentum, t, s::Sign, pulse::Pulse, ionization_potential::AbstractIonizationPotential, ::Type2)
#     S_V = volkov_phase(momentum, t, pulse, ionization_potential)
#     return -1im * ionization_potential.value * t + 1im * S_V
# end

# # Double Derivative of Volkov Phase
# function S_double_prime(E::ElectricField, momentum::Abstractmomentum, A::VectorPotential)
#     result = -(E.x * (tilde.x + A.x) + E.y * (tilde.y + A.y) + E.z * (tilde.z + A.z))
#     return result
# end

# # Calculation for saddle point solutions
# function find_saddle_points(pulse::Pulse, momentum::Abstractmomentum, ionization_potential::AbstractIonizationPotential, sign::Sign, type::SaddleType)
#     len = pulse.num_cycles > 4.0 ? pulse.num_cycles * 10 : 20
#     real_grid = range(0, pulse.T_p, length=len)
#     imag_grid = range(0, pulse.T, length=len)
#     time_grid = [Complex(re, im) for im in imag_grid, re in real_grid]

#     fs(t_s) = saddle_point_eq(t_s, pulse, momentum, ionization_potential, sign, type)

#     roots = []
#     for guess in time_grid
#         result = nlsolve(x -> [fs(x[1])], [guess])
#         if result.f_converged
#             t_s_root = round(real(result.zero[1]), digits=2) + round(imag(result.zero[1]), digits=2) * im
#             if (-pulse.T_p < real(t_s_root) < pulse.T_p) && (0 < imag(t_s_root))
#                 push!(roots, t_s_root)
#             end
#         end
#     end

#     isempty(roots) && push!(roots, 0 + 0im)
#     return sort(unique(roots), by=x -> real(x))
# end

# # Total contributions from the Pulse Shape Integrals
# function Integrals(pulse::Pulse, momentum::Abstractmomentum, ionization_potential::AbstractIonizationPotential, sign::Sign, envelope, type::SaddleType)
#     saddles = find_saddle_points(pulse, momentum, ionization_potential, sign, type)
#     total = 0.0 + 0.0im
#     for ts in saddles
#         A = vector_potential(ts, pulse)
#         E = electric_field(ts, pulse)
#         fval = isa(type, Type1) ? envelope.sin2(ts) : (A.x^2 + A.y^2 + A.z^2)
#         Spp = S_double_prime(E, momentum, A)
#         prefac = sqrt(2 * pi / (1im * Spp))
#         contrib = fval * exp(phase(momentum, ts, sign, pulse, ionization_potential, type)) * prefac
#         total += contrib
#     end
#     return isa(type, Type1) ? pulse.A0 * exp(-1im * signval(sign) * pulse.phi_CEP) * total : total
# end


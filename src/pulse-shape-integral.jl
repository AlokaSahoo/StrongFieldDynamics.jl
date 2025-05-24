using QuadGK

export sin2Sv


@doc raw"""
The pulse shape integral
"""
function F1_integral(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ ; sign=1)

    integrand(t) = pulse.f(t) * exp(-im*(a_electron.ε + sign*pulse.ω)*t + im*pulse.Sv(t, θ, pulse, a_electron, p_electron))
    return pulse.A₀ * exp(-im * sign * pulse.ϕ) * quadgk(integrand, 0.0, pulse.Tp, rtol=1e-8)[1]
end


@doc raw"""
The pulse shape integral
"""
function F2_integral(pulse::Pulse, a_electron::AtomicElectron)
    return zero(ComplexF64)
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
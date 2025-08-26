using ApproxFun
using LinearAlgebra
using QuadGK
using SpecialFunctions

export sin2Sv


"""
    F1_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F1_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

    # g(t)  = StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * t
    gp(t) = StrongFieldDynamics.Sv_prime_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω)

    ga = StrongFieldDynamics.Sv_general(pulse.duration[begin], θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * pulse.duration[begin]
    gb = StrongFieldDynamics.Sv_general(pulse.duration[end], θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * pulse.duration[end]

    res = levin_integrate_approxfun(pulse.f, gp, pulse.duration[begin], pulse.duration[end], ga, gb)

    return pulse.A₀ * exp(-im * sign * pulse.cep) * res
end


"""
    F2_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F2_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    a2(t) = ( pulse.A₀^2 / (1 + pulse.ϵ^2) ) * pulse.f(t)^2 * ( cos(pulse.ω*t + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*t + pulse.cep)^2 )
    # g(t)  = StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε * t
    gp(t) = StrongFieldDynamics.Sv_prime_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε 

    ga = StrongFieldDynamics.Sv_general(pulse.duration[begin], θ, ϕ, pulse, p_electron) - a_electron.ε * pulse.duration[begin]
    gb = StrongFieldDynamics.Sv_general(pulse.duration[end], θ, ϕ, pulse, p_electron) - a_electron.ε * pulse.duration[end]

    return levin_integrate_approxfun(a2, gp, pulse.duration[begin], pulse.duration[end], ga, gb)
end



"""
    Sv_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

Computes the Volkov phase for a sin² pulse envelope at time `t`.
"""
function Sv_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    term1 = p_electron.ε * t

    integrand2(τ)  = pulse.f(τ) * ( cos(pulse.ω*τ + pulse.cep) * cos(ϕ) + pulse.helicity * pulse.ϵ * sin(pulse.ω*τ + pulse.cep) * sin(ϕ) )
    int2, _ = quadgk(integrand2, pulse.duration[begin], t, maxevals=1e5) 
    term2 = int2 * pulse.A₀ * p_electron.p * sin(θ) / sqrt(1+pulse.ϵ^2)

    integrand3(τ)  = pulse.f(τ)^2 * ( cos(pulse.ω*τ + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*τ + pulse.cep)^2 )
    int3, _ = quadgk(integrand3, pulse.duration[begin], t, maxevals=1e5) 
    term3 = int3 * (pulse.A₀)^2 / (1+pulse.ϵ^2) / 2.0

    return term1 + term2 + term3

end


"""
    Sv_prime_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

Computes the time derivative of the Volkov phase for a sin² pulse envelope at time `t`.
"""
function Sv_prime_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    term1 = p_electron.ε

    term2 = ( pulse.A₀ * p_electron.p * sin(θ) / sqrt(1+pulse.ϵ^2) ) * 
            pulse.f(t) * ( cos(pulse.ω*t + pulse.cep) * cos(ϕ) + pulse.helicity*pulse.ϵ * sin(pulse.ω*t + pulse.cep) * sin(ϕ) )

    term3 = ( (pulse.A₀^2) / (1+pulse.ϵ^2) / 2.0 ) * pulse.f(t)^2 * ( cos(pulse.ω*t + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*t + pulse.cep)^2 )

    return term1 + term2 + term3

end

"""
    F1_integral_levin_approxfun(color::Pulse, pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F1_integral_levin_approxfun(color::Pulse, pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

    # g(t)  = StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * t
    gp(t) = StrongFieldDynamics.Sv_prime_general_cartesian(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * color.ω)

    ga = StrongFieldDynamics.Sv_general_cartesian(pulse.duration[begin], θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * color.ω) * pulse.duration[begin]
    gb = StrongFieldDynamics.Sv_general_cartesian(pulse.duration[end], θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * color.ω) * pulse.duration[end]

    res = levin_integrate_approxfun(color.f, gp, pulse.duration[begin], pulse.duration[end], ga, gb)

    return color.A₀ * exp(-im * sign * color.cep) * res
end


"""
    F2_integral_levin_approxfun(pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F2_integral_levin_approxfun(pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    a2(t)   = pulse.Ax(t)^2 + pulse.Ay(t)^2
    # a2(t) = ( pulse.A₀^2 / (1 + pulse.ϵ^2) ) * pulse.f(t)^2 * ( cos(pulse.ω*t + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*t + pulse.cep)^2 )
    # g(t)  = StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε * t
    gp(t) = StrongFieldDynamics.Sv_prime_general_cartesian(t, θ, ϕ, pulse, p_electron) - a_electron.ε 

    ga = StrongFieldDynamics.Sv_general_cartesian(pulse.duration[begin], θ, ϕ, pulse, p_electron) - a_electron.ε * pulse.duration[begin]
    gb = StrongFieldDynamics.Sv_general_cartesian(pulse.duration[end], θ, ϕ, pulse, p_electron) - a_electron.ε * pulse.duration[end]

    return levin_integrate_approxfun(a2, gp, pulse.duration[begin], pulse.duration[end], ga, gb)
end


"""
    F1_integral_levin_approxfun(pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F1_integral_levin_approxfun(pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    # g(t)  = StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * t
    gp(t) = StrongFieldDynamics.Sv_prime_general_cartesian(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * color.ω)

    ga = StrongFieldDynamics.Sv_general_cartesian(pulse.duration[begin], θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * color.ω) * pulse.duration[begin]
    gb = StrongFieldDynamics.Sv_general_cartesian(pulse.duration[end], θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * color.ω) * pulse.duration[end]

    res = levin_integrate_approxfun(color.f, gp, pulse.duration[begin], pulse.duration[end], ga, gb)

    return color.A₀ * exp(-im * sign * color.cep) * res
end



"""
    Sv_general_cartesian(t::Float64, θ::Float64, ϕ::Float64, pulse::Union{Pulse, MultiColor}, p_electron::ContinuumElectron)

Computes the Volkov phase for a general pulse envelope at time `t`.
"""
function Sv_general_cartesian(t::Float64, θ::Float64, ϕ::Float64, pulse::Union{Pulse, MultiColor}, p_electron::ContinuumElectron)

    p = spherical2cartesian(p_electron.p, θ, ϕ)

    # First term
    term1 = p_electron.ε * t

    int_pxAx, _ = quadgk(t -> p.x * pulse.Ax(t), pulse.duration[begin], t, maxevals=1e5)
    int_pyAy, _ = quadgk(t -> p.y * pulse.Ay(t), pulse.duration[begin], t, maxevals=1e5)

    int_Ax2 , _ = quadgk(t -> pulse.Ax(t)^2, pulse.duration[begin], t, maxevals=1e5)
    int_Ay2 , _ = quadgk(t -> pulse.Ay(t)^2, pulse.duration[begin], t, maxevals=1e5)

    return ( term1 + int_pxAx + int_pyAy + 0.5 * int_Ax2 + 0.5 * int_Ay2 )

end


"""
    Sv_prime_general_cartesian(t::Float64, θ::Float64, ϕ::Float64, pulse::Union{Pulse, MultiColor}, p_electron::ContinuumElectron)

Computes the Volkov phase for a sin² pulse envelope at time `t`.
"""
function Sv_prime_general_cartesian(t::Float64, θ::Float64, ϕ::Float64, pulse::Union{Pulse, MultiColor}, p_electron::ContinuumElectron)

    p = spherical2cartesian(p_electron.p, θ, ϕ)

    return ( p_electron.ε + p.x * pulse.Ax(t) + p.y * pulse.Ay(t) + 0.5 * pulse.Ax(t)^2 + 0.5 * pulse.Ay(t)^2 )

end


function pA(t::Float64, θ::Float64, ϕ::Float64, pulse::StrongFieldDynamics.MultiColor, p_electron::ContinuumElectron)
    p_electron.p * sin(θ) / sqrt(2) * ( pulse.colors[1].A₀ * pulse.colors[1].f(t) * cos(-ϕ + pulse.colors[1].ω * t + pulse.colors[1].cep) +
                                        pulse.colors[2].A₀ * pulse.colors[2].f(t) * cos( ϕ + pulse.colors[2].ω * t + pulse.colors[2].cep) )
end


function squareA(t::Float64, pulse::MultiColor)
    term1 = 0.5 * (pulse.colors[1].A₀^2 * pulse.colors[1].f(t)^2 + pulse.colors[2].A₀^2 * pulse.colors[2].f(t)^2)
    costerm = (pulse.colors[1].ω + pulse.colors[2].ω) * t + (pulse.colors[1].cep + pulse.colors[2].cep)
    term2 = pulse.colors[1].A₀ * pulse.colors[1].f(t) * pulse.colors[2].A₀ * pulse.colors[2].f(t) * cos(costerm)

    return term1 + term2
end


"""
    levin_integrate_approxfun(f::Function, gp::Function, a::Float64, b::Float64, ga::Float64, gb::Float64)

Compute the highly oscillatory integral ∫[a,b] f(x) exp(ig(x)) dx using Levin's method with ApproxFun.

This implementation uses ApproxFun's automatic differentiation and adaptive function approximation
to solve the differential equation ψ'(x) + ig'(x)ψ(x) = f(x) with ψ(a) = 0.

Arguments:
- f: amplitude function
- gp: derivative of phase function g'(x)
- a, b: integration bounds
- ga, gb: g(a) and g(b) at the integration limits

Returns:
- ComplexF64: Levin integral approximation ψ(b)e^{ig(b)} - ψ(a)e^{ig(a)}
"""
function levin_integrate_approxfun(f::Function, gp::Function, a::Float64, b::Float64, ga::Float64, gb::Float64)
    domain = a..b
    space = Chebyshev(domain)
    
    # Create functions on the domain
    F = Fun(f, space)
    Gp = Fun(gp, space)
    
    # Set up the differential operator: d/dx + i*g'(x)
    D = Derivative(space)
    A = D + im * Gp
    
    # Could be controversial
    # Solve the boundary value problem: A*ψ = F with ψ(a) = 0
    # Using the constraint ψ(a) = 0
    B = Evaluation(space, a)
    psi = [B; A] \ [0; F]
    
    # println("psi(0) $(psi(0))")
    
    # Compute Levin integral
    result = psi(b) * exp(im * gb) - psi(a) * exp(im * ga)
    
    return result

end

#=======================================================================================================================

=======================================================================================================================#

"""
    F1_integral_quadgk(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)


"""
function F1_integral_quadgk(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

    # integrand(t) = pulse.f(t) * exp(-im*(a_electron.ε + sign*pulse.ω)*t + im*sin2Sv_quad(t, θ, ϕ, pulse, p_electron))
    # res, _ = quadgk(integrand, 0.0, pulse.Tp, maxevals=10^1)
    integrand_real(t) = pulse.f(t) * cos(-(a_electron.ε + sign*pulse.ω)*t + sin2Sv_quadgk(t, θ, ϕ, pulse, p_electron))
    integrand_imag(t) = pulse.f(t) * sin(-(a_electron.ε + sign*pulse.ω)*t + sin2Sv_quadgk(t, θ, ϕ, pulse, p_electron))
    real_part, _ = quadgk(integrand_real, 0.0, pulse.Tp, maxevals=10^5)
    imag_part, _ = quadgk(integrand_imag, 0.0, pulse.Tp, maxevals=10^5)

    return pulse.A₀ * exp(-im * sign * pulse.cep) * (real_part + im * imag_part)
    # return pulse.A₀ * exp(-im * sign * pulse.cep) * res
end


"""
    F2_integral_quadgk(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)


"""
function F2_integral_quadgk(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    function A2(t::Float64)
        a = ( ( pulse.A₀ / sqrt(1 + pulse.ϵ^2) )  * pulse.f(t) * exp(-im * (pulse.ω * t + pulse.cep)) ) .* pulse.helicity
        return abs2(a)
    end

    integrand(t) = A2(t) * exp( (-im * a_electron.ε * t) + im * sin2Sv_quadgk(t, θ, ϕ, pulse, p_electron) ) 

    return quadgk(integrand, 0.0, pulse.Tp, maxevals=10^5)[1]
end


"""
    F1_integral_quadgk(color::Pulse, pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)


"""
function F1_integral_quadgk(color::Pulse, pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

    # integrand(t) = pulse.f(t) * exp(-im*(a_electron.ε + sign*pulse.ω)*t + im*sin2Sv_quad(t, θ, ϕ, pulse, p_electron))
    # res, _ = quadgk(integrand, 0.0, pulse.Tp, maxevals=10^1)
    integrand_real(t) = color.f(t) * cos(-(a_electron.ε + sign*color.ω)*t + StrongFieldDynamics.Sv_general_cartesian(t, θ, ϕ, pulse, p_electron))
    integrand_imag(t) = color.f(t) * sin(-(a_electron.ε + sign*color.ω)*t + StrongFieldDynamics.Sv_general_cartesian(t, θ, ϕ, pulse, p_electron))
    real_part, _ = quadgk(integrand_real, 0.0, color.Tp, maxevals=10^5)
    imag_part, _ = quadgk(integrand_imag, 0.0, color.Tp, maxevals=10^5)

    return color.A₀ * exp(-im * sign * color.cep) * (real_part + im * imag_part)
    # return pulse.A₀ * exp(-im * sign * pulse.cep) * res
end


#=============================================== From JAC ==========================================================#
"""
    F1_integral_jac(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, thetap::Float64, phip::Float64 ; sign=1)
    
Returns the integration result of type `ComplexF64`.
"""
function F1_integral_jac(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, thetap::Float64, phip::Float64 ; sign=1)
    wa = zero(ComplexF64)   
    np = pulse.np;   Tp = 2pi * np / pulse.ω
    omega = pulse.ω
    phiCep = pulse.cep;
    sinSqrArg = 0.5 * omega / np;
    lambda = pulse.helicity
    
    p = sqrt(2.0*p_electron.ε)
    px = p*sin(thetap)*cos(phip)
    py = p*sin(thetap)*sin(phip)
    
    epsilon = pulse.ϵ
    
    A0eps = pulse.A₀/sqrt(1.0 + epsilon^2)
    
    #Define Gauss-Legendre grid, convergence is typically good for orderGL = 100 * np (time consuming for np > 10); tested up to np = 20
    if  np <= 10     orderGL = 200*np
    else             orderGL = 2000
    end
    tgrid, weights = QuadGK.gauss(orderGL,0.0,Tp)
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        t = tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        cosIntegral = 0.25 / (omega * (np^2-1)) * (  2*sin(phiCep) + 2 * (np^2-1) * sin(phiCep + omega*t) - np * ( (1+np)*sin( phiCep + (np-1)/np * omega*t ) + (np-1) * sin( phiCep + (np+1)/np * omega*t )  ) )
        
        sinIntegral = 0.25 / (omega * (np^2-1)) * ( -2*cos(phiCep) - 2 * (np^2-1) * cos(phiCep + omega*t) + np * ( (1+np)*cos( phiCep + (np-1)/np * omega*t ) + (np-1) * cos( phiCep + (np+1)/np * omega*t )  ) )
        
        cos2Integral = sin(2*phiCep)/omega * ( -6 - np/(np-1) - np/(np+1) + 8*np/(2*np-1) + 8*np/(2*np+1) )
                        + 12*t + 6/omega * cos(2*omega*t) * sin(2*phiCep) + 6/omega * cos(2*phiCep) * sin(2*omega*t) - 16/omega*np*sin(omega*t/np) + 2/omega*np*sin(2*omega*t/np)
                        - 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t) + np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t))
                        + np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np * omega*t)) - 8*np/(omega*(2*np-1))*sin(2*phiCep + (2*np-1)/np * omega*t )
        cos2Integral = cos2Integral / 64
        
        sin2Integral = 12*t + 6/omega * sin(2*phiCep) * ( 1/(1-5*np^2+4*np^4) - cos(2*omega*t) ) - 6/omega * cos(2*phiCep)*sin(2*omega*t)
                        - 16/omega * np * sin(omega*t/np) + 2/omega * np * sin(2*omega*t/np) + 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t)
                        - np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t )) - np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np*omega*t)) + 8/omega * np/(2*np-1) * sin(2*phiCep + (2*np-1)/np*omega*t)
        sin2Integral = sin2Integral / 64
        
        SVolkov = p_electron.ε*t + A0eps*px*cosIntegral + A0eps*lambda*epsilon*py*sinIntegral + 0.5 * A0eps^2 * ( cos2Integral + epsilon^2 * sin2Integral )
        
        #Compute integrand at gridpoint t
        # if  plus    integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( a_electron.ε + beam.omega ) * t - SVolkov ) )
        # else        integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( a_electron.ε - beam.omega ) * t - SVolkov ) )
        # end

        integrand = sin( sinSqrArg * t )^2 * exp( -im * ( ( a_electron.ε + sign * omega ) * t - SVolkov ) )
        
        #Gauss-Legendre sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    # if  plus    wa = beam.A0 * exp(-im * phiCep) * wa
    # else        wa = beam.A0 * exp(im * phiCep) * wa
    # end

    wa = pulse.A₀ * exp(-sign * im * phiCep) * wa

    return wa
end

"""
    F2_integral_jac(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, thetap::Float64, phip::Float64)

Returns the integration result of type ComplexF64.
"""
function F2_integral_jac(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, thetap::Float64, phip::Float64)

    wa = zero(ComplexF64)   
    np = pulse.np;   Tp = 2pi * np / pulse.ω
    omega = pulse.ω
    phiCep = pulse.cep;
    sinSqrArg = 0.5 * omega / np;
    lambda = pulse.helicity
    
    p  = sqrt(2.0*p_electron.ε)
    px = p*sin(thetap)*cos(phip)
    py = p*sin(thetap)*sin(phip)
    
    epsilon = pulse.ϵ
    
    A0eps = pulse.A₀/sqrt(1.0 + epsilon^2)
    
    #Define Gauss-Legendre grid, convergence is typically good for orderGL = 100 * np (time consuming for np > 10); tested up to np = 20
    if  np <= 10     orderGL = 500*np
    else             orderGL = 1000
    end
    tgrid, weights = QuadGK.gauss(orderGL,0.0,Tp)
    
    
    #Sum over grid and compute Gauss-Legendre sum
    for    j = 1:orderGL
        t = tgrid[j]
        
        #Compute Volkov phase at gridpoint t
        cosIntegral = 0.25 / (omega * (np^2-1)) * (  2*sin(phiCep) + 2 * (np^2-1) * sin(phiCep + omega*t) - np * ( (1+np)*sin( phiCep + (np-1)/np * omega*t ) + (np-1) * sin( phiCep + (np+1)/np * omega*t )  ) )
        
        sinIntegral = 0.25 / (omega * (np^2-1)) * ( -2*cos(phiCep) - 2 * (np^2-1) * cos(phiCep + omega*t) + np * ( (1+np)*cos( phiCep + (np-1)/np * omega*t ) + (np-1) * cos( phiCep + (np+1)/np * omega*t )  ) )
        
        cos2Integral = sin(2*phiCep)/omega * ( -6 - np/(np-1) - np/(np+1) + 8*np/(2*np-1) + 8*np/(2*np+1) )
                        + 12*t + 6/omega * cos(2*omega*t) * sin(2*phiCep) + 6/omega * cos(2*phiCep) * sin(2*omega*t) - 16/omega*np*sin(omega*t/np) + 2/omega*np*sin(2*omega*t/np)
                        - 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t) + np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t))
                        + np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np * omega*t)) - 8*np/(omega*(2*np-1))*sin(2*phiCep + (2*np-1)/np * omega*t )
        cos2Integral = cos2Integral / 64
        
        sin2Integral = 12*t + 6/omega * sin(2*phiCep) * ( 1/(1-5*np^2+4*np^4) - cos(2*omega*t) ) - 6/omega * cos(2*phiCep)*sin(2*omega*t)
                        - 16/omega * np * sin(omega*t/np) + 2/omega * np * sin(2*omega*t/np) + 8*np/(omega*(1+2*np)) * sin(2*phiCep + (2+1/np)*omega*t)
                        - np/(omega*(np-1)) * sin(2*(phiCep + (np-1)/np *omega*t )) - np/(omega*(1+np)) * sin(2*(phiCep + (np+1)/np*omega*t)) + 8/omega * np/(2*np-1) * sin(2*phiCep + (2*np-1)/np*omega*t)
        sin2Integral = sin2Integral / 64
        
        SVolkov = p_electron.ε*t + A0eps*px*cosIntegral + A0eps*lambda*epsilon*py*sinIntegral + 0.5 * A0eps^2 * ( cos2Integral + epsilon^2 * sin2Integral )
        
        #Compute integrand at gridpoint t
        integrand = sin( sinSqrArg * t )^4 * ( cos(omega*t+phiCep)^2 + epsilon^2 * sin(omega*t+phiCep)^2 ) * exp( -im * ( a_electron.ε * t - SVolkov ) )
        
        #Gauss-Legendre sum
        wa = wa + weights[j] * integrand
    end
    
    #Multiply with global factor
    wa = A0eps^2 * wa
    
    return( wa )
end


"""
    sin2Sv(t::Float64, θ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

Defines the Volkov phase integral over time for a sin² envelope. Only for circular polarization from the Birger PRA 2020.
"""
function sin2Sv(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    # Phtotelectron energy
    εp = p_electron.ε
    # Laser pulse 
    ω  = pulse.ω ;   Up = pulse.Up ;    np = pulse.np;     ξ  = pulse.cep - pulse.helicity * ϕ
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

# """
# sin2Sv_quad(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

# """
# function sin2Sv_quad(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

#     # Phtotelectron energy
#     εp = p_electron.ε
#     # Laser pulse 
#     ω  = pulse.ω ;   Up = pulse.Up ;    np = pulse.np;     ξ  = pulse.cep - pulse.helicity * ϕ
#     # a  =  pulse.A₀ * p_electron.p * sin(θ) / (sqrt(2) * ω)

#     function sin2env(x) 
#         if x < 0 || x >= pulse.Tp
#             return 0.
#         end
#         sin(ω*x/2.0/np)^2
#     end

#     term1 = εp * t

#     int2, _ = quadgk(τ -> sin2env(τ) * cos(ω*τ + ξ), 0, t, maxevals=10^5) 
#     term2 = int2 * pulse.A₀ * p_electron.p * sin(θ) / sqrt(2.0)

#     int3, _ = quadgk(τ -> (sin2env(τ))^2, 0, t, maxevals=10^5) 
#     term3 = int3 * (pulse.A₀)^2 / 4.0

#     return term1 + term2 + term3
# end


"""
    sin2Sv_quadgk(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

"""
function sin2Sv_quadgk(t::Float64, θ::Float64, ϕ::Float64, pulse::StrongFieldDynamics.Pulse, p_electron::StrongFieldDynamics.ContinuumElectron)

    # Phtotelectron energy
    # εp = p_electron.ε
    # Laser pulse 
    # ω  = pulse.ω ;   Up = pulse.Up ;    np = pulse.np;     ξ  = pulse.cep - pulse.helicity * ϕ
    # a  =  pulse.A₀ * p_electron.p * sin(θ) / (sqrt(2) * ω)

    term1 = p_electron.ε * t

    int2, _ = quadgk(τ -> pulse.f(τ) * ( cos(pulse.ω * τ + pulse.cep) * cos(ϕ) + pulse.helicity * pulse.ϵ * sin(pulse.ω*τ + pulse.cep) * sin(ϕ) ), 0, t, maxevals=10^5) 
    term2 = int2 * pulse.A₀ * p_electron.p * sin(θ) / sqrt(1.0 + pulse.ϵ^2)

    int3, _ = quadgk(τ -> ( pulse.f(τ)^2 * (cos(pulse.ω * τ + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω * τ + pulse.cep)^2) ), 0, t, maxevals=10^5) 
    term3 = int3 * (pulse.A₀)^2 / 2.0 / (1.0 + pulse.ϵ^2)

    return term1 + term2 + term3
end


"""
sin2Sv_prime(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)


"""
function sin2Sv_prime(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    return p_electron.ε + 
            ( pulse.A₀ * p_electron.p * sin(θ) / sqrt(2.0) ) * pulse.f(t) * cos(pulse.ω * t + pulse.cep - pulse.helicity * ϕ) + 
            (pulse.A₀)^2 * (pulse.f(t))^2 / 4.0

end

"""
    sin2Sv_prime_b(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

Only for circular polarizaton,
Takes the analytical derivative from the expression in Appendix B, Birger 2020 for sin2 pulse
"""
function sin2Sv_prime_b(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    # Phtotelectron energy
    εp = p_electron.ε
    # Laser pulse 
    ω  = pulse.ω ;   Up = pulse.Up ;    np = pulse.np;     ξ  = pulse.cep - pulse.helicity * ϕ
    a  =  pulse.A₀ * p_electron.p * sin(θ) / (sqrt(2) * ω)

    ck = [-1/4, 1/2, -1/4]

    if t < 0 || t ≥ pulse.Tp
        return εp
    else
        # Compute the additional terms evaluated at t
        term1 = (3/8) * Up
        term2 = - (Up * np) / (2 * ω) * (ω / np) * cos((ω * t) / np)
        term3 = (Up * np) / (16 * ω) * (2 * ω /np) * cos(2 * (ω * t) / np)
        sum_term = sum((ck[k+2] / (1 + k/np)) * ( ω * (1 + k/np) * cos(ω * (1 + k/np) * t + ξ) ) for k in -1:1)

        return εp + term1 + term2 + term3 + a * sum_term
    end

end


"""
    gaussianSv(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)


"""
function gaussianSv(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    ξ  = pulse.cep - pulse.helicity * ϕ

    term1 = p_electron.ε * t

    term2 = pulse.Up * sqrt(π/log(4)) * (pulse.Tp / t) * erf( 2 * sqrt(log(4)) * (t/pulse.Tp) )

    term3 = pulse.A₀ * p_electron.p * sin(θ) / sqrt(2.0) * (pulse.Tp / 8) * sqrt(π/log(2)) * exp(-1 * pulse.ω^2 * pulse.Tp^2 / log(65536)) *
                        ( exp(-1 * im * ξ) * erf( (im*pulse.Tp^2 * pulse.ω + 8 * log(2) * t) / (4 * pulse.Tp * sqrt(log(2))) ) - 
                               exp(im * ξ) * erf( (im*pulse.Tp^2 * pulse.ω - 8 * log(2) * t) / (4 * pulse.Tp * sqrt(log(2))) ))

    return term1 + term2 + term3
end


############################################################# Numerical Phase Danish ################################################################
# Concrete types for Vector Potential components
# struct VectorPotential
#     x::Number
#     y::Number
#     z::Number
# end

# struct Observable
#     x::Number
#     y::Number
#     z::Number
# end


# Envelope function
function envelope_function(time::Number, pulse::Pulse)
    argument = pulse.ω * time / (2 * pulse.np)
    return sin(argument)^2
end

# Vector potential function returning a concrete VectorPotential object
function vector_potential(time::Number, pulse::Pulse)::Cartesian
    sqrt_factor = sqrt(1 + pulse.ϵ^2)
    f = envelope_function(time, pulse)
    x = pulse.A₀ * f / sqrt_factor * cos(pulse.ω * time + pulse.cep)
    y = pulse.ϵ * pulse.helicity * pulse.A₀ * f / sqrt_factor * sin(pulse.ω * time + pulse.cep)
    return Cartesian(x, y, 0.0)
end

function A2_danish(t, pulse)
    return integrate_A2(t, :x, pulse) + integrate_A2(t, :y, pulse) + integrate_A2(t, :z, pulse)
end

# Numerical integration for A² and p⋅A components
function integrate_A2(ts, component::Symbol, pulse::Pulse)
    integrand(t) = (component == :x ? vector_potential(t, pulse).x^2 :
                     component == :y ? vector_potential(t, pulse).y^2 :
                     vector_potential(t, pulse).z^2)
    result, _ = quadgk(integrand, 0.0, ts, maxevals=1e5)
    return result
end

function integrate_A(ts, component::Symbol, pulse::Pulse)
    integrand(t) =  (component == :x ? vector_potential(t, pulse).x :
                              component == :y ? vector_potential(t, pulse).y :
                              vector_potential(t, pulse).z)
    result, _ = quadgk(integrand, 0.0, ts, maxevals=1e5)
    return result
end

# Calculate S(t) numerically
function S(t, momentum::Cartesian, pulse::Pulse)
    kinetic_term = 0.5 * (momentum.x^2 + momentum.y^2 + momentum.z^2) * t
    potential_term = integrate_A2(t, :x, pulse) + integrate_A2(t, :y, pulse) + integrate_A2(t, :z, pulse)
    interaction_term = momentum.x * integrate_A(t, :x, pulse) + momentum.y * integrate_A(t, :y, pulse) + momentum.z * integrate_A(t, :z, pulse)
    # ip_term = ionization_potential.value * t

    return kinetic_term + 0.5 * potential_term + interaction_term #+ ip_term
end

function sin2Sv_danish(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)
    # momentum = Cartesian(p_electron.p * sin(θ) * cos(ϕ), p_electron.p * sin(θ) * sin(ϕ), p_electron.p * cos(θ))
    momentum = spherical2cartesian(p_electron.p, θ, ϕ)

    return S(t, momentum, pulse)

end

function sin2Sv_prime_danish(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)
    momentum = Cartesian(p_electron.p * sin(θ) * cos(ϕ), p_electron.p * sin(θ) * sin(ϕ), p_electron.p * cos(θ))
    A = vector_potential(t, pulse)
    E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
    return E
end

#======================================Obsolate===============================================#
# """
#     sin2Sv_danish(t::Float64, theta::Float64, phi::Float64, pulse::Pulse, p_electron::ContinuumElectron)


# """
# function sin2Sv_danish(t::Float64, theta::Float64, phi::Float64, pulse::Pulse, p_electron::ContinuumElectron)

#     # --- Convert Spherical Momentum to Cartesian ---
#     px = p_electron.p * sin(theta) * cos(phi)
#     py = p_electron.p * sin(theta) * sin(phi)
#     # The full squared magnitude of the momentum is simply p_electron.p^2
#     p_squared = p_electron.p^2

#     # Extract other parameters from the Pulse object
#     A₀ = pulse.A₀
#     ϵ = pulse.ϵ
#     λ_h = pulse.helicity
#     ω = pulse.ω
#     np = pulse.np
#     cep = pulse.cep

#     # --- Calculation based on Eq. A.16 ---
#     ω_vec = [ω * (1 - 1 / np), ω, ω * (1 + 1 / np)]
#     A_vec = [-0.25 * A₀, 0.5 * A₀, -0.25 * A₀]

#     # Term 1: Free-electron kinetic energy (uses the full p²) 
#     term1 = 0.5 * p_squared * t

#     # Term 2: Ponderomotive energy contribution
#     term2 = (t / 4.0) * sum(A_vec.^2)

#     # Term 3: Oscillatory diagonal term
#     term3_sum = sum((A_vec[i]^2 / (8.0 * ω_vec[i])) * sin(2.0 * ω_vec[i] * t + 2.0 * cep) for i in 1:3)
#     term3 = ((1.0 - ϵ^2) / (1.0 + ϵ^2)) * term3_sum

#     # Terms 4 & 5: Off-diagonal contributions
#     term4_sum = 0.0
#     term5_sum = 0.0
#     for i in 1:2
#         for j in (i+1):3
#             term4_sum += (A_vec[i] * A_vec[j] / (2.0 * (ω_vec[i] - ω_vec[j]))) * sin((ω_vec[i] - ω_vec[j]) * t)
#             term5_sum += (A_vec[i] * A_vec[j] / (2.0 * (ω_vec[i] + ω_vec[j]))) * sin((ω_vec[i] + ω_vec[j]) * t + 2.0 * cep)
#         end
#     end
#     term4 = term4_sum
#     term5 = ((1.0 - ϵ^2) / (1.0 + ϵ^2)) * term5_sum

#     # Terms 6 & 7: Momentum-field coupling terms (use pₓ and pᵧ) [cite: 95, 96]
#     term6 = (px / sqrt(1.0 + ϵ^2)) * sum((A_vec[i] / ω_vec[i]) * sin(ω_vec[i] * t + cep) for i in 1:3)
#     term7 = -ϵ * λ_h * (py / sqrt(1.0 + ϵ^2)) * sum((A_vec[i] / ω_vec[i]) * cos(ω_vec[i] * t + cep) for i in 1:3)

#     # Combine all terms for the final Volkov phase
#     S_p_t = term1 + term2 + term3 + term4 + term5 + term6 + term7

#     return S_p_t
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

# abstract type AbstractCoordinate end

# struct Cartesian <: AbstractCoordinate
#     x::Float64
#     y::Float64
#     z::Float64
# end

# struct Spherical <: AbstractCoordinate
#     r::Float64
#     θ::Float64
#     ϕ::Float64
# end

# # Saddle Point Equation for two pulse shape integrals
# function saddle_point_eq(t, pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, s::Sign, ::Type1)
#     A = vector_potential(t, pulse)
#     E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
#     return E - (ionization_potential + signval(s) * pulse.omega)
# end

# function saddle_point_eq(t, pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, s::Sign, ::Type2)
#     A = vector_potential(t, pulse)
#     E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
#     return E - ionization_potential
# end

# # Corresponding phase for the pulse shape integrals
# function phase(momentum::Cartesian, t, s::Sign, pulse::Pulse, ionization_potential::Float64, ::Type1)
#     S_V = volkov_phase(momentum, t, pulse, ionization_potential)
#     return -1im * (ionization_potential + signval(s) * pulse.omega) * t + 1im * S_V
# end

# function phase(momentum::Cartesian, t, s::Sign, pulse::Pulse, ionization_potential::Float64, ::Type2)
#     S_V = volkov_phase(momentum, t, pulse, ionization_potential)
#     return -1im * ionization_potential * t + 1im * S_V
# end

# # Double Derivative of Volkov Phase
# function S_double_prime(E::ElectricField, momentum::Cartesian, A::VectorPotential)
#     result = -(E.x * (tilde.x + A.x) + E.y * (tilde.y + A.y) + E.z * (tilde.z + A.z))
#     return result
# end

# # Calculation for saddle point solutions
# function find_saddle_points(pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, sign::Sign, type::SaddleType)
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
# function Integrals(pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, sign::Sign, envelope, type::SaddleType)
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

################################################################## (Modified) Danish's Part #################################################################

# abstract type SaddleType end

# struct Type1 <: SaddleType end  # For F₁ (with ±ω term)
# struct Type2 <: SaddleType end  # For F₂ (no ω term)

# abstract type AbstractCoordinate end

# struct Cartesian <: AbstractCoordinate
#     x::Float64
#     y::Float64
#     z::Float64
# end

# struct Spherical <: AbstractCoordinate
#     r::Float64
#     θ::Float64
#     ϕ::Float64
# end

# # Saddle Point Equation for two pulse shape integrals
# function saddle_point_eq(t, pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, s::Sign, ::Type1)
#     A = vector_potential(t, pulse)
#     E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
#     return E - (ionization_potential + sign * pulse.omega)
# end

# function saddle_point_eq(t, pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, s::Sign, ::Type2)
#     A = vector_potential(t, pulse)
#     E = 0.5 * ((momentum.x + A.x)^2 + (momentum.y + A.y)^2 + (momentum.z + A.z)^2)
#     return E - ionization_potential
# end

# # Corresponding phase for the pulse shape integrals
# function phase(momentum::Cartesian, t, s::Sign, pulse::Pulse, ionization_potential::Float64, ::Type1)
#     S_V = volkov_phase(momentum, t, pulse, ionization_potential)
#     return -1im * (ionization_potential + sign * pulse.omega) * t + 1im * S_V
# end

# function phase(momentum::Cartesian, t, sign::Int64, pulse::Pulse, ionization_potential::Float64, ::Type2)
#     S_V = volkov_phase(momentum, t, pulse, ionization_potential)
#     return -1im * ionization_potential * t + 1im * S_V
# end

# # Double Derivative of Volkov Phase
# function S_double_prime(E::ElectricField, momentum::Cartesian, A::VectorPotential)
#     result = -(E.x * (tilde.x + A.x) + E.y * (tilde.y + A.y) + E.z * (tilde.z + A.z))
#     return result
# end

# # Calculation for saddle point solutions
# function find_saddle_points(pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, sign::Sign, type::SaddleType)
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
# function Integrals(pulse::Pulse, momentum::Cartesian, ionization_potential::Float64, sign::Sign, envelope, type::SaddleType)
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
#     return isa(type, Type1) ? pulse.A0 * exp(-1im * sign * pulse.phi_CEP) * total : total
# end

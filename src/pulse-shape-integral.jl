using QuadGK

export sin2Sv

include("levin-integration.jl")

"""
    F1_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F1_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

    # g(t)  = StrongFieldDynamics.sin2Sv_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * t
    gp(t) = StrongFieldDynamics.sin2Sv_prime_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω)

    ga = 0.0
    gb = StrongFieldDynamics.sin2Sv_general(pulse.Tp, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * pulse.Tp

    res = levin_integrate_approxfun(pulse.f, gp, 0.0, pulse.Tp, ga, gb)

    return pulse.A₀ * exp(-im * sign * pulse.cep) * res
end


"""
    F2_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F2_integral_levin_approxfun(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    a2(t) = ( pulse.A₀^2 / (1 + pulse.ϵ^2) ) * pulse.f(t)^2 * ( cos(pulse.ω*t + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*t + pulse.cep)^2 )
    # g(t)  = StrongFieldDynamics.sin2Sv_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε * t
    gp(t) = StrongFieldDynamics.sin2Sv_prime_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε 

    ga = 0.0        # for sin2 pulse
    gb = StrongFieldDynamics.sin2Sv_general(pulse.Tp, θ, ϕ, pulse, p_electron) - a_electron.ε * pulse.Tp

    return levin_integrate_approxfun(a2, gp, 0.0, pulse.Tp, ga, gb)
end


"""
    F1_integral_levin(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F1_integral_levin(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64 ; sign=1)

    # g(t)  = StrongFieldDynamics.sin2Sv_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * t
    gp(t) = StrongFieldDynamics.sin2Sv_prime_general(t, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω)

    ga = 0.0
    gb = StrongFieldDynamics.sin2Sv_general(pulse.Tp, θ, ϕ, pulse, p_electron) - (a_electron.ε + sign * pulse.ω) * pulse.Tp

    res = levin_integrate(pulse.f, gp, 0.0, pulse.Tp, ga, gb, 100)

    return pulse.A₀ * exp(-im * sign * pulse.cep) * res
end


"""
    F2_integral_levin(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

Computes the pulse shape integral of the form:
Returns the integration result of type ComplexF64.
"""
function F2_integral_levin(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    a2(t) = ( pulse.A₀^2 / (1 + pulse.ϵ^2) ) * pulse.f(t)^2 * ( cos(pulse.ω*t + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*t + pulse.cep)^2 )
    # g(t)  = StrongFieldDynamics.sin2Sv_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε * t
    gp(t) = StrongFieldDynamics.sin2Sv_prime_general(t, θ, ϕ, pulse, p_electron) - a_electron.ε 

    ga = 0.0        # for sin2 pulse
    gb = StrongFieldDynamics.sin2Sv_general(pulse.Tp, θ, ϕ, pulse, p_electron) - a_electron.ε * pulse.Tp

    return levin_integrate(a2, gp, 0.0, pulse.Tp, ga, gb, 100)
end


"""
    F1_integral_quadgk(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64 ; sign=1)


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



#=============================================== From JAC ==========================================================#
@doc raw"""

    F1_integral_jac(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, thetap::Float64, phip::Float64 ; sign=1)
    
Computes the pulse shape integral of the form:
\mathcal{F}_1\left[\pm \omega ; f ; \mathbf{p}\right] = A_0 e^{\mp i \phi_{\text{CEP}}} \int_{-\infty}^{\infty} d\tau f(\tau) e^{-i(\varepsilon_i \pm \omega)\tau + i S_p(\tau)}, \\
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

@doc raw"""

    F2_integral_jac(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, thetap::Float64, phip::Float64)

Computes the pulse shape integral of the form:
\mathcal{F}_2\left[f ; \mathbf{p}\right] &= \int_{-\infty}^{\infty} d\tau A^2(\tau) e^{-i \varepsilon_p \tau + i S_p(\tau)}
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
    sin2Sv_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)


"""
function sin2Sv_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    term1 = p_electron.ε * t

    integrand2(τ)  = pulse.f(τ) * ( cos(pulse.ω*τ + pulse.cep) * cos(ϕ) + pulse.helicity * pulse.ϵ * sin(pulse.ω*τ + pulse.cep) * sin(ϕ) )
    int2, _ = quadgk(integrand2, 0, t, maxevals=1e5) 
    term2 = int2 * pulse.A₀ * p_electron.p * sin(θ) / sqrt(1+pulse.ϵ^2)

    integrand3(τ)  = pulse.f(τ)^2 * ( cos(pulse.ω*τ + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*τ + pulse.cep)^2 )
    int3, _ = quadgk(integrand3, 0, t, maxevals=1e5) 
    term3 = int3 * (pulse.A₀)^2 / (1+pulse.ϵ^2) / 2.0

    return term1 + term2 + term3

end


"""
    sin2Sv_prime_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)


"""
function sin2Sv_prime_general(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    term1 = p_electron.ε

    term2 = ( pulse.A₀ * p_electron.p * sin(θ) / sqrt(1+pulse.ϵ^2) ) * 
            pulse.f(t) * ( cos(pulse.ω*t + pulse.cep) * cos(ϕ) + pulse.helicity*pulse.ϵ * sin(pulse.ω*t + pulse.cep) * sin(ϕ) )

    term3 = ( (pulse.A₀^2) / (1+pulse.ϵ^2) / 2.0 ) * pulse.f(t)^2 * ( cos(pulse.ω*t + pulse.cep)^2 + pulse.ϵ^2 * sin(pulse.ω*t + pulse.cep)^2 )

    return term1 + term2 + term3

end


"""
    sin2Sv(t::Float64, θ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

Defines the Volkov phase integral over time for a sin² envelope. Only for circular polarization 
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

"""
sin2Sv_quad(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

"""
function sin2Sv_quad(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    # Phtotelectron energy
    εp = p_electron.ε
    # Laser pulse 
    ω  = pulse.ω ;   Up = pulse.Up ;    np = pulse.np;     ξ  = pulse.cep - pulse.helicity * ϕ
    # a  =  pulse.A₀ * p_electron.p * sin(θ) / (sqrt(2) * ω)

    function sin2env(x) 
        if x < 0 || x >= pulse.Tp
            return 0.
        end
        sin(ω*x/2.0/np)^2
    end

    term1 = εp * t

    int2, _ = quadgk(τ -> sin2env(τ) * cos(ω*τ + ξ), 0, t, maxevals=10^5) 
    term2 = int2 * pulse.A₀ * p_electron.p * sin(θ) / sqrt(2.0)

    int3, _ = quadgk(τ -> (sin2env(τ))^2, 0, t, maxevals=10^5) 
    term3 = int3 * (pulse.A₀)^2 / 4.0

    return term1 + term2 + term3
end


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

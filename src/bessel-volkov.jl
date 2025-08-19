struct Cartesian 
    x::Number
    y::Number
    z::Number
end

function polar2cartesian(r, θ, ϕ)
    x = r * sin(θ) * cos(ϕ)
    y = r * sin(θ) * sin(ϕ)
    z = r * cos(θ)

    return Cartesian(x, y, z)
end

# # Define the alpha functions
# function alpha_neg1(factor::Factor, pulse::Pulse)
#     return sqrt(factor.kappa / (4 * pi)) * factor.cn1 * besselj(pulse.mγ + 1, factor.kappa * pulse.r)
# end

# function alpha_0(factor::Factor, pulse::Pulse)
#     return sqrt(factor.kappa / (2 * pi)) * factor.c0 * besselj(pulse.mγ, factor.kappa * pulse.r)
# end

# function alpha_pos1(factor::Factor, pulse::Pulse)
#     return sqrt(factor.kappa / (4 * pi)) * factor.cp1 * besselj(pulse.mγ - 1, factor.kappa * pulse.r)
# end

# # Define the beta functions
# function beta_neg1(factor::Factor, pulse::Pulse)
#     return (pulse.mγ - 1) * pulse.ϕr + factor.k_z * factor.z
# end

# function beta_0(factor::Factor, pulse::Pulse)
#     return pulse.mγ * pulse.ϕr + factor.k_z * factor.z
# end

# function beta_pos1(factor::Factor, pulse::Pulse)
#     return (pulse.mγ + 1) * pulse.ϕr + factor.k_z * factor.z
# end


# Function to calculate the Fourier coefficients: beta, gamma, sigma
function CalculateFourierCoefficients(tilde::Cartesian, pulse::Pulse)
    
    # Calculate alpha
    alphaMinus1 = sqrt(pulse.χ/4/pi) * pulse.cms[-1] * besselj(pulse.mγ + 1, pulse.χ * pulse.rp)
    alpha0      = sqrt(pulse.χ/2/pi) * pulse.cms[0]  * besselj(pulse.mγ,     pulse.χ * pulse.rp)
    alphaPlus1  = sqrt(pulse.χ/4/pi) * pulse.cms[1]  * besselj(pulse.mγ - 1, pulse.χ * pulse.rp)

    # Calculate beta
    beta = (tilde.x^2 + tilde.y^2 + tilde.z^2) / 2 + (3 / 32) * pulse.A₀^2 * (2 * alphaMinus1^2 + 2 * alphaPlus1^2 + alpha0^2)
    
    # Calculate gamma
    gamma = (pulse.A₀ / (2 * pulse.ω)) * [
        -tilde.y * alphaPlus1, 
        tilde.z * alpha0, 
        tilde.y * alphaMinus1,
        (pulse.np * tilde.y * alphaPlus1) / (2 * (pulse.np - 1)), 
        -(pulse.np * tilde.z * alpha0) / (2 * (pulse.np - 1)), 
        -(pulse.np * tilde.y * alphaMinus1) / (2 * (pulse.np - 1)),
        (pulse.np * tilde.y * alphaPlus1) / (2 * (pulse.np + 1)), 
        -(pulse.np * tilde.z * alpha0) / (2 * (pulse.np + 1)), 
        -(pulse.np * tilde.y * alphaMinus1) / (2 * (pulse.np + 1))
    ]
    
    # Calculate sigma
    sigma = (pulse.A₀ / (2 * pulse.ω)) * [
        -pulse.A₀ * pulse.np / 4 * (2 * alphaPlus1^2 + 2 * alphaMinus1^2 + alpha0^2), 
        pulse.A₀ * pulse.np / 32 * (2 * alphaPlus1^2 + 2 * alphaMinus1^2 + alpha0^2), 
        -tilde.x * alphaPlus1, 
        3 * pulse.A₀ / 32 * (alpha0^2 - 4 * alphaMinus1 * alphaPlus1),
        -tilde.x * alphaMinus1, 
        pulse.A₀ * pulse.np / (8 * (2 * pulse.np - 1)) * (4 * alphaMinus1 * alphaPlus1 - alpha0^2),
        (pulse.np * tilde.x * alphaPlus1) / (2 * (pulse.np - 1)), 
        pulse.A₀ * pulse.np / (64 * (pulse.np - 1)) * (alpha0^2 - 4 * alphaMinus1 * alphaPlus1),
        (pulse.np * tilde.x * alphaMinus1) / (2 * (pulse.np - 1)),
        pulse.A₀ * pulse.np / (64 * (pulse.np + 1)) * (alpha0^2 - 4 * alphaMinus1 * alphaPlus1), 
        (pulse.np * tilde.x * alphaPlus1) / (2 * (pulse.np + 1)), 
        (pulse.np * tilde.x * alphaMinus1) / (2 * (pulse.np + 1)), 
        pulse.A₀ * pulse.np / (8 * (2 * pulse.np + 1)) * (4 * alphaMinus1 * alphaPlus1 - alpha0^2)
    ]

    return beta, gamma, sigma
end

# Function to calculate frequencies and phases: omegaC, omegaS, phiC, phiS
function CalculateFrequenciesAndPhases(pulse::Pulse)
    # Calculate omegaC
    omegaC = pulse.ω * [1, 1, 1, (pulse.np - 1) / pulse.np, (pulse.np - 1) / pulse.np, (pulse.np - 1) / pulse.np, (pulse.np + 1) / pulse.np, (pulse.np + 1) / pulse.np, (pulse.np + 1) / pulse.np]
    
    # Calculate omegaS
    omegaS = pulse.ω * [
        -1 / pulse.np, -2 / pulse.np, 1, 2, 1, 
        (2 * pulse.np - 1) / pulse.np, (pulse.np - 1) / pulse.np, 2 * (pulse.np - 1) / pulse.np, (pulse.np - 1) / pulse.np,
        2 * (pulse.np + 1) / pulse.np, (pulse.np + 1) / pulse.np, (pulse.np + 1) / pulse.np, (2 * pulse.np + 1) / pulse.np
    ]
    
    # Calculate phiC
    phiC = pulse.ϕr * [pulse.mγ - 1, pulse.mγ, pulse.mγ + 1, pulse.mγ - 1, pulse.mγ, pulse.mγ + 1, pulse.mγ - 1, pulse.mγ, pulse.mγ + 1]
    
    # Calculate phiS
    phiS = pulse.ϕr * [
        0, 0, pulse.mγ - 1, 2 * pulse.mγ, pulse.mγ + 1, 
        2 * pulse.mγ, pulse.mγ - 1, 2 * pulse.mγ, pulse.mγ + 1, 
        2 * pulse.mγ, pulse.mγ - 1, pulse.mγ + 1, 2 * pulse.mγ
    ]

    return omegaC, omegaS, phiC, phiS
end

# Function to calculate the Volkov phase
function VolkovPhase(t::Float64, θ::Float64, ϕ::Float64, pulse::Pulse, p_electron::ContinuumElectron)

    tilde = polar2cartesian(p_electron.p, θ, ϕ)

    beta, gamma, sigma = CalculateFourierCoefficients(tilde,  pulse)
    omegaC, omegaS, phiC, phiS = CalculateFrequenciesAndPhases(pulse)
    sumCos = sum(gamma[j] * cos(phiC[j] - omegaC[j] * t) for j in 1:9)
    sumSin = sum(sigma[l] * sin(phiS[l] - omegaS[l] * t) for l in 1:13)

    return beta * t + sumCos + sumSin
end
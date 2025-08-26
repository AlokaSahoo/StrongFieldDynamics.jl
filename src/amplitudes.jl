using Dierckx
using QuadGK

"""
    probality(pulse::AbstractPulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

Calculates the ionization probability
"""
function probability(pulse::AbstractPulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    prob = zero(Float64)

    for mj = -a_electron.j : a_electron.j
        for msp = -1//2:1//2
            amplitude = StrongFieldDynamics.T0( pulse, a_electron, p_electron, mj, msp, θ, ϕ)
            prob += abs2(amplitude)
        end
    end

    prob = prob * p_electron.p / (2.0*a_electron.j + 1)

end

"""
    T0(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, mj::Rational{Int64}, msp::Rational{Int64}, θ::Float64, ϕ::Float64)

Computes the direct scattering amplitude 
"""
function T0(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, mj::Rational{Int64}, msp::Rational{Int64}, θ::Float64, ϕ::Float64)
    # Maximum number of partial waves
    lp_max = 3

    l = a_electron.l
    j = a_electron.j
    r = a_electron.r
    u = pulse.u

    # Intitialize term1, term2, term3
    term1 = zero(ComplexF64) ; term2 = zero(ComplexF64) ; term3 = zero(ComplexF64) ;

    # Evaluation for term1 and term2
    for lp in 0:lp_max
        for jp in abs(lp - 1//2):(lp + 1//2)
            # reduced matrix element
            p_partialwave = compute_partial_wave(lp, jp, p_electron, a_electron)
            matrix_elem12 = matrix_element(p_partialwave, a_electron, r)
                
            for q in -1:1

                if matrix_elem12 == 0.0 continue end

                # Term 1 contribution
                factor1 = (-1)^q * u[q] * Ylm(lp, (mj - msp - q), θ, ϕ) *
                            ClebschGordan(lp, (mj - msp - q), 1//2, msp, jp, (mj - q) ) * 
                            ClebschGordan(j, mj, 1, -q, jp, (mj - q) )
                term1 += factor1 * matrix_elem12 

                # Term 2 contribution 
                factor2 = -u[-q] * Ylm(lp, (mj - msp + q), θ, ϕ) *
                            ClebschGordan(lp, (mj - msp + q), 1//2, msp, jp, (mj + q)) * 
                            ClebschGordan(j, mj, 1, q, jp, mj + q)
                term2 += factor2 * matrix_elem12
            end

            # Term 3 A^2 Only compute term3 if lp == l
            if lp == l
                factor3 = Ylm(l, mj - msp, θ, ϕ) *
                        ClebschGordan(l, mj - msp, 1//2, msp, j, mj) *
                        inner_product(p_partialwave, a_electron, r)
                term3 = factor3 * (-im / sqrt(2 * pi)) * StrongFieldDynamics.F2_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ)
            end
        end
    end
    term1 = term1 * (-im * sqrt(2 / pi) ) * StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ ; sign=1)
    term2 = term2 * (-im * sqrt(2 / pi) ) * StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ ; sign=-1)

    # Total result
    return term1 + term2 + term3
end


"""
    T0(pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, mj::Rational{Int64}, msp::Rational{Int64}, θ::Float64, ϕ::Float64)

Computes the direct scattering amplitude 
"""
function T0(pulse::MultiColor, a_electron::AtomicElectron, p_electron::ContinuumElectron, mj::Rational{Int64}, msp::Rational{Int64}, θ::Float64, ϕ::Float64)
    # Maximum number of partial waves
    lp_max = 5

    l = a_electron.l
    j = a_electron.j
    r = a_electron.r
    # u = pulse.u

    # Intitialize term1, term2, term3
    term1 = zero(ComplexF64) ; term2 = zero(ComplexF64) ; term3 = zero(ComplexF64) ;

    # Calculating the pulse shape before hand
    pulse_shape_integral = [ [0.0 + 0.0*im, 0.0 + 0.0*im] , [0.0 + 0.0*im, 0.0 + 0.0*im] ]

    pulse_shape_integral[1][1] = StrongFieldDynamics.F1_integral_levin_approxfun(pulse.colors[1], pulse, a_electron, p_electron, θ, ϕ ; sign=1)
    pulse_shape_integral[1][2] = StrongFieldDynamics.F1_integral_levin_approxfun(pulse.colors[1], pulse, a_electron, p_electron, θ, ϕ ; sign=-1)

    pulse_shape_integral[2][1] = StrongFieldDynamics.F1_integral_levin_approxfun(pulse.colors[2], pulse, a_electron, p_electron, θ, ϕ ; sign=1)
    pulse_shape_integral[2][2] = StrongFieldDynamics.F1_integral_levin_approxfun(pulse.colors[2], pulse, a_electron, p_electron, θ, ϕ ; sign=-1)

    # Evaluation for term1 and term2
    for lp in 0:lp_max
        for jp in abs(lp - 1//2):(lp + 1//2)
            # reduced matrix element
            p_partialwave = compute_partial_wave(lp, jp, p_electron, a_electron)
            matrix_elem12 = matrix_element(p_partialwave, a_electron, r)
            
            if matrix_elem12 == 0.0 continue end

            for (i, color) in enumerate(pulse.colors)
                for q in -1:1
                    
                    # Term 1 contribution
                    factor1 = (-1)^q * color.u[q] * Ylm(lp, (mj - msp - q), θ, ϕ) *
                                ClebschGordan(lp, (mj - msp - q), 1//2, msp, jp, (mj - q) ) * 
                                ClebschGordan(j, mj, 1, -q, jp, (mj - q) )
                    term1 += factor1 * (-im * sqrt(2 / pi) ) * matrix_elem12 * pulse_shape_integral[i][1]

                    # Term 2 contribution 
                    factor2 = -color.u[-q] * Ylm(lp, (mj - msp + q), θ, ϕ) *
                                ClebschGordan(lp, (mj - msp + q), 1//2, msp, jp, (mj + q)) * 
                                ClebschGordan(j, mj, 1, q, jp, mj + q)
                    term2 += factor2 * (-im * sqrt(2 / pi) ) * matrix_elem12 * pulse_shape_integral[i][2]
                end

            end

            # Term 3 A^2 Only compute term3 if lp == l
            if lp == l
                factor3 = Ylm(l, mj - msp, θ, ϕ) *
                        ClebschGordan(l, mj - msp, 1//2, msp, j, mj) *
                        inner_product(p_partialwave, a_electron, r)
                term3 = factor3 * (-im / sqrt(2 * pi)) * StrongFieldDynamics.F2_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ)
            end

        end
    end
    # term1 = term1 * (-im * sqrt(2 / pi) ) # * StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ ; sign=1)
    # term2 = term2 * (-im * sqrt(2 / pi) ) # * StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ ; sign=-1)

    # Total result
    return term1 + term2 + term3
end


"""
    matrix_element(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

Computes the reduced matrix element `<εp lp jp || p || n l j>` for the momentum operator 
between a continuum electron partial wave and an atomic (bound) electron state.

# Arguments
- `p_partialwave::PartialWave`: Continuum electron partial wave with quantum numbers lp, jp
- `a_electron::AtomicElectron`: Atomic (bound) electron with quantum numbers l, j  
- `r::Vector{Float64}`: Radial grid points for numerical integration

# Returns
- `ComplexF64`: The reduced matrix element value

# References
- Equation (16) from PRA 2023 paper for the radial integral formulation
"""
function matrix_element(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

    lp = p_partialwave.l ;    jp = p_partialwave.j
    l  = a_electron.l    ;    j  = a_electron.j

    # Placeholder for the inner matrix element <εp lp || p || n l>
    prefactor = (-im)^(lp+1) * (- sqrt(2*lp + 1) ) * ClebschGordan(lp, 0, 1, 0, l, 0) 

    # To further skip calculations if prefactor is zero
    if prefactor == zero(ComplexF64) return zero(ComplexF64) end

    # Continuum electron partial waves radial part
    cP = Dierckx.Spline1D(r, p_partialwave.P)

    # atomic (bound) electron radial part
    aP = Dierckx.Spline1D(r, a_electron.P)

    # Equation (16) PRA 2023
    integral =  quadgk(r -> ( cP(r) / r ) * ( r * Dierckx.derivative(aP, r) - 0.5 * (lp - l) * (lp + l + 1) * aP(r) ), r[1], r[end])[1]
    
    reduced_matrix_element = prefactor * integral
    
    # Initialize the sum
    result_sum = 0.0
    
    for m in -l:l
        for mp in -lp:lp
            for ms in -1//2:1//2
                
                result_sum += ClebschGordan(lp, mp, 1//2, ms, jp, (mp + ms)) * ClebschGordan(l, m, 1//2, ms, j, (m + ms))
                # Debug output for Clebsch-Gordan coefficients (commented)
                # println("ClebschGordan($lp, $mp, 1//2, $ms, $jp, $(mp+ms)) ClebschGordan($l, $m, 1//2, $ms, $j, $(m + ms))")
                # println(ClebschGordan(lp, mp, 1//2, ms, jp, (mp + ms)), " ", ClebschGordan(l, m, 1//2, ms, j, (m + ms)))
            end
        end
    end

    result = result_sum * reduced_matrix_element

    # Debug output for intermediate results (commented)
    # println("matrix inner $(reduced_matrix_element) and result $(result)")
    
    return result
end


"""
    matrix_element_2(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

Computes the reduced matrix element `<εp lp jp || p || n l j>` for the momentum operator 
between a continuum electron partial wave and an atomic (bound) electron state.

# Arguments
- `p_partialwave::PartialWave`: Continuum electron partial wave with quantum numbers lp, jp
- `a_electron::AtomicElectron`: Atomic (bound) electron with quantum numbers l, j  
- `r::Vector{Float64}`: Radial grid points for numerical integration

# Returns
- `ComplexF64`: The reduced matrix element value

# References
- Refer page 169 - 172 Atomic Structure theory by W R Johnson
"""
function matrix_element_2(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

    lp = p_partialwave.l ;    jp = p_partialwave.j
    l  = a_electron.l    ;    j  = a_electron.j

    C_reduced = if l == (lp - 1)
        (-sqrt(lp)) 
    elseif l == (lp + 1)
        sqrt(lp + 1)
    else
        0.0
    end
    # C_reduced = (-1)^(lp) * sqrt((2*lp+1)*(2*l+1)) * wigner3j(lp, 1, l, 0, 0, 0)

    prefactor = im * (p_partialwave.ε - a_electron.ε) * C_reduced

    # Continuum electron partial waves radial part
    cP = Dierckx.Spline1D(r, p_partialwave.P)

    # atomic (bound) electron radial part
    aP = Dierckx.Spline1D(r, a_electron.P)

    # Equation (16) PRA 2023
    integral =  quadgk(r -> ( r * cP(r) * aP(r) ), r[1], r[end])[1]
    
    reduced_matrix_element = prefactor * integral
    
    # Initialize the sum
    result_sum = 0.0
    
    for m in -l:l
        for mp in -lp:lp
            for ms in -1//2:1//2
                
                result_sum += ClebschGordan(lp, mp, 1//2, ms, jp, (mp + ms)) * ClebschGordan(l, m, 1//2, ms, j, (m + ms))
                # Debug output for Clebsch-Gordan coefficients (commented)
                # println("ClebschGordan($lp, $mp, 1//2, $ms, $jp, $(mp+ms)) ClebschGordan($l, $m, 1//2, $ms, $j, $(m + ms))")
                # println(ClebschGordan(lp, mp, 1//2, ms, jp, (mp + ms)), " ", ClebschGordan(l, m, 1//2, ms, j, (m + ms)))
            end
        end
    end

    result = result_sum * reduced_matrix_element

    # Debug output for intermediate results (commented)
    # println("matrix inner $(reduced_matrix_element) and result $(result)")
    
    return result
end


"""
    inner_product(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

Calculates the overlap integral between continuum and bound electron radial wavefunctions.

This function computes the inner product ⟨ψ_continuum|ψ_bound⟩ of the radial parts of
the electron wavefunctions, which appears in the third term of the scattering amplitude
calculation and represents the direct overlap between initial and final states.

# Arguments
- `p_partialwave::PartialWave`: Continuum electron partial wave containing radial wavefunction P
- `a_electron::AtomicElectron`: Atomic electron containing bound state radial wavefunction P  
- `r::Vector{Float64}`: Radial grid points used for numerical integration

# Returns
- `ComplexF64`: The overlap integral multiplied by the phase factor (-i)^l

# Implementation Details
- Uses cubic spline interpolation (Dierckx.jl) for smooth integration
- Applies numerical quadrature (QuadGK.jl) for accurate integral evaluation
- Includes the angular momentum phase factor (-i)^l_p for proper quantum mechanical normalization
"""
function inner_product(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})
    # Interpolate the radial functions using cubic splines for smooth integration
    cP = Dierckx.Spline1D(r, p_partialwave.P)
    aP = Dierckx.Spline1D(r, a_electron.P)

    # Compute the radial overlap integral ∫ P_continuum(r) * P_bound(r) dr
    result = quadgk(r -> cP(r) * aP(r), r[1], r[end])[1]

    # Apply the angular momentum phase factor (-i)^l_p
    result = result * (-im)^(p_partialwave.l)

    return result
end


"""
    probality_uncoupled(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

Calculates the ionization probability
"""
function probability_uncoupled(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, θ::Float64, ϕ::Float64)

    prob = zero(Float64)

    for m = -a_electron.l//1 : a_electron.l//1
        amplitude = StrongFieldDynamics.T0_uncoupled( pulse, a_electron, p_electron, m, θ, ϕ)
        prob += abs2(amplitude)
    end

    prob = prob * p_electron.p / (2.0*a_electron.j + 1)

end


"""
    T0_uncoupled(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, m::Rational{Int64}, θ::Float64, ϕ::Float64)

Computes the direct scattering amplitude 
"""
function T0_uncoupled(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron, m::Rational{Int64}, θ::Float64, ϕ::Float64)
    # Maximum number of partial waves
    lp_max = 3

    l = a_electron.l
    j = a_electron.j
    r = a_electron.r
    u = pulse.u

    # Intitialize term1, term2, term3
    term1 = zero(ComplexF64) ; term2 = zero(ComplexF64) ; term3 = zero(ComplexF64) ;

    # Evaluation for term1 and term2
    for lp in 0:lp_max
            # reduced matrix element
            p_partialwave = StrongFieldDynamics.compute_partial_wave(lp, 0//1, p_electron, a_electron)
            matrix_elem12 = reduced_matrix_element_uncoupled(p_partialwave, a_electron, r)
                
            for q in -1:1

                if matrix_elem12 == 0.0 continue end

                # Term 1 contribution
                factor1 = (-1)^q * u[q] * Ylm(lp, (m - q), θ, ϕ) *
                            ClebschGordan(l, m, 1, -q, lp, (m - q) )
                term1 += factor1 * matrix_elem12 

                # Term 2 contribution 
                factor2 = -u[-q] * Ylm(lp, (m + q), θ, ϕ) *
                            ClebschGordan(l, m, 1, q, lp, (m + q) )
                term2 += factor2 * matrix_elem12
            end

            # Only compute term3 if lp == l
            if lp == l
                factor3 = Ylm(l, m, θ, ϕ) *
                        inner_product(p_partialwave, a_electron, r)
                term3 = factor3 * (-im / sqrt(2 * pi)) * StrongFieldDynamics.F2_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ)
            end
    end
    term1 = term1 * (-im * sqrt(2 * pi) ) * StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ ; sign=1)
    term2 = term2 * (-im * sqrt(2 * pi) ) * StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, θ, ϕ ; sign=-1)

    # Total result
    return term1 + term2 + term3
end


"""
    reduced_matrix_element_uncoupled(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

Computes the reduced matrix element `<εp lp jp || p || n l j>` for the momentum operator 
between a continuum electron partial wave and an atomic (bound) electron state.

# Arguments
- `p_partialwave::PartialWave`: Continuum electron partial wave with quantum numbers lp, jp
- `a_electron::AtomicElectron`: Atomic (bound) electron with quantum numbers l, j  
- `r::Vector{Float64}`: Radial grid points for numerical integration

# Returns
- `ComplexF64`: The reduced matrix element value

# References
- Equation (16) from PRA 2023 paper for the radial integral formulation
"""
function reduced_matrix_element_uncoupled(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

    lp = p_partialwave.l ;    jp = p_partialwave.j
    l  = a_electron.l    ;    j  = a_electron.j

    # Placeholder for the inner matrix element <εp lp || p || n l>
    prefactor = (-im)^(lp+1) * (- sqrt(2*lp + 1) ) * ClebschGordan(lp, 0, 1, 0, l, 0) 

    # To further skip calculations if prefactor is zero
    if prefactor == zero(ComplexF64) return zero(ComplexF64) end

    # Continuum electron partial waves radial part
    cP = Dierckx.Spline1D(r, p_partialwave.P)

    # atomic (bound) electron radial part
    aP = Dierckx.Spline1D(r, a_electron.P)

    # Equation (16) PRA 2023
    integral =  quadgk(r -> ( cP(r) / r ) * ( r * Dierckx.derivative(aP, r) - 0.5 * (lp - l) * (lp + l + 1) * aP(r) ), r[1], r[end])[1]
    
    result = prefactor * integral
    
    return result
end
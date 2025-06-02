using Dierckx
using QuadGK

"""
    probality(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

Calculates the ionization probability
"""
function probability(pulse::Pulse, a_electron::AtomicElectron, p_electron::ContinuumElectron)

    prob = zero(Float64)

    for mj = -a_electron.j : a_electron.j
        for msp = -1//2:1//2
            amplitude = StrongFieldDynamics.T0( pulse,a_electron, p_electron, mj, msp, deg2rad(90), 0.0)
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
    lp_max = 10

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
            r, P, δ = bessel_electron(p_electron.ε, lp, r)
            p_partialwave = StrongFieldDynamics.PartialWave(p_electron.ε, lp, jp, P, δ)
            matrix_elem12 = reduced_matrix_element(p_partialwave, a_electron, r)
                
            for q in -1:1
                if matrix_elem12 == 0.0 continue end

                # Term 1 contribution
                factor1 = (-1)^q * u[q] * Ylm(lp, (mj - msp - q), θ, ϕ) *
                            ClebschGordan(lp, (mj - msp - q), 1//2, msp, jp, (mj - q) ) * 
                            ClebschGordan(j, mj, 1, -q, jp, (mj - q) )
                term1 += factor1 * matrix_elem12 

                # Term 2 contribution 
                factor2 = conj(u[q]) * Ylm(lp, (mj - msp + q), θ, ϕ) *
                            ClebschGordan(lp, (mj - msp + q), 1//2, msp, jp, (mj + q)) * 
                            ClebschGordan(j, mj, 1, q, jp, mj + q)
                term2 += factor2 * matrix_elem12
                println("factor1 $factor1, factor2 $factor2")
            end

            # Only compute term3 if lp == l
            if lp == l
                factor3 = Ylm(l, mj - msp, θ, ϕ) *
                        ClebschGordan(l, mj - msp, 1//2, msp, j, mj) *
                        inner_product(p_partialwave, a_electron, r)
                term3 = factor3 * (-im / sqrt(2 * pi)) * F2_integral(pulse, a_electron, p_electron, θ)
            end
        end
    end
    term1 = term1 * (-im * sqrt(2 / pi) ) * F1_integral(pulse, a_electron, p_electron, θ ; sign=1)
    term2 = term2 * (-im * sqrt(2 / pi) ) * F1_integral(pulse, a_electron, p_electron, θ ; sign=-1)

    println("Term 1 $term1, Term 2 $term2 and Term3 $term3")

    # Total result
    return term1 + term2 + term3
end


"""
    reduced_matrix_element(εp, lp, jp, n, l, j)

Computes the reduced matrix element `<εp lp jp || p || n l j>`
"""
function reduced_matrix_element(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

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
    
    matrix_element_inner = prefactor * integral
    
    # Initialize the sum
    result_sum = 0.0
    
    for m in -l:l
        for mp in -lp:lp
            for ms in -1//2:1//2
                
                result_sum += ClebschGordan(lp, mp, 1//2, ms, jp, (mp + ms)) * ClebschGordan(l, m, 1//2, ms, j, (m + ms))
                # println("ClebschGordan($lp, $mp, 1//2, $ms, $jp, $(mp+ms)) ClebschGordan($l, $m, 1//2, $ms, $j, $(m + ms))")
                # println(ClebschGordan(lp, mp, 1//2, ms, jp, (mp + ms)), " ", ClebschGordan(l, m, 1//2, ms, j, (m + ms)))
            end
        end
    end

    result = result_sum * matrix_element_inner

    # println("matrix inner $(matrix_element_inner) and result $(result)")
    
    return result
end


"""
    inner_product(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})

Calculates the inner product of continuum and bound (atomic) electron radial wavefunction
"""
function inner_product(p_partialwave::PartialWave, a_electron::AtomicElectron, r::Vector{Float64})
    # Interpolated the radial function to Spline for inegration
    cP = Dierckx.Spline1D(r, p_partialwave.P)
    aP = Dierckx.Spline1D(r, a_electron.P)

    result = quadgk(r -> cP(r) * aP(r), r[1], r[end])[1]

    # Multiply im^(-l)
    result = result * (-im)^(p_partialwave.l)

    return result
end
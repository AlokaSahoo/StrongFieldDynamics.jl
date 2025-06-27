using LinearAlgebra
using QuadGK
using ApproxFun


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
    
    # Solve the boundary value problem: A*ψ = F with ψ(a) = 0
    # Using the constraint ψ(a) = 0
    B = Evaluation(space, a)
    psi = [B; A] \ [0; F]
    
    # println("psi(0) $(psi(0))")
    
    # Compute Levin integral
    result = psi(b) * exp(im * gb) - psi(a) * exp(im * ga)
    
    return result

end

# # Old Version
# function levin_integrate_approxfun(f::Function, gp::Function, a::Float64, b::Float64, ga::Float64, gb::Float64)
#     domain = a..b
#     x = Fun(identity, domain)
#     F = Fun(f, x.space)
#     Gp = Fun(gp, x.space)
#     A = Derivative() + im * Gp
#     p = [Evaluation(0); A] \ [0; F]
#     return p(b) * exp(im * gb) - p(a) * exp(im * ga)
# end



"""
    BarycentricInterpolator

A Julia implementation of barycentric interpolation for polynomial interpolation.
"""
struct BarycentricInterpolator
    nodes::Vector{Float64}
    weights::Vector{Float64}
    values::Vector{ComplexF64}
    
    function BarycentricInterpolator(nodes::Vector{Float64})
        n = length(nodes)
        weights = zeros(Float64, n)
        
        # Compute barycentric weights
        for i in 1:n
            weights[i] = 1.0
            for j in 1:n
                if i != j
                    weights[i] /= (nodes[i] - nodes[j])
                end
            end
        end
        
        new(nodes, weights, zeros(ComplexF64, n))
    end
end

"""
    set_values!(interp, values)

Set the function values at the interpolation nodes.
"""
function set_values!(interp::BarycentricInterpolator, values::Vector{<:Number})
    interp.values .= ComplexF64.(values)
end

"""
    evaluate(interp, x)

Evaluate the interpolating polynomial at point x using barycentric formula.
"""
function evaluate(interp::BarycentricInterpolator, x::Float64)
    # Check if x is one of the nodes
    for i in 1:length(interp.nodes)
        if abs(x - interp.nodes[i]) < 1e-14
            return interp.values[i]
        end
    end
    
    # Barycentric interpolation formula
    numerator = 0.0 + 0.0im
    denominator = 0.0
    
    for i in 1:length(interp.nodes)
        term = interp.weights[i] / (x - interp.nodes[i])
        numerator += term * interp.values[i]
        denominator += term
    end
    
    return numerator / denominator
end

"""
    evaluate_vector(interp, x_vec)

Evaluate the interpolating polynomial at multiple points.
"""
function evaluate_vector(interp::BarycentricInterpolator, x_vec::Vector{Float64})
    return [evaluate(interp, x) for x in x_vec]
end

"""
    finite_difference_derivative(interp, nodes, dx=1e-8)

Compute derivatives of Lagrange basis functions using finite differences.
"""
function finite_difference_derivative(interp::BarycentricInterpolator, nodes::Vector{Float64}, dx::Float64=1e-8)
    n = length(nodes)
    L_derivs = zeros(ComplexF64, n, n)
    
    # Identity matrix for basis function values
    I_matrix = Matrix{Float64}(I, n, n)
    
    for i in 1:n
        # Create interpolator for i-th Lagrange basis function
        basis_interp = BarycentricInterpolator(nodes)
        set_values!(basis_interp, I_matrix[i, :])
        
        # Compute derivative using central difference
        nodes_plus = nodes .+ dx
        nodes_minus = nodes .- dx
        
        vals_plus = evaluate_vector(basis_interp, nodes_plus)
        vals_minus = evaluate_vector(basis_interp, nodes_minus)
        
        L_derivs[i, :] = (vals_plus - vals_minus) / (2 * dx)
    end
    
    return L_derivs
end

"""
    chebyshev_nodes(a, b, n)

Generate n Chebyshev nodes in the interval [a, b].
"""
function chebyshev_nodes(a::Float64, b::Float64, n::Int)
    k = 0:(n-1)
    nodes = 0.5 * (a + b) .+ 0.5 * (b - a) .* cos.(π .* (k .+ 0.5) ./ n)
    return sort(nodes)
end

"""
    levin_integrate(f, g, gp, a, b, n=20)

Compute the highly oscillatory integral ∫[a,b] f(x) exp(ig(x)) dx using Levin's method.

Arguments:
- f: amplitude function
- g: phase function  
- gp: derivative of phase function
- a, b: integration bounds
- n: number of collocation points (default: 20)

Returns:
- Levin integral approximation
"""
function levin_integrate(f::Function, g::Function, gp::Function, a::Float64, b::Float64, n::Int)
    # Generate Chebyshev nodes
    nodes = chebyshev_nodes(a, b, n)
    
    # Create barycentric interpolator and compute derivatives
    interp = BarycentricInterpolator(nodes)
    L_derivs = finite_difference_derivative(interp, nodes)
    
    # Assemble the linear system: A * c = f(nodes)
    A = zeros(ComplexF64, n, n)
    rhs = zeros(ComplexF64, n)
    
    # Identity matrix for basis function values at nodes
    I_matrix = Matrix{Float64}(I, n, n)
    
    for k in 1:n
        for j in 1:n
            A[k, j] = L_derivs[j, k] + 1im * gp(nodes[k]) * I_matrix[j, k]
        end
        rhs[k] = f(nodes[k])
    end
    
    # Solve for coefficients
    c = A \ rhs
    
    # Reconstruct ψ(x) and evaluate at boundaries
    psi_interp = BarycentricInterpolator(nodes)
    set_values!(psi_interp, c)
    
    psi_a = evaluate(psi_interp, a)
    psi_b = evaluate(psi_interp, b)
    
    # Compute Levin integral: ψ(b)e^{ig(b)} - ψ(a)e^{ig(a)}
    I_levin = psi_b * exp(1im * g(b)) - psi_a * exp(1im * g(a))
    
    return I_levin
end

"""
    levin_integrate(f, g, gp, a, b, n=20)

Compute the highly oscillatory integral ∫[a,b] f(x) exp(ig(x)) dx using Levin's method.

Arguments:
- f: amplitude function
- g: phase function  
- gp: derivative of phase function
- ga, gb: integration bounds
- n: number of collocation points (default: 20)

Returns:
- ComplexF64, Levin integral approximation
"""
function levin_integrate(f::Function, gp::Function, a::Float64, b::Float64, ga::Float64, gb::Float64, n::Int)
    # Generate Chebyshev nodes
    nodes = chebyshev_nodes(a, b, n)
    
    # Create barycentric interpolator and compute derivatives
    interp = BarycentricInterpolator(nodes)
    L_derivs = finite_difference_derivative(interp, nodes)
    
    # Assemble the linear system: A * c = f(nodes)
    A = zeros(ComplexF64, n, n)
    rhs = zeros(ComplexF64, n)
    
    # Identity matrix for basis function values at nodes
    I_matrix = Matrix{Float64}(I, n, n)
    
    for k in 1:n
        for j in 1:n
            A[k, j] = L_derivs[j, k] + 1im * gp(nodes[k]) * I_matrix[j, k]
        end
        rhs[k] = f(nodes[k])
    end
    
    # Solve for coefficients
    c = A \ rhs
    
    # Reconstruct ψ(x) and evaluate at boundaries
    psi_interp = BarycentricInterpolator(nodes)
    set_values!(psi_interp, c)
    
    psi_a = evaluate(psi_interp, a)
    psi_b = evaluate(psi_interp, b)
    
    # Compute Levin integral: ψ(b)e^{ig(b)} - ψ(a)e^{ig(a)}
    I_levin = psi_b * exp(1im * gb) - psi_a * exp(1im * ga)
    
    return I_levin
end


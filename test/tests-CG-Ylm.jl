using StrongFieldDynamics
using Test
using LinearAlgebra

@testset "Clebsch-Gordan Coefficients" begin
    
    @testset "Basic Properties and Selection Rules" begin
        # Test triangle inequality: |j1 - j2| ≤ J ≤ j1 + j2
        @test ClebschGordan(1, 0, 1, 0, 3, 0) ≈ 0.0  # J > j1 + j2
        @test ClebschGordan(2, 1, 1, 0, 0, 1) ≈ 0.0  # J < |j1 - j2|
        
        # Test magnetic quantum number conservation: m1 + m2 = M
        @test ClebschGordan(1, 1, 1, 0, 2, 0) ≈ 0.0  # m1 + m2 ≠ M
        @test ClebschGordan(1, 1, 1, -1, 2, 1) ≈ 0.0  # m1 + m2 ≠ M
        
        # Test when |m| > j
        @test ClebschGordan(1, 2, 1, 0, 2, 2) ≈ 0.0  # |m1| > j1
        @test ClebschGordan(1, 0, 1, 2, 2, 2) ≈ 0.0  # |m2| > j2
        @test ClebschGordan(1, 1, 1, 1, 2, 3) ≈ 0.0  # |M| > J
    end
    
    @testset "Known Exact Values" begin
        # Simple cases with j = 1/2
        @test ClebschGordan(1//2, 1//2, 1//2, 1//2, 1, 1) ≈ 1.0
        @test ClebschGordan(1//2, 1//2, 1//2, -1//2, 1, 0) ≈ 1/√2
        @test ClebschGordan(1//2, -1//2, 1//2, 1//2, 1, 0) ≈ 1/√2
        @test ClebschGordan(1//2, -1//2, 1//2, -1//2, 1, -1) ≈ 1.0
        
        @test ClebschGordan(1//2, 1//2, 1//2, -1//2, 0, 0) ≈ 1/√2
        @test ClebschGordan(1//2, -1//2, 1//2, 1//2, 0, 0) ≈ -1/√2
        
        # Cases with j = 1
        @test ClebschGordan(1, 1, 1, 1, 2, 2) ≈ 1.0
        @test ClebschGordan(1, 1, 1, 0, 2, 1) ≈ 1/√2
        @test ClebschGordan(1, 0, 1, 1, 2, 1) ≈ 1/√2
        @test ClebschGordan(1, 1, 1, -1, 2, 0) ≈ 1/√6
        @test ClebschGordan(1, 0, 1, 0, 2, 0) ≈ sqrt(2/3)
        @test ClebschGordan(1, -1, 1, 1, 2, 0) ≈ 1/√6
        
        # Test J = 0 cases (singlet states)
        @test ClebschGordan(1, 1, 1, -1, 0, 0) ≈ 1/√3
        @test ClebschGordan(1, 0, 1, 0, 0, 0) ≈ -1/√3
        @test ClebschGordan(1, -1, 1, 1, 0, 0) ≈ 1/√3
    end
    
    @testset "Symmetry Relations" begin
        # Test symmetry: C(j1,m1,j2,m2;J,M) = (-1)^(j1+j2-J) C(j2,m2,j1,m1;J,M)
        j1, m1, j2, m2, J, M = 1, 1, 1//2, 1//2, 3//2, 3//2
        cg1 = ClebschGordan(j1, m1, j2, m2, J, M)
        cg2 = ClebschGordan(j2, m2, j1, m1, J, M)
        phase = (-1)^(j1 + j2 - J)
        @test cg1 ≈ phase * cg2
        
        # Test another symmetry case
        j1, m1, j2, m2, J, M = 1, 0, 1, 1, 1, 1
        cg1 = ClebschGordan(j1, m1, j2, m2, J, M)
        cg2 = ClebschGordan(j2, m2, j1, m1, J, M)
        phase = (-1)^(j1 + j2 - J)
        @test cg1 ≈ phase * cg2
    end
    
    @testset "Orthogonality Relations" begin
        # Test orthogonality over J,M: ∑_{J,M} C(j1,m1,j2,m2;J,M) * C(j1,m1',j2,m2';J,M) = δ_{m1,m1'} δ_{m2,m2'}
        j1, j2 = 1//2, 1//2
        m1, m2 = 1//2, -1//2
        m1_prime, m2_prime = -1//2, 1//2
        
        sum_orthog = 0.0
        # Sum over all allowed J and M values
        for J in abs(j1-j2):(j1+j2)
            for M in -J:J
                cg1 = ClebschGordan(j1, m1, j2, m2, J, M)
                cg2 = ClebschGordan(j1, m1_prime, j2, m2_prime, J, M)
                sum_orthog += cg1 * cg2
            end
        end
        @test abs(sum_orthog) < 1e-12  # Should be zero since (m1,m2) ≠ (m1',m2')
        
        # Test normalization for same quantum numbers
        m1, m2 = 1//2, -1//2
        sum_norm = 0.0
        for J in abs(j1-j2):(j1+j2)
            for M in -J:J
                cg = ClebschGordan(j1, m1, j2, m2, J, M)
                sum_norm += cg^2
            end
        end
        @test sum_norm ≈ 1.0 atol=1e-12  # Should be normalized
        
        # Test completeness relation: ∑_{j1,m1,j2,m2} C(j1,m1,j2,m2;J,M) * C(j1,m1,j2,m2;J',M') = δ_{J,J'} δ_{M,M'}
        # This is more complex and computationally intensive, so we test a simple case
        J, M = 1, 0
        J_prime, M_prime = 1, 1
        j1, j2 = 1//2, 1//2
        
        sum_completeness = 0.0
        for m1 in [-j1, j1], m2 in [-j2, j2]
            cg1 = ClebschGordan(j1, m1, j2, m2, J, M)
            cg2 = ClebschGordan(j1, m1, j2, m2, J_prime, M_prime)
            sum_completeness += cg1 * cg2
        end
        @test abs(sum_completeness) < 1e-12  # Should be zero since (J,M) ≠ (J',M')
    end
    
    @testset "Error Handling" begin
        # Test that function handles invalid inputs gracefully
        @test ClebschGordan(1, 2, 1, 0, 2, 2) ≈ 0.0  # Invalid m1
        @test ClebschGordan(-1, 0, 1, 0, 2, 0) ≈ 0.0  # Negative j1
    end
end

@testset "Spherical Harmonics (Ylm)" begin
    
    @testset "Basic Properties and Edge Cases" begin
        # Test Y_0^0 (constant function)
        @test Ylm(0, 0, 0.0, 0.0) ≈ 1/√(4π)
        @test Ylm(0, 0, π/2, π/4) ≈ 1/√(4π)
        @test Ylm(0, 0, π, 2π) ≈ 1/√(4π)
        
        # Test invalid m values
        @test Ylm(1, 2, π/2, 0.0) ≈ 0.0  # |m| > l
        @test Ylm(2, 3, π/4, π/2) ≈ 0.0  # |m| > l
        @test Ylm(1, -2, 0.0, 0.0) ≈ 0.0  # |m| > l
    end
    
    @testset "Known Exact Values" begin
        # l = 1 spherical harmonics
        # Y_1^{-1} = √(3/8π) * sin(θ) * e^{-iφ}
        θ, φ = π/2, 0.0  # θ = π/2 → sin(θ) = 1, φ = 0 → e^{-iφ} = 1
        @test real(Ylm(1, -1, θ, φ)) ≈ √(3/(8π))
        @test imag(Ylm(1, -1, θ, φ)) ≈ 0.0
        
        # Y_1^0 = √(3/4π) * cos(θ)
        θ = 0.0  # cos(0) = 1
        @test real(Ylm(1, 0, θ, 0.0)) ≈ √(3/(4π))
        @test imag(Ylm(1, 0, θ, 0.0)) ≈ 0.0
        
        θ = π  # cos(π) = -1
        @test real(Ylm(1, 0, θ, 0.0)) ≈ -√(3/(4π))
        @test imag(Ylm(1, 0, θ, 0.0)) ≈ 0.0
        
        # Y_1^1 = -√(3/8π) * sin(θ) * e^{iφ}
        θ, φ = π/2, 0.0  # sin(π/2) = 1, e^{i*0} = 1
        @test real(Ylm(1, 1, θ, φ)) ≈ -√(3/(8π))
        @test imag(Ylm(1, 1, θ, φ)) ≈ 0.0
        
        # Test with φ = π/2
        θ, φ = π/2, π/2  # e^{i*π/2} = i
        @test real(Ylm(1, 1, θ, φ)) ≈ 0.0 atol=1e-15
        @test imag(Ylm(1, 1, θ, φ)) ≈ -√(3/(8π))
    end
    
    @testset "Orthogonality Relations" begin
        # Test orthogonality: ∫ Y_l^m(θ,φ)* Y_l'^m'(θ,φ) dΩ = δ_{l,l'} δ_{m,m'}
        # Use simple numerical integration for a few test cases
        
        function numerical_orthogonality_test(l1, m1, l2, m2, n_θ=20, n_φ=40)
            # Gauss-Legendre quadrature would be better, but this is sufficient for testing
            θ_points = range(0, π, length=n_θ)
            φ_points = range(0, 2π, length=n_φ)
            
            integral = 0.0 + 0.0im
            for θ in θ_points
                for φ in φ_points
                    y1 = Ylm(l1, m1, θ, φ)
                    y2 = Ylm(l2, m2, θ, φ)
                    # dΩ = sin(θ) dθ dφ, with weights for numerical integration
                    weight = sin(θ) * (π/n_θ) * (2π/n_φ)
                    integral += conj(y1) * y2 * weight
                end
            end
            return integral
        end
        
        # Test orthogonality for different (l,m) pairs
        @test abs(numerical_orthogonality_test(1, 0, 1, 1)) < 0.1  # Should be ≈ 0
        @test abs(numerical_orthogonality_test(1, -1, 1, 1)) < 0.1  # Should be ≈ 0
        @test abs(numerical_orthogonality_test(1, 0, 2, 0)) < 0.1  # Should be ≈ 0
        
        # Test normalization (same l,m)
        @test abs(real(numerical_orthogonality_test(0, 0, 0, 0)) - 1.0) < 0.1  # Should be ≈ 1
        @test abs(real(numerical_orthogonality_test(1, 0, 1, 0)) - 1.0) < 0.1  # Should be ≈ 1
        @test abs(real(numerical_orthogonality_test(1, 1, 1, 1)) - 1.0) < 0.1  # Should be ≈ 1
    end
    
    @testset "Complex Conjugate Relations" begin
        # Test: Y_l^{-m} = (-1)^m * (Y_l^m)*
        
        l, m = 2, 1
        θ, φ = π/3, π/4
        
        y_pos_m = Ylm(l, m, θ, φ)
        y_neg_m = Ylm(l, -m, θ, φ)
        expected = (-1)^m * conj(y_pos_m)
        
        @test real(y_neg_m) ≈ real(expected) atol=1e-12
        @test imag(y_neg_m) ≈ imag(expected) atol=1e-12
        
        # Test another case
        l, m = 3, 2
        θ, φ = π/6, 3π/4
        
        y_pos_m = Ylm(l, m, θ, φ)
        y_neg_m = Ylm(l, -m, θ, φ)
        expected = (-1)^m * conj(y_pos_m)
        
        @test real(y_neg_m) ≈ real(expected) atol=1e-12
        @test imag(y_neg_m) ≈ imag(expected) atol=1e-12
    end
    
    @testset "Symmetry Properties" begin
        # Test behavior at special angles
        
        # At θ = 0 (north pole), only m = 0 terms survive
        θ = 0.0
        φ = π/3  # φ should not matter at poles
        
        @test abs(Ylm(1, 1, θ, φ)) < 1e-12  # Should be ≈ 0
        @test abs(Ylm(1, -1, θ, φ)) < 1e-12  # Should be ≈ 0
        @test abs(Ylm(2, 1, θ, φ)) < 1e-12  # Should be ≈ 0
        @test abs(Ylm(2, -2, θ, φ)) < 1e-12  # Should be ≈ 0
        
        # At θ = π (south pole), only m = 0 terms survive
        θ = π
        @test abs(Ylm(1, 1, θ, φ)) < 1e-12  # Should be ≈ 0
        @test abs(Ylm(1, -1, θ, φ)) < 1e-12  # Should be ≈ 0
        @test abs(Ylm(3, 2, θ, φ)) < 1e-12  # Should be ≈ 0
    end
    
    @testset "Numerical Stability" begin
        # Test at boundary values and potential numerical trouble spots
        
        # Very small angles
        θ, φ = 1e-10, 0.0
        @test isfinite(Ylm(2, 1, θ, φ))
        @test isfinite(Ylm(3, -2, θ, φ))
        
        # Angles very close to π
        θ, φ = π - 1e-10, π
        @test isfinite(Ylm(2, 1, θ, φ))
        @test isfinite(Ylm(4, 3, θ, φ))
        
        # Large l values
        θ, φ = π/2, π/4
        @test isfinite(Ylm(10, 5, θ, φ))
        @test isfinite(Ylm(15, -7, θ, φ))
    end
    
    @testset "Consistency with Real Spherical Harmonics" begin
        # For m = 0, spherical harmonics should be real
        l = 3
        m = 0
        θ, φ = π/3, π/5
        
        y = Ylm(l, m, θ, φ)
        @test abs(imag(y)) < 1e-12
        
        # Test another case
        l = 5
        m = 0
        θ, φ = 2π/3, 4π/7
        
        y = Ylm(l, m, θ, φ)
        @test abs(imag(y)) < 1e-12
    end
end

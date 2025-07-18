using StrongFieldDynamics
using OffsetArrays
using Test

@testset "StrongFieldDynamics.jl" begin
    
    # Include specialized test files
    @testset "Clebsch-Gordan and Spherical Harmonics" begin
        include("tests-CG-Ylm.jl")
    end
    
    # Add other test categories here as the package grows
    @testset "Core Types and Structures" begin
        # TODO: Add tests for Pulse, AtomicElectron, ContinuumElectron, etc.

        @testset "Pulse" begin
            pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=0.5)
            
            @test pulse.I₀ ≈ 0.0014247266774422734
            @test pulse.A₀ ≈ 0.2838197988315132
            @test pulse.ω  ≈ 0.05695419066173491
            @test pulse.λ  ≈ 15117.809007262262
            @test pulse.Tp ≈ 220.63996465148577
            @test pulse.u  ≈ OffsetVector([0.31622776601683794, 0.0, -0.9486832980505138], -1:1)
        end

    end
    
    @testset "Physical Calculations" begin
        # TODO: Add tests for SFA calculations, photoionization, etc.

        settings = Settings()

        a_electron = StrongFieldDynamics.compute_atomic_electron(36, settings.ionization_scheme) ;
        p_electron = StrongFieldDynamics.ContinuumElectron(0.5, sqrt(2*0.5), settings.continuum_solution) ;

        @testset "Volkov Phase Circular" begin
            @testset "Circular" begin
                pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)

                @test StrongFieldDynamics.sin2Sv(pulse.Tp-1, float(pi/4), float(pi/4), pulse, p_electron) ≈ 111.4862413771776
                @test StrongFieldDynamics.sin2Sv_general(pulse.Tp-1, float(pi/4), float(pi/4), pulse, p_electron) ≈ 111.4862413771776
            end

            @testset "Elliptical" begin
                pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=0.5)

                @test StrongFieldDynamics.sin2Sv_general(pulse.Tp-1, float(pi/4), float(pi/4), pulse, p_electron) ≈ 111.48624327809965
            end
        end

        @testset "Pulse Shape Integration" begin

            pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)
            @test isapprox(StrongFieldDynamics.F1_integral_levin_approxfun(pulse, a_electron, p_electron, float(pi)/2, float(0), sign=1), -8.389285477361603e-5 + 3.036140632685656e-5im, rtol = 1e-5)
            
        end

        @testset "Distributions" begin
            pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)
            @test isapprox(compute_energy_distribution(36, pulse; settings = settings, energy_range=(5*pulse.ω, 5*pulse.ω), n_points=1).distribution[1], 7.345370827208923e-7, rtol=1e-5 )
        end

    end
    
    @testset "Utility Functions" begin
        # TODO: Add tests for helper functions, conversions, etc.
    end
end

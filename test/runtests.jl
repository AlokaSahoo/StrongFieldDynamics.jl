using StrongFieldDynamics
using Test

@testset "StrongFieldDynamics.jl" begin
    
    # Include specialized test files
    @testset "Clebsch-Gordan and Spherical Harmonics" begin
        include("tests-CG-Ylm.jl")
    end
    
    # Add other test categories here as the package grows
    @testset "Core Types and Structures" begin
        # TODO: Add tests for Pulse, AtomicElectron, ContinuumElectron, etc.
    end
    
    @testset "Physical Calculations" begin
        # TODO: Add tests for SFA calculations, photoionization, etc.
    end
    
    @testset "Utility Functions" begin
        # TODO: Add tests for helper functions, conversions, etc.
    end
end

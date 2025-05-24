module StrongFieldDynamics

# using Base.Threads: @threads
using PhysicalConstants
# using Distributed: pmap
# using JAC
using LinearAlgebra
using QuasiArrays

include("types-fucntions.jl")

include("amplitudes.jl")
include("electron-wavefunction.jl")
include("laser-pulse.jl")
include("plots.jl")
include("pulse-shape-integral.jl")
include("units.jl")

end
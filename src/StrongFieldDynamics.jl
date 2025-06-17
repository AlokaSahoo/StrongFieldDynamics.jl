module StrongFieldDynamics

using Base.Threads: @threads
using Distributed: pmap
# using JenaAtomicCalculator
using LinearAlgebra
using QuasiArrays

include("types-fucntions.jl")

include("amplitudes.jl")
include("electron-wavefunction.jl")
include("plots.jl")
include("pulse-shape-integral.jl")
include("units.jl")

end
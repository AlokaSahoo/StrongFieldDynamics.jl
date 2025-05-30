# using Revise
using StrongFieldDynamics

using PhysicalConstants.CODATA2022: α, a_0
using QuasiArrays

const SOL = 1/α
const I_AU = 3509.44758
const lamda_AU = 18.89725988 # cm^-1
const Ip_AU = 0.036749405469679 
const T_AU = 2.418884326509 * 10e-17


function fn_pulse()
    # Parameters
    λ = 800 #nm
    I = 5.0    # 10^13 W/cm^2
    ncycles = 2

    # IP = 14 #eV

    λ *= lamda_AU
    ω = 2 * pi * SOL / λ
    momentum = sqrt(2*5*ω) # Why 5 
    I /= I_AU
    # IP *= Ip_AU
    E_0 = sqrt(8 * pi * I / SOL) / ω
    lambda_val = 1.0
    cep = 0.0
    T  = 2pi / ω
    Tp = ncycles * 2pi / ω      # ncycles * T
    ϕ = pi

    # u0 = 0
    up(ϵ, pol) = -1/sqrt(2*(1 + ϵ^2)) * (1 + pol * ϵ)     # u+ = − 1/√(2(1 + ε²)) * (1 + λε)
    um(ϵ, pol) =  1/sqrt(2*(1 + ϵ^2)) * (1 - pol * ϵ)     # u− = 1/√(2(1 + ε²)) * (1 − λε)
    u0 = 0
    ϵ = 1
    pol = 1
    u = QuasiVector([um(ϵ, pol), u0, up(ϵ, pol)], -1:1)
    A₀ = sqrt( 8pi * I / ω^2 / SOL )
    Up = A₀^2 / 4.0 
    f(t) = sin(ω*t/2/ncycles)^2

    A(t) = A₀ * f(t) * exp(-im * (ω*t + ϕ))     # u will be multiplied later

    return StrongFieldDynamics.Pulse(I, A₀, λ, ω, ncycles, Tp, Up, f, ϕ, pol, ϵ, u, sin2Sv)
end

pulse = fn_pulse() ;
r, aP = StrongFieldDynamics.atomic_electron(3,0) ;
IP = 14/27.21138 #eV
nP(r, IP) = 2^2.5 * IP^1.5 * r * exp(-sqrt(2*IP)*r)
aP = nP.(r, IP)
a_electron = StrongFieldDynamics.AtomicElectron(-IP, 3, 0, 1//2, r, aP) ;

εₚ = 5.0
p = sqrt(2*εₚ)
p_electron = StrongFieldDynamics.ContinuumElectron(εₚ, p, :bessel) ;
probability = StrongFieldDynamics.probability(pulse, a_electron, p_electron)

sin2Sv(0.5, 1.5, pulse, a_electron, p_electron)

StrongFieldDynamics.T0( pulse,a_electron, p_electron, 1//2, 1//2, deg2rad(90), 0.0)
StrongFieldDynamics.T0( pulse,a_electron, p_electron, 1//2, -1//2, deg2rad(90), 0.0)
StrongFieldDynamics.T0( pulse,a_electron, p_electron, -1//2, 1//2, deg2rad(90), 0.0)
StrongFieldDynamics.T0( pulse,a_electron, p_electron, -1//2, -1//2, deg2rad(90), 0.0)

energies = 0.01:0.02:10.0

p_electrons = ContinuumElectron[]

for ep in energies
    p_electron = StrongFieldDynamics.ContinuumElectron(ep, sqrt(2*ep), :bessel)
    push!(p_electrons, p_electron)
end


probabilities = Float64[]

for p_electron in p_electrons
    probability = StrongFieldDynamics.probability(pulse, a_electron, p_electron)
    push!(probabilities, probability)
end

using CairoMakie
CairoMakie.plot(energies, probabilities)

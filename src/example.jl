using Revise
using StrongFieldDynamics

using PhysicalConstants.CODATA2022: α, a_0
# using UnitfulAtomic
using QuasiArrays

const SOL = 1/α
const I_AU = 3509.44758
const lamda_AU = 18.89725988 # cm^-1
const Ip_AU = 0.036749405469679 
const T_AU = 2.418884326509 * 10e-17


function fn_pulse()
    # Parameters
    λ = 800 #nm
    I = 90.0    # 10^13 W/cm^2
    ncycles = 8

    IP = 15.7596 #eV

    λ *= lamda_AU
    ω = 2 * pi * SOL / λ
    momentum = sqrt(2*5*ω) # Why 5 
    I /= I_AU
    IP *= Ip_AU
    E_0 = sqrt(8 * pi * I / SOL) / ω
    ϵ = 1.0 # 0.56
    lambda_val = 1.0
    cep = 0.0
    # T  = 
    Tp = 2pi * ncycles / ω
    ϕ = 0.0 

    # u0 = 0
    up(ϵ, pol) = -1/sqrt(2*(1 + ϵ^2)) * (1 + pol * ϵ)     # u+ = − 1/√(2(1 + ε²)) * (1 + λε)
    um(ϵ, pol) =  1/sqrt(2*(1 + ϵ^2)) * (1 - pol * ϵ)     # u− = 1/√(2(1 + ε²)) * (1 − λε)
    u0 = 0
    ϵ = 0.0
    pol = 1
    u = QuasiVector([um(ϵ, pol), 0, up(ϵ, pol)], -1:1)
    A₀ = sqrt( 8pi * I / ω^2 / SOL )
    Up = A₀^2 / 4.0 
    f(t) = sin(ω*t/2/ncycles)

    return StrongFieldDynamics.Pulse(I, A₀, λ, ω, ncycles, Tp, Up, f, ϕ, pol, ϵ, u, sin2Sv)
end

pulse = fn_pulse()
r, aP = StrongFieldDynamics.atomic_electron(3,0)
a_electron = StrongFieldDynamics.AtomicElectron(0.198, 3, 0, 1//2, r, aP)
εₚ = 5.0
p = sqrt(2*εₚ)
p_electron = StrongFieldDynamics.ContinuumElectron(εₚ, p, :bessel)

# sin2Sv(0.5, 1.5, pulse, a_electron, p_electron)

StrongFieldDynamics.T0( pulse,a_electron, p_electron, 1//2, 1//2, 0.0, deg2rad(90))


"""
# Electron Wavefunction Module

This module provides functions for computing atomic and continuum electron wavefunctions
used in strong field dynamics calculations. It supports:

- Atomic bound state wavefunctions for various elements (Li, Ne, Ar, Kr, Xe)
- Continuum electron wavefunctions (free and distorted waves)
- Partial wave computations for photoionization processes

The module interfaces with Fortran libraries for distorted wave calculations and uses
pre-computed atomic data files for accurate bound state representations.
"""

using Dierckx
using SpecialFunctions
using DelimitedFiles

# Path to the Fortran shared object library for distorted partial waves
dir = @__DIR__
freeSchrodinger = dir*"/../deps/mod_sfree.so"


"""
    compute_atomic_electron(Z::Int64, n::Int64, l::Int64; ip::Float64=0.0) → AtomicElectron

Computes the atomic electron wavefunction for specified quantum numbers.

# Arguments
- `Z::Int64`: Atomic number (nuclear charge)
- `n::Int64`: Principal quantum number 
- `l::Int64`: Orbital angular momentum quantum number
- `ip::Float64=0.0`: Ionization potential in atomic units (Hartree). If 0.0, uses default values.

# Returns
- `AtomicElectron`: Structure containing the radial wavefunction and quantum numbers

# Supported Atoms and Orbitals
- **Lithium (Z=3)**: 1s, 2s orbitals
- **Neon (Z=10)**: 1s, 2p orbitals  
- **Argon (Z=18)**: 1s, 3p orbitals
- **Krypton (Z=36)**: 1s, 4p orbitals
- **Xenon (Z=54)**: 1s, 5p orbitals

# Physics
The function loads pre-computed Hartree-Fock or Dirac-Fock wavefunctions from data files.
These represent accurate many-electron atomic wavefunctions that account for electron 
correlation and relativistic effects where applicable.

# Examples
```julia
# Lithium 2s electron
li_2s = compute_atomic_electron(3, 2, 0)

# Neon 2p electron with custom ionization potential
ne_2p = compute_atomic_electron(10, 2, 1; ip=0.8)  # 0.8 Hartree

# Argon 3p electron
ar_3p = compute_atomic_electron(18, 3, 1)
```

# Notes
- Ionization potentials are converted from eV to atomic units (÷27.21138)
- Data files contain radial positions and corresponding wavefunction values
- Interpolation is performed to ensure smooth wavefunction representation
"""
function compute_atomic_electron(Z::Int64, scheme::IonizationScheme; ip::Float64=0.0)

    hydrogenic(r, IP) = 2^2.5 * IP^1.5 * r * exp(-sqrt(2*IP)*r)

    if     (Z == 1)         # Hydrogen 1s ionization
        n = 1 ; l = 0 ;
        j = 1//2    ;  
        if ip == 0.0  ip = 0.5 end  # Default IP: 13.6 eV

        data = readdlm(dir * "/../deps/H-I.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,2]  # r in bohr, P(r) radial wavefunction

    elseif     (Z == 3)         # Lithium 2s ionization
        n = 2 ; l = 0 ;
        j = 1//2    ;  
        if ip == 0.0  ip = 5.3917/27.21138 end  # Default IP: 5.3917 eV

        data = readdlm(dir * "/../deps/Li-I.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,4]  # r in bohr, P(r) radial wavefunction

    elseif (Z == 10)        # Neon 2p+ ionization
        n = 2 ; l = 1 ;
        j = 3//2    ;  
        if ip == 0.0  ip = 21.5645/27.21138 end  # Default IP: 21.5645 eV

        data = readdlm(dir * "/../deps/Ne-I.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,8]

    elseif (Z == 18)        # Argon 3p+ ionization
        n = 3 ; l = 1 ;
        j = 3//2    ;  
        if ip == 0.0  ip = 15.7596/27.21138 end  # Default IP: 15.7596 eV

        data = readdlm(dir * "/../deps/Ar-I.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,14]

    elseif (Z == 36)        # Krypton 4p+ ionization
        n = 4 ; l = 1 ;
        j = 3//2    ;  
        if ip == 0.0  ip = 13.9996/27.21138 end  # Default IP: 13.9996 eV

        data = readdlm(dir * "/../deps/Kr-I.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,24]

    elseif (Z == 54)        # Xenon 5p+ ionization
        n = 5 ; l = 1 ;
        j = 3//2    ;  
        if ip == 0.0  ip = 12.1298/27.21138 end  # Default IP: 12.1298 eV

        data = readdlm(dir * "/../deps/Xe-I.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,34]

    else
        error("Bound-state electron wavefunction not found for Z=$Z, n=$n, l=$l. " *
              "Supported combinations: Li(Z=3): 1s,2s; Ne(Z=10): 1s,2p; " *
              "Ar(Z=18): 1s,3p; Kr(Z=36): 1s,4p; Xe(Z=54): 1s,5p")
    end

    # Interpolate the wavefunction data for smooth representation
    if scheme == Atomic
        itp = Dierckx.Spline1D(r_, P_)
        P = itp.(r_)

    elseif scheme == Hydrogenic
        n = 1 ; l = 0 ;
        j  = 1//2 ;
        P = hydrogenic.(r_, ip)
    end

    return AtomicElectron(Z, n, l, j, (-ip), r_, P)
end


"""
    compute_potential(a_electron::AtomicElectron) → Vector{Float64}

Computes the effective ionic potential for continuum electron calculations.

# Arguments
- `a_electron::AtomicElectron`: Atomic electron structure containing atomic number and radial grid

# Returns
- `Vector{Float64}`: Effective potential V(r) values on the atomic electron's radial grid (atomic units)

# Physics
The function loads pre-computed effective ionic potentials that represent the interaction 
between a continuum electron and the remaining ion after photoionization. 
The effective potential has the asymptotic behavior:
```
V(r) → -Z_eff/r    as r → ∞
```
where Z_eff is the effective ionic charge seen by the continuum electron.

# Supported Atoms
- **Lithium (Z=3)**: Li⁺ ion potential for Li → Li⁺ + e⁻ ionization
- **Neon (Z=10)**: Ne⁺ ion potential for Ne → Ne⁺ + e⁻ ionization
- **Argon (Z=18)**: Ar⁺ ion potential for Ar → Ar⁺ + e⁻ ionization
- **Krypton (Z=36)**: Kr⁺ ion potential for Kr → Kr⁺ + e⁻ ionization
- **Xenon (Z=54)**: Xe⁺ ion potential for Xe → Xe⁺ + e⁻ ionization

# Usage
```julia
# Compute atomic electron wavefunction
atom = compute_atomic_electron(18, Atomic)  # Ar 3p electron

# Get the corresponding ionic potential
V_ion = compute_potential(atom)

# Use for distorted wave calculations
energy = 1.0  # 1 Hartree
l = 1         # p-wave
r, P, δ = distorted_electron(energy, l, atom.r, V_ion)
```

# Data Source
Potentials are loaded from pre-computed data files:
- `Li-II-rV.dat`: Li⁺ ionic potential
- `Ne-II-rV.dat`: Ne⁺ ionic potential  
- `Ar-II-rV.dat`: Ar⁺ ionic potential
- `Kr-II-rV.dat`: Kr⁺ ionic potential
- `Xe-II-rV.dat`: Xe⁺ ionic potential

# Notes
- The potential is interpolated onto the atomic electron's radial grid using spline interpolation
- Data files contain r*V(r) values which are converted to V(r) by dividing by r
- Units are atomic units (Hartree for energy, Bohr radius for distance)
- The potential includes relativistic effects for heavier atoms

# Errors
Throws an error if the atomic number Z is not supported. Currently supports
Z = 3, 10, 18, 36, 54 corresponding to Li, Ne, Ar, Kr, Xe respectively.

# See Also
- [`compute_atomic_electron`](@ref): Computes the initial bound state
- [`distorted_electron`](@ref): Uses this potential for continuum wavefunctions
- [`compute_partial_wave`](@ref): Higher-level interface for partial wave calculations
"""
function compute_potential(a_electron::AtomicElectron)
    if a_electron.Z == 3
        data = readdlm(dir * "/../deps/Li-II-rV.dat")
        r_, rV_ = data[:,1], data[:,2]
    elseif a_electron.Z == 10
        data = readdlm(dir * "/../deps/Ne-II-rV.dat")
        r_, rV_ = data[:,1], data[:,2]
    elseif a_electron.Z == 18
        data = readdlm(dir * "/../deps/Ar-II-rV.dat")
        r_, rV_ = data[:,1], data[:,2]
    elseif a_electron.Z == 36
        data = readdlm(dir * "/../deps/Kr-II-rV.dat")
        r_, rV_ = data[:,1], data[:,2]
    elseif a_electron.Z == 54
        data = readdlm(dir * "/../deps/Xe-II-rV.dat")
        r_, rV_ = data[:,1], data[:,2]
    else
        error("Supported combinations: Li(Z=3); Ne(Z=10); " *
              "Ar(Z=18); Kr(Z=36); Xe(Z=54)")
    end

    # Interpolating to the radial grid
    itp = Dierckx.Spline1D(r_, rV_)
    rV = itp.(a_electron.r)

    return rV
end


"""
    bessel_electron(ε::Float64, l::Int64, r::Vector{Float64}) → (Vector{Float64}, Vector{Float64}, Float64)

Computes the free electron continuum wavefunction using spherical Bessel functions.

# Arguments
- `ε::Float64`: Kinetic energy of the free electron (in atomic units)
- `l::Int64`: Orbital angular momentum quantum number
- `r::Vector{Float64}`: Radial grid points (in bohr)

# Returns
- `r::Vector{Float64}`: Input radial grid
- `P::Vector{Float64}`: Radial wavefunction P(r) = r * jₗ(pr)
- `δ::Float64`: Phase shift (always 0.0 for free electrons)

# Physics
For a free electron in the absence of any potential, the radial wavefunction is:
```
P(r) = r * jₗ(pr)
```
where jₗ is the spherical Bessel function of order l, and p = √(2ε) is the momentum.

This represents the exact solution to the Schrödinger equation with V(r) = 0.

# Example
```julia
r_grid = range(0.1, 50.0, 500)
energy = 1.0  # 1 Hartree kinetic energy
l = 1  # p-wave

r, P, phase = bessel_electron(energy, l, r_grid)
```
"""
function bessel_electron(ε::Float64, l::Int64, r::Vector{Float64})
    # Linear momentum of the electron: p = √(2mε) with m=1 in atomic units
    p = sqrt(2ε)    

    # Radial part of the continuum electron: P(r) = r * jₗ(pr)
    # The factor of r converts from reduced radial function u(r) to P(r)
    P = @. r * SpecialFunctions.sphericalbesselj(l, p * r)

    return r, P, 0.0  # Phase shift is zero for free electrons
end


"""
    distorted_electron(ε::Float64, l::Int64, r::Vector{Float64}, rV::Vector{Float64}) → (Vector{Float64}, Vector{Float64}, Float64)

Computes the distorted wave continuum electron wavefunction in an atomic potential.

# Arguments
- `ε::Float64`: Kinetic energy of the electron (in atomic units)
- `l::Int64`: Orbital angular momentum quantum number  
- `r::Vector{Float64}`: Radial grid points (in bohr)
- `rV::Vector{Float64}`: Potential energy array V(r) at grid points (in atomic units)

# Returns
- `r::Vector{Float64}`: Input radial grid
- `P::Vector{Float64}`: Normalized radial wavefunction P(r)
- `δ::Float64`: Total phase shift (inner + Coulomb contributions)

# Physics
Solves the radial Schrödinger equation with the given potential:
```
d²P/dr² + [2ε - 2V(r) - l(l+1)/r²]P = 0
```

The wavefunction asymptotically behaves as:
```
P(r) → sin(pr - lπ/2 + δₗ) / p    as r → ∞
```

This accounts for scattering from the atomic potential, providing more accurate
wavefunctions for photoionization calculations than free electron approximations.

# Implementation
Uses a Fortran library (mod_sfree.so) for numerical integration of the 
radial Schrödinger equation with proper boundary conditions.

# Example  
```julia
r_grid = range(0.1, 50.0, 500)
V_coulomb = -2.0 ./ r_grid  # Hydrogen-like potential
energy = 0.5  # 0.5 Hartree kinetic energy

r, P, delta = distorted_electron(energy, 0, r_grid, V_coulomb)
println("s-wave phase shift: ", delta, " radians")
```
"""
function distorted_electron(ε::Float64, l::Int64, r::Vector{Float64}, rV::Vector{Float64})

    # Prepare arrays for Fortran interface
    npts=length(r)   # Number of radial points
    r0=zeros(npts+1); r0[2:npts+1]=r      # Fortran 1-based indexing
    rV0=zeros(npts+1); rV0[2:npts+1]=rV   # Potential array with 1-based indexing
    
    # Phase shift references for Fortran output
    iPhase=Ref{Float64}(0.0)  # Inner phase shift
    cPhase=Ref{Float64}(0.0)  # Coulomb phase shift
    
    # Output arrays (oversized for safety)
    r_=zeros(Float64,25000)
    P=zeros(Float64,25000)   # Large component (radial wavefunction)
    Q=zeros(Float64,25000)   # Small component (for relativistic case)

    # Call Fortran subroutine for numerical integration
    ccall((:mysfree, freeSchrodinger), Cvoid, 
            (Ref{Int64},Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int64}, 
             Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}), 
            npts, r0, rV0, ε, l, r_, P, Q, iPhase, cPhase)
    
    println("Inner phase shift = ", iPhase.x, "  Coulomb phase shift = ", cPhase.x)

    # Linear momentum for normalization
    p = sqrt(2ε)

    return r, (P[2:npts+1] ./ p), iPhase.x + cPhase.x  # Return normalized P(r) and total phase
end


"""
    compute_partial_wave(l::Int64, j::Rational{Int64}, p_electron::ContinuumElectron, a_electron::AtomicElectron) → PartialWave

Computes a partial wave for photoionization calculations.

# Arguments
- `l::Int64`: Orbital angular momentum quantum number
- `j::Rational{Int64}`: Total angular momentum quantum number (j = l ± 1/2)
- `p_electron::ContinuumElectron`: Continuum electron specification (energy, solution type)
- `a_electron::AtomicElectron`: Bound atomic electron data

# Returns
- `PartialWave`: Structure containing energy, quantum numbers, wavefunction, and phase shift

# Physics
Computes the continuum electron partial wave corresponding to the photoionization
transition from the bound atomic orbital. The choice of solution method (Bessel vs 
Distorted) affects the accuracy:

- **Bessel**: Free electron approximation, ignores atomic potential
- **Distorted**: Includes scattering from atomic potential (more accurate)

# Example
```julia
# Set up bound and continuum electrons
atom = compute_atomic_electron(18, 3, 1)  # Ar 3p
continuum = ContinuumElectron(1.0, :Distorted)  # 1 Hartree, distorted wave

# Compute p-wave (l=1) partial wave
partial_wave = compute_partial_wave(1, 3//2, continuum, atom)
```

# Notes
- The radial grid from the atomic electron is used for the continuum calculation
- Phase shifts are important for interference effects in photoionization
- Distorted waves require the atomic potential (currently needs rV to be defined)
"""
function compute_partial_wave(l::Int64, j::Rational{Int64}, p_electron::ContinuumElectron, a_electron::AtomicElectron)
    P = zeros(Float64, length(a_electron.r))  ; δ = 0.0

    if p_electron.solution == Bessel
        r, P, δ = bessel_electron(p_electron.ε, l, a_electron.r)
    elseif p_electron.solution == Distorted
        # Note: rV needs to be defined in the calling scope or passed as parameter
        rV = compute_potential(a_electron)
        r, P, δ = distorted_electron(p_electron.ε, l, a_electron.r, rV)
    else
        @error "$(p_electron.solution) - not implemented yet"
    end

    return PartialWave(p_electron.ε, l, j, P, δ)
end

# """
#     continuumElectron(E::Float64, l::Int64, pot::Vector{Float64})

# Generates the continuum (photo) electron radial wavefunctions for energy `E` and orbital angular momentum `l`
# """
# function continuumElectron(E::Float64, l::Int64)

#     α = 0.0072973525643  ;      wc = 1/α

#     if l == 0 κ = -1 else κ = l end

#     # grid = pot.grid
#     r = range(1e-6, 50.0, 500)

#     # Vitp = Dierckx.Spline1D(grid.r, - pot.Zr ./ grid.r)
#     # V(r) = Vitp(r)
#     V(r) = -2/r

#     # Define the system of ODEs
#     function dirac_system!(du, u, p, r)
#         P, Q = u  # u[1] = P(r), u[2] = Q(r)
#         # κ, wc = p

#         # Equations
#         du[1] = -(κ / r) * P + ((E + 2.0 * wc^2 - V(r)) / wc) * Q  # dP/dr
#         du[2] = (κ / r) * Q + ((-E + V(r)) / wc) * P   # dQ/dr
#     end

#     # if ( κ < 0 && abs(κ) > 30 ) κ = ( abs(κ) - 1 ) end

#     # # Define the radial range to solve over
#     # if κ < 30
#     #     rspan = (1e-6, grid.r[grid.NoPoints]) 
#     # elseif κ < 60
#     #     rspan = (1e-2, grid.r[grid.NoPoints])
#     # elseif κ < 80
#     #     rspan = (0.5, grid.r[grid.NoPoints])
#     # elseif κ < 150
#     #     rspan = (2.0, grid.r[grid.NoPoints])
#     # else
#     #     rspan = (5.0, grid.r[grid.NoPoints])
#     # end

#     rspan = (1e-6, 50.0)

#     # Initial conditions: 
#     u0 = computeInitialCondition(rspan[2], E, l)
#     # u0 = big.(computeInitialCond(rspan[2], E, sh, pot))
#     # u0 =[1e-30, 1e-30]

#     # Define the ODE problem
#     prob = OrdinaryDiffEq.ODEProblem(dirac_system!, u0, rspan)

#     # Solve the ODE problem
#     sol = OrdinaryDiffEq.solve(prob, Vern9(), abstol=1e-15, reltol=1e-15);

#     # P = sol.(r)[1,:]         ;    Q = sol(r)[2,:] 
#     # Pprime = sol.(r, Val{1})[1,:]     ;    Qprime = sol(r, Val{1})[2,:] 

#     P       = [sol(r_val)[1] for r_val in r]
#     Pprime  = [sol(r, Val{1})[1] for r_val in r]

#     return( P, Pprime )
# end



# """
#     Continuum.computeInitialCondition(r::Float64, energy::Float64, sh::Subshell, pot::Radial.Potential)

# Assymptotic condition as the Intial condition to solve the Dirac radial euation.
# """
# function computeInitialCondition(r::Float64, energy::Float64, l::Int64)
#     Zbar = 1    #-Radial.determineZbar(pot)  ; 
#     # mtp = size( cOrbital.P, 1)  ;             #energy = cOrbital.energy;      
#     # kappa = sh.kappa  ;                       l = Basics.subshell_l(sh) 
#     # α = Defaults.getDefaults("alpha")  ;      wc = 1/α

#     α = 0.0072973525643  ;      wc = 1/α

#     if l == 0 kappa = -1 else kappa = l end

#     q  = sqrt( energy * (energy + 2 * wc^2) ) / wc  ;    x  = q * r

#     ## println("Normalization with Coulomb functions")
#     # sa = "Normalization with Coulomb functions"
#     λ  = sqrt(kappa^2 - Zbar^2 / wc^2)  ; λm1 = λ - 1.0
#     η  = Zbar * α * (energy + wc^2) / sqrt( energy * (energy + 2 * wc^2) ) 

#     xTP = η + sqrt(η^2 + λ*(λ+1.0))
#     if x < xTP println("The kr is less than the Coulomb turning point") end

#     Δ = angle(SpecialFunctions.gamma(λm1 + 1 + im * η))
#     if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end

#     θ  = x - λm1*pi/2 - η*log(2x) + Δ                               #Eq 7.3
#     if θ > 1e4   θ = mod(θ, 2pi) end

#     hgfλ = twoFzero(im*η - λm1, im*η + λm1 + 1, im*2*x)
#     hgfλm1 = twoFzero(im*η - λm1 + 1, im*η + λm1 + 2, im*2*x)
#     hgfλm1 = im * hgfλm1 * (im*η - λm1) * (im*η + λm1 + 1) / ( 2.0 * x^2 )

#     GiFλm1 = hgfλ * exp(im*θ)  ;   GPiFPλm1 = ( hgfλm1 + im *(1.0 - η/x) * hgfλ ) * exp(im*θ)
#     gm_1 = GiFλm1.re ; fm_1 = GiFλm1.im   ;   gpm_1 = GPiFPλm1.re ; fpm_1 = GPiFPλm1.im

#     if abs(gm_1*fpm_1 - fm_1*gpm_1 - 1.0) > 1e-15 
#         ef= 0.0 ; eg = 0.0
#         f, fp, g, gp = GSL.sf_coulomb_wave_FG_e(η, x, λm1, 0, ef, eg)
#         gm_1 = g.val ; fm_1 = f.val ; gpm_1 = gp.val ; fpm_1 = fp.val
#     end

#     f = λ * ((λ/x + η/λ)*fm_1 - fpm_1) / sqrt(λ^2 + η^2)
#     g = λ * ((λ/x + η/λ)*gm_1 - gpm_1) / sqrt(λ^2 + η^2)

#     N  =  ( α^2 * Zbar^2 * (energy + 2 * wc^2)^2 + (kappa + λ)^2 * wc^2 * q^2 )^(-0.5) / λ

#     # Coulomb phase shift for λ
#     rnu = angle(Zbar * α * (energy + 2 * wc^2) - im * (kappa+λ) * sqrt(energy*(energy+2*wc^2)))
#     Δ =  rnu - (λ-l-1.0)*pi/2 + Δ
#     if ( Zbar < 0.0 && kappa < 0)  Δ = Δ - pi ; N = -N end
#     if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end

#     fu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * f + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
#     gu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * g + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )

#     fl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * f + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
#     gl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * g + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )

#     return [fu, fl]
# end 



# """
# `function twoFzero(CA::ComplexF64, CB::ComplexF64, CZ::ComplexF64)`
#     ... Calculates the Hypergeometric function 2F0(CA,CB;1/CZ) hypergeometric asymptotic series.
#         Taken from Radial package by Salvat et al.
#         A ComplexF64 value is returned.  
# """
# function twoFzero(CA::ComplexF64, CB::ComplexF64, CZ::ComplexF64)

#     EPS=1.0E-16; ACCUR=0.5E-15; NTERM=75

#     RRP=1.0
#     RRN=0.0
#     RIP=0.0
#     RIN=0.0
#     CDF=1.0 + 0.0im
#     ERR2=0.0
#     ERR3=1.0
#     AR=0.0
#     AF=0.0
#     CF=0.0 + 0.0im

#     for I = 1 : NTERM
#         J=I-1
#         CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
#         ERR1=ERR2
#         ERR2=ERR3
#         ERR3=abs(CDF)
#         if (ERR1 > ERR2 && ERR2 < ERR3) break end
#         AR=CDF.re
#         if(AR > 0.0) 
#         RRP=RRP+AR
#         else
#         RRN=RRN+AR
#         end
#         AI=(-im*CDF).re
#         if(AI > 0.0) 
#         RIP=RIP+AI
#         else
#         RIN=RIN+AI
#         end
#         CF=complex(RRP+RRN,RIP+RIN)
#         AF=abs(CF)
#         if(AF > 1.0e25) 
#         CF=0.0 + 0im
#         ERR=1.0
#         break
#         end
#         if(ERR3 < 1.0e-25*AF || ERR3 < EPS) 
#         ERR=EPS
#         break
#         end    
#     end

#     # ****  Round off error.

#     TR=abs(RRP+RRN)
#     if(TR > 1.0e-25) 
#     ERRR=(RRP-RRN)*ACCUR/TR
#     else
#     ERRR=1.0e0
#     end
#     TI=abs(RIP+RIN)
#     if(TI > 1.0e-25) 
#     ERRI=(RIP-RIN)*ACCUR/TI
#     else
#     ERRI=1.0
#     end

#     #  ****  ... and truncation error.
#     if(AF > 1.0e-25) 
#     ERR=max(ERRR,ERRI)+ERR2/AF
#     else
#     ERR=max(ERRR,ERRI)
#     end

#     return CF

# end


# """
#     continuumElectron(E::Float64, l::Int64, pot::Vector{Float64})

# Generates the continuum (photo) electron radial wavefunctions for energy `E` and orbital angular momentum `l`
# """
# function continuumElectron(E::Float64, l::Int64, pot::Vector{Float64})
#     k = sqrt(2E)
#     Z = 1  # Nuclear charge (Hydrogen atom)

#     # Define the potential V(r)
#     V(r) = -Z / r
    
#     # Define the ODE system
#     function radial_se!(ddu, du, u, p, r)
#         E, l = p  # Parameters: Energy E and angular momentum l
#         # du[1] = u[2]  # u[1] = u(r), u[2] = du/dr
#         ddu[1] = ( l * (l + 1) / r^2 - 2 * (E - V(r)) ) * u[1]
#     end

#     # Initial conditions (regular solution near r=0)
#     r0 = 1e-6               # Avoid r=0 singularity
#     u0 = r0^(l + 1)         # u(r) ~ r^{l+1} near r=0
#     v0 = (l + 1) * r0^l     # du/dr ~ (l+1) r^l near r=0

#     # Solve the ODE
#     rspan = (r0, 10.0)  # Radial range
#     params = (E, l)  # Pass E and l as parameters
#     prob = SecondOrderODEProblem(radial_se!, [v0], [u0], rspan, params)
#     sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

#     # Extract the solution
#     r = range(rspan[1], rspan[2], 500)
#     u = [sol(r_val)[1] for r_val in r]
#     R = u ./ r  # Radial wavefunction R(r) = u(r)/r  

#     uprime = [sol(r_val, Val{1})[1] for r_val in r]
#     Rprime = - uprime ./ r^2

#     return R, Rprime
# end



# """
#     generateOrbitalSciML(energy::Float64, sh::Subshell, pot::Radial.Potential)

# Gernerates the continuum elctron radial wavefunctions.
# """
# function generateOrbitalSciML(energy::Float64, sh::Subshell, pot::Radial.Potential)

#     α = Defaults.getDefaults("alpha")  ;      wc = 1/α
#     E = energy  ;   κ = sh.kappa
#     grid = pot.grid

#     Vitp = Dierckx.Spline1D(grid.r, - pot.Zr ./ grid.r)
#     V(r) = Vitp(r)

#     # Define the system of ODEs
#     function dirac_system!(du, u, p, r)
#         P, Q = u  # u[1] = P(r), u[2] = Q(r)
#         # κ, wc = p

#         # Equations
#         du[1] = -(κ / r) * P + ((E + 2.0 * wc^2 - V(r)) / wc) * Q  # dP/dr
#         du[2] = (κ / r) * Q + ((-E + V(r)) / wc) * P   # dQ/dr
#     end

#     if ( κ < 0 && abs(κ) > 30 ) κ = ( abs(κ) - 1 ) end

#     # Define the radial range to solve over
#     if κ < 30
#         rspan = (1e-6, grid.r[grid.NoPoints]) 
#     elseif κ < 60
#         rspan = (1e-2, grid.r[grid.NoPoints])
#     elseif κ < 80
#         rspan = (0.5, grid.r[grid.NoPoints])
#     elseif κ < 150
#         rspan = (2.0, grid.r[grid.NoPoints])
#     else
#         rspan = (5.0, grid.r[grid.NoPoints])
#     end

#     # Initial conditions: 
#     u0 = Continuum.computeInitialCondition(rspan[2], energy, sh, pot)
#     # u0 = big.(computeInitialCond(rspan[2], energy, sh, pot))
#     # u0 =[1e-30, 1e-30]

#     # Define the ODE problem
#     prob = OrdinaryDiffEq.ODEProblem(dirac_system!, u0, rspan)

#     # Solve the ODE problem
#     sol = OrdinaryDiffEq.solve(prob, Vern9(), abstol=1e-15, reltol=1e-15);

#     P = sol(grid.r)[1,:]         ;    Q = sol(grid.r)[2,:] 
#     Pprime = sol(grid.r, Val{1})[1,:]     ;    Qprime = sol(grid.r, Val{1})[2,:] 

#     cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, grid)

#     println("energy $energy kappa $(sh.kappa)")

#     return( cOrbital )
# end


# """
#     Continuum.computeInitialCondition(r::Float64, energy::Float64, sh::Subshell, pot::Radial.Potential)

# Assymptotic condition as the Intial condition to solve the Dirac radial euation.
# """
# function computeInitialCondition(r::Float64, energy::Float64, sh::Subshell, pot::Radial.Potential)
#     Zbar = -Radial.determineZbar(pot)  ; 
#     # mtp = size( cOrbital.P, 1)  ;             #energy = cOrbital.energy;      
#     kappa = sh.kappa  ;                       l = Basics.subshell_l(sh) 
#     α = Defaults.getDefaults("alpha")  ;      wc = 1/α
#     q  = sqrt( energy * (energy + 2 * wc^2) ) / wc  ;    x  = q * r

#     ## println("Normalization with Coulomb functions")
#     # sa = "Normalization with Coulomb functions"
#     λ  = sqrt(kappa^2 - Zbar^2 / wc^2)  ; λm1 = λ - 1.0
#     η  = Zbar * α * (energy + wc^2) / sqrt( energy * (energy + 2 * wc^2) ) 

#     xTP = η + sqrt(η^2 + λ*(λ+1.0))
#     if x < xTP println("The kr is less than the Coulomb turning point") end

#     Δ = angle(SpecialFunctions.gamma(λm1 + 1 + im * η))
#     if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end

#     θ  = x - λm1*pi/2 - η*log(2x) + Δ                               #Eq 7.3
#     if θ > 1e4   θ = mod(θ, 2pi) end

#     hgfλ = twoFzero(im*η - λm1, im*η + λm1 + 1, im*2*x)
#     hgfλm1 = twoFzero(im*η - λm1 + 1, im*η + λm1 + 2, im*2*x)
#     hgfλm1 = im * hgfλm1 * (im*η - λm1) * (im*η + λm1 + 1) / ( 2.0 * x^2 )

#     GiFλm1 = hgfλ * exp(im*θ)  ;   GPiFPλm1 = ( hgfλm1 + im *(1.0 - η/x) * hgfλ ) * exp(im*θ)
#     gm_1 = GiFλm1.re ; fm_1 = GiFλm1.im   ;   gpm_1 = GPiFPλm1.re ; fpm_1 = GPiFPλm1.im

#     if abs(gm_1*fpm_1 - fm_1*gpm_1 - 1.0) > 1e-15 
#         ef= 0.0 ; eg = 0.0
#         f, fp, g, gp = GSL.sf_coulomb_wave_FG_e(η, x, λm1, 0, ef, eg)
#         gm_1 = g.val ; fm_1 = f.val ; gpm_1 = gp.val ; fpm_1 = fp.val
#     end

#     f = λ * ((λ/x + η/λ)*fm_1 - fpm_1) / sqrt(λ^2 + η^2)
#     g = λ * ((λ/x + η/λ)*gm_1 - gpm_1) / sqrt(λ^2 + η^2)

#     N  =  ( α^2 * Zbar^2 * (energy + 2 * wc^2)^2 + (kappa + λ)^2 * wc^2 * q^2 )^(-0.5) / λ

#     # Coulomb phase shift for λ
#     rnu = angle(Zbar * α * (energy + 2 * wc^2) - im * (kappa+λ) * sqrt(energy*(energy+2*wc^2)))
#     Δ =  rnu - (λ-l-1.0)*pi/2 + Δ
#     if ( Zbar < 0.0 && kappa < 0)  Δ = Δ - pi ; N = -N end
#     if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end

#     fu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * f + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
#     gu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * g + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )

#     fl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * f + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
#     gl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * g + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )

#     return [fu, fl]
# end 


#=======================================================================================================
using DifferentialEquations

using Plots
using SpecialFunctions # For loggamma (complex gamma function)
Plots.GRBackend()

# --- 1. Define Potentials ---

# Coulomb Potential (attractive for electron, Z is positive nuclear charge)
function coulomb_potential(r, Z_nucleus)
    if r < 1e-9 # Avoid division by zero, though r_min in solver should prevent this.
        return Z_nucleus > 0 ? -1e9 : 0.0 # Large negative for attractive
    end
    return -Z_nucleus / r
end

# Woods-Saxon Potential (example short-range)
function woods_saxon_potential(r, V0_ws, R0_ws, a_ws)
    if r < 1e-9 && R0_ws < 1e-3 # Only if R0 is also tiny, otherwise r=0 is fine
         return -V0_ws / (1 + exp(-R0_ws / a_ws))
    end
    return -V0_ws / (1 + exp((r - R0_ws) / a_ws))
end

# Combined potential function to be used in the Schrodinger equation
function total_potential(r, Z_nucleus, short_range_func, sr_params)
    V_c = coulomb_potential(r, Z_nucleus)
    V_sr = short_range_func(r, sr_params...)
    return V_c + V_sr
end


# --- 2. Set up the ODE system ---
# d²u/dr² = [l(l+1)/r² + 2(V(r) - E)] u(r)
function schrodinger_radial_rhs!(du, u_vec, p, r)
    # u_vec[1] is u(r)
    # u_vec[2] is u'(r)
    # p: (E, l, Z_nucleus, short_range_func, sr_params)

    E, l, Z_nucleus, short_range_func, sr_params = p

    V_r = total_potential(r, Z_nucleus, short_range_func, sr_params)
    
    centrifugal_term = if r > 1e-9
        l * (l + 1) / (r^2)
    else
        l > 0 ? 1e18 : 0.0 # Should not be hit if r_min > 0
    end
    
    du[1] = u_vec[2] # du/dr
    du[2] = (centrifugal_term + 2.0 * (V_r - E)) * u_vec[1] # d²u/dr²
    
    return nothing
end

# --- Parameters ---
E_continuum = 0.5  # Energy in Hartrees (>0)
l_quantum = 0      # Orbital angular momentum
Z_nucleus = 2.0    # Nuclear charge for Coulomb potential (e.g., 1.0 for Hydrogen-like)
                   # Set Z_nucleus = 0.0 to have only short-range potential

# Woods-Saxon potential parameters
V0_ws = 0.0       # Depth in Hartrees. Negative for attractive. (Original WS is often defined with V0 positive, V = -V0/...)
R0_ws = 0.0        # Radius in Bohr radii
a_ws  = 0.5        # Diffuseness in Bohr radii
sr_potential_params = (V0_ws, R0_ws, a_ws)
chosen_sr_potential_func = woods_saxon_potential

# To test pure Coulomb, set V0_ws = 0.0
# sr_potential_params_nul = (0.0, R0_ws, a_ws)


# Integration range
r_min = 1e-5       # Start integration slightly away from r=0
r_max = 40.0       # Max radius. Increase if wavelength is long or decay of V_sr is slow.
r_span = (r_min, r_max)

# --- 3. Initial Conditions u(r_min) and u'(r_min) ---
# u_l(r) ~ C * r^(l+1) near origin. Set C=1.
u_initial = r_min^(l_quantum + 1)
du_initial = (l_quantum + 1) * r_min^l_quantum
u0 = [u_initial, du_initial]

# --- 4. Solve the ODE ---
params_ode = (E_continuum, l_quantum, Z_nucleus, chosen_sr_potential_func, sr_potential_params)
# For pure Coulomb test:
# params_ode = (E_continuum, l_quantum, Z_nucleus, chosen_sr_potential_func, sr_potential_params_nul)


prob = ODEProblem(schrodinger_radial_rhs!, u0, r_span, params_ode)
# alg = Tsit5()
alg = Vern7() # Higher accuracy solver
sol = solve(prob, alg, reltol=1e-9, abstol=1e-9, dense=true, saveat=0.1)

# --- 5. Plot Results ---
r_plot = sol.t
u_plot = [s[1] for s in sol.u]

plot_u = plot(r_plot, u_plot,
              xlabel="r (Bohr radii)",
              ylabel="u_l(r)",
              title="Reduced Radial Wavefunction u_l(r)\nE=$(E_continuum) H, l=$l_quantum, Z=$Z_nucleus, V0_sr=$(sr_potential_params[1])",
              legend=false,
              linewidth=2)
display(plot_u)

# --- 6. Phase Shift Calculation ---
if length(sol.t) > 10 && E_continuum > 0
    # Choose r_match towards the end of the integration range
    # Ensure it's far enough for V_short_range to be negligible
    # And for asymptotic Coulomb forms to be valid
    num_points_to_skip = min(10, length(sol.t) ÷ 2) # Skip some points from the end to avoid boundary effects if any
    r_match_idx = length(sol.t) - num_points_to_skip 
    r_match = sol.t[r_match_idx]
    
    u_val_match = sol.u[r_match_idx][1]
    du_val_match = sol.u[r_match_idx][2]

    # Check if short-range potential is negligible at r_match
    V_sr_at_match = chosen_sr_potential_func(r_match, sr_potential_params...)
    V_coulomb_at_match = coulomb_potential(r_match, Z_nucleus)
    
    println("\n--- Phase Shift Calculation ---")
    println("Matching at r_match = $r_match")
    println("  V_short_range(r_match) = $V_sr_at_match")
    println("  V_Coulomb(r_match) = $V_coulomb_at_match")
    println("  Relative V_sr to E: $(abs(V_sr_at_match / E_continuum))")

    if abs(V_sr_at_match / E_continuum) < 0.005 || sr_potential_params[1] == 0.0 # Threshold or V0_sr is zero
        k_asym = sqrt(2 * E_continuum) # Asymptotic wave number k = sqrt(2mE/ħ²), m=1, ħ=1
        
        eta = 0.0 # Sommerfeld parameter
        sigma_l = 0.0 # Coulomb phase shift

        if Z_nucleus != 0.0
            eta = -Z_nucleus / k_asym # For attractive potential V = -Z/r
            # sigma_l = arg( gamma(l_quantum + 1 + im*eta) )
            # Using loggamma for better stability: arg(z) = imag(log(z))
            sigma_l = imag(loggamma(complex(l_quantum + 1, eta)))
        end
        
        # Asymptotic form: u(r) ~ C * sin(k*r - η*ln(2kr) - lπ/2 + σ_l + δ_l^short)
        # Derivative: u'(r) ~ C * (k - η/r) * cos(k*r - η*ln(2kr) - lπ/2 + σ_l + δ_l^short)
        
        term_kr = k_asym * r_match
        term_eta_log = (Z_nucleus != 0.0) ? eta * log(2 * k_asym * r_match) : 0.0
        term_l_pi = l_quantum * π / 2.0
        
        # Argument of tan: X = kr - η*ln(2kr) - lπ/2 + σ_l + δ_l^short
        # tan(X) = u(r_match) * (k - η/r_match) / u'(r_match)
        
        # Effective k for derivative part in Coulomb field: k_eff_prime = (k - η/r_match)
        # This is d/dr (kr - η*ln(2kr))
        k_eff_prime = k_asym - ( (Z_nucleus != 0.0) ? eta / r_match : 0.0 )

        # Use atan(y,x) for correct quadrant. tan(X) = Numerator / Denominator
        # Numerator: u_val_match * k_eff_prime
        # Denominator: du_val_match
        X_numeric = atan(u_val_match * k_eff_prime, du_val_match)
        
        # δ_l^short = X_numeric - (kr_match - η*ln(2kr_match) - lπ/2 + σ_l)
        delta_l_short_rad_raw = X_numeric - (term_kr - term_eta_log - term_l_pi + sigma_l)
        
        # Normalize phase shift to [0, π) or (-π/2, π/2] or (-π, π]
        # Standard for scattering is often [0, π) or applying mod π
        delta_l_short_rad = mod(delta_l_short_rad_raw, π) # Julia's mod(x,p) gives result in [0,p) if x>0, or (p,0] if x<0 and p > 0
                                                         # For p=pi, gives [0,pi) or (-pi,0]. If negative, add pi.
        if delta_l_short_rad < -1e-9 # Check for slightly negative from mod operation.
            delta_l_short_rad += π
        end
        
        println("Asymptotic k = $k_asym")
        if Z_nucleus != 0.0
            println("Sommerfeld parameter η = $eta")
            println("Coulomb phase shift σ_$l_quantum = $(rad2deg(sigma_l)) degrees = $(sigma_l) radian")
        end
        println("Short-range phase shift δ_$(l_quantum)^short = $(rad2deg(delta_l_short_rad)) degrees = $(delta_l_short_rad) radian")
        if Z_nucleus != 0.0
             total_phase_contribution = rad2deg(mod(sigma_l + delta_l_short_rad_raw, π)) # or mod2pi
             println("Total phase shift (σ_l + δ_l^short) mod π = $total_phase_contribution degrees = $(mod(sigma_l + delta_l_short_rad_raw, π)) radian")
        end

    else
        println("Short-range potential V_sr($r_match) = $V_sr_at_match is not negligible compared to E = $E_continuum.")
        println("Phase shift calculation might be inaccurate. Increase r_max or check parameters.")
    end
else
    println("Solution has too few points or E_continuum <= 0, cannot calculate phase shift.")
end


println("\nNote: Atomic units used (ħ=1, m_e=1, e=1). Energy in Hartrees, distances in Bohr radii.")
println("Z_nucleus is the charge of the nucleus (e.g., 1 for Hydrogen).")
println("V0_ws for Woods-Saxon is the depth; negative for attractive potential.")

# --- Test Cases Suggestion ---
# 1. Pure Coulomb: Set V0_ws = 0.0. Then δ_l^short should be very close to 0.
#    sr_potential_params_test = (0.0, R0_ws, a_ws)
#    Then re-run with:
#    params_ode = (E_continuum, l_quantum, Z_nucleus, chosen_sr_potential_func, sr_potential_params_test)
#    prob = ODEProblem(schrodinger_radial_rhs!, u0, r_span, params_ode)
#    sol = solve(prob, alg, reltol=1e-9, abstol=1e-9, dense=true, saveat=0.1)
#    ... and then phase shift part.

# 2. No Coulomb (Z_nucleus = 0.0), only short-range.
#    Then η=0, σ_l=0. δ_l^short becomes the total phase shift δ_l.
#    This should match the results from the previous simpler script.

=======================================================================================================#
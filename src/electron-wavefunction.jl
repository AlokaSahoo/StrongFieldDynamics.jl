using Dierckx
# using OrdinaryDiffEq: ODEProblem, solve, Vern9, Tsit5
using SpecialFunctions
using DelimitedFiles

# calling the fortran shared object for distorted patial waves
dir = @__DIR__
freeSchrodinger = dir*"/../deps/mod_sfree.so"


"""
    atomic_electron(n::Int64, l::Int64)

Computes the atomic electron wavefunction for the corresponding `n` & `l` Principal and Orbital quantum numbers, repectively. 
"""
function atomic_electron(n::Int64, l::Int64)
    if (n == 3 && l == 0)
        data = readdlm(dir * "/../deps/Li-III.dat", skipstart = 2)
        r_, P_ = data[:,1], data[:,2]
    else
        error("Bound-state electron wavefunction not found.")
    end

    # Interpolating the 
    itp = Dierckx.Spline1D(r_, P_)
    P = itp.(r_)

    return r_, P
end


"""
    bessel_electron(ε::Flost64, l::Int64, r::Vector{Float64})

Computes the radial part of the continuum (free) electron with Bessel function 
                        `P (r) = r jₗ(pr)`
Returns the radial wavefunction of type `Vector{Float64}`
"""
function bessel_electron(ε::Float64, l::Int64, r::Vector{Float64})
    # linear momentum of the electron
    p = sqrt(2ε)    

    # Radial part of the continuum electron P = r * jₗ(pr)
    P = @. r * SpecialFunctions.sphericalbesselj(l, p * r)

    return r, P, 0.0
end


"""
    distorted_electron(ε::Float64, l::Int64, r::Vector{Float64}, pot::Vector{Float64})

Computes the distorted wave function of continuum photo electron with Bessel function
Returns the radial wavefunction and phase-shift of type `(Vector{Float64}, Float64)`
"""
function distorted_electron(ε::Float64, l::Int64, r::Vector{Float64}, rV::Vector{Float64})

    # Radial part and phase shift of the continuum electron
    npts=length(r)   # No of radila points
    r0=zeros(npts+1); r0[2:npts+1]=r
    rV0=zeros(npts+1); rV0[2:npts+1]=rV
    iPhase=Ref{Float64}(0.0)
    cPhase=Ref{Float64}(0.0)
    r_=zeros(Float64,25000)
    P=zeros(Float64,25000)
    Q=zeros(Float64,25000)

    #radial so file name defined at the toP

    ccall((:mysfree, freeSchrodinger), Cvoid, 
            (Ref{Int64},Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}), 
            npts, r0, rV0, ε, l, r_, P, Q, iPhase, cPhase)
    
    print("Inner phase shift = ", iPhase.x, "  coulomb phase shift = ", cPhase.x)

    # liner momentum
    p = sqrt(2ε)

    return r, (P[2:npts+1] ./ p), iPhase.x + cPhase.x
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
# function continuumElectron-old(E::Float64, l::Int64, pot::Vector{Float64})
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
#using DifferentialEquations

# using Plots
# using SpecialFunctions # For loggamma (complex gamma function)
# Plots.GRBackend()

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
using StrongFieldDynamics

#=----------------------------------------------------------------------------------------------
                Energy Distrbution
-----------------------------------------------------------------------------------------------=#
pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)

a_electron = StrongFieldDynamics.compute_atomic_electron(36, 4, 1) ;
r  = a_electron.r
IP = -a_electron.ε   #14/27.21138 #eV
nP(r, IP) = 2^2.5 * IP^1.5 * r * exp(-sqrt(2*IP)*r)
aP = nP.(r, IP)
a_electron = StrongFieldDynamics.AtomicElectron(36, 1, 0, 1//2, a_electron.ε, r, aP) ;

ed = compute_energy_distribution(a_electron, pulse; energy_range=(1e-8, 10*pulse.ω), n_points=400, coupled=true);
StrongFieldDynamics.plot_energy_distribution(ed)

#=----------------------------------------------------------------------------------------------
                Angular Distrbution
-----------------------------------------------------------------------------------------------=#

ad = compute_angular_distribution(a_electron, pulse; energy=4.5*pulse.ω, n_ϕ =200, θ = pi/2)
StrongFieldDynamics.plot_angular_distribution(ad)


#=----------------------------------------------------------------------------------------------
                Momentum Distrbution
-----------------------------------------------------------------------------------------------=#

pulse = StrongFieldDynamics.Pulse(I₀ = 1.5e14, λ=800, cycles=10, cep=float(0.), helicity=1, ϵ=1.0)

a_electron = StrongFieldDynamics.compute_atomic_electron(36, 4, 1) ;
r  = a_electron.r
IP = -a_electron.ε   #14/27.21138 #eV
nP(r, IP) = 2^2.5 * IP^1.5 * r * exp(-sqrt(2*IP)*r) ;
aP = nP.(r, IP) ;
a_electron = StrongFieldDynamics.AtomicElectron(36, 1, 0, 1//2, a_electron.ε, r, aP) ;

md = compute_momentum_distribution(a_electron, pulse, energy_range=(0.0, 10*pulse.ω), n_p=150, n_phi=300) ;
StrongFieldDynamics.plot_momentum_distribution(md)



# ps = md.p
# phis = md.φ

# dist = [ md.distribution[p, 1, phi] for phi in eachindex(phis), p in eachindex(ps)]

# using CairoMakie
# f = Figure(size = (800, 500))

# ax = PolarAxis(f[1, 1], title = "Surface")
# # rs = 0:10
# # phis = range(0, 2pi, 37)
# # cs = [r+cos(4phi) for phi in phis, r in rs]
# p = surface!(ax, phis, ps, zeros(size(dist)), color = dist, shading = NoShading, colormap = :coolwarm)
# ax.gridz = 100
# tightlimits!(ax) # surface plots include padding by default
# Colorbar(f[2, 1], p, vertical = false, flipaxis = false)
# display(f)

# f = Figure(size = (800, 500))

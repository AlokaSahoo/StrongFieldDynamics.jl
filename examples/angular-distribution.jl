using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)

a_electron = StrongFieldDynamics.compute_atomic_electron(36, 4, 1) ;
r  = a_electron.r
IP = -a_electron.ε   #14/27.21138 #eV
nP(r, IP) = 2^2.5 * IP^1.5 * r * exp(-sqrt(2*IP)*r)
aP = nP.(r, IP)
a_electron = StrongFieldDynamics.AtomicElectron(36, 1, 0, 1//2, a_electron.ε, r, aP) ;

ad = compute_angular_distribution(a_electron, pulse; energy=4.5*pulse.ω, n_ϕ =200, θ = pi/2)
StrongFieldDynamics.plot_angular_distribution(ad)

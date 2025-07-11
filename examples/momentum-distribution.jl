using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 1.5e14, λ=800, cycles=10, cep=float(0.), helicity=1, ϵ=1.0)

a_electron = StrongFieldDynamics.compute_atomic_electron(36, 4, 1) ;
r  = a_electron.r
IP = -a_electron.ε   #14/27.21138 #eV
nP(r, IP) = 2^2.5 * IP^1.5 * r * exp(-sqrt(2*IP)*r) ;
aP = nP.(r, IP) ;
a_electron = StrongFieldDynamics.AtomicElectron(36, 1, 0, 1//2, a_electron.ε, r, aP) ;

md = compute_momentum_distribution(a_electron, pulse, energy_range=(0.0, 10*pulse.ω), n_p=150, n_phi=300) ;
StrongFieldDynamics.plot_momentum_distribution(md)


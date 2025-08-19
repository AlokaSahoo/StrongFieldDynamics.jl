using Distributed

addprocs(48)

@everywhere using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 1.5e14, λ=800, cycles=10, cep=float(pi), helicity=1, ϵ=1.0)

settings = Settings(ionization_scheme=Atomic) #Settings(ionization_scheme=Atomic)

md = compute_momentum_distribution(36, pulse; settings = settings, energy_range=(0.0, 12*pulse.ω), n_p=150, n_phi=200) ;
# StrongFieldDynamics.plot_momentum_distribution(md, save_path="/home/alok/Softwares/StrongFieldDynamics.jl/examples/md-14Aug25.png")
StrongFieldDynamics.plot_momentum_distribution(md)

for worker in workers()
	rmprocs(worker)
end


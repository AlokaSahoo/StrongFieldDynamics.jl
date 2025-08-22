using Distributed

addprocs(48)

@everywhere using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)

settings = Settings() #Settings(ionization_scheme=Atomic)

md = compute_momentum_distribution(18, pulse; settings = settings, energy_range=(0.0, 12*pulse.ω), n_p=150, n_phi=200) ;
# StrongFieldDynamics.plot_momentum_distribution(md, save_path="/home/alok/Softwares/StrongFieldDynamics.jl/examples/md-14Aug25.png")
StrongFieldDynamics.plot_momentum_distribution(md, yScale=log10)
StrongFieldDynamics.plot_momentum_distribution(md, save_path="./examples/md-danish-8np-0.png")

for worker in workers()
	rmprocs(worker)
end
using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0)

settings = Settings()
        
ed = compute_energy_distribution(36, pulse; settings = settings, energy_range=(1e-8, 10*pulse.ω), n_points=400);
StrongFieldDynamics.plot_energy_distribution(ed)

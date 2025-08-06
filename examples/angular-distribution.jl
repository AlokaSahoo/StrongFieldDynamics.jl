using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 1.5e14, λ=800, cycles=10, cep=float(pi), helicity=1, ϵ=0.25)

settings = Settings()
        
ad = compute_angular_distribution(36, pulse; settings = settings, energy=4.5*pulse.ω, θ_range = (pi/2, pi/2), n_θ=1, n_ϕ =200)
StrongFieldDynamics.plot_angular_distribution(ad)

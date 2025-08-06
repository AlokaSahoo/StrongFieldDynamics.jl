using StrongFieldDynamics

pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0, envelope=Sin2)
pulse = StrongFieldDynamics.Pulse(I₀ = 5e13, λ=800, cycles=2, cep=float(pi), helicity=1, ϵ=1.0, envelope=Gaussian,
                                                                                         duration=(-pulse.Tp/2, pulse.Tp/2))


settings = Settings()
        
ed = compute_energy_distribution(36, pulse; settings = settings, energy_range=(1e-8, 10*pulse.ω), n_points=300);
StrongFieldDynamics.plot_energy_distribution(ed)

p_electron = StrongFieldDynamics.ContinuumElectron(0.5, sqrt(2*0.5), Bessel) ;

theta = pi/3
phi   = pi/2

StrongFieldDynamics.gaussianSv(pulse.Tp, float(theta), float(phi), pulse, p_electron)
StrongFieldDynamics.Sv_general(pulse.Tp, float(theta), float(phi), pulse, p_electron)
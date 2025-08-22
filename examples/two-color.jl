using StrongFieldDynamics

using Distributed, SlurmClusterManager

addprocs(SlurmManager())

@everywhere using StrongFieldDynamics

@everywhere unit(t) = 1.0

pulse1 = StrongFieldDynamics.Pulse(I₀ = 8.5e13, λ=780, cycles=15, cep=float(0), helicity=1, ϵ=1.0, envelope = unit) ;
pulse2 = StrongFieldDynamics.Pulse(I₀ = 5.6e13, λ=390, cycles=30, cep=float(0), helicity=-1, ϵ=1.0, envelope = unit) ;

pulse = StrongFieldDynamics.add_pulses(pulse1, pulse2)

# Plot the Vector potential
t = -pulse.Tp:0.01:pulse.Tp
StrongFieldDynamics.plot_vector_potential_trajectory(pulse, t, title="two-color Unit", save_path="./examples/vector-potential-unit.png")


settings = Settings()

# Momentum distribution
md = compute_momentum_distribution(18, pulse; settings = settings, energy_range=(1e-8, 12*pulse.ω[1]), n_p=150, n_phi=200) ;
StrongFieldDynamics.plot_energy_distribution(ed, title="two-color", save_path="./examples/two-colour-Unit-22Aug25.png")

# Energy distribution
ed = compute_energy_distribution(18, pulse; settings = settings, energy_range=(1e-8, 12*pulse.ω[1]), n_points=300);
StrongFieldDynamics.plot_energy_distribution(ed, title="two-color")



t = -pulse.Tp:0.01:pulse.Tp
StrongFieldDynamics.plot_vector_potential_trajectory(pulse, t, title="two-color")
StrongFieldDynamics.plot_vector_potential_trajectory(pulse1, t, title="monochromatic Sin2")

p_electron = StrongFieldDynamics.ContinuumElectron(0.5, sqrt(2*0.5), Bessel) ;

t = 0.5
θ = pi/3
ϕ = pi/3

@time StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse1, p_electron)
@time StrongFieldDynamics.Sv_general_cartesian(t, θ, ϕ, pulse1, p_electron)
@time StrongFieldDynamics.sin2Sv(t, θ, ϕ, pulse1, p_electron)

@time StrongFieldDynamics.Sv_general_prime_cartesian(t, θ, ϕ, pulse1, p_electron)
@time StrongFieldDynamics.Sv_prime_general(t, θ, ϕ, pulse1, p_electron)

@time StrongFieldDynamics.Sv_general_cartesian(t, θ, ϕ, pulse, p_electron)


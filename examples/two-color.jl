using StrongFieldDynamics

# using Distributed, SlurmClusterManager

# addprocs(SlurmManager())

# @everywhere using StrongFieldDynamics

# @everywhere unit(t) = 1.0

unit(t) = 1.0

pulse1 = StrongFieldDynamics.Pulse(I₀ = 8.5e13, λ=780, cycles=8, cep=float(0), helicity=1, ϵ=1.0, envelope = unit) ;
pulse2 = StrongFieldDynamics.Pulse(I₀ = 5.6e13, λ=390, cycles=16, cep=float(0), helicity=-1, ϵ=1.0, envelope = unit) ;

pulse = StrongFieldDynamics.add_pulses(pulse1, pulse2)

# Plot the Vector potential
t = -pulse.Tp:0.01:pulse.Tp
StrongFieldDynamics.plot_vector_potential_trajectory(pulse, t, title="two-color Unit")


settings = Settings()

# Energy distribution
ed = compute_energy_distribution(18, pulse; settings = settings, energy_range=(1e-8, 12*pulse.ω[1]), n_points=300);
StrongFieldDynamics.plot_energy_distribution(ed, title="two-color")

# Momentum distribution
md = compute_momentum_distribution(18, pulse; settings = settings, energy_range=(1e-8, 12*pulse.ω[1]), n_p=150, n_phi=200) ;
StrongFieldDynamics.plot_energy_distribution(ed, title="two-color", save_path="./examples/two-colour-Unit-22Aug25.png")

#================================= Testing the Volkov phase ============================================#
t = -pulse.Tp:0.01:pulse.Tp
StrongFieldDynamics.plot_vector_potential_trajectory(pulse, t, title="two-color")
StrongFieldDynamics.plot_vector_potential_trajectory(pulse1, t, title="monochromatic Sin2")

p_electron = StrongFieldDynamics.ContinuumElectron(0.5, sqrt(2*0.5), Bessel) ;

t = 10.5
θ = pi/4
ϕ = pi/4

@time StrongFieldDynamics.Sv_general(t, θ, ϕ, pulse1, p_electron)
@time StrongFieldDynamics.Sv_general_cartesian(t, θ, ϕ, pulse1, p_electron)
@time StrongFieldDynamics.sin2Sv(t, θ, ϕ, pulse1, p_electron)

@time StrongFieldDynamics.Sv_prime_general_cartesian(t, θ, ϕ, pulse1, p_electron)
@time StrongFieldDynamics.Sv_prime_general(t, θ, ϕ, pulse1, p_electron)

@time StrongFieldDynamics.Sv_general_cartesian(t, θ, ϕ, pulse, p_electron)

#================================= Testing the Pulse Shape Integral ====================================#
a_electron = StrongFieldDynamics.compute_atomic_electron(18, Atomic)

# monochromatic first
@time StrongFieldDynamics.F1_integral_quadgk(pulse1, a_electron, p_electron, θ, ϕ ; sign=1)
@time StrongFieldDynamics.F1_integral_levin_approxfun(pulse1, a_electron, p_electron, θ, ϕ ; sign=1)

@time StrongFieldDynamics.F1_integral_levin_approxfun(pulse1, pulse, a_electron, p_electron, θ, ϕ ; sign=1)
@time StrongFieldDynamics.F1_integral_quadgk(pulse1, pulse, a_electron, p_electron, θ, ϕ ; sign=1)


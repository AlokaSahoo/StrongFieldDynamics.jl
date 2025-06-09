# Theory

This page provides an overview of the theoretical foundations behind **StrongFieldDynamics.jl**.

## Strong Field Approximation (SFA)

The SFA is a widely used approach for modeling ionization of atoms in intense laser fields. It treats the interaction of the electron with the laser field non-perturbatively, while the atomic potential is included approximately.

- **Direct Ionization:** The electron is ionized directly by the laser field.
- **Rescattering:** After ionization, the electron can be driven back to the parent ion by the laser field, leading to rescattering phenomena.

## Partial Wave Expansion

The electron continuum in the laser field is given by Volkov-type states which, similarly, can be expanded into (distorted) partial waves

```math
\left| \chi_{\mathbf{p}}(t) \right\rangle = \sqrt{\frac{2}{\pi}} e^{-i S_V(t)} \sum_{\ell_p=0}^{\infty} \sum_{m_p=-\ell_p}^{\ell_p} Y_{\ell_p m_p}^*\left(\vartheta_p, \varphi_p\right) \left| \varepsilon_p \ell_p m_p \right\rangle \otimes \left| \chi_{f, m_s^{\prime}} \right\rangle,
```

In the velocity gauge, the interaction between the continuum electron and laser ﬁeld is thereby accounted for by the so-called Volkov phase

```math
S_V(t) = \frac{1}{2} \int^t d\tau \left[ \mathbf{p} + \mathbf{A}(\tau) \right]^2
```

This expansion is crucial for capturing the angular structure of the ionization process.

## Atomic Potential Effects

The package provides several methods to include the influence of the atomic potential on the ionized electron, improving upon the plain SFA.

## Direct Ionization Amplitude

The direct ionization amplitude is computed using the partial wave expansion of the Volkov states. The amplitude for direct ionization can be expressed as:

```math
\begin{aligned}
T_0\left(\mathbf{p}, m_j, m_s^{\prime}\right) = & -\frac{i}{\sqrt{2 \pi}} \mathcal{F}_1[\omega ; f ; \mathbf{p}] \left( \sum_{\ell_p=0}^{\infty} \sum_{j_p \geq 1/2} \sum_{q=0, \pm 1} (-1)^q u_q Y_{\ell_p, m_j - m_s^{\prime} - q}\left(\vartheta_p, \varphi_p\right) \right. \\
& \quad \left. \times \left\langle \ell_p \left(m_j - m_s^{\prime} - q\right), \frac{1}{2} m_s^{\prime} \middle| j_p \left(m_j - q\right) \right\rangle \left\langle j m_j, 1 (-q) \middle| j_p \left(m_j - q\right) \right\rangle \left\langle \varepsilon_p \ell_p j_p \middle\| \mathbf{p} \middle\| n \ell j \right\rangle \right) \\
& -\frac{i}{\sqrt{2 \pi}} \mathcal{F}_1[-\omega ; f ; \mathbf{p}] \left( \sum_{\ell_p=0}^{\infty} \sum_{j_p \geq 1/2} \sum_{q=0, \pm 1} u_q^* Y_{\ell_p, m_j - m_s^{\prime} + q}\left(\vartheta_p, \varphi_p\right) \right. \\
& \quad \left. \times \left\langle \ell_p \left(m_j - m_s^{\prime} + q\right), \frac{1}{2} m_s^{\prime} \middle| j_p \left(m_j + q\right) \right\rangle \left\langle j m_j, 1 q \middle| j_p \left(m_j + q\right) \right\rangle \left\langle \varepsilon_p \ell_p j_p \middle\| \mathbf{p} \middle\| n \ell j \right\rangle \right) \\
& -\frac{i}{\sqrt{2 \pi}} \mathcal{F}_2[f ; \mathbf{p}] Y_{\ell, m_j - m_s^{\prime}}\left(\vartheta_p, \varphi_p\right) \left\langle \ell \left(m_j - m_s^{\prime}\right), \frac{1}{2} m_s^{\prime} \middle| j m_j \right\rangle \left\langle \varepsilon_p \ell j m_j \middle| n \ell j m_j \right\rangle,
\end{aligned}
```

And the Ionization Probability is given by

```math
\mathrm{P}(\mathbf{p}) = = \frac{p}{2 j + 1} \sum_{m_j = -j}^{j} \sum_{m_s^{\prime} = \pm 1/2} \left| T\left(\mathbf{p}, m_j, m_s^{\prime}\right) \right|^2
```

## References

- Birger Böning, Stephan Fritzsche, "Partial-wave representation of the strong-field approximation", Phys. Rev. A 102, 053108 (2020).
- Birger Böning, Stephan Fritzsche, "Partial-wave representation of the strong-field approximation. II. Coulomb asymmetry in the photoelectron angular distribution of many-electron atoms", Phys. Rev. A 107, 023108 (2023).


For more details, see the source code and comments in the package.

# Theory

This page provides an overview of the theoretical foundations behind **StrongFieldDynamics.jl**.

## Strong Field Approximation (SFA)

The SFA is a widely used approach for modeling ionization of atoms in intense laser fields. It treats the interaction of the electron with the laser field non-perturbatively, while the atomic potential is included approximately.

- **Direct Ionization:** The electron is ionized directly by the laser field.
- **Rescattering:** After ionization, the electron can be driven back to the parent ion by the laser field, leading to rescattering phenomena.

## Partial Wave Expansion

The electron wavefunction is expanded in terms of spherical harmonics, allowing for the analysis of angular momentum components and angular distributions of photoelectrons.

$$ \Psi(\mathbf{r}, t) = \sum_{l,m} R_{l}(r, t) Y_{l}^{m}(\theta, \phi) $$

This expansion is crucial for capturing the angular structure of the ionization process.

## Atomic Potential Effects

The package provides several methods to include the influence of the atomic potential on the ionized electron, improving upon the plain SFA.

## References

- Birger Böning, Stephan Fritzsche, "Partial-wave representation of the strong-field approximation", Phys. Rev. A 102, 053108 (2020).
- Birger Böning, Stephan Fritzsche, "Partial-wave representation of the strong-field approximation. II. Coulomb asymmetry in the photoelectron angular distribution of many-electron atoms", Phys. Rev. A 107, 023108 (2023).


For more details, see the source code and comments in the package.

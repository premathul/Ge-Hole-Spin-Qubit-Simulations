# BdG strain-cycle seed: audited scope

This is an **independent toy calculation**, not a realistic Ge hole-spin simulation, and not a demonstration of fermion-parity pumping. It is isolated to avoid touching recovered simulation data.

## Run

```bash
python -m pip install numpy pytest
python -m pytest -q topological_phonon_seed/tests
python -m topological_phonon_seed.kitaev_strain
```

The chain BdG matrix is checked for Hermiticity and particle-hole symmetry. A known static class-D Kitaev-chain invariant is benchmarked. Two phenomenological controls trace a gapped loop. **Neither the static invariant nor loop closure diagnoses a parity pump.** There is no cycle invariant in this checkpoint; its status is explicitly NOT COMPUTED.

## Research-critical next steps

1. Construct confined four-band Ge Luttinger–Kohn + Bir–Pikus states with documented units, strain tensor conventions and boundary conditions.
2. Obtain two *physical* independent elastic strain tensor eigenfields, then calculate their matrix elements in the confined hole basis.
3. Introduce physically antisymmetric superconducting pairing, including interface and magnetic-field assumptions.
4. Define and independently benchmark the correct topological **cycle** invariant (including gauge and boundary-condition subtleties), and test parity transport and edge spectral flow.
5. Evaluate the full gap throughout the cycle, finite size, disorder and adiabatic driving constraints.
6. Only after these gates pass, estimate realistic mode amplitudes and possible photon/spin readout.

Negative results must be retained. No originality or experimental feasibility is claimed.

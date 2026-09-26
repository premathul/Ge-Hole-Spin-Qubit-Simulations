# Recovered simulation status

Last recovered state: 2026-09-25.

## Zeeman / g observable

The production campaign used strict low-field solves with overlap tracking of the qubit Kramers subspace.

Recovered state:
- coarse: 30/30 per dot
- baseline: 30/30 per dot
- fine-left: 18/30
- fine-right: queued

Baseline left-dot principal effective g values:
- 13.82824
- 0.229149
- 0.123098

Baseline right-dot principal effective g values:
- 13.74079
- 0.212792
- 0.178042

These were provisional pending mesh convergence.

## Tunnel coupling

Recovered fine-mesh resonance estimate:

`t_c = 7.1834 ± 0.6895 µeV`

or

`t_c / h = 1.7369 ± 0.1667 GHz`

The successive mesh correction was still about 9.6%, so this was approaching convergence rather than fully mesh converged.

## Charge stability

Recovered conditional quantities:
- `U_L = 4.17197 meV`
- `U_R = 4.14020 meV`
- `U_m = 0.907526 meV`

Charge cells from `(0,0)` through `(2,2)` were resolved in the many-body calculation.

## T2* and T1

At the recovered checkpoint, the production workflow had not yet promoted a final publication-ready `T2*` or phonon-limited `T1`.

The intended dependency chain was:

mesh-converged `G` → electrical susceptibility → charge-noise/sweet-spot analysis → conditional `T2*` → finite-q spin-phonon `T1`.

Any earlier lifetime numbers should be checked against the final production outputs before being treated as validated.

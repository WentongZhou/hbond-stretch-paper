# Target H-bond force constant across environments (WP4)

k_X = d²X/dR² at R_min of the target O···O (kcal/mol/Å²), f_X = k_X/k_total. CT = electrons lost by the acceptor fragment of the target bond at R_min.

| environment | R_min Å | E_int | k_total | ω(H₂O) | k_Elec | k_ExRep | k_OrbRel | k_CorrDisp | f_Elec | f_ExRep | f_OrbRel | f_CorrDisp | CT Mull | CT IAO |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| dimer | 2.9660 | -5.088 | 19.424 | 159.5 | -32.63 | 81.65 | -24.93 | -4.67 | -1.680 | +4.204 | -1.283 | -0.240 | +0.0257 | +0.0448 |
| trimer acid-ctrl | 2.8750 | -7.008 | 26.071 | 184.8 | -41.34 | 109.99 | -36.33 | -6.25 | -1.585 | +4.219 | -1.394 | -0.240 | +0.0315 | +0.0564 |
| trimer acid-ctrl pair | 2.9673 | -5.098 | 19.393 | 159.4 | -32.62 | 82.10 | -25.30 | -4.78 | -1.682 | +4.233 | -1.305 | -0.247 | +0.0238 | +0.0445 |
| trimer base-ctrl | 2.8945 | -6.957 | 25.226 | 181.7 | -42.18 | 108.06 | -34.99 | -5.66 | -1.672 | +4.284 | -1.387 | -0.225 | +0.0334 | +0.0565 |
| trimer base-ctrl pair | 2.9718 | -5.099 | 19.278 | 158.9 | -32.24 | 81.30 | -25.01 | -4.77 | -1.672 | +4.217 | -1.297 | -0.247 | +0.0232 | +0.0441 |

## Where does Δk come from?

Δk_X = k_X(env) − k_X(ref); share = Δk_X / Δk_total. A share > 1 means the channel over-explains the change and is compensated by the others.

| environment | reference | Δk_total | Δk/k | Δω cm⁻¹ | ΔR_min Å | Δk_Elec (share) | Δk_ExRep (share) | Δk_OrbRel (share) | Δk_CorrDisp (share) | ΔCT IAO |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| trimer acid-ctrl | dimer | +6.647 | +0.342 | +25.3 | -0.0911 | -8.70 (-1.31) | +28.34 (+4.26) | -11.41 (-1.72) | -1.58 (-0.24) | +0.0117 |
| trimer acid-ctrl pair | dimer | -0.031 | -0.002 | -0.1 | +0.0012 | +0.01 (-0.38) | +0.45 (-14.63) | -0.37 (+12.18) | -0.12 (+3.83) | -0.0003 |
| trimer base-ctrl | dimer | +5.802 | +0.299 | +22.3 | -0.0715 | -9.55 (-1.65) | +26.41 (+4.55) | -10.06 (-1.73) | -1.00 (-0.17) | +0.0118 |
| trimer base-ctrl pair | dimer | -0.146 | -0.008 | -0.6 | +0.0058 | +0.39 (-2.70) | -0.35 (+2.42) | -0.09 (+0.58) | -0.10 (+0.69) | -0.0007 |

## Apparent exponent n in k ∝ δq^n across environments

| proxy | n | R² | points |
|---|---:|---:|---:|
| Mulliken CT | 0.89 | 0.920 | 5 |
| IAO CT | 1.17 | 0.993 | 5 |

Attractive channels at R_min (share of the restoring force −dX/dR among channels with dX/dR > 0):

| environment | Elec | OrbRel | CorrDisp |
|---|---:|---:|---:|
| dimer | 0.595 | 0.306 | 0.099 |
| trimer acid-ctrl | 0.566 | 0.336 | 0.099 |
| trimer acid-ctrl pair | 0.592 | 0.309 | 0.099 |
| trimer base-ctrl | 0.580 | 0.329 | 0.091 |
| trimer base-ctrl pair | 0.592 | 0.309 | 0.099 |

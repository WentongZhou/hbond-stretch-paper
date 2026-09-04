# Force-constant partition: dfcheck_acid_pair

R_min = 2.9713 Å, E_int = -4.747 kcal/mol, k_total = 18.384 kcal/mol/Å² = 12.77 N/m, ω(H2O) = 155.2 cm⁻¹, ω(D2O) = 147.2 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -4.747 | 0.000 | 18.384 | +1.000 | -87.71 |
| Elec | -7.299 | 12.450 | -29.558 | -1.608 | 87.06 |
| Exch | -8.038 | 25.881 | -82.224 | -4.473 | 261.54 |
| Rep | 13.973 | -47.596 | 160.421 | +8.726 | -543.54 |
| ExRep | 5.934 | -21.715 | 78.197 | +4.253 | -281.99 |
| OrbRel | -2.302 | 7.019 | -24.940 | -1.357 | 93.14 |
| Corr | -0.365 | 1.622 | -5.153 | -0.280 | 14.75 |
| Disp | -0.715 | 0.623 | -0.162 | -0.009 | -0.67 |
| CorrDisp | -1.080 | 2.245 | -5.315 | -0.289 | 14.08 |
| Steric | 4.854 | -19.470 | 72.882 | +3.964 | -267.91 |

Sum of primary slopes at R_min (should be ~0): 0.0000 kcal/mol/Å
Sum of primary k minus k_total (closure): -6.93e-08 kcal/mol/Å²
Electrons lost by fragment 'acceptor_water' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.016493198565530928, 'slope_e_per_A': -0.025684001725742093}, IAO {'value_e': 0.04039262905575654, 'slope_e_per_A': -0.08467353243297689}
Max |closure| over scan: 1.14e-12 Eh

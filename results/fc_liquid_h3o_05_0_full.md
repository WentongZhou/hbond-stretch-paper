# Force-constant partition: liquid_h3o_05_0_full

R_min = 3.0025 Å, E_int = -13.222 kcal/mol, k_total = 100.963 kcal/mol/Å² = 70.15 N/m, ω(H2O) = 363.6 cm⁻¹, ω(D2O) = 344.8 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -13.222 | -0.000 | 100.963 | +1.000 | 90.16 |
| Elec | -39.256 | 0.488 | -40.867 | -0.405 | -68.69 |
| Exch | -66.001 | 7.366 | -196.906 | -1.950 | -188.84 |
| Rep | 118.088 | -13.211 | 407.363 | +4.035 | 434.42 |
| ExRep | 52.087 | -5.845 | 210.457 | +2.085 | 245.58 |
| OrbRel | -16.517 | 2.707 | -55.832 | -0.553 | -76.96 |
| Corr | -4.257 | 2.130 | -13.697 | -0.136 | -10.07 |
| Disp | -5.279 | 0.519 | 0.901 | +0.009 | 0.30 |
| CorrDisp | -9.536 | 2.650 | -12.796 | -0.127 | -9.77 |
| Steric | 42.551 | -3.195 | 197.661 | +1.958 | 235.81 |

Sum of primary slopes at R_min (should be ~0): -0.0002 kcal/mol/Å
Sum of primary k minus k_total (closure): 1.15e-07 kcal/mol/Å²
Electrons lost by fragment 'water_2nd' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.03023898811159412, 'slope_e_per_A': -0.000581871385413842}, IAO {'value_e': 0.037768983543099786, 'slope_e_per_A': 0.06863618302934735}
Max |closure| over scan: 3.21e-12 Eh

At the reference geometry R0 = 2.9338 Å (residual force -6.775 kcal/mol/Å): k = 97.137, ω = 356.6 cm⁻¹, f_Elec = -0.382, f_ExRep = +2.045, f_OrbRel = -0.538, f_CorrDisp = -0.125

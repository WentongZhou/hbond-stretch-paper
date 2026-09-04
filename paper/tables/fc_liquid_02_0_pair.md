# Force-constant partition: liquid_02_0_pair

R_min = 2.9992 Å, E_int = -4.068 kcal/mol, k_total = 16.192 kcal/mol/Å² = 11.25 N/m, ω(H2O) = 145.6 cm⁻¹, ω(D2O) = 138.1 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -4.068 | -0.001 | 16.192 | +1.000 | -76.20 |
| Elec | -6.266 | 10.998 | -26.251 | -1.621 | 75.91 |
| Exch | -6.875 | 21.946 | -69.153 | -4.271 | 217.85 |
| Rep | 11.832 | -39.940 | 133.565 | +8.249 | -447.46 |
| ExRep | 4.957 | -17.994 | 64.412 | +3.978 | -229.61 |
| OrbRel | -1.629 | 5.017 | -17.766 | -1.097 | 65.93 |
| Corr | -0.382 | 1.306 | -4.004 | -0.247 | 12.27 |
| Disp | -0.748 | 0.673 | -0.199 | -0.012 | -0.70 |
| CorrDisp | -1.130 | 1.979 | -4.202 | -0.260 | 11.57 |
| Steric | 3.827 | -16.015 | 60.210 | +3.718 | -218.04 |

Sum of primary slopes at R_min (should be ~0): -0.0007 kcal/mol/Å
Sum of primary k minus k_total (closure): -7.02e-09 kcal/mol/Å²
Electrons lost by fragment 'acceptor_water' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.013887342958149286, 'slope_e_per_A': -0.013494389988685648}, IAO {'value_e': 0.03383056126772435, 'slope_e_per_A': -0.06857191607819127}
Max |closure| over scan: 1.20e-12 Eh

At the reference geometry R0 = 2.9117 Å (residual force -1.741 kcal/mol/Å): k = 23.999, ω = 177.3 cm⁻¹, f_Elec = -1.410, f_ExRep = +3.654, f_OrbRel = -1.020, f_CorrDisp = -0.223

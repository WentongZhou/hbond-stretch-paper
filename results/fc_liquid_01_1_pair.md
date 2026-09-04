# Force-constant partition: liquid_01_1_pair

R_min = 3.2428 Å, E_int = -2.146 kcal/mol, k_total = 7.525 kcal/mol/Å² = 5.23 N/m, ω(H2O) = 99.3 cm⁻¹, ω(D2O) = 94.1 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -2.146 | -0.000 | 7.525 | +1.000 | -34.47 |
| Elec | -2.698 | 4.021 | -9.217 | -1.225 | 27.47 |
| Exch | -2.895 | 9.228 | -28.796 | -3.827 | 89.11 |
| Rep | 4.698 | -15.732 | 52.039 | +6.915 | -172.10 |
| ExRep | 1.803 | -6.504 | 23.243 | +3.089 | -82.99 |
| OrbRel | -0.513 | 1.319 | -4.403 | -0.585 | 16.08 |
| Corr | -0.144 | 0.529 | -1.759 | -0.234 | 5.43 |
| Disp | -0.595 | 0.634 | -0.338 | -0.045 | -0.45 |
| CorrDisp | -0.739 | 1.163 | -2.097 | -0.279 | 4.98 |
| Steric | 1.064 | -5.341 | 21.145 | +2.810 | -78.01 |

Sum of primary slopes at R_min (should be ~0): -0.0002 kcal/mol/Å
Sum of primary k minus k_total (closure): -1.51e-07 kcal/mol/Å²
Electrons lost by fragment 'acceptor_water' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.01010384035637911, 'slope_e_per_A': -0.012167958151412835}, IAO {'value_e': 0.01491396458676556, 'slope_e_per_A': -0.029678364629009548}
Max |closure| over scan: 1.17e-12 Eh

At the reference geometry R0 = 2.9235 Å (residual force -5.030 kcal/mol/Å): k = 27.682, ω = 190.4 cm⁻¹, f_Elec = -0.883, f_ExRep = +2.563, f_OrbRel = -0.507, f_CorrDisp = -0.173

# Force-constant partition: liquid_06_1_embed

R_min = 3.1785 Å, E_int = -2.463 kcal/mol, k_total = 9.230 kcal/mol/Å² = 6.41 N/m, ω(H2O) = 109.9 cm⁻¹, ω(D2O) = 104.3 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -2.463 | -0.000 | 9.230 | +1.000 | -43.73 |
| Elec | -3.386 | 5.710 | -13.457 | -1.458 | 37.47 |
| Exch | -4.082 | 12.878 | -39.505 | -4.280 | 118.98 |
| Rep | 6.747 | -22.283 | 72.246 | +7.827 | -232.42 |
| ExRep | 2.665 | -9.405 | 32.741 | +3.547 | -113.44 |
| OrbRel | -0.699 | 2.175 | -7.518 | -0.814 | 26.79 |
| Corr | -0.357 | 0.832 | -2.231 | -0.242 | 6.02 |
| Disp | -0.687 | 0.688 | -0.304 | -0.033 | -0.57 |
| CorrDisp | -1.044 | 1.520 | -2.536 | -0.275 | 5.45 |
| Steric | 1.621 | -7.885 | 30.205 | +3.272 | -107.99 |

Sum of primary slopes at R_min (should be ~0): -0.0002 kcal/mol/Å
Sum of primary k minus k_total (closure): 5.60e-08 kcal/mol/Å²
Electrons lost by fragment 'acceptor_water' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.015450198022294774, 'slope_e_per_A': -0.012064608171329943}, IAO {'value_e': 0.0223496596118674, 'slope_e_per_A': -0.04200323559340561}
Max |closure| over scan: 1.73e-12 Eh

At the reference geometry R0 = 3.1047 Å (residual force -0.810 kcal/mol/Å): k = 12.886, ω = 129.9 cm⁻¹, f_Elec = -1.282, f_ExRep = +3.270, f_OrbRel = -0.756, f_CorrDisp = -0.231

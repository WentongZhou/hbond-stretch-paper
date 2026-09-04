# Force-constant partition: liquid_01_1_embed

R_min = 3.2129 Å, E_int = -2.546 kcal/mol, k_total = 8.330 kcal/mol/Å² = 5.79 N/m, ω(H2O) = 104.4 cm⁻¹, ω(D2O) = 99.1 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -2.546 | -0.000 | 8.330 | +1.000 | -36.47 |
| Elec | -3.117 | 4.165 | -8.962 | -1.076 | 28.05 |
| Exch | -3.130 | 9.752 | -29.895 | -3.589 | 92.40 |
| Rep | 5.094 | -16.670 | 54.291 | +6.517 | -179.75 |
| ExRep | 1.964 | -6.918 | 24.396 | +2.929 | -87.35 |
| OrbRel | -0.595 | 1.555 | -4.929 | -0.592 | 17.37 |
| Corr | -0.184 | 0.554 | -1.851 | -0.222 | 5.95 |
| Disp | -0.614 | 0.644 | -0.324 | -0.039 | -0.48 |
| CorrDisp | -0.798 | 1.198 | -2.175 | -0.261 | 5.47 |
| Steric | 1.166 | -5.720 | 22.221 | +2.668 | -81.89 |

Sum of primary slopes at R_min (should be ~0): -0.0003 kcal/mol/Å
Sum of primary k minus k_total (closure): 8.10e-07 kcal/mol/Å²
Electrons lost by fragment 'acceptor_water' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.010239603776941768, 'slope_e_per_A': -0.008953061922039854}, IAO {'value_e': 0.016978929387384795, 'slope_e_per_A': -0.0312890765806992}
Max |closure| over scan: 1.55e-12 Eh

At the reference geometry R0 = 2.9235 Å (residual force -4.603 kcal/mol/Å): k = 26.538, ω = 186.4 cm⁻¹, f_Elec = -0.852, f_ExRep = +2.550, f_OrbRel = -0.521, f_CorrDisp = -0.178

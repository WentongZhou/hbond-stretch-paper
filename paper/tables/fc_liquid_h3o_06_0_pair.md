# Force-constant partition: liquid_h3o_06_0_pair

R_min = 3.0528 Å, E_int = -3.905 kcal/mol, k_total = 14.285 kcal/mol/Å² = 9.92 N/m, ω(H2O) = 136.8 cm⁻¹, ω(D2O) = 129.7 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -3.905 | -0.000 | 14.285 | +1.000 | -65.58 |
| Elec | -5.617 | 8.975 | -20.782 | -1.455 | 60.35 |
| Exch | -5.684 | 18.148 | -56.946 | -3.986 | 178.08 |
| Rep | 9.616 | -32.443 | 108.199 | +7.574 | -360.73 |
| ExRep | 3.932 | -14.295 | 51.253 | +3.588 | -182.65 |
| OrbRel | -1.370 | 3.746 | -12.760 | -0.893 | 47.12 |
| Corr | -0.140 | 0.890 | -3.154 | -0.221 | 10.31 |
| Disp | -0.709 | 0.683 | -0.272 | -0.019 | -0.71 |
| CorrDisp | -0.849 | 1.573 | -3.426 | -0.240 | 9.60 |
| Steric | 3.083 | -12.721 | 47.827 | +3.348 | -173.05 |

Sum of primary slopes at R_min (should be ~0): -0.0003 kcal/mol/Å
Sum of primary k minus k_total (closure): 1.26e-07 kcal/mol/Å²
Electrons lost by fragment 'acceptor_water' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.008326742485662952, 'slope_e_per_A': -0.0040517306427056456}, IAO {'value_e': 0.026986868197425525, 'slope_e_per_A': -0.055459881220623085}
Max |closure| over scan: 9.10e-13 Eh

At the reference geometry R0 = 2.9215 Å (residual force -2.531 kcal/mol/Å): k = 25.129, ω = 181.4 cm⁻¹, f_Elec = -1.215, f_ExRep = +3.233, f_OrbRel = -0.820, f_CorrDisp = -0.198

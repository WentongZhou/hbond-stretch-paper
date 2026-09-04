# Force-constant partition: liquid_h3o_02_0_full

R_min = 2.9196 Å (no minimum in the scanned range), E_int = -13.016 kcal/mol, k_total = 43.488 kcal/mol/Å² = 30.21 N/m, ω(H2O) = 238.6 cm⁻¹, ω(D2O) = 226.3 cm⁻¹ (mass effect only)

| component | X(R_min) kcal/mol | dX/dR kcal/mol/Å | k_X kcal/mol/Å² | k_X / k_total | cubic kcal/mol/Å³ |
|---|---:|---:|---:|---:|---:|
| Total | -13.016 | -14.015 | 43.488 | +1.000 | -88.85 |
| Elec | -26.313 | 10.834 | -4.282 | -0.098 | 59.49 |
| Exch | -35.895 | 43.057 | -60.838 | -1.399 | 201.39 |
| Rep | 63.039 | -78.945 | 120.843 | +2.779 | -424.39 |
| ExRep | 27.144 | -35.888 | 60.004 | +1.380 | -223.00 |
| OrbRel | -8.657 | 6.499 | -7.919 | -0.182 | 60.34 |
| Corr | -1.842 | 3.170 | -4.634 | -0.107 | 15.45 |
| Disp | -3.348 | 1.370 | 0.319 | +0.007 | -1.13 |
| CorrDisp | -5.190 | 4.540 | -4.315 | -0.099 | 14.32 |
| Steric | 21.954 | -31.348 | 55.690 | +1.281 | -208.68 |

Sum of primary slopes at R_min (should be ~0): -14.0148 kcal/mol/Å
Sum of primary k minus k_total (closure): -1.94e-07 kcal/mol/Å²
Electrons lost by fragment 'water_2nd' (acceptor of the target bond) at R_min: Mulliken {'value_e': 0.020450171019046465, 'slope_e_per_A': -0.007200920776387645}, IAO {'value_e': 0.043895498064372415, 'slope_e_per_A': -0.05096211274909126}
Max |closure| over scan: 2.25e-12 Eh

At the reference geometry R0 = 2.9196 Å (residual force -14.015 kcal/mol/Å): k = 43.488, ω = 238.6 cm⁻¹, f_Elec = -0.098, f_ExRep = +1.380, f_OrbRel = -0.182, f_CorrDisp = -0.099

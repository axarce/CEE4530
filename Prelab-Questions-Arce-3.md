## Laboratory 3 Prelab Questions
### Anthony Arce

1. How many grams of NaHCO3 would be required to keep the ANC levels in a lake above 50 µeq/L for 3 hydraulic residence times given an influent pH of 3.0 and a lake volume of 4 L, if the current lake ANC is 0 µeq/L?


```python
from aide_design.play import *
LakeVol = 4*u.L
NaHCO3MolWt = 84 * u.g / u.mol
pH_Influent = 3.0
H_Conc = 10 ** -pH_Influent * u.mol / u.L
ResidenceTimes = 3  # unitless
ANC_In = H_Conc
print(ANC_In)
ANC_Out = 50*(10 ** -6) * u.mol/ u.L
print(ANC_Out)

ANC_0 = (ANC_Out - ANC_In * (1 - np.exp(ResidenceTimes))) * np.exp(ResidenceTimes)
print(ANC_0)

NaHCO3_0 = ANC_0
NaHCO3_AddWt = LakeVol * NaHCO3_0 * NaHCO3MolWt
print(NaHCO3_AddWt)
```

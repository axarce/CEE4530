### Laboratory 6: Gas Transfer Prelaboratory Questions
## Anthony Arce
```python
from aide_design.play import*
```

## 1)	Calculate the mass of sodium sulfite needed to reduce all the dissolved oxygen in 4 L of pure water in equilibrium with the atmosphere and at 30°C.

```Python
Volume = 4 * u.L
Temp = 30 * u.degC
dO2 = 7.6 * u.mg / u.L
dO2_molWt = 32000 * u.mg / u.mol
Na2SO3_molWt = 126000 * u.mg / u.mol
mass_dO2 = dO2 * Volume
mol_ratio = 2 # 2 moles sulfite for 1 mol oxygen gas
mass_Na2SO3 = (mass_dO2/dO2_molWt) * mol_ratio * Na2SO3_molWt
print(mass_Na2SO3)
```
The mass of sodium sulfite needed to reduce all the dissolved oxygen in 4 liters of water at 30 degrees Celsius in equilibrium with the atmosphere is 239.4 milligrams.

## 2)	Describe your expectations for dissolved oxygen concentration as a function of time during a reaeration experiment.  Assume you have added enough sodium sulfite to consume all of the oxygen at the start of the experiment. What would the shape of the curve look like?

The shape of the curve would approximately be a parabola starting from the initial DO concentration going to nearly zero that is concave up with a tail end that will reach the saturation DO point eventually.

## 3)	Why is kv,i  not zero when the gas flow rate is zero? How can oxygen transfer into the reactor even when no air is pumped into the diffuser?

The volumetric gas transfer coefficient can be non-zero due to the exposure of the reactor to the atmosphere, allowing for gas transfer from the aqueous phase to the gaseous phase.


## 4)	Describe your expectations for k,vi  as a function of gas flow rate. Do you expect a straight line? Why?

Expectations for kv,i as a function of gas flow rate would decrease as the gas flow rate increases as it becomes more difficult to supersaturate the water beyond its saturation point and excess oxygen will simply bubble out.

## 5)	A dissolved oxygen probe was placed in a small vial in such a way that the vial was sealed. The water in the vial was sterile. Over a period of several hours the dissolved oxygen concentration gradually decreased to zero. Why? (You need to know how dissolved oxygen probes work to answer this!)

The dissolved oxygen concentration will tend to zero in a sterile vial as the dissolved oxygen probe will reduce the dissolved oxygen in the closed vial to H20, using up the available O2. Since the system is closed to the atmosphere, there is no aeration of the water and therefore cannot replenish as the probe consumes O2.

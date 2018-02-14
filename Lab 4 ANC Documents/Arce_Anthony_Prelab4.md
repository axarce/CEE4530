## Laboratory 4 PreLab Questions
##### Ajara Cobourne

1. The ANC of
```python
from aide_design.play import*
Kw = 10**(-14) * (u.mole/u.L)**2
K1_carbonate = 10**(-6.37)*u.mol/u.L
K2_carbonate = 10**(-10.25)*u.mol/u.L
K_Henry_CO2 = 10**(-1.5) * u.mole/(u.L*u.atm) # ** is ^
P_CO2 = 10**(-3.5) * u.atm

pH = 3.5
ANC_Cayuga = .0016*u.mol/u.L #of carbonates
ANC_Wolf = .000070*u.mol/u.L

def invpH(pH):
  return 10**(-pH)*u.mol/u.L

def alpha0_carbonate(pH):
   alpha0_carbonate = 1/(1+(K1_carbonate/invpH(pH))*(1+(K2_carbonate/invpH(pH))))
   return alpha0_carbonate

def alpha1_carbonate(pH):
  alpha1_carbonate = 1/((invpH(pH)/K1_carbonate) + 1 + (K2_carbonate/invpH(pH)))
  return alpha1_carbonate

def alpha2_carbonate(pH):
  alpha2_carbonate = 1/(1+(invpH(pH)/K2_carbonate)*(1+(invpH(pH)/K1_carbonate)))
  return alpha2_carbonate

# Alpha plotting
# Create a uniform spaced array from 3 to 12
pH_graph = np.linspace(3,12,500)
plt.plot(pH_graph, alpha0_carbonate(pH_graph),'r', pH_graph, alpha1_carbonate(pH_graph),'b',pH_graph, alpha2_carbonate(pH_graph),'g')
plt.xlabel('pH')
plt.ylabel('Fraction of total carbonates')
plt.legend(['alpha_0', 'alpha_1', 'alpha_2'], loc='best')
#plt.savefig('images/alphagraph.png')
#plt.show()

def ANC_closed(pH,Total_Carbonates):
  return Total_Carbonates*(alpha1_carbonate(pH)+2*alpha2_carbonate(pH)) + Kw/invpH(pH) - invpH(pH)

def ANC_open(pH):
  return ANC_closed(pH,P_CO2*K_Henry_CO2/alpha0_carbonate(pH))

import scipy
from scipy import optimize
optimize.brentq(ANC_open(pH),5, 7)
```
2. For a pure strong acid at pH 3.2 we can simplify the equation for
Acid Neutralizing Capacity (ANC)
$${\text{ANC}} = [HCO_3^ - {\text{] + 2[CO}}_3^{ - 2}{\text{] + [O}}{{\text{H}}^{\text{ - }}}{\text{] - [}}{{\text{H}}^{\text{ + }}}{\text{]}}$$
to   (ANC)
$${\text{ANC}} = [HCO_3^ - {\text{] + 2[CO}}_3^{ - 2}{\text{] + [O}}{{\text{H}}^{\text{ - }}}{\text{] - [}}{{\text{H}}^{\text{ + }}}{\text{]}}$$
because XYZ
```Python
from aide_design.play import*
import scipy
from scipy import optimize

ANC_in = -10**-3.5*u.mol/u.L
print(ANC_in)
ANC_out_Wolf = ANC_in*(1-np.exp(-0.2)) + (ANC_Wolf*np.exp(-0.2))
ANC_out_Cayuga =  ANC_in*(1-np.exp(-0.2)) + (ANC_Cayuga*np.exp(-0.2))
print(ANC_Wolf)
print(ANC_Cayuga)
print(ANC_out_Wolf, ANC_out_Cayuga)
ANC_o_Wolf = -1.122e-08#*u.mol/u.L
ANC_o_Cayuga = 0.001253#*u.mol/u.L

Kw = 10**(-14) * (u.mole/u.L)**2
K1_carbonate = 10**(-6.37)*u.mol/u.L
K2_carbonate = 10**(-10.25)*u.mol/u.L
K_Henry_CO2 = 10**(-1.5) * u.mole/(u.L*u.atm)
P_CO2 = 10**(-3.5) * u.atm

def ANC_closed(pH,Total_Carbonates):
  return Total_Carbonates*(alpha1_carbonate(pH)+2*alpha2_carbonate(pH)) + Kw/invpH(pH) - invpH(pH)

def ANC_open(pH):
  return ANC_closed(pH,P_CO2*K_Henry_CO2/alpha0_carbonate(pH))

#pH = 3.2

# Strip the units off of the ANC function so that scipy can calculate the root.
def ANC_open_unitless(pH):
  return (ANC_open(pH)).magnitude

def pH_open():
  return optimize.brentq(ANC_open_unitless, 1, 12)
#help(optimize.brentq)
pH_open() #ph_open is ok 5.68
print('The pH of is', pH_open())
#plt ANC graph
plt.plot(pH_graph, ANC_open(pH_graph),'r')
plt.xlabel('pH')
plt.ylabel('ANC')
plt.yscale('log')
plt.savefig('images/ANCgraph.png')
plt.show()

def wavenumber(T, h):
  """Return the wavenumber of wave using period and water height from bed."""
  q = 10  # this is a guess to find pH
  diff = (P_CO2.magnitude*K_Henry_CO2/alpha0_carbonate(q).magnitude)*(alpha1_carbonate(q).magnitude + 2*alpha2_carbonate(q).magnitude) + (Kw.magnitude/ 10**-q) - 10**-q
  while diff<0:
      LHS = ANC_o_Cayuga
      RHS = (0.0003162*0.03162/alpha0_carbonate.magnitude)*(alpha1_carbonate.magmagnitude + 2*alpha2_carbonate.magnitude) + (1e-14/ 10**-q.magnitude) - 10**-q.magnitude
      diff = LHS - RHS
      q = q - 0.0001 #iterate to solve
  return q
  print(q)
```

2. ANC influent is 0.00063095734 mol/L This is a strong acid so the hydrogen ion concentration is much higher than the carbonate concentration. For a pure strong acid at pH 3.2 we can simplify the equation for
Acid Neutralizing Capacity (ANC)
$${\text{ANC}} = [HCO_3^ - {\text{] + 2[CO}}_3^{ - 2}{\text{] + [O}}{{\text{H}}^{\text{ - }}}{\text{] - [}}{{\text{H}}^{\text{ + }}}{\text{]}}$$
to   (ANC)
$${\text{ANC}} = - {{\text{[H]}}^{\text{ + }}}{\text{}}$$
```Python
print(10**-3.2)

```

3. [H+] is not a conserved species because there are carbonates for them to react with
$${H_2}CO_3^* \overset {K_1} \longleftrightarrow {H^+} + HCO_3^- $$


$$HCO_3^ - \overset {{K_2}} \longleftrightarrow {H^ + } + CO_3^{ - 2}$$

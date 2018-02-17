
```Python
from aide_design.play import*
```
###### Anthony Arce
###### Elle Blake
###### Ajara Cobourne

#### **Data Analysis for Week 3 Lab: Acid Rain**
## **Q1. Experiment 1**
```Python
# import your file
file = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain_Good Data.csv') # the second experiment, I'm calling it #1 since we're treating it as #1 for analysis
array1 = np.array(file)
Flow_Rate = 4.499 * u.milliliter / u.s
Volume = 4 * u.L
ResidenceTime = (Volume / Flow_Rate).to(u.min)
Time_Min = array1[:,1] * u.min
pH1 = array1[:,2]
# plotting
plt.figure('ax',(10,8))
plt.plot(Time_Min, pH1)
# put in your x and y variables
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('pH', fontsize=15)
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\exp1.jpg')
plt.show()

```
The results are shown in Figure 1.
 ![graph](C:\Users\Anthony\github\CEE4530_axa2\images\exp1.jpg)

## Q2 Experiment 1. Calculate alpha 0, 1, and 2 based on the pH measures at each time point
```Python

Kw = 10**(-14) * (u.mole/u.L)**2
K1_carbonate = 10**(-6.37)*u.mol/u.L
K2_carbonate = 10**(-10.25)*u.mol/u.L
K_Henry_CO2 = 10**(-1.5) * u.mole/(u.L*u.atm)
P_CO2 = 10**(-3.5) * u.atm


ph_exp1 = np.array([7.91356, 7.242599, 6.663552, 4.388395, 3.649247])
ph_exp2 = np.array([7.966759, 8.002065,  7.61718, 7.44151, 6.447749])


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

alpha0_exp1 = alpha0_carbonate(ph_exp1)
alpha1_exp1 = alpha1_carbonate(ph_exp1)
alpha2_exp1 = alpha2_carbonate(ph_exp1)

alphatotal_exp1 = alpha0_exp1 + alpha1_exp1 + alpha2_exp1

# Alpha total array is unable to be displayed without removing units of dimensionless, so .magnitude is necessary to obtain the desired array.
print(alphatotal_exp1.magnitude)
```


## **Q3 Experiment 1.**
```python
pH_0min = 7.91356
H_0min = 10**-7.91356
ANC_in = -10 ** -3 *  u.mol / u.L
ANC_0 =  (.6235 * u.g / (84.0661 * (u.g / u.mol)) ) / (4*u.L)

ANC_Out = np.zeros(len(Time_Min)) * u.mol/u.L
print(len(Time_Min))
print(len(ANC_Out))
for i in range(0,len(Time_Min)):
  ANC_Out[i] = ANC_in * (1-np.exp(-Time_Min[i]/ResidenceTime)) + ANC_0 * np.exp(-Time_Min[i]/ResidenceTime)
  i = i + 1

plt.figure('ax',(10,8))
plt.plot(Time_Min, ANC_Out)
# put in your x and y variables
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('ANC_Effluent (mg/L)', fontsize=15)
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\Question3.jpg')
plt.show()

```



## **Q4 Experiment 1**

Calculate ANC using the equation (1.11):
$ANC = C_{T} * (alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$

```python
ANC_Out_Closed = np.zeros(len(pH1)) * u.mol / u.L
for i in range(0,len(pH1)):
  ANC_Out_Closed[i] = ANC_in * (alpha1_carbonate(pH1[i]) + alpha2_carbonate(pH1[i])) + Kw / invpH(pH1[i]) - invpH(pH1[i])
  i = i + 1




```

## Q5 Experiment 1
Calculate ANC using the equation (1.15):
  $ANC = \frac{P_{CO_{2}}K_{H}}{a_{0}} * (alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$



```python
ANC_Out_Open = np.zeros(len(pH1)) * u.mol / u.L
for i in range(0, len(pH1)):
  ANC_Out_Open[i] = (P_CO2  * K_Henry_CO2) * (alpha1_carbonate(pH1[i]) + alpha2_carbonate(pH1[i])) / (alpha0_carbonate(pH1[i])) + (Kw / invpH(pH1[i])) - invpH(pH1[i])
  i = i + 1

plt.figure('ax',(10,8))
ANC_Effluent = plt.plot(Time_Min, ANC_Out, Time_Min, ANC_Out_Closed, Time_Min, ANC_Out_Open)
# put in your x and y variables
plt.legend(['Conservative', 'Closed', 'Open'], loc = 'best')
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('ANC (mg/L)', fontsize=15)
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\ANC_Effluent.jpg')
plt.show()

```
## Q6 Experiment 1
Equation 1.21 gives: $ANC_{out} = \left [ ANC_{in}\cdot \left ( 1-e ^{\frac{-t}{\Theta }} \right )\right ]+   ANC_{0}\cdot e^{\frac{-t}{\Theta }}$
Modeling the system as a completely mixed flow reaction, we can determine $C_{T}$ as a function of time:
$C = C_{in}(1-e^{\frac{-t}{\theta }})+C_{0}e^{\frac{-t}{\theta}}$

For a conservative species $C_{T}$ becomes:

$C_{T}=C_{T_{0}}\cdot e^{\frac{-t}{\theta}}$


```python
C_0 = 0
NaHCO3_in = (0.6235 / 84.007) * u.mol
C_in = NaHCO3_in/Volume
C_in
array1 = np.array(file)
Flow_Rate = 4.499 * (u.milliliter / u.s)
Volume = 4 * u.L
ResidenceTime = (Volume / Flow_Rate).to(u.min)
Time_Min = array1[:,1] * u.min

C_T_conservative = np.zeros(len(Time_Min)) * u.mol / u.L
for i in range(0,len(Time_Min)):
  C_T_conservative[i] = C_in * np.exp(-Time_Min[i]/ResidenceTime)
  i = i + 1

C_T_conservative
# plotting
plt.figure('Conservative Ct',(10,8))
plt.plot(Time_Min, C_T_conservative)
# put in your x and y variables
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('C_T', fontsize=15)
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\C_T.jpg')
plt.show()

```
# Q7 Experiment 1
Derive an equation for $C_{T}$ as a function of ANC and pH using the following equation:
$ANC = C_{T}(alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$

Solving for $C_{T}$ (assuming $[OH]^{-}$ is negligible) yields:
$C_T = \frac{ANC_{t} + [H^{+}]}{(a_{1} + 2 a_{2})}$

```python
C_T_closed = np.zeros(len(Time_Min)) * u.mol / u.L

for i in range(0, len(Time_Min)):
  C_T_closed[i] = (ANC_Out[i] + invpH(pH1[i])) / (alpha1_carbonate(pH1[i]) + alpha2_carbonate(pH1[i]))
  i = i + 1



```
Equilibrium $C_T = \frac{P_{CO2} + K_{H}}{a_{0}}$


# Q8 Experiment 1
```python
C_T_equil = np.zeros(len(pH1)) * u.mol / u.L

for i in range(0, len(pH1)):
  C_T_equil[i] = (P_CO2 * K_Henry_CO2)  / alpha0_carbonate(pH1[i])
  i = i + 1

plt.figure('Ct',(10,8))
plt.plot(Time_Min, C_T_conservative, Time_Min, C_T_closed, Time_Min, C_T_equil)
# put in your x and y variables
plt.legend(['Conservative', 'Closed', 'Equilibrium'], loc = 'best')
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('C_T (mol/L)', fontsize=15)
plt.yscale('log')
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\C_T.jpg')
plt.show()
```


## Q9 Experiment 1
The lake is best modeled after a closed system, it behaves as a non-volatile system. To change the lake to a volatile system, there would need to be mixing due to a wind current to drive equilibrium with the atmosphere.



## Q10. Analyze the data from the 2nd experiment and graph the data appropriately what did you learn from the 2nd experiment


```python
# Create a uniform spaced array from 3 to 12
pH_graph = np.linspace(3,12,50)
plt.plot(pH_graph, alpha0_carbonate(pH_graph),'r', pH_graph, alpha1_carbonate(pH_graph),'b')
plt.xlabel('pH')
plt.ylabel('Fraction of total carbonates')
plt.legend(['alpha_0', 'alpha_1'], loc='best')

plt.savefig('images/AcidRain_Comparison.png')
plt.show()
```

The comparison of the first and second experiment are shown in Figure 1.
 ![graph](images/AcidRain_Comparison.png)

** Figure 1. Lake pH as a function of hydraulic residence time(?).**

## Q10. Analyze the data from the 2nd experiment and graph the data appropriately what did you learn from the 2nd experiment
## Experiment 2
```python
Flow_Rate2 = 5.0779 * u.milliliter / u.s
Volume = 4 * u.L
ResidenceTime = (Volume / Flow_Rate).to(u.min)
exp2 = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain _2.csv')
exp2_array = np.array(exp2)
Time_Min2 = exp2_array[:,2]
ph_exp2 = exp2_array[:,3]
plt.figure('ax',(10,8))
plt.plot(Time_Min2, ph_exp2)
# put in your x and y variables
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('pH', fontsize=15)
plt.savefig()
plt.show()
```

# Additional Questions

1) If we assume that the NaCO3 does not get mixed in with the water ANC_out = (ANC_in*(1-exp^-(time/theta))+ANC_0*exp^-(time/theta)) #will ask at OH

2) Some complicating factors are:

  -The extent of mixing could affect the results due to the residence time not being consistent throughout the experiment. Further, the solution may not be uniform if there were reduced mixing, causing the pH probe to get inaccurate readings that are not representative of the "lake" as a whole.

  -The solubility of CaCO3 is less soluble than NaCO3 in water; however, it becomes more soluble in more acidic water.

  -Density of a CaCO3 slurry is 2.71 g/cm³

  -Too much calcium ions could increase the hardness of the water

Q11. Don’t forget to write a brief paragraph on conclusions and on suggestions using Markdown.

​Conclusions:

The pH probe turned out to be very sensitive since it is really measuring voltage and not pH. In our second experiment we did not ground the apparatus and our graph reflects that. Although the general shape of a titration curve is evident there are many fluctuations.

Suggestions:
Try to change something that does not affect the pH probe since it is too sensitive, The da
Q12. Verify that your report and graphs meet the requirements as outlined in the course materials.

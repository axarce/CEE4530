```python
from aide_design.play import*
```
###### Anthony Arce
###### Elle Blake
###### Ajara Cobourne

###### All laboratory reports (all 3 flavors!) must include conclusions, suggestions/comments and well done graphs!

#### **Data Analysis for Week 3 Lab: Acid Rain**
# **Q1. Experiment 1**
```python
# import your file
file = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain_Good Data.csv') # the second experiment, I'm calling it #1 since we're treating it as #1 for analysis
array1 = np.array(file)
Flow_Rate = 4.499 * u.milliliter / u.s
Volume = 4 * u.L
ResidenceTime = Volume / Flow_Rate
Time_Min = array1[:,1]
pH1 = array1[:,2]
# plotting
plt.figure('ax',(10,8))
plt.plot(Time_Min, pH1)
# put in your x and y variables
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('pH', fontsize=15)
plt.savefig('images/AcidRain_1.png')
plt.show()

```
The alpha terms representing the carbonate species are shown in Figure 1.
 ![graph](images/AcidRain_1.png)

# ** Q2 Experiment 1 Calculate alpha1, 2, and 3 based on the pH measures at each time point**
```python
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
alphatotal_exp1
```

# **Q3 Experiment 1**
```python
pH_0min = 7.91356
H_0min = -10^7.91356
ANC_0 = (.6235/84.0661)/4

def ANC_out(ANC_in, ANC_0, time, theta)
  ANC_out = (ANC_in*(1-exp^-(time/theta))+ANC_0*exp^-(time/theta))
  return ANC_out

  plt.figure('ax',(10,8))
  plt.plot(Time_Min, pH1)
  # put in your x and y variables
  plt.xlabel('Time (min)', fontsize=15)
  plt.ylabel('pH', fontsize=15)
  plt.savefig('images/AcidRain_1_conservative.png')
  plt.show()

```
# **Q4 Experiment 1**

Calculate ANC using the equation (1.11):
$ANC = C_{T}*(alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$

```python
C_T*ANC_0*(alpha1_carbonate) + 2*alpha2_carbonate

def ANC_out(ANC_in, ANC_0, time, theta)
  ANC_out = (ANC_in*(1-exp^-(time/theta))+ANC_0*exp^-(time/theta)) return ANC_out

```



Calculate ANC using the equation (1.15):
  $ANC = \frac{P_{CO_{2}}K_{H}}{a_{0}}*(alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$
```python
ANC=P_CO2*K_Henry_CO2/alpha0_exp1*(alpha1_exp1+2*alpha2_exp1)+K_Henry_CO2/H_0min-H_0min
return ANC
```
###Q6 Experiment 1
Equation 1.21 gives: $ANC_{out} = -\left [ ANC_{in}\cdot \left ( 1-e ^{\frac{-t}{\Theta }} \right )\right ]+   ANC_{0}\cdot e^{\frac{-t}{\Theta }}$


** Q10. Analyze the data from the 2nd experiment and graph the data appropriately what did you learn from the 2nd experiment **


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

## Q10. Analyze the data from the 2nd experiment and graph the data appropriately what did you learn from the 2nd experiment **
```python
Flow_Rate2 = 5.0779 * u.milliliter / u.s
Volume = 4 * u.L
ResidenceTime = Volume / Flow_Rate
exp2 = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain _2.csv')
exp2_array = np.array(exp2)
Time_Min2 = exp2_array[:,2]
ph_exp2 = exp2_array[:,3]
plt.figure('ax',(10,8))
plt.plot(Time_Min2, ph_exp2)
# put in your x and y variables
plt.xlabel('Time (min)', fontsize=15)
plt.ylabel('pH', fontsize=15)
plt.show()
```

Q11. Don’t forget to write a brief paragraph on conclusions and on suggestions using Markdown.
Q12. Verify that your report and graphs meet the requirements as outlined in the course materials.

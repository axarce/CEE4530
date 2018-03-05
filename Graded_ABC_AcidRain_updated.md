```Python
from aide_design.play import*
from scipy import stats
```
###### Anthony Arce
###### Elle Blake
###### Ajara Cobourne

#### **Laboratory Assignments 3 and 4: Acid Rain and Acid Neutralizing Capacity**


### Introduction and Objectives
We are doing this experiment to test the boundaries of the acid neutralizing capacity (ANC) of our lake and examining various methods of determining ANC. To meet our goal we have to test the ANC at various points in the process of remediating a lake with NaCO3. We hoped to observe hydrogen ion conservation, learn how to use the Gran plot to calculate the maximum amount of acid the lake can take before no longer being able to neutralize added acid. Finally we wish to examine whether the ANC of our experiment matched the theoretical ANCs. We expect the results to be applicable to future remediation projects we conduct as Environmental Engineers.This would be particularly useful in projects comparing the qualities of different bodies of water, as ANC is an important characteristic to consider. We will be able to understand more about ANC through this experiment by using the following derivations:

$ANC_{out} = \left [ ANC_{in}\cdot \left ( 1-e ^{\frac{-t}{\Theta }} \right )\right ]+   ANC_{0}\cdot e^{\frac{-t}{\Theta }}$

Modeling the system as a completely mixed flow reaction, we can determine $C_{T}$ as a function of time:

$C = C_{in}(1-e^{\frac{-t}{\theta }})+C_{0}e^{\frac{-t}{\theta}}$

For a conservative species $C_{T}$ becomes:

$C_{T}=C_{T_{0}}\cdot e^{\frac{-t}{\theta}}$

There are many different ways to calculate the ANC when modeling the lake using different methods. Assuming that the lake can be modeled as a completely mixed flow reactor, we calculate the "conservative ANC" using the following derivation.

Equation 1.21 gives: $ANC_{out} = \left [ ANC_{in}\cdot \left ( 1-e ^{\frac{-t}{\Theta }} \right )\right ]+   ANC_{0}\cdot e^{\frac{-t}{\Theta }}$
Modeling the system as a completely mixed flow reaction, we can determine $C_{T}$ as a function of time:
$C = C_{in}(1-e^{\frac{-t}{\theta }})+C_{0}e^{\frac{-t}{\theta}}$

For a conservative species $C_{T}$ becomes:

$C_{T}=C_{T_{0}}\cdot e^{\frac{-t}{\theta}}$
However, if we were to assume that no carbonates are exchanged with the atmosphere throughout the experiment, we would describe ANC as that of a closed system and calculate it using the following equation (1.13):
$ANC = C_{T} * (alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$
In order to get a full understanding of ANC in different conditions, we can also calculate the "open ANC" under the assumption that carbonates are at equilibrium with the atmosphere using the following equation (1.15):
  $ANC = \frac{P_{CO_{2}}K_{H}}{a_{0}} * (alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$


### Procedure(s)

## Laboratory 3: Acid Rain
For the first trial of this experiment we followed the instructions as given in the CEE 4530 Spring 2018 lab Manual and used a tube size of 17. For the second trial of the experiment we continuously monitored pH using a calibrated pH probe. We pumped in acid rain while draining out lake water at a constant flow rate and using a tube size of 17. We added bromocresol green to the lake and added of NaHCO3  while continuously stirring. Lake water samples were taken every 5 minutes from time zero minutes to 20 minutes. Finally we measured the effluent flow rate. We deviated from the lab manual and did not ground our equipment, which allowed electricity to flow through the lake throughout the second trail.

## Laboratory 4: Acid Neutralizing Capacity
The procedure was followed as stated in the CEE 4530 Spring 2018 Lab Manual with a few modifications. The team did not measure out 50 ml for each sample, but measured out volumes between 40 ml and 50 ml and recorded the mass and volume of “lake water” for each trail. The titrant used was 0.1 M HCl instead of the concentration stated in the lab manual. Additionally, when adding titrant to lake samples that had a pH higher than approximately 6.5, the titrant was added in increments larger than 0.1 mL.


### Data Analysis


```Python
# import your file
# Q1 Plot measured pH of the lake versus dimensionless hydraulic residence time (t/q).
file = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain_Good Data.csv') # the second experiment, I'm calling it #1 since we're treating it as #1 for analysis
array1 = np.array(file)
Flow_Rate = 4.499 * u.milliliter / u.s
Volume = 4 * u.L
ResidenceTime = (Volume / Flow_Rate).to(u.min)
ResidenceTime_DimLess = array1[:,1] * u.min / ResidenceTime
Time_Min = array1[:,1]*u.min
pH1 = array1[:,2]
# plotting
plt.figure()
plt.plot(ResidenceTime_DimLess, pH1)
# put in your x and y variables
plt.xlabel('Residence Time')
plt.ylabel('pH')
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\exp1.jpg')
plt.show()

```

```python
#Q2 Calculate ANC assuming lake is a completely mixed flow reactor. Expected theoretical ANC. labeled as "convervative ANC"
Kw = 10**(-14) * (u.mole/u.L)**2
K1_carbonate = 10**(-6.37)*u.mol/u.L
K2_carbonate = 10**(-10.25)*u.mol/u.L
K_Henry_CO2 = 10**(-1.5) * u.mole/(u.L*u.atm)
P_CO2 = 10**(-3.5) * u.atm
pH_0min = 7.91356
H_0min = 10**-7.91356
ANC_in = -10 ** -3 *  u.mol / u.L
ANC_0 =  (.6235 * u.g / (84.0661 * (u.g / u.mol)) ) / (4*u.L)
C_0 = 0
NaHCO3_in = (0.6235 / 84.007) * u.mol
C_in = NaHCO3_in/Volume


array1 = np.array(file)
Flow_Rate = 4.499 * (u.milliliter / u.s)
Volume = 4 * u.L
ResidenceTime = (Volume / Flow_Rate).to(u.min)
print(ResidenceTime)
Time_Min = array1[:,1] * u.min
pH_0min = 7.91356
H_0min = 10**-7.91356
ANC_in = -10 ** -3 *  u.mol / u.L
ANC_0 =  (.6235 * u.g / (84.0661 * (u.g / u.mol)) ) / (4*u.L)
ANC_Out = np.zeros(len(Time_Min)) * u.mol/u.L

for i in range(0,len(Time_Min)):
  ANC_Out[i] = ANC_in * (1-np.exp(-Time_Min[i]/ResidenceTime)) + ANC_0 * np.exp(-Time_Min[i]/ResidenceTime)
  i = i + 1


# plotting
plt.plot(Time_Min/ResidenceTime, ANC_Out)
plt.xlabel('Residence Time(s) (unitless)')
plt.ylabel('ANC (mg/L)')
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain Q2 image.jpg')
plt.show()

```

```python
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

# Q3 Using eq 1.11 ANC effluent Closed system
ANC_Out_Closed = np.zeros(len(pH1)) * u.mol / u.L

for i in range(0,len(pH1)):
  ANC_Out_Closed[i] = ANC_in * (alpha1_carbonate(pH1[i]) + alpha2_carbonate(pH1[i])) + Kw / invpH(pH1[i]) - invpH(pH1[i])
  i = i + 1

```

```python
#Q4 calculate ANC assuming an open system, and that carbonates are at equilibrium with the atmosphere (equation 1.15)
ANC_Out_Open = np.zeros(len(pH1)) * u.mol / u.L
for i in range(0, len(pH1)):
  ANC_Out_Open[i] = (P_CO2  * K_Henry_CO2) * (alpha1_carbonate(pH1[i]) + alpha2_carbonate(pH1[i])) / (alpha0_carbonate(pH1[i])) + (Kw / invpH(pH1[i])) - invpH(pH1[i])
  i = i + 1


ANCMeasured = np.array([0.00139,0.000916,0.000282,0.000035,-0.0001853])
Time = np.array([0,5,10,15,20])*u.min


plt.figure()
ANC_Effluent = plt.plot(Time_Min/ResidenceTime, ANC_Out, Time_Min/ResidenceTime, ANC_Out_Closed, Time_Min/ResidenceTime, ANC_Out_Open,Time/ResidenceTime, ANCMeasured)


plt.legend(['Conservative', 'Closed', 'Open','Measured'], loc = 'best')
plt.xlabel('Residence Times')
plt.ylabel('ANC (mg/L)')
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain Q3 image.jpg')

plt.show()

```

```python
#Q5 analyze data from second experiment and graph
Flow_Rate2 = 5.0779 * u.milliliter / u.s
Volume = 4 * u.L
ResidenceTime = (Volume / Flow_Rate).to(u.min)
exp2 = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain _2.csv')
exp2_array = np.array(exp2)
Time_Min2 = exp2_array[:,2]
ph_exp2 = exp2_array[:,3]
plt.figure()
plt.plot(Time_Min2, ph_exp2)
# put in your x and y variables
plt.xlabel('Time (min)')
plt.ylabel('pH')
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\Acid Rain Q5 image.jpg')
plt.show()
```
```python
#Q1, plot the titration curve of the t=0 sample with 0.05N HCl
titrPHfile = pd.read_csv(r'C:\Users\Anthony\github\CEE4530_axa2\Time0TitrantVolpH.csv')
titrPHarray = np.array(titrPHfile)


VTitrant = titrPHarray[:,0] * u.mL
pHtime0 = titrPHarray[:,1]
EquivVol = np.array([0.6915]) * u.mL
# pHregresion is the value of pH associated with the equivalent volume calculated by using linear regression on the first three points for the time 0.
pHregression = np.array([5.752641])

plt.plot(VTitrant, pHtime0, 'r-', EquivVol, pHregression, 'bo')

plt.xlabel('Titrant Vol (mL)')
plt.ylabel('pH')
plt.legend(['Titration Curve', 'Equivalent Vol'])
plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\titrantpHtime0.jpg')
plt.show()
```

```python
#prepare a Gran plot using titration curve data of the t=0 sample
def invpH(pH):
  return 10**(-pH)
#Define the gran function.

def F1(V_sample,V_titrant,pH):
  return (V_sample + V_titrant)/V_sample * invpH(pH)
VSample = 50*u.mL
#Create an array of the F1 values.
F1_data = F1(VSample,VTitrant,pHtime0)
F1_data
#By inspection I guess that there are 4 good data points in the linear region.
N_good_points = 3
#use scipy linear regression. Note that the we can extract the last n points from an array using the notation [-N:]
slope, intercept, r_value, p_value, std_err = stats.linregress(VTitrant[-N_good_points:],F1_data[-N_good_points:])
#reattach the correct units to the slope and intercept.
intercept = intercept*u.mole/u.L
slope = slope*(u.mole/u.L)/u.mL
V_eq = -intercept/slope

print(V_eq)
#The equivilent volume agrees well with the value calculated by ProCoDA.
#create an array of points to draw the linear regression line
x=[V_eq.magnitude,VTitrant[-1].magnitude ]
y=[0,(VTitrant[-1]*slope+intercept).magnitude]
#Now plot the data and the linear regression
plt.plot(VTitrant, F1_data,'o',x,y,'r')
plt.xlabel('Titrant Volume (mL)')
plt.ylabel('Gran function (mole/L)')
plt.legend(['data'])

plt.savefig(r'C:\Users\Anthony\github\CEE4530_axa2\images\Gram.jpg')
plt.show()

```

```python
#plot the measured ANC of the lake on the same graph as was used to plot the conservative, volatile, and nonvolatile ANC models.
# ANCMeasured is just an array of the calculated ANC concentrations from the Gran analysis excel files. The laboratory group ABC wrote down the measured ANC values.

ANCMeasured = np.array([0.00139,0.000916,0.000282,0.000035,-0.0001853])
RelResTime = np.array([0,5,10,15,20])*u.min
#Conservative, Volatile, and nonvolatile ANC models observed
ANC_Conserv = 0.001851 mole / liter
```

### Results

The pH of the lake versus hydraulic residence are in Figure 1.
 ![graph](C:\Users\Anthony\github\CEE4530_axa2\images\exp1.jpg)

Figure 1. Measured pH of lake versus dimensionless hydraulic residence time


The XYZ representing the blah blah are shown in Figure 2.
![graph](C:\Users\Anthony\github\CEE4530_axa2\Acid Rain Q2 image.jpg)

Figure 2. Effluent Lake ANC modeled as a CMFR with conservative ANC


The XYZ representing the blah blah are shown in Figure 3.
 ![graph](C:\Users\Anthony\github\CEE4530_axa2\Acid Rain Q3 image.jpg)

Figure 3. Effluent Lake ANC with a closed system


The XYZ representing the blah blah are shown in Figure 5.
 ![graph](C:\Users\Anthony\github\CEE4530_axa2\Acid Rain Q5 image.jpg)

Figure 4. Effluent Lake ANC with electricity flowing throughout.



### Discussion

In order to get a full understanding of ANC in different conditions, as a completely mixed flow reactor, as a closed system, as well as an open system.

In modeling the system as a completely mixed flow reactor, we can determine ANC using the following equation (1.21) to determine the "conservative ANC", see figure 2 from the results section:
$ANC_{out} = \left [ ANC_{in}\cdot \left ( 1-e ^{\frac{-t}{\Theta }} \right )\right ]+   ANC_{0}\cdot e^{\frac{-t}{\Theta }}$



Alternatively, if we were to assume that no carbonates are exchanged with the atmosphere throughout the experiment, we would describe ANC as that of a closed system and calculate it using the following equation (1.13), see figure 3 from the results section:
$ANC = C_{T} * (alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$


Finally, we can also calculate the "open ANC" under the assumption that carbonates are at equilibrium with the atmosphere using the following equation (1.15), see figure 4 from the results section:
  $ANC = \frac{P_{CO_{2}}K_{H}}{a_{0}} * (alpha_{1} + 2alpha_{2}) + \frac{K_{w}}{\left [ H^{+} \right ]}-\left [ H^{+} \right ]$

The differences amongst these plot demonstrate that


### Conclusions
The conclusions section should not include any new observation. It is the place to summarize the results in a few sentences. Make sure you connect your conclusions to your objectives for doing the research.

The pH probe turned out to be very sensitive since it is really measuring voltage and not pH. In our second experiment we did not ground the apparatus and our graph reflects that. Although the general shape of a titration curve is evident there are many fluctuations.

#### Suggestions and Comments

The data acquisition software inserts comments directly into the data files and makes it difficult to extract data during analysis. However, using the code Monroe shared with the class this problem has been addressed. The magnetic stirrer has a short cord and getting the pH probe to stretch from the computer to the lake was tricky. The setup of the lake proved to have some issues as water leaked out of the pipe before reaching the beaker as we measured flow rate. Additionally since the sinks drain on opposite sides of the bench and are shallow and small it was difficult to measure the effluent flow rate by measuring the exit volume of “lake water” over a short period of time. The experiment could be easier to understand if there were more discussion regarding the trend that we will likely to see from our ProCoda data. By failing to ground our experiment in one of the trials, it took our group a bit to recognize that there was an error in the ProCoda data that we were collecting. I believe having more data sharing from Trail 1 between groups would provide valuable insight into ANC properties since each group does not have time to do every experiment. Another way to analyze the data would be using the Gran plot which involves a titration where the equivalence volume is estimated based on pH and titrant volume. Overall the lab did provide valuable insight into managing acid rain lake remediation and the difference between theoretical and experimental methods.

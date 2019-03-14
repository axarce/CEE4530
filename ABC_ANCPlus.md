
```Python
from aide_design.play import*
from scipy import stats
import importlib
```
###### Anthony Arce
###### Elle Blake
###### Ajara Cobourne

#### **Laboratory Assignments 3 and 4: Acid Rain and Acid Neutralizing Capacity**

##Objectives

We are doing this experiment to test the boundaries of the acid neutralizing capacity (ANC) of our lake and examining various methods of determining ANC. To meet our goal we have to test the ANC at various points in the process of remediating a lake with NaCO3. We hoped to observe hydrogen ion conservation, learn how to use the Gran plot to calculate the maximum amount of acid the lake can take before no longer being able to neutralize added acid. Finally we wish to examine whether the ANC of our experiment matched the theoretical ANCs. We expect the results to be applicable to future remediation projects we conduct as Environmental Engineers.This would be particularly useful in projects comparing the qualities of different bodies of water, as ANC is an important characteristic to consider.

##Procedures

The procedure was followed as stated in the CEE 4530 Spring 2018 Lab Manual with a few modifications. The team did not measure out 50 ml for each sample, but measured out volumes between 40 ml and 50 ml and recorded the mass and volume of “lake water” for each trail. The titrant used was 0.1 M HCl instead of the concentration stated in the lab manual. Additionally, when adding titrant to lake samples that had a pH higher than approximately 6.5, the titrant was added in increments larger than 0.1 mL.

##Results

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
slope, intercept, r_value, p_value, std_err = stats.linregress(titrantVolTime0[-N_good_points:],F1_data[-N_good_points:])
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

##Discussion



We defined the First Gran function as:
$F_{1}=\frac{V_{s}+V_{t}}{V_{s}}\left[ H^{+} \right ]$ The $F_{1}$ is then plotted as a function of $V_{t}$ and as a result, we gain a straight line slope = $\frac{N_{t}}{V_{s}}$.
We then are able to find the ANC value using the equation $V_{s}ANC = V_{e}N_{t}$ and solving for ANC = $\frac{V_{e}N_{t}}{V_{s}}$. This resulted in our being able to produce a Gran plot (see Results section, figure X).




## Suggestions and Comments

It would be beneficial for the data to automatically merge into one succinct Excel file.  It was tedious to keep track and organize 10 separate files from 10 trials. If we were analyzing the lake at more points in time it would become difficult to manage the data files.  We believe that it would have been beneficial to have a more effective way to clean the pH probe before trials. Though we spent a significant amount of time cleaning off the probe with a water bottle, we were never sure that the probe was completely clean. We believe that there could have been a more detailed explanation/walk through of the gran plot analysis as well as a simpler explanation of the chemistry behind ANC values and the titration process. An interesting modification could involve looking at how the titration works for a polyprotic acid and using a graph of the derivative of our titration plot to find maxima and determine equivalence points.

# Verify that your report and graphs meet the requirements as outlined in the course materials.

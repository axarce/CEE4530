#### Import statements to run code for Nutrient Removal Project.

```Python
from aide_design.play import*
from scipy import stats
import importlib
import scipy
from scipy import special
from scipy.optimize import curve_fit
import collections
import os

def Column_of_data(data_file_path,start,end,column,units):
    """ This function extracts a column of data from a ProCoDA data file.
    The file must be the original tab delimited file.
    Parameters
    ----------
    data_file_path : string of the file name or file path.
    If the file is in the working directory, then the file name is sufficient.
    Example data_file_path = 'Reactor_data.txt'
    start: index of first row of data to extract from the data file
    end: index of last row of data to extract from the data
    If the goal is to extract the data up to the end of the file use -1
    column: index of the column that you want to extract. Column 0 is time.
    The first data column is column 1.
    units: string of the units you want to apply to the data.
    Example 'mg/L'S
    If an empty string, '', is passed then no units are applied.
    Returns
    -------
    numpy array of experimental data with the units applied.
    """
    df = pd.read_csv(data_file_path,delimiter='\t')
    if units == '':
        data = np.array(pd.to_numeric(df.iloc[start:end,column]))
    else:
        data = np.array(pd.to_numeric(df.iloc[start:end,column]))*u(units)
    return data;
```

###### Anthony Arce
###### Elle Blake
###### Ajara Cobourne

#### Nutrient Removal: Final Report

## Introduction

Humans have been, and continue, to interact with the environment in harmful ways. One example is when effluent discharge with particulate solids and high levels of nutrients are mixed into bodies of water. Water that has too much nutrients can lead to several harmful outcomes such as low dissolved oxygen, dead marine animals, murky water, increased turbidity, algal blooms, and diminished populations of plants and animals that are essential to a healthy ecosystem, and the effects of BOD are affected by seasonal conditions (EPA 2001). In order to better understand the impacts that low water quality may have and the potential solutions to this problem, our group modeled a wastewater treatment plant (WWTP) through a sequencing batch reactor. The goal of a WWTP is to reduce BOD (Biological Oxygen Demand) and provide harmless effluent water to avoid reduce the impact that such wastewater has. The EPA standard is that the 30 day average 5 day BOD (BOD5) of the effluent be less than 30 mg/L. Further, the typical BOD test takes 5 days; however, given the time constraints of this experiment period, measuring the oxygen uptake rate is sufficient for our team to infer the BOD of the effluent in this experiment. Ultimately, in this experiment we will look at how to efficiently reduce BOD in effluent water by exploring the optimal oxygen flow rate during aeration. A study done in 2000 explored the relationship between airflow rates and oxygen transfer efficiency and concluded that "an increase in the air flow rate leads to a decrease in the oxygen transfer efficiency", with similar results under both wastewater treatment plant process conditions and in clean water (Gillot & Héduit, 2000).

In this experiment we model a WWTP in a small beaker containing aerobic bacteria whose substrate use can be models by the Monod relationship: $\frac{dL}{dt} = \frac{-kLX}{K_s+L}$

Assuming Ks, the half velocity constant is relatively large we can simplify the equation and integrate to obtain a new equations and say that the oxygen use rate is the same as the substrate use rate.

$L = L_o e^{-k_{ox}t}$

$\frac{dL}{dt} = -k_{ox}L = -k_{ox}L_oe^{-k_{ox}t}$

When aeration is the only thing affecting [O] we can utilize this equation.

$ln \frac{C^* - C}{C^* - C_o}$

However, the aerobic bacteria change in the oxygen deficit like this:
$\frac{\partial D}{\partial t} = k_e + k_{ox}L_oe^{-k_{ox}t}-\hat{k}_{v,l}D$
This equation holds under the assumption that we have enough oxygen to maintain microbial consumption.
As we seek to maximize OTE
We modified the oxygen transfer efficiency equation:
OTE = $\frac{\hat{k}_{v,l}(C^{*}-C)VRT}{Q_{air}P_{air}f_{O_{2}}MW_{O_{2}}}$

to account for a controlled molar airflow rate, resulting in the following equation :

$OTE = \frac{\dot{n}_{aq O_{2}}}{f_{O_{2}}\dot{n}_{air}} = \frac{V\hat{k}_{v,l}(C^{*}-C)}{f_{O_{2}}\dot{n}_{air}MW_{O_{2}}}$.

We used these concepts to investigate the relationship between OTE and air flow rates.


## Objectives
As observed in the Gas Transfer lab, as water becomes more saturated with air, the slower gas transfer occurs. An ideal reactor would maintain a very low level of oxygen in the tank, just enough to allow the microbe community to thrive, but not enough to decrease the reactor oxygen transfer efficiency. By varying the oxygen flow rate during the aeration sequence of the wastewater treatment in our sequencing batch reactor, we will attempt to determine the optimal flow rate for wastewater treatment. Ultimately we hope that our results will provide insight for larger scale wastewater treatment plants so that they are better able to improve the quality of effluent in an economically and energy efficient manner. We expect that our results will mirror of those Gillot and Héduit's 2000 study and hypothesize that as we increase the airflow rate, the oxygen transfer efficiency will decrease as well.

## Procedures
A schematic of our experimental setup is found in Figure 1.
![figure1](images\schematic.jpg)

Figure 1. Schematic of the laboratory setup for this experiment

We chose to model the WWTP as a small beaker and used ProCoDa II, a process control and data acquisition system to automate the process.
We programmed ProCoDA to do the following:
1. Fill batch reactor with wastewater (to heat depth)
2. Heat the wastewater
3. Drain excess wastewater to reach appropriate aeration water depth
4. Aerate
5. Settle
6. Drain excess wastewater so that mostly settled sludge remains for next cycle
7. (Back to #1)

Process #2 was required due to the fact that we were using synthetic wastewater to operate our WWTP. If left at room temperature the wastewater would degrade, but we wanted room temperature wastewater to mimic the conditions of a real WWTP. The water did not heat up to room temperature within a reasonable time period so we added a heater to ProCoDA. We added a temperature probe and proportional–integral–derivative controller to achieve room temperature wastewater and because the heater required a minimum water level to function well, we added sequence #3, “drain to aeration depth”, so that during aeration the wastewater would not bubble over onto the lab bench.

To initialize the first trial of the experiment, 400 ml of aerobic sludge was added to the sequencing batch reactor. The sludge was then aerated and then let it settle and continue from that step in the cycle. Cycle aeration times were programmed for 4 hours and the following air flow rates were tested: 250, 500, 850, 1000 and 3000 $\frac{\mu L}{s}$


## Results and Discussion
The code written in the Appendix of this report provided the following results from our experimental data.

In obtaining the oxygen transfer efficiencies (OTE), a representative subset of each dissolved oxygen timeseries data is obtained by examining the dataset and determining an aeration period from minimal dissolved oxygen at 0 mg/L to about 4 mg/L. Looking at the subset of data, linear regression is done on a log transformation of the dataset relative to the saturation concentration of dissolved oxygen to determine the slope of aeration to obtain k,vi values in which to calculate OTE values.

The log transformation equation used to determine kv,i when applying linear regression with respect to time:
$ln \frac{C^* - C}{C^* - C_o}$




| Oxygen Transfer Efficiency | $\frac{\mu L}{s}$ |
| ---------------- | ----- |
| $OTE_{250}$    | 2.15e-03   |
| $OTE_{500}$    | 5.651e-03    |
| $OTE_{850}$    | 1.874e-03  |
| $OTE_{1000}$   | 1.77e-03   |
| $OTE_{3000}$   | 4.141e-05  |
| $OTE_{10000}$  | 1.348e-05      |


The results provide evidence that match the expectations of the theory; as the oxygen deficit decreases the rate of gas transfer decreases and becomes less efficient. The OTE exponentially decreased with airflow rate which implies one should try to keep a low airflow rate and sufficiently large oxygen deficit to maintain a high functioning wastewater treatment plant. The 250 micromolar flow rate used as the original testing value did not match the expected trend of oxygen transfer efficiency, but both samples of dissolved oxygen for the 250 micromolar flow rate presented difficult challenges in that negative values and readings far above saturation were present throughout both trials. The trend overall is consistent with gas transfer theory.


![figure2](images\OTE.jpg)

Figure 2. Oxygen transfer efficiencies for flow rates of 250, 500, 850, 1000, 3000 and 10000


## Suggestions and Comments
If we had to do this over again, in terms of equipment we would like to have more reliable/consistent equipment. Some of the equipment proved to be unreliable, particularly the DO probe. At times we had readings of 240 mg O2/L so we had to frequently adjust.

A low-cost solution to any problems with dissolved oxygen probes without necessarily having to upgrade could be as simple as attaching multiple, 2 or 3, DO probes to the simulated wastewater treatment reactor at different points in the reactor in order to reduce the chances of failure of dissolved oxygen readings overnight. By examining the DO readings for each of the probes, groups can make rational judgments as to what to interpret from the multiple readings, whether one was giving numbers beyond or under physical limits an

Perhaps because our tank was so small we lost quite a bit of sludge during each trial, despite a 45 minute settling period. We did not have to stir our tank because it was so small, which is a good thing in terms of energy, however having to refrigerate and heat the WW was an energy intensive process). The purpose of this project is to address real world concerns and what was effective in our 1L tank may not be effective in treating the 6.5 million gallons of sewage per day, coming from Ithaca every day.

The final presentation of four Nutrient Removal Projects was redundant. I think we could have achieved more variety if we had a CEE Lab material Walkthrough Day about a week before the project proposal is due (and the AguaClara Lab-if it seems necessary). It is daunting to not know what tools are available, especially for people who's project team is not Lab-based. We ended up not fleshing out the creative ideas that might have yielded more interesting results.

### References/Bibliography
1.https://www.epa.gov/sites/production/files/2014-08/documents/nutrient-memo-nov142001.pdf  
2.https://www.sciencedirect.com/science/article/pii/S0043135499003231
3.http://www.ingentaconnect.com/contentone/wef/wer/1999/00000071/00000004/art00012
4.https://www.sciencedirect.com/science/article/pii/S0960852417320825


### Appendix: Python Code & Original Proposal


(1) Python code used to perform analysis on our experimental data:
```Python
cd C:\Users\Anthony\github\CEE4530_axa2\Final Project Files


# Extracting Dissolved Oxygen Data and Timeseries Data from the ProCoDa files
|DO_250Day1 = Column_of_data('250flowDay1.xls', 0, -1, 5, 'mg/L')
DO_250Day2 = Column_of_data('250FlowDay2.xls', 0, -1, 5, 'mg/L')
DO_500Day3 = Column_of_data('500flow.xls', 0, -1, 5, 'mg/L')
DO_850Day4 = Column_of_data('850Flow.xls', 0, -1, 5, 'mg/L')
DO_MaxDay5 = Column_of_data('fullflow.xls', 0, -1, 5, 'mg/L')
DO_1000Day6 = Column_of_data('1000flow.xls', 0, -1, 5, 'mg/L')
DO_3000Day7 = Column_of_data('3000flow.xls', 0, -1, 5, 'mg/L')

Time_Day1 = Column_of_data('250flowDay1.xls', 0, -1, 1, 's')
Time_Day2 = Column_of_data('250FlowDay2.xls', 0, -1, 1, 's')
Time_Day3 = Column_of_data('500flow.xls', 0, -1, 1, 's')
Time_Day4 = Column_of_data('850Flow.xls', 0, -1, 1, 's')
Time_Day5 = Column_of_data('fullflow.xls', 0, -1, 1, 's')
Time_Day6 = Column_of_data('1000flow.xls', 0, -1, 1, 's')
Time_Day7 = Column_of_data('3000flow.xls', 0, -1, 1, 's')

# Plotting of the particularly poor DO data for a 250 umol day to illustrate error.
DO_plot_1 = plt.plot(Time_Day1.to(u.hr), DO_250Day1.to(u.mg/u.L), 'ro')
plt.xlabel(r'$time (hr)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.title('250 umol Airflow')
plt.savefig(r'images\250Trail1.jpg')
plt.show()

# The second trial of a 250
DO_plot_2 = plt.plot(Time_Day2.to(u.hr), DO_250Day2.to(u.mg/u.L), 'ro')
plt.xlabel(r'$time (hr)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.title('250 umol Airflow')
plt.savefig(r'images\250Trail2.jpg')
plt.show()

temp = 293
P_O2 = 0.21
C_star = (P_O2*np.exp((1727/temp)-2.105))*u.mg/u.L
print(C_star)

def C_ratio(C_Star, DO, Time):
  C_ratio = np.zeros(len(Time))
  for i in range(0, len(Time)):
    if DO[i] > C_Star:
      DO[i] = C_star - 0.01 * u.mg/u.L
    if DO[i] < 0 * u.mg/u.L:
      DO[i] = 0 * u.mg/u.L
    C_ratio[i] = np.log((C_Star - DO[i])/(C_Star))
    i = i + 1
  return C_ratio

# Upon examination of datasets for any point where dissolved oxygen dropped to 0, we examine aeration from a 0 point to about 4 mg/L to get a fair representation of kv,i values for each flow rate.
C_ratio_250_1 = C_ratio(C_star, DO_250Day1, Time_Day1)
C_ratio_250_2 = C_ratio(C_star, DO_250Day2[625:795], Time_Day2[625:795])
C_ratio_500 = C_ratio(C_star, DO_500Day3[599:655], Time_Day3[599:655])
C_ratio_850 = C_ratio(C_star, DO_850Day4[742:806], Time_Day4[742:806])
C_ratio_max = C_ratio(C_star, DO_MaxDay5[3429:3720], Time_Day5[3429:3720])
C_ratio_1000 = C_ratio(C_star, DO_1000Day6[645:691], Time_Day6[645:691])
C_ratio_3000 = C_ratio(C_star, DO_3000Day7[319:421], Time_Day7[319:421])

kv_250, intercept_250, r_value_250, p_value, std_err = stats.linregress(Time_Day2[625:795], C_ratio_250_2)
kv_500, intercept_500, r_value_30, p_value, std_err = stats.linregress(Time_Day3[599:655], C_ratio_500)
kv_850, intercept_850, r_value_70, p_value, std_err = stats.linregress(Time_Day4[742:806], C_ratio_850)
kv_max, intercept_max, r_value_160, p_value, std_err = stats.linregress(Time_Day5[3429:3720], C_ratio_max)
kv_1000, intercept_1000, r_value_370, p_value, std_err = stats.linregress(Time_Day6[645:691], C_ratio_1000)
kv_3000, intercept_3000, r_value_850, p_value, std_err = stats.linregress(Time_Day7[319:421], C_ratio_3000)


k_values = np.array([kv_250, kv_500, kv_850, kv_1000, kv_3000, kv_max])*-1/u.s
k_values
flow_values = np.array([250, 500, 850, 1000, 3000, 10000])*u.umol/u.s

oxygen_deficit = 6*u.mg/u.L
Volume = (11/16) * u.L
f_O2 = 0.21
MW_O2 = 32 * u.g/u.mol

OTE = np.zeros(len(flow_values))
for i in range(0, len(k_values)):
  OTE[i] = ((Volume * (k_values[i]) * oxygen_deficit) / (f_O2 * (flow_values[i]) * MW_O2))
  i = i + 1


plt.plot(flow_values, OTE, 'ro')
plt.xlabel('Air Flow (umol / s)')
plt.ylabel('OTE')
plt.title('Oxygen transfer efficiency as a function of Air Flow Rate')
plt.savefig(r'images\OTE.jpg')
plt.show()

```


(2) Original Research proposal

#### Introduction

>In the anthropocene humans have been, and continue to, interact with the environment in harmful ways. One example is when effluent discharge with particulate solids and high levels of nutrients mixed into bodies of water. Water that has too much nutrients can lead to several harmful outcomes such as low dissolved oxygen, dead fish, murky water, increased turbidity, algal blooms, and diminished populations of plants and animals that are essential to a healthy ecosystem, and the effects of BOD are affected by seasonal conditions (EPA 2001). A good way to examine potential solutions to this problem is through a model wastewater treatment plant (WWTP). The goal of a WWTP is to reduce BOD (Biological Oxygen Demand) and provide harmless effluent water to avoid harming the environment. The EPA standard is that the 30 day average 5 day BOD (BOD5) of the effluent be less than 30 mg/L. The typical BOD test takes 5 days, but measuring the oxygen uptake rate is sufficient for us to infer the BOD of the effluent in this experiment. In this experiment we will look at how to reduce BOD in effluent water.


#### Objectives

>As observed in the Gas Transfer lab, as water becomes more saturated with ~~water~~, the slower gas transfer occurs. An ideal reactor would maintain a very low level on oxygen in the tank, just enough to allow the microbe community to thrive, but not enough to decrease the reactor efficiency. We will try to run our reactor at constant high efficiency and will vary temperature and observe how that affects the BOD effluent. We hypothesize that BOD and temperature will vary directly. We will estimate BOD removal by measuring the oxygen uptake rate. If we increase the concentration of oxygen in the reactor tank, cut off the airflow and monitor how much the oxygen concentration in our batch reactor changes we can obtain the oxygen uptake rate.


>•Key design parameters

>-Flow rates -  0, 160, 850 μM/s.

>-Volumes - 600ml

>-Concentrations - diluted stock concentration of synthetic wastewater

>-Range of parameters that you are varying - Temperature: 275, 295, 315 K

#### Timeline

>4/11/18: Goals - Set up experimental apparatus, including ~~ProCoDa~~ ProCoDA (Process Control and Data Acquisition) and a test run with tap water. If time allows run 1-2 trials with wastewater

>4/18/18: Goals - Run 6 trials with synthetic wastewater

>4/25/18: Goals - Gather additional data as needed and begin analysis

>5/02/18: Goals - Complete analysis

>5/09/18: Goal: Present final results in class

>One challenge we might face when doing this experiment might be getting a successful reactor setup. Maintaining a constant temperature might also be difficult. There will not be a TA bench model and we will have to do more trial and error than we are used to. Monroe said that this was a notoriously difficult set up in the past, so we may not be able to collect data quickly. Additionally we will have to modify some of the equipment in the lab. We do not have any small beakers that have enough spouts to correctly model a batch reactor so modifying what we do have may prove to be a challenge.

#### Expectations
>A 1999 study by Griffin, Bhattarai and Xiang found that there was a significant difference between BOD levels in effluent water at different temperatures. They observed a seasonal variation where at temperatures less than 293 K, BOD removal had more variation and was less effective. However at temperatures above 293 K more BOD removal occurred. In February 2018 Waki, Yasuda, Fukumoto, Béline, and Magrí published the results of a study where wastewater was treated at either a continuously low DO, or a continuously high DO level at temperatures ranging from 283 to 303 K. They found that nutrient and BOD removal were highest at low DO levels and high temperatures. We expect to get similar results to these studies.


>•Resources needed to conduct experiments – What tools will you use?

>-1 L modified reactor

>-ProCoDa II

>-Sensors - (pressure, temperature, dissolved oxygen, flow rate)

>-Synthetic Wastewater

>-Activated Sludge

>-A filter

>-A scale

>-A warmer

>-Ice

>-Glass beaker

>-Magnetic Stirrer

>-Accumulator

>-Needle Valve

>-Solenoid Valve

>-Air Supply

#### References/Bibliography

>1.https://www.epa.gov/sites/production/files/2014-08/documents/nutrient-memo-nov142001.pdf


>2.http://www.ingentaconnect.com/contentone/wef/wer/1999/00000071/00000004/art00012


>3.https://www.sciencedirect.com/science/article/pii/S0960852417320825

#### Some Relevant Equations:
>$\frac{dL}{dt} = \frac{-kLX}{K_s+L}$  the Monod relationship, substrate use by bacteria
Assuming Ks, the half velocity constant is relatively large we can simplify the equation and integrate to obtain the following:
$L = L_o e^{-k_{ox}t}$
We can say that the oxygen use rate is the same as the substrate use rate so

>$\frac{dL}{dt} = -k_{ox}L = -k_{ox}L_oe^{-k_{ox}t}$

>$ln \frac{C^* - C}{C^* - C_o}$ when aeration is the only thing affecting [O]

>The change in the oxygen deficit: $\frac{\partial D}{\partial t} = k_e + k_{ox}L_oe^{-k_{ox}t}-\hat{k}_{v,l}D$
. This equation holds under the assumption that we have enough oxygen to maintain microbial consumption.

```Python
from aide_design.play import*
from scipy import stats
import importlib
import scipy
from scipy import special
from scipy.optimize import curve_fit
import collections
import os
import tkinter as tk
from tkinter import filedialog
root = tk.Tk()
root.withdraw()

def aeration_data(DO_column):
    """ This function extracts the data from folder containing tab delimited files of aeration data.
    The file must be the original tab delimited file.
    All text strings below the header must be removed from these files.
    The file names must be the air flow rates with units of micromoles/s.
    An example file name would be "300.xls" where 300 is the flowr ate in micromoles/s
    The function opens a file dialog for the user to select the directory containing the data.
    Parameters
    ----------
    DO_column: index of the column that contains the dissolved oxygen concentration data.
    Returns
    -------
    filepaths: list of all file paths in the directory sorted by flow rate
    airflows: sorted numpy array of air flow rates with units of micromole/s attached
    DO_data: sorted list of numpy arrays. Thus each of the numpy data arrays can have different lengths to accommodate short and long experiments
    time_data: sorted list of numpy arrays containing the times with units of seconds. Each
    """

    dirpath = filedialog.askdirectory()
    #return the list of files in the directory
    filenames = os.listdir(dirpath)
    #extract the flowrates from the filenames and apply units
    airflows=((np.array([i.split('.', 1)[0] for i in filenames])).astype(np.float32))
    #sort airflows and filenames so that they are in ascending order of flow rates
    idx   = np.argsort(airflows)
    airflows = (np.array(airflows)[idx])*u.umole/u.s
    filenames = np.array(filenames)[idx]

    filepaths = [os.path.join(dirpath, i) for i in filenames]
    #DO_data is a list of numpy arrays. Thus each of the numpy data arrays can have different lengths to accommodate short and long experiments
    # cycle through all of the files and extract the column of data with oxygen concentrations and the times
    DO_data=[Column_of_data(i,0,-1,DO_column,'mg/L') for i in filepaths]
    time_data=[(ftime(i,0,-1)).to(u.s) for i in filepaths]
    aeration_collection = collections.namedtuple('aeration_results','filepaths airflows DO_data time_data')
    aeration_results = aeration_collection(filepaths, airflows, DO_data, time_data)
    return aeration_results

def O2_sat(Pressure_air,Temperature):
    """
    This equation is valid for 278 K < T < 318 K
    Parameters
    ----------
    Pressure_air: air pressure with appropriate units.
    Temperature: water temperature with appropriate units
    Returns
    -------
    Saturated oxygen concentration in mg/L
    """
    fraction_O2 = 0.21
    Pressure_O2 = Pressure_air *fraction_O2
    return (Pressure_O2.to(u.atm).magnitude)*u.mg/u.L*np.exp(1727/Temperature.to(u.K).magnitude - 2.105)

def Gran(data_file_path):
    """ This function extracts the data from a ProCoDA Gran plot file.
    The file must be the original tab delimited file.
    Parameters
    ----------
    data_file_path : string of the file name or file path.
    If the file is in the working directory, then the file name is sufficient.
    Example data_file_path = 'Reactor_data.txt'
    Returns
    -------
    V_titrant (mL) as numpy array
    ph_data as numpy array (no units)
    V_sample (mL) volume of the original sample that was titrated
    Normality_titrant (mole/L) normality of the acid used to titrate the sample
    V_equivalent (mL) volume of acid required to consume all of the ANC
    ANC (mole/L) Acid Neutralizing Capacity of the sample
    """
    df = pd.read_csv(data_file_path,delimiter='\t',header=5)
    V_t = np.array(pd.to_numeric(df.iloc[0:,0]))*u.mL
    pH = np.array(pd.to_numeric(df.iloc[0:,1]))
    df = pd.read_csv(data_file_path,delimiter='\t',header=-1,nrows=5)
    V_S = pd.to_numeric(df.iloc[0,1])*u.mL
    N_t = pd.to_numeric(df.iloc[1,1])*u.mole/u.L
    V_eq = pd.to_numeric(df.iloc[2,1])*u.mL
    ANC_sample = pd.to_numeric(df.iloc[3,1])*u.mole/u.L
    Gran_collection = collections.namedtuple('Gran_results','V_titrant ph_data V_sample Normality_titrant V_equivalent ANC')
    Gran = Gran_collection(V_titrant=V_t, ph_data=pH,V_sample=V_S, Normality_titrant=N_t, V_equivalent=V_eq, ANC=ANC_sample )
    return Gran;


def ftime(data_file_path,start,end):
    """ This function extracts the column of times from a ProCoDA data file.
    The file must be the original tab delimited file.
    Parameters
    ----------
    data_file_path : string of the file name or file path.
    If the file is in the working directory, then the file name is sufficient.
    Example data_file_path = 'Reactor_data.txt'
    start: index of first row of data to extract from the data file
    end: index of last row of data to extract from the data
    If the goal is to extract the data up to the end of the file use -1
    Returns
    -------
    numpy array of experimental times starting at 0 day with units of days.
    """
    df = pd.read_csv(data_file_path,delimiter='\t')
    start_time = pd.to_numeric(df.iloc[start,0])*u.day
    day_times = pd.to_numeric(df.iloc[start:end,0])
    time_data = np.subtract((np.array(day_times)*u.day),start_time)
    return time_data;

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

def notes(data_file_path):
    """This function extracts any experimental notes from a ProCoDA data file.
    The file must be the original tab delimited file.
    Parameters
    ----------
    data_file_path : string of the file name or file path.
    If the file is in the working directory, then the file name is sufficient.
    Example data_file_path = 'Reactor_data.txt'
    Returns
    -------
    dataframe showing the rows of the data file that contain text notes
    inserted during the experiment.
    Use this to identify the section of the data file that you want to extract.
    """
    df = pd.read_csv(data_file_path,delimiter='\t')
    text_row = df.iloc[0:-1,0].str.contains('[a-z]','[A-Z]')
    text_row_index = text_row.index[text_row == True].tolist()
    notes = df.loc[text_row_index]
    return notes


#carbonates
#The following code defines the carbonate system and provides functions for calculating Acid Neutralizing Capacity.
Kw = 10**(-14) * (u.mole/u.L)**2
K1_carbonate = 10**(-6.37)*u.mol/u.L
K2_carbonate = 10**(-10.25)*u.mol/u.L
K_Henry_CO2 = 10**(-1.5) * u.mole/(u.L*u.atm)
P_CO2 = 10**(-3.5) * u.atm

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

def ANC_closed(pH,Total_Carbonates):
    return Total_Carbonates*(alpha1_carbonate(pH)+2*alpha2_carbonate(pH)) + Kw/invpH(pH) - invpH(pH)

def ANC_open(pH):
    return ANC_closed(pH,P_CO2*K_Henry_CO2/alpha0_carbonate(pH))

# Reactors
# The following code is for reactor responses to tracer inputs.
def CMFR(t,C_initial,C_influent):
    """ This function calculates the output concentration of a completely mixed flow reactor given an influent and initial concentration.
    Parameters
    ----------
    C_initial : The concentration in the CMFR at time zero.
    C_influent : The concentration entering the CMFR.
    t: time made dimensionless by dividing by the residence time of the CMFR. t can be a single value or a numpy array.
    Returns
    -------
    Effluent concentration
    """
    return C_influent * (1-np.exp(-t)) + C_initial*np.exp(-t)

def E_CMFR_N(t, N):
    """ This function calculates a dimensionless measure of the output tracer concentration from a spike input to a series of completely mixed flow reactors.
    Parameters
    ----------
    t: time made dimensionless by dividing by the residence time of the reactor. t can be a single value or a numpy array.
    N : The number of completely mixed flow reactors (CMFR) in series. This would logically be constrained to real numbers greater than 1.
    Returns
    -------
    (Concentration * volume of 1 CMFR) / (mass oftracer)
    """
    #make sure that time is dimensionless and not a mixed time unit
    if hasattr(t, 'magnitude'):
      t.ito(u.dimensionless)
    return (N**N)/special.gamma(N) * (t**(N-1))*np.exp(-N*t)

def E_Advective_Dispersion(t, Pe):
    """ This function calculates a dimensionless measure of the output tracer concentration from a spike input to reactor with advection and dispersion.
    Parameters
    ----------
    t: time made dimensionless by dividing by the reactor residence time. t can be a single value or a numpy array.
    Pe : The ratio of advection to dispersion ((mean fluid velocity)/(Dispersion*flow path length))
    Returns
    -------
    (Concentration * volume of reactor) / (mass of tracer)
    """
    #make sure that time is dimensionless and not a mixed time unit
    if hasattr(t, 'magnitude'):
      t.ito(u.dimensionless)
    #replace any times at zero with a number VERY close to zero to avoid divide by zero errors

    t[t==0]=10**(-50)
    return (Pe/(4*np.pi*t))**(0.5)*np.exp((-Pe*((1-t)**2))/(4*t))

def Tracer_CMFR_N(t_seconds, t_bar, C_bar, N):
    """ Used by Solver_CMFR_N. All inputs and outputs are unitless.
    This is The model function, f(x, ...). It takes the independent variable as the first argument and the parameters to fit as separate remaining arguments.
    Parameters
    ----------
    t_seconds : Array of times (units of seconds, but unitless)
    t_bar : Average time spent in the total reactor (units of seconds, but unitless).
    C_bar : (Mass of tracer)/(volume of the total reactor) unitless.
    N : The number of completely mixed flow reactors (CMFR) in series. This would logically be constrained to real numbers greater than 1.
    Returns
    -------
    (C_bar*E_CMFR_N(t_seconds/t_bar, N))
    The model concentration as a function of time
    """

    return C_bar*E_CMFR_N(t_seconds/t_bar, N)

def Solver_CMFR_N(t_data, C_data, theta_guess, C_bar_guess):
    """ Use non-linear least squares to fit the function, Tracer_CMFR_N(t_seconds, t_bar, C_bar, N), to reactor data.
    Parameters
    ----------
    t_data : Array of times with units
    C_data : Array of tracer concentration data with units
    theta_guess : Estimate of time spent in the total reactor with units.
    C_bar_guess : Estimate of (Mass of tracer)/(volume of the total reactor) with units.
    Returns
    -------
    a tuple with theta (units of s), C_bar (same units as C_bar_guess), and N as the best fit to the data.
    """

    C_unitless = C_data.magnitude
    C_units = str(C_bar_guess.units)
    t_seconds = (t_data.to(u.s)).magnitude
    # assume that a guess of 1 reactor in series is close enough to get a solution
    p0 = [theta_guess.to(u.s).magnitude, C_bar_guess.magnitude,1]
    popt, pcov = curve_fit(Tracer_CMFR_N, t_seconds, C_unitless, p0)
    Solver_theta = popt[0]*u.s
    Solver_C_bar = popt[1]*u(C_units)
    Solver_N = popt[2]
    Reactor_results = collections.namedtuple('Reactor_results','theta C_bar N')
    CMFR = Reactor_results(theta=Solver_theta, C_bar = Solver_C_bar, N = Solver_N)
    return CMFR


def Tracer_AD_Pe(t_seconds, t_bar, C_bar, Pe):
    """ Used by Solver_AD_Pe. All inputs and outputs are unitless.
    This is The model function, f(x, ...). It takes the independent variable as the first argument and the parameters to fit as separate remaining arguments.
    Parameters
    ----------
    t_seconds : Array of times (units of seconds, but unitless)
    t_bar : Average time spent in the reactor (units of seconds, but unitless).
    C_bar : (Mass of tracer)/(volume of the reactor) unitless.
    Pe : The Peclet number for the reactor.
    Returns
    -------
    C_bar*E_Advective_Dispersion(t_seconds/t_bar, Pe)
    The model concentration as a function of time
    """

    return C_bar*E_Advective_Dispersion(t_seconds/t_bar, Pe)

def Solver_AD_Pe(t_data, C_data, theta_guess, C_bar_guess):
    """ Use non-linear least squares to fit the function, Tracer_AD_Pe(t_seconds, t_bar, C_bar, Pe), to reactor data.
    Parameters
    ----------
    t_data : Array of times with units
    C_data : Array of tracer concentration data with units
    theta_guess : Estimate of time spent in one CMFR with units.
    C_bar_guess : Estimate of (Mass of tracer)/(volume of one CMFR) with units.
    Returns
    -------
    a tuple with theta (units of s), C_bar (same units as C_bar_guess), and Pe as the best fit to the data.
    """

    C_unitless = C_data.magnitude
    C_units = str(C_bar_guess.units)
    t_seconds = (t_data.to(u.s)).magnitude
    # assume that a guess of 1 reactor in series is close enough to get a solution
    p0 = [theta_guess.to(u.s).magnitude, C_bar_guess.magnitude,5]
    popt, pcov = curve_fit(Tracer_AD_Pe, t_seconds, C_unitless, p0, bounds=(0,inf))
    Solver_theta = popt[0]*u.s
    Solver_C_bar = popt[1]*u(C_units)
    Solver_Pe = popt[2]
    Reactor_results = collections.namedtuple('Reactor_results','theta C_bar Pe')
    AD = Reactor_results(theta=Solver_theta, C_bar = Solver_C_bar, Pe = Solver_Pe)
    return AD
```
###### Anthony Arce
###### Elle Blake
###### Ajara Cobourne

#### Data Analysis for Lab 6: Gas Transfer

## Introduction and Objectives

We decided to conduct this experiment to examine and understand the relationship between gas transfer and flow rates. Traditionally the oxygen transfer efficiency is quite low. We hoped to gain a clearer understanding of gas exchange through this experiment as it is useful in wastewater treatment plant optimization. Wastewater treatment plants require consistent aerobic degradation, which can be maintained through oxygen transfer into activated sludge tanks. By doing this experiment, we are better able to understand and make decisions regarding the most useful technologies used to enhance such gas transfer by creating a high interface surface area.

## Procedures
We followed the procedures as stated in the CEE 4530 Spring 2018 lab manual with a few modifications. Instead of letting the accumulator pressure get to 75%, we set the needle valve to approximately 50% (20,000 Pa). To deoxygenate the 600 ml solution between trials we added 0.4 ml of 100 mg/ml of sodium sulfite and we conducted trials at air flow rates of approximately 0, 30, 70, 160, 370, 850 and 2000 μM/s.


## Data Analysis
```python
cd C:\Users\Anthony\github\CEE4530_axa2\Gas Transfer Lab\datafiles
#1 Air Flow rate (test @ 200microM/s)
R = 8314 * u.L * u.Pa / (u.mol * u.K)#for Pa, not kPa
T = 295*u.K
V = 0.6*u.L
Pressure_200 = Column_of_data('200.txt', 0, -1, 2, 'Pa')
Time_200 = ftime('200.txt', 0, -1).to(u.s)
slope, intercept, r_value, p_value, std_err = stats.linregress(Time_200, Pressure_200)
PressureSlope = slope * u.Pa / u.s

N = PressureSlope * V / (R * T)


#3 Plot the representative data set showing dissolved oxygen vs. time
#Import the dissolved oxygen data from each respective trial
cd C:\Users\Anthony\github\CEE4530_axa2\Gas Transfer Lab\datafiles



Flow0_DO = Column_of_data('0.txt',0,-1,2,'mg/L')
Flow30_DO = Column_of_data('30.txt',0,-1,2,'mg/L')
Flow70_DO = Column_of_data('70.txt',0,-1,2,'mg/L')
Flow160_DO = Column_of_data('160.txt',0,-1,2,'mg/L')
Flow370_DO = Column_of_data('370.txt',0,-1,2,'mg/L')
Flow850_DO = Column_of_data('850.txt',0,-1,2,'mg/L')
Flow2000_DO = Column_of_data('2000.txt',0,-1,2,'mg/L')

#Import the time data from each respective trial
Flow0_time = ftime('0.txt',0,-1).to(u.s)
Flow30_time = ftime('30.txt',0,-1).to(u.s)
Flow70_time = ftime('70.txt',0,-1).to(u.s)
Flow160_time = ftime('160.txt',0,-1).to(u.s)
Flow370_time = ftime('370.txt',0,-1).to(u.s)
Flow850_time = ftime('850.txt',0,-1).to(u.s)
Flow2000_time = ftime('2000.txt',0,-1).to(u.s)

#3 and 7: Plot the representative data set W/reaeration model
cd  C:\Users\Anthony\github\CEE4530_axa2\Gas Transfer Lab
DO_plot = plt.plot(Flow160_time.to(u.min), Flow160_DO.to(u.mg/u.L), 'ro')
Time_Min = Flow160_time.to(u.min)

plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.title('Dissolved Oxygen Vs. Time of 160 Flow Trial')
plt.savefig(r'images\160Trial.jpg')
plt.show()

# 4: Calculate C*
temp = 295
P_O2 = 0.21
C_star = (P_O2*np.exp((1727/temp)-2.105))*u.mg/u.L
                      print(C_star)

def C_ratio(C_Star, DO, Time):
  C_ratio = np.zeros(len(Time))
  for i in range(0, len(Time)):
    C_ratio[i] = np.log((C_Star - DO[i])/(C_Star - DO[0]))
    i = i + 1
  return C_ratio

C_ratio_0 = C_ratio(C_star, Flow0_DO, Flow0_time)
C_ratio_30 = C_ratio(C_star, Flow30_DO, Flow30_time)
C_ratio_70 = C_ratio(C_star, Flow70_DO, Flow70_time)
C_ratio_160 = C_ratio(C_star, Flow160_DO, Flow160_time)
C_ratio_370 = C_ratio(C_star, Flow370_DO, Flow370_time)
C_ratio_850 = C_ratio(C_star, Flow850_DO, Flow850_time)
C_ratio_2000 = C_ratio(C_star, Flow2000_DO, Flow2000_time)


# 5: Estimate K,vi using linear regression
kv_0, intercept_0, r_value_0, p_value, std_err = stats.linregress(Flow0_time, C_ratio_0)
kv_30, intercept_30, r_value_30, p_value, std_err = stats.linregress(Flow30_time, C_ratio_30)
kv_70, intercept_70, r_value_70, p_value, std_err = stats.linregress(Flow70_time, C_ratio_70)
kv_160, intercept_160, r_value_160, p_value, std_err = stats.linregress(Flow160_time, C_ratio_160)
kv_370, intercept_370, r_value_370, p_value, std_err = stats.linregress(Flow370_time, C_ratio_370)
kv_850, intercept_850, r_value_850, p_value, std_err = stats.linregress(Flow850_time, C_ratio_850)
kv_2000, intercept_2000, r_value_2000, p_value, std_err = stats.linregress(Flow2000_time, C_ratio_2000)


#6: Create graph of C_ratio vs time with lin regression

line = np.zeros(len(Flow70_time))
for i in range(0,len(Flow70_time)):
  line[i] = kv_70 * (Flow70_time[i]).magnitude + intercept_70
  i = i + 1

C_ratio_Plot = plt.plot(Flow70_time.to(u.min), C_ratio_70, Flow70_time.to(u.min), line)
plt.xlabel('Time (min)')
plt.title('ln[(C* - C)/(C* - C0) of 70umol/s flow')
plt.legend(['C_ratio of Data', 'Linear Fit'])
plt.savefig(r'images\Cratio_70.jpg')
plt.show()

#7: C(t) can be modeled as C(t) = C_star - (C_star - C_0)exp(-kt)
# For the representative dataset of 160mmol/s:
C_160 = np.zeros(len(Flow160_time)) * u.mg / u.L
C0 = 0.5 * u.mg / u.L
for i in range(0, len(Flow160_time)):
  C_160[i] = C_star - (C_star - C0)*np.exp(kv_160/u.s*Flow160_time[i])
  i = i + 1

DO_plot = plt.plot(Flow160_time.to(u.min), Flow160_DO.to(u.mg/u.L), 'ro', Flow160_time.to(u.min), C_160, 'b-')


plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.title('Dissolved Oxygen Vs. Time of 160 Flow Trial')
plt.legend(['DO Data', 'Reaeration Model'])
plt.savefig(r'images\160Trial.jpg')
plt.show()


#8 Plot Kv,i as a function of airflow rate

k_values = np.array([kv_0, kv_30, kv_70, kv_160, kv_370, kv_850, kv_2000])*-1/u.s
k_values
flow_values = np.array([0, 30, 70, 160, 370, 850, 2000])*u.umol/u.s

k_plot = plt.plot(flow_values[1:], k_values[1:])
plt.xlabel('Air Flow (umol / s)')
plt.ylabel('K,vi')
plt.title('K,vi as a function of Air Flow')
plt.savefig(r'images\kvalues.jpg')
plt.show()

#9
Flow160 = np.array(Flow160_DO)
Flow_edit=np.array(Flow160)
Flow_edit=Flow_edit[12:18]

Time0 = np.array(Flow160_time)
Time_edit=np.array(Time0)
Time_edit=Time_edit[12:18]

plt.plot(Time_edit, Flow_edit)
plt.xlabel('Time (sec)')
plt.ylabel('DO mg/L')
plt.title('Subset of aeration data to get a better linear growth')
plt.savefig(r'images\Linearize.jpg')
plt.show()

#P10. lot OTE as a function of airlow rate
oxygen_deficit = 6*u.mg/u.L
Volume = 0.5 * u.L
f_O2 = 0.21
MW_O2 = 32 * u.g/u.mol

OTE = np.zeros(len(flow_values))
for i in range(1, len(k_values)):
  OTE[i] = ((Volume * (k_values[i]) * oxygen_deficit) / (f_O2 * (flow_values[i]) * MW_O2))
  i = i + 1

plt.plot(flow_values[1:], OTE[1:])
plt.xlabel('Air Flow (umol / s)')
plt.ylabel('OTE')
plt.title('Oxygen transfer efficiency as a function of Air Flow Rate')
plt.savefig(r'images\OTE.jpg')
plt.show()


#plot the molar rate of oxygen dissolution into the aqueous phase

molar_rate = np.zeros(len(flow_values))*u.umol/u.s
for i in range(1, len(k_values)):
  molar_rate[i] = ((Volume / MW_O2) * k_values[i] * oxygen_deficit)
  i = i + 1
molar_rate

flow_values = flow_values[1:]
molar_rate = molar_rate[1:]
molar_rate
plt.plot(flow_values, molar_rate)
plt.xlabel('Air Flow (umol / s)')
plt.ylabel('Molar Transfer umol / s')
plt.title('Molar O2  as a function of Air Flow')
plt.savefig(r'images\O2transfer.jpg')
plt.show()

```

## Results and Discussion

In order to calculate C* with our data, we used the formula:

$P_{O_{2}}*e^{(\frac{1727}{T}-2.105)}$ =  $C^{*}$

We used linear regression to evaluate equation 1.5:
$ln\frac{C^{*}-C}{C^{*}-C_{0}} = -{k}_{v,l}(t-t_{0})$ in order to get ${k}_{v,l}$ as the slope of the line. Given the apparatus of our experiment, we determined that this simple gas transfer model is appropriate as the gas transfer coefficient is independent of the dissolved oxygen concentration.

!(https://github.com/hispanicberniebro/CEE4530/blob/master/Gas%20Transfer%20Lab/images/160Trial.jpg)
Figure 1: Dissolved Oxygen v. Time and Reaeration Model v. time

Figure 1b: Dissolved Oxygen v. Time (data subset to obtain more linear plot)

The representative of the reaeration model is slightly concave while both the experimental data and the model demonstrate positive slopes. The oxygen deficit is decreasing logarithmically (Figure 1). The slope of the linearized data line is the ${k}_{v,l}$ value. A more linear subset of this data is shown in Figure 1b where residual sulfite may have prevented oxygen level from changing at their true rate.

!(https://github.com/hispanicberniebro/CEE4530/blob/master/Gas%20Transfer%20Lab/images/Cratio_70.jpg)
Figure 2: Linearized Data v. Time

The representative plot of the linearized data versus time demonstrates a downward sloping line, as such the oxygen deficit is decreasing logarithmically (Figure 2). The slope of the linearized data line is the ${k}_{v,l}$ value (Table 1).

!(https://github.com/hispanicberniebro/CEE4530/blob/master/Gas%20Transfer%20Lab/images/kvalues.jpg)
Figure 3: k as a function of airflow rate

The plot of the k value as a function of time demonstrates a sharps peak near the beginning of the trial and a smaller peak around an airflow of 350 micro moles per second ( Figure 3). Afterwards our k values demonstrated a declining k value trend as the flow rate increased. This demonstrates that the efficiency of the aeration decreases as we increase the flowrate.


!(https://github.com/hispanicberniebro/CEE4530/blob/master/Gas%20Transfer%20Lab/images/OTE.jpg)
Figure 4: OTE as a function of airflow rate

The oxygen transfer efficiency (OTE) displays an exponentially declining slope and approaches zero. This supports the idea that gas transfer efficiency decreases as airflow increases (Figure 4).

We modified the oxygen transfer efficiency  equation 1.9:
OTE = $\frac{\hat{k}_{v,l}(C^{*}-C)VRT}{Q_{air}P_{air}f_{O_{2}}MW_{O_{2}}}$

to consider when the molar airflow rate is controlled, resulting in us using  equation 1.10:

$OTE = \frac{\dot{n}_{aq O_{2}}}{f_{O_{2}}\dot{n}_{air}} = \frac{V\hat{k}_{v,l}(C^{*}-C)}{f_{O_{2}}\dot{n}_{air}MW_{O_{2}}}$

!(https://github.com/hispanicberniebro/CEE4530/blob/master/Gas%20Transfer%20Lab/images/O2transfer.jpg)
Figure 5: Molar rate of oxygen dissolution as a function of time

The molar rate of oxygen dissolution has a shape similar to Figure 3. There is a sharp peak near the beginning of the trial and another smaller peak around the 300 flow rate, after which the plots slopes down linearly (Figure 5). The plot of the molar rate of oxygen dissolution decreased with increasing airflow rate.  

The results from our efficiency plot versus flowrate demonstrate that the efficiency experiences exponential decay with increasing flow rate.
Furthermore, in order to determine the molar rate of oxygen dissolution into the aqueous phase as a function of airflow rate, we used equation 1.7:  
  $\dot{n} = \frac{V}{MW_{O_{2}}}\frac{dC}{dt}$

The air flow rate we obtained was 193 μM/s which is within 20% of the target flow of 200 μM/s (Table 1).

Table 1: Experiment Parameters
| Parameters             | Values    |
| :-------------         | :-------- |
| Test flow rate (μM/s)  | 193       |
| C*                     | 8.923 mg/L          |
| $\hat{k}_{0}$         |  2.495 E -4         |
| $\hat{k}_{30}$         |    2.7 E - 2       |
| $\hat{k}_{70}$         |    2.95 E -2        |
| $\hat{k}_{160}$        |     2.35 E -2      |
| $\hat{k}_{370}$        |      2.39  E -2  |
|$\hat{k}_{850}$         |   2.31 E -2        |
|$\hat{k}_{2000}$         |   2.08 E -2        |


## Conclusions

The results of our experiment gave us a C* value of 8.925 mg/L. The results provide evidence that match the expectations of the theory; as the oxygen deficit decreases the rate of gas transfer decreases and becomes less efficient. The OTE exponentially decreased with airflow rate which implies one should try to keep a low airflow rate and sufficiently large oxygen deficit to maintain a high functioning waste water treatment plant. With this experimental understanding of gas transfer rates in relation to efficiency we will be able to apply this knowledge to future wastewater treatment plant optimization projects.  

## Suggestions
While conducting our experiment some large air bubbles stuck to our DO probe and caused some faulty reading and we had to redo one trial. Additionally the probe seemed to lose it’s calibration over the course of the week as everyone had to recalibrate the apparatus  at the beginning of lab. The experimental apparatus was difficult to initially set up; however, it was useful that we were given a whole lab period to do so. Over the course of the week, some water from our apparatus likely evaporated, but we were not aware of this at the time so we could not correct it. A way in which we could modify the experimental apparatus would be to get more reliable dissolved oxygen probes that did not require as much recalibration throughout the experiment. Due to so much of our analysis relying on our dissolved oxygen data, it is important that we have accurate data, which is a direct result of the dissolved oxygen probe quality. The needle valve was also bit difficult to adjust when we were trying to get to 50% pressure in the accumulator. Adding some type of activated sludge would be an interesting way to further explore the idea of oxygen transfer.

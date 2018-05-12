### Import statements to run code for Nutrient Removal Project.

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


## Objectives


## Data Analysis

```Python
cd C:\Users\Anthony\github\CEE4530_axa2\Final Project Files

COD_Feed = 325/100 * u.mg/u.L

DO_250Day1 = Column_of_data('250flowDay1.xls', 0, -1, 5, 'mg/L')
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


DO_plot_1 = plt.plot(Time_Day1.to(u.hr), DO_250Day1.to(u.mg/u.L), 'ro')
plt.xlabel(r'$time (hr)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.title('250 umol Airflow')
plt.savefig(r'images\250Trail1.jpg')
plt.show()


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
for i in range(1, len(k_values)):
  OTE[i] = ((Volume * (k_values[i]) * oxygen_deficit) / (f_O2 * (flow_values[i]) * MW_O2))
  i = i + 1

plt.plot(flow_values, OTE)
plt.xlabel('Air Flow (umol / s)')
plt.ylabel('OTE')
plt.title('Oxygen transfer efficiency as a function of Air Flow Rate')
plt.savefig(r'images\OTE.jpg')
plt.show()

```

## Laboratory 2
## Anthony Arce
2. Create a graph of absorbance at 660 nm vs. concentration of methylene blue in Atom using the exported data file. Does absorbance at 660 nm increase linearly with concentration of methylene blue?



```python
from aide_design.play import *

help(np.array)
concentration = np.array([0,1,2,3,4,5]*u.mg/u.L)
absorbance660 = np.array([.0005,0.2216,0.3686,0.5812,0.7560,0.9999])

plt.plot(concentration,absorbance660)
plt.ylabel('Absorbance')
plt.xlabel('Concentration Methylene Blue (mg/L)')
plt.show()
```
3. 3)	Plot epsilon as a function of wavelength for each of the standards on a single graph. Note that the path length is 1 cm. Make sure you include units and axis labels on your graph. If Beer’s law is obeyed what should the graph look like?

```python

File = r'C:\Users\Anthony\github\CEE4530_axa2\CEE4530_Lab1data_ohe3.csv'
data = pd.read_csv(File)
info = np.array(data)

plt.plot(info[:,0], info[:,1:7])
plt.xlabel('Wavelengh (nm)')
plt.ylabel('Absorbance')
plt.show()

```

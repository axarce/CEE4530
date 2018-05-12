
###### Anthony Arce

###### Elle Blake

###### Ajara Cobourne

<span style="color:blue">
I interpret your proposal to be measuring the rate of BOD degradation as a function of temperature. This is a reasonable relationship to explore. In order to do this, you will need a method to control the temperature of the sequencing batch reactor. ProCoDA can use PID (look it up), a temperature probe, and a heating element to control temperature. You could conceivable use a glass beaker on a hot plant or we could get a small immersion heater.
There are some complexities here because measuring BOD is such a pain. Temperature also influences the equilibrium oxygen concentration, the first order reaeration rate constant, and the diffusion of oxygen into the oxygen probe. This will make it challenging to separate these 3 influences as you are interpreting the results.
Another critique is that changing the water temperature of wastewater treatment isn't economically viable.
I suggest you continue to consider hypotheses to test.

</span>

<span style="color:blue">alternative hypotheses: 1) The volume of air required to achieve a target BOD removal efficiency can be reduced by operating at a lower DO concentration in the SBR. 2) The solids concentration (and hence BOD removal efficiency) can be increased by using a tube settler to withdraw settled waste during the time normally used for settling. (It should be possible to continue aerating while draining through the tube settler if the inlet for the tube settler is far from the aeration source or if a bubble removal system is added to the tube settler.) 3) The final BOD can be estimated from the initial BOD and continuously measured oxygen dissolution rate. 4) more ideas here  </span>

### Final Project Nutrient Removal

### Introduction

In the anthropocene humans have been, and continue to, interact with the environment in harmful ways. One example is when effluent discharge with particulate solids and high levels of nutrients mixed into bodies of water. Water that has too much nutrients can lead to several harmful outcomes such as low dissolved oxygen, dead fish, murky water, increased turbidity, algal blooms, and diminished populations of plants and animals that are essential to a healthy ecosystem, and the effects of BOD are affected by seasonal conditions (EPA 2001). A good way to examine potential solutions to this problem is through a model wastewater treatment plant (WWTP). The goal of a WWTP is to reduce BOD (Biological Oxygen Demand) and provide harmless effluent water to avoid harming the environment. The EPA standard is that the 30 day average 5 day BOD (BOD5) of the effluent be less than 30 mg/L. The typical BOD test takes 5 days, but measuring the oxygen uptake rate is sufficient for us to infer the BOD of the effluent in this experiment. In this experiment we will look at how to reduce BOD in effluent water.


### Objectives

As observed in the Gas Transfer lab, as water becomes more saturated with ~~water~~, the slower gas transfer occurs. An ideal reactor would maintain a very low level on oxygen in the tank, just enough to allow the microbe community to thrive, but not enough to decrease the reactor efficiency. We will try to run our reactor at constant high efficiency and will vary temperature and observe how that affects the BOD effluent. We hypothesize that BOD and temperature will vary directly. <span style="color:red">(State this more precisely. What is your hypothesis?) </span>We will estimate BOD removal by measuring the oxygen uptake rate. If we increase the concentration of oxygen in the reactor tank, cut off the airflow and monitor how much the oxygen concentration in our batch reactor changes we can obtain the oxygen uptake rate.


•Key design parameters

-Flow rates -  0, 160, 850 μM/s.

-Volumes - 600ml

-Concentrations - diluted stock concentration of synthetic wastewater

-Range of parameters that you are varying - Temperature: 275, 295, 315 K

### Timeline

4/11/18: Goals - Set up experimental apparatus, including ~~ProCoDa~~ ProCoDA (Process Control and Data Acquisition) and a test run with tap water. If time allows run 1-2 trials with wastewater

4/18/18: Goals - Run 6 trials with synthetic wastewater

4/25/18: Goals - Gather additional data as needed and begin analysis

5/02/18: Goals - Complete analysis

5/09/18: Goal: Present final results in class

One challenge we might face when doing this experiment might be getting a successful reactor setup. Maintaining a constant temperature might also be difficult. There will not be a TA bench model and we will have to do more trial and error than we are used to. Monroe said that this was a notoriously difficult set up in the past, so we may not be able to collect data quickly. Additionally we will have to modify some of the equipment in the lab. We do not have any small beakers that have enough spouts to correctly model a batch reactor so modifying what we do have may prove to be a challenge.

### Expectations
A 1999 study by Griffin, Bhattarai and Xiang found that there was a significant difference between BOD levels in effluent water at different temperatures. They observed a seasonal variation where at temperatures less than 293 K, BOD removal had more variation and was less effective. However at temperatures above 293 K more BOD removal occurred. In February 2018 Waki, Yasuda, Fukumoto, Béline, and Magrí published the results of a study where wastewater was treated at either a continuously low DO, or a continuously high DO level at temperatures ranging from 283 to 303 K. They found that nutrient and BOD removal were highest at low DO levels and high temperatures. We expect to get similar results to these studies.


•Resources needed to conduct experiments – What tools will you use?
-1 L modified reactor

-ProCoDa II

-Sensors - (pressure, temperature, dissolved oxygen, flow rate)

-Synthetic Wastewater

-Activated Sludge

-A filter

-A scale

-A warmer

-Ice

-Glass beaker

-Magnetic Stirrer

-Accumulator

-Needle Valve

-Solenoid Valve

-Air Supply



### References/Bibliography


1.https://www.epa.gov/sites/production/files/2014-08/documents/nutrient-memo-nov142001.pdf


2.http://www.ingentaconnect.com/contentone/wef/wer/1999/00000071/00000004/art00012


3.https://www.sciencedirect.com/science/article/pii/S0960852417320825




#### Some Relevant Equations:
<span style="color:red">(I fixed the equations below so that they display correctly. You sometimes forgot the $ and the equations had extra "\[". ) </span>

$\frac{dL}{dt} = \frac{-kLX}{K_s+L}$  the Monod relationship, substrate use by bacteria

Assuming Ks, the half velocity constant is relatively large we can simplify the equation and integrate to obtain the following:

$L = L_o e^{-k_{ox}t}$

We can say that the oxygen use rate is the same as the substrate use rate so

$\frac{dL}{dt} = -k_{ox}L = -k_{ox}L_oe^{-k_{ox}t}$

$ln \frac{C^* - C}{C^* - C_o}$ when aeration is the only thing affecting [O]

The change in the oxygen deficit: $\frac{\partial D}{\partial t} = k_e + k_{ox}L_oe^{-k_{ox}t}-\hat{k}_{v,l}D$

This equation holds under the assumption that we have enough oxygen to maintain microbial consumption.

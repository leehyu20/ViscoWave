# VISCOWAVE VERSION 1.0.0

## INTRODUCTION

**ViscoWave** is a dynamic Finite Layer Method (FLM) which is capable of modeling the dynamics (i.e., transient stress wave propagation) of a variety of flexible pavement structures with elastic/viscoelastic layers, with/without bedrock (Lee, 2013 & 2014).  As shown in the figure below, the FLM is similar to the Finite Element Method (FEM) in the sense that the entire pavement structure is broken down into smaller elements for the analysis.  
<img src="Release Notes/Figures/FEM_vs_FLM.jpg" ></img>

The FEM requires each layer of the pavement be divided into a large number of smaller elements to ensure accurate and reliable results.  Due to the large number of elements, FEM usually is computationally expensive, which may be undesirable for backcalculation.  

On the other hand, the fundamental concept behind the FLM is to define an entire layer (or sublayer) as a single element; thus for a 3-layer pavement system, FLM only needs 3 elements.  This is made possible by taking advantage of the simplified geometry of each pavement layer (see above figure on the right) and solving the wave equations analytically (or semi-analytically). In other words, the phenomenon of stress wave propagation within a given layer is solved in closed-form or semi closed-form.  The interaction between multiple layers (or elements in this case) and the boundary conditions are then accounted for in the same manner as the FEM.

A unique characteristic of ViscoWave which separates it from other available FLM solutions (e.g., SAPSI-M by Chatti and Yun, 1996; 3D-Move by Siddharthan et. al., 2000; LAMDA by Al-Khoury et. al., 2001) is that ViscoWave is based on impulse loading and responses. The following figure shows the conceptual schematics on how the impulse loading and responses are used in ViscoWave for obtaining the overall pavement response under a time-dependent, Falling Weight Deflectometer (FWD) type of loading. The details regarding the formulation of ViscoWave is documented in Lee (2013, 2014).

<img src="Release Notes/Figures/Impulse_Load_Response.jpg " ></img>


## INSTALLING VISCOWAVE

Installation of ViscoWave is simple. Download the "ViscoWave_Installer.zip" and extract the files to a local directory. Then double click the Setup.exe file. Follow the on-screen instructions to complete installation.  

**Notes on ViscoWave installation:**

**1.	Installation of ViscoWave requires that the user has Administrator Access to the computer.**

**2.	ViscoWave (and its installer) has only been tested on 64-bit Windows 10 Environment, with 32-bit or 64-bit Microsot Excel.**

## GETTING STARTED

Once ViscoWave is installed, you should see a ViscoWave shortcut created on your Windows Desktop. Double-click on the shortcut to open the ViscoWave interface built in the Microsoft Excel environment. 

**Notes on ViscoWave interface:**

**1.	This is an interim interface that may be improved/changed with future release of ViscoWave.**

**2.	This is a “Read-Only” file. Please save the template with a different file name if you do not want to lose your changes.**

**3.	You need to have some basic knowledge of Microsoft Excel.**

**ViscoWave Excel interface** has the following 4 tabs:

1.	**VW** – Tab for ViscoWave forward simulation and backcalculation. 
2.	**Default Inputs** – Tab for storing the default (or seed) pavement parameters (Modulus, Poisson’s Ratio, Density, Thickness, and Damping).
3.	**Dynamic_Modulus_Calc** – Tab for converting the Relaxation Modulus (in time-domain) to Dynamic Modulus (in frequency-domain), or vice versa. 
4.	**AGPL License V3** – Tab showing the license agreement (AGPL Version 3.0) for ViscoWave. 

A little more description on the VW and Dynamic_Modulus_Calc tabs are provided in the following. 

### VW Tab for ViscoWave Simulation

Running ViscoWave is easy!  Simply fill in the pavement structure information as well as sensor offsets, FWD plate radius, and the FWD load time history as shown in the figures below. The simulation results will automatically be updated. 

**Notes on ViscoWave inputs:**

**1.	To simulate a viscoelastic layer, the modulus of that particular layer must be assigned a value of ZERO in the General Pavement Structure table, and the sigmoidal function coefficients should be provided in the table below it.**

**2.	To simulate a half-space (i.e., a semi-infinite subgrade), the thickness of the last layer must be specified to ZERO. However, It is generally recommended NOT to use a half-space for backcalculation purposes – use a very thick subgrade instead (e.g., 300 in. to 500 in. subgrade thickness).**

**3.	Current ViscoWave interface allows up to 6 pavement layers (including 3 viscoelastic layers) and 9 deflection sensors. This is simply due to the interface, not the limitation of ViscoWave code.**

<img src="Release Notes/Figures/ViscoWave_Input_Screen.jpg " ></img>

The VW tab also allows you to generate a few synthetic load time histories: Half-Sine, Haver-Sine, and Normal (Gaussian) shaped loads. Simply fill in the necessary inputs for the load you want then click on the corresponding button to generate the load. 

**Notes on ViscoWave time histories:**

**1.	Time Increment (dt) must be 0.2 milliseconds (Hard Coded in C++ Code).**

**2.	Maximum time should NOT exceed 0.07 seconds. It is generally recommended not to modify the time column. In other words, use the time from 0 to 0.0598 sec. at an interval of 0.002 sec.**

<img src="Release Notes/Figures/ViscoWave_TimeHistories.jpg " ></img>


To run backcalculation, first import a CSV file including the FWD load and deflection time histories. Fill in the pavement inputs. Then, try changing some of the pavement parameters (i.e., modulus) until you see a reasonable seed value to start with. Then “Run Backcalculation” button to let the Excel Solver find you the optimum values. 

**Notes on CSV file and Backcalculation:**

**1.	A couple of sample CSV files are provided in the folder named "Sample FWD Time Histories".**

**2.	Future release of ViscoWave will include algorithms for creating CSV files from FWD output files (for Dynatest, JILS, and Kuab FWDs).**

**3.	Current version of ViscoWave only allows for unconstrained optimization (i.e., backcalculation without any upper or lower bounds specified to the variables). As such, it is strongly recommended that you start with a good Seed value for backcalculation.  A future release of ViscoWave may allow for constrained optimization.**

<img src="Release Notes/Figures/ViscoWave_Defl_Plots.jpg " ></img>

### Dynamic_Modulus_Calc Tab for Viscoelastic Modulus

As discussed above, the Sigmoidal Function for ViscoWave (in the VW tab) corresponds to the Relaxation Modulus in time-domain.  The Dynamic_Modulus_Calc tab allows for converting the Relaxation Modulus (in time-domain) to Dynamic modulus (in frequency-domain). If you have Dynamic Modulus Sigmoidal Coefficients (say from laboratory testing of asphalt mixtures), and want to use them in ViscoWave, this tab also allows you to convert the Dynamic Modulus to Relaxation Modulus. 


1.	To convert from Relaxation Modulus (E(t)) to Dynamic Modulus (|E*|), simply fill in the Relaxation coefficients (Top Left table) and click on the button “Convert E(t) to |E*|”.
2.	To convert from Dynamic Modulus (|E*|) to Relaxation Modulus (E(t)), simply fill in the Relaxation coefficients (Top Right table) and click on the button “Convert |E*| to E(t)”.

The tab also provides the Prony series coefficients that can be used to calculate both the Relaxation and Dynamic moduli. 

<img src="Release Notes/Figures/Sigmoidal_and_Prony_Coefficients.jpg" ></img>

The tab also shows the plots for the Relaxation and Dynamic moduli from both the Sigmoidal functions and Prony series. 

<img src="Release Notes/Figures/Modulus_Plots.jpg " ></img>


## COMING SOON…

ViscoWave will continue to be updated for the next few months.  Stay tuned for the following updates to be released!

1.	To include algorithms for creating CSV files from FWD output files (for Dynatest, JILS, and Kuab FWDs). 
2.	To implement constrained optimization for backcalculation with lower and upper bounds for layer moduli. 
3.	A detailed user manual. 

## REFERENCES

Al-Khoury, R., Skarpas, A., Kasbergen, C., and Blaauwendraad, J., Spectral element technique for efficient parameter identification of layered media. I. Forward calculation. International Journal of Solids and Structures, Vol. 38, No. 9, 2001, pp. 1605-1623.

Chatti, K., and Yun, K.K., SAPSI-M: Computer Program for Analyzing Asphalt Concrete Pavements under Moving Arbitrary Loads. In Transportation Research Record 1539, TRB, National Research Council, Washington, D.C., 1996, pp. 88–95. 

Lee, H.S. 2013. Development of a New Solution for Viscoelastic Wave Propagation of Pavement Structures and Its Use in Dynamic Backcalculation. Ph. D. Dissertation. Michigan State University. 

Lee, H.S. 2014. ViscoWave – A New Solution for Viscoelastic Wave Propagation of Layered Structures Subjected to an Impact Load. International Journal of Pavement Engineering, 15(6): 542-557.

Siddharthan, R.V., Krishnamenon, N., and Sebaaly, P.E., Finite-Layer Approach to Pavement Response Evaluation. In Transportation Research Record 1706, TRB, National Research Council, Washington, D.C., 2000, pp. 43–49. 

# MI-CL-Final-Project
Code to implement my final project for MGTECON634 Spring '21.

The three main sections of my short paper each have their own accompanying script. 

"Propensity Score Prediction" corresponds to "pscores.R", which uses several methods to predict propensity scores on my data and estimates the DGT draft specification with these alternative propensity scores. 

"Event Study Design" corresponds to "event_studies.R", which estimates several pooled event studies with the simple estimators, constructed from out-of-the-box ML methods, described in the paper. 

"Heterogeneous Treatment Effects" corresponds to "causal_forests.R", which runs causal forests on my key outcome variables and implements some calibration tests.  

Some of the estimation steps require more memory than I have available on a local machine, and run on the Sherlock HPC cluster; as such, some estimation code is commented out and followed by reading in the results from the same code run on Sherlock.  

The code runs on pre-publication data from:

*Diamond, R., A. Guren, and R. Tan (2020, June).  The Effect of Foreclosures on Homeowners,Tenants, and Landlords.  Technical Report w27358, National Bureau of Economic Research,Cambridge, MA*

I do not have permission from the authors to share this data. 

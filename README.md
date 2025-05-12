# SCUBIDO
The **SCUBIDO** (Simulating Climate Using Bayesian Inference from proxy Data Observations) R package has been created to produce Bayesian reconstructions of climate variability based on XRF core scanning data. 
The statistics behind this approach can be found in Boyall et al (in prep) and should be cited alongside this package if used for future work. 
The R package can be downloaded using the following:

```r
devtools::install_github("LauraBoyall/SCUBIDO")
library(SCUBIDO)
```
You can find detailed instructions on how to run this package here (https://lauraboyall.github.io/SCUBIDO/) and on this Youtube tutorial here: https://www.youtube.com/watch?v=3xsaj7CTqPk. 

**NOTE** that this package uses JAGS (Just Another Gibbs Sampler) and therefore requires the software to be downloaded to your device before the code can run. You can access the latest version for download here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/ 

**Stuff to remember**
The model files are saved to your working directory with the system data in the title. This allows the model to be re-uploaded rather than having to be re-run, however we would recomend that 
if you run the model again on another day you go to the original saved file and change the date to the date you upload, otherwise the model will re-run and not upload. 
For example if your file is saved as *cal_model_100924.RDs* for the 10th of September 2024 but you are now wanting to explore the model again in October just manually change to a different 
date of that day e.g. *cal_model_031024.RDs*

Any additional help needed please contact Laura.Boyall.2016@live.rhul.ac.uk or SCUBIDO.info@gmail.com

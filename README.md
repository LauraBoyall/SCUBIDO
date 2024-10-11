# SCUBIDO
<img src="man/SCUBIDO logo .png" align="right" width="120" />
The **SCUBIDO** (Simulating Climate Using Bayesian Inference from proxy Data Observations) R package has been created to produce Bayesian reconstructions of climate variability based on XRF core scanning data. 
The statistics behind this approach can be found in Boyall et al (in prep) and should be cited alongside this package if used for future work. 
The R package can be downloaded using the following:
```r
devtools::install_github("LauraBoyall/SCUBIDO")
library(SCUBIDO)
```
## Instructions
This package has three functions: `SCUBIDO_input` which is used to sort the modern calibration data (includes Age BP, an instrumental 
climate record and then overlapping XRF data) and the fossil/proxy data (containing Age BP and the rest of the XRF data which does not overlap with the climate data. 
This function sorts the data and should be saved as an object to be used in the next function. 
The next function is `SCUBIDO_cal` which runs the calibration model which looks at the direct and covariant relationships between the XRF data and climate. This should also be saved into an 
object to be used in the following function. 
The final function is `SCUBIDO_reconstruct` which applies a what was learnt in the calibration period to the rest of the proxy data. **Note that this function can take a very long time to run**. 

## Stuff to note
The `modern_data` data frame should be set out as age BP in the first column, the climate timeseries in the second and then the rest of the XRF data in the remaining columns.

The `modern_data` data frame should be set out as age BP in the first column, and then then the rest of the XRF data in the remaining columns, thus one less column than `modern_data`

The model files are saved to your working directory with the system data in the title. This allows the model to be re-uploaded rather than having to be re-run, however we would recomend that 
if you run the model again on another day you go to the original saved file and change the date to the date you upload, otherwise the model will re-run and not upload. 
For example if your file is saved as *cal_model_100924.RDs* for the 10th of September 2024 but you are now wanting to explore the model again in October just manually change to a different 
date of that day e.g. *cal_model_031024.RDs*

Any additional help needed please contact Laura.Boyall.2016@live.rhul.ac.uk

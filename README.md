# SCUBIDO - Simulating Climate Using Bayesian Inference with proxy Data Observations 

`SCUBIDO` is an R package which has been designed to transform multivariate XRF (X-Ray fluorescence) core scanning data into a quantitative climate timeseries. 
This reconstruction uses Bayesian Inference alongside additional statistical approaches to reconstruct the climate providing a posterior distribution of climate
with quantified uncertainties. 

This package relies on JAGS (Just Another Gibbs Sampler) for the Markov Chain Monte Carlo algorithms. This software should be downloaded prior to running `SCUBIDO` by clicking 
[here](https://sourceforge.net/projects/mcmc-jags/). 

A detailed summary of the modelling approach can be found in Parnell et al. (2015, Appl. Statist.), and a more specific overview in Boyall et al. (In prep). 

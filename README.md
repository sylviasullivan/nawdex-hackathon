# Scripts for NAWDEX hacakthon 

A collection of analysis scripts for ICON simulations for the NAWDEX field campaign. The Voigt group performed analysis at KIT in 2020 and 2021. Write-up done by S Sullivan in 2022.

The directory structure is as follows:

* bashscripts: scripts to pull data from the DKRZ archive
* ccfsampling: calculation of cloud controlling factor sampling
* classsampling: calculation of cloud cover, cloud fraction, and cloud condensate profiles filtered by cloud class
* domain-mean: calculation of domain-mean time-mean cloud radiative heating rates and other temperature tendencies
* globalmodels: extraction of cloud radiative heating rates from coarse resolution AMIP-like models
* means: calculation of mean quantities, where mean can be over time and/or space
* postprocessing: scripts and notebooks to calculate postprocessed data, such as ocean masks, re-diagnosed radiaive heating rates, CRE from CERES, ...
* sandbox: a playground of scripts that were written but not polished, yet still might be useful to keep
* shared: a collection of utility dictionaries and functions, e.g., the simulation dictionary and a function to drop the first day of a dataset 
* utilities: utility scripts to, for example, select analysis days, plot domain-mean values, or clean up axes

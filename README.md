# MSS_NPS

To get started, download and unzip US-level ACS data for people and households from

https://www2.census.gov/programs-surveys/acs/data/pums/2023/5-Year/csv_pca.zip

and 

https://www2.census.gov/programs-surveys/acs/data/pums/2023/5-Year/csv_hca.zip

Place these in the `data` directory, then run `process_population.R` to generate the population file. The simulation logic is contained in `main.R`.

The code for each model resides in the `model/` directory and helper functions are in `utils.R`.

*Any dependencies should be loaded only in `main.R`*

# ACS Data

For an overview of the different variables in the ACS, consult the data dictionary at 

https://www2.census.gov/programs-surveys/acs/tech_docs/pums/data_dict/PUMS_Data_Dictionary_2023.pdf

# Simulation study structure

We conduct an empirical simulation study treating the ACS data from California as our population. Results are reported at the PUMA level.

We take as a response the binary variable HICOV (health insurance coverage) with the following covariates
* age (categorized as 1. under 18, 2. 18-33, 3. 34-64, and 4. 65+)
* race (Hispanic, white, Asian, Black, other), and
* sex (male, female)
* we may also want to explore using "ACCESSINET" (access to internet) as a covariate that may be predictive of one's inclusion into the non-probability sample, but do not do this in the main code yet. 

The output for each model should contains the following fields
* `model`: the name of the method
* `PUMA`: PUMA (area) ID
* `point_est`: a point estimate of the proportion of HICOV by PUMA
* `lower_CI`: the lower CI bound at confidence level `alpha` (this parameter is set in main.R)
* `upper_CI`: the upper CI bound at confidence level `alpha`

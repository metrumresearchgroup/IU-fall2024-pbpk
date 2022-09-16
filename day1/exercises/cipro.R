
# Cipro PK: Case-Study -------------------------------------------

#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3555057/pdf/bcp0075-0180.pdf

library(mrgsolve)
library(dplyr)
library(PKPDmisc)
library(purrr)

# Cipro PK: warm-up ------------------------------------



# Cipro PK: infusion doses ------------------------------------

#' - Simulate 400 mg as IV infusion over 1 hour
#' - Plot ciprofloxacin concentration over 12 hours

#' 400 mg IV over 1 hour



#' - Simulate 400 mg q12h infusion doses for one day
#' - Plot over 48 hours of dosing

# Cipro PK: covariates ----------------------------------

#' - Weight: 40 to 140 kg
#' - Creatinine clearance: 90, 60,45, 30, 15 mL/min
#' - Sex: male / female
#' - Inputs: idata_set, event object
#' - Outputs: ciprofloxacin concentrations & covariate vallues



# Cipro PK: AUC ------------------------------

# Cipro PK: organ-specific exposure -------------------------------



# Cipro PK: simulation with uncertainty ---------------------------

#' - ciprofloxacin 400 mg IV over 2 hours
#' - What is the distribution of AUC after the last dose on day 3 of 
#'   q12h dosing?




# Cipro PK: AUC:MIC -----------------------------------------------

#' First, make all of the dosing regimens

#' Test it out

#' Pull them together into a list

#' The posterior

#' For regimen (x) simulate from 24 to 48 hours

#' Add a factor

#' mic 0.25

#' mic 1


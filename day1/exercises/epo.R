
# EPO route of administration --------------------------

source("src/global.R")
library(dplyr)
library(mrgsolve)
rm(list = ls())

#' Load the epo model out of the model directory
#' 
#' In this model
#' - The SC dosing compartment is 1
#' - The IV dosing compartment is 2
#' - The hemoglobin output variable is HGBi
#' 
#' For a flat 40,000 IU/week dose, simulate hemoglobin changes
#' over 1 month for IV and SC administration
#' 
#' Use the ev_days() function to help construct a simulation 
#' comparing 40,000 IU/week flat dose and 100 IU/kg TIW (M,W,F)
#' over 1 month
#' 




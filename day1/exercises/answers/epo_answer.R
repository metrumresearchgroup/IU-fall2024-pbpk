
# EPO route of administration answer ------------

source(here("day1/src/global.R"))
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

mod <- mread("epo", here("day1/model")) %>% zero_re()

# Compare IV and SC ====================

data <- expand.ev(amt = 40000, cmt = c(1,2), ii = 168, addl = 4)

mod %>% 
  data_set(data) %>%
  mrgsim(end = 700, delta = 0.1) %>% 
  plot(HGBi~time)


# Compare QW and TIW ====================

flat <- ev_days(
  ev(amt = 40000, cmt = 1), 
  days="m", 
  addl = 3
)

wt <- ev_days(
  ev(amt = 70*100, cmt = 1), 
  days="m,w,f", 
  addl = 3
)
 
data <- as_data_set(flat,wt)

mod %>% 
  data_set(data) %>%
  mrgsim(end = 700, delta = 0.1) %>% 
  plot(HGBi~time)



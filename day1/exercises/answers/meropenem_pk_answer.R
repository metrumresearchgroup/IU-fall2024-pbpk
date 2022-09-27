
# Meropenem PK answer ----------------------------------

source(here("day1/src/global.R"))
library(dplyr)
library(mrgsolve)
rm(list = ls())

#' - Load the `meropenem` model from the model directory
#' - Simulate the following scenarios:
#'   - 100 mg IV bolus q8h x 3
#'   - 100 mg IV over 3 hours q8h x3
#'   
#' Look at the `CC` output

mod <- mread("meropenem", here("day1/model"))

mod %>% 
  ev(amt = 100, ii = 8, addl = 2) %>%
  mrgsim(end = 24, delta = 0.1) %>% 
  plot(CC~time)

mod %>% 
  ev(amt = 100, ii = 8, addl = 2, tinf = 3) %>%
  mrgsim(end = 24, delta = 0.1) %>% 
  plot(CC ~ time)






# GCSF dosing answer --------------------------------------

source(here("day1/src/global.R"))
library(dplyr)
library(mrgsolve)
rm(list = ls())


#' - Model: model/gcsf.cpp
#' 
#' - Simulate 2.5, 5, and 10 mcg/kg assuming 50, 70, 90 kg individual
#' 
#' - Do daily SC dosing x 7d

mod <- mread_cache("gcsf", here("day1","model")) %>% zero_re()

data <- 
  expand.ev(
    dose = c(2.5, 5, 10), 
    WT = c(50,70,90), 
    ii = 24, 
    total = 7
  ) %>% mutate(amt = WT * dose)

out <- mrgsim_d(mod, data, carry_out = "dose,WT")

plot(out, RESP~time|factor(dose), groups = WT, scales = "same")

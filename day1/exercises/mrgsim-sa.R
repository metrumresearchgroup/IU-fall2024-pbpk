source(here("day1/src/global.R"))
library(dplyr)
library(mrgsolve)
library(mrgsim.sa)
rm(list = ls())

mod <- mread("yoshikado", here("day1/model")) %>% 
  update(delta = 0.1, end = 12)

#' A single pitavastatin dose
pit <- ev(amt = 30, cmt = 1)

#' A single pitavastatin dose 30 min after CsA
csa <- ev(amt = 2000, cmt = 2)
ddi <- seq(csa, wait = 0.5, pit)


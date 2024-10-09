library(here)
source(here("day1/src/global.R"))
library(data.table)
library(mrgsolve)

#' 
#' Simulate Cyclosporine PK from PBPK model
#' 
#' - load the `csa` model from the `day1/model` directory
#'   - interrogate the model compartments and parameters
#' - prepare a dose of 2000 ug/kg x1 administered to `gut`
#' - plot `CSA` and `CSAliv` outputs versus time
#' 
#' - Next, administer 2000 ug/kg BID for 1 week
#' 

mod <- mread("csa", project = here("day1","model"))

mod
see(mod)
param(mod)
init(mod)


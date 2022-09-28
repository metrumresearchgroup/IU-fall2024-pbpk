library(here)
source(here("day1/src/global.R"))


#' 
#' # Rifampin DDI with midazolam
#' 
#' - model file: `day1/model/rifampicin_midazolam.cpp`
#' - look at the model, parameters and compartments
#' - Dosing regimen 1
#'   - rifampin 600 mg po daily x10 `then`
#'   - midazolam 3 mg x1 12 hours after last rifampin dose
#' - Dosing regimen 2
#'   - midazolam 3 mg x1 administered at the same time as 
#'     we administered in regimen 1
#'     
#' - Rifampin goes in `Xgutlumen`, midazolam goes in `mgutlumen`
#'     
#' - Plot rifampin concentration over time for regimen 1
#' - Plot midazolam concentration over time for regimen 1 and 2
#'
#' - Outputs
#'  - Ccentral - rifampin concentration
#'  - mCcentral - midazolam concentration
#' 

mod <- mread("rifampicin_midazolam", here("day1/model"))






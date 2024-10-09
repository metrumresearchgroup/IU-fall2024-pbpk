library(here)
source(here("day1/src/global.R"))


#' 
#' # Rifampin DDI with midazolam
#' 
#' - model file: `day1/model/rifampicin_midazolam.cpp`
#' - look at the model, parameters and compartments
#' - Dosing regimen
#'   - rifampin 600 mg po daily x10 `then`
#'   - midazolam 3 mg x1 24 hours after last rifampin dose
#'     
#' - Rifampin goes in `Xgutlumen`, midazolam goes in `Mgutlumen`
#'     
#' - Plot rifampin concentration over time
#' - Plot midazolam concentration over time
#'
#' - Outputs
#'  - Ccentral - rifampin concentration
#'  - mCcentral - midazolam concentration
#' 

mod <- mread("rifampicin_midazolam", here("day1/model"))

e1 <- ev(amt = 600, cmt = "Xgutlumen", ii = 24, addl = 9)
e2 <- ev(amt = 3, cmt = "Mgutlumen")
e <- seq(e1, e2)
mod %>% ev(e) %>% mrgsim(end = 260) %>% plot(Ccentral + mCcentral ~ time/24)




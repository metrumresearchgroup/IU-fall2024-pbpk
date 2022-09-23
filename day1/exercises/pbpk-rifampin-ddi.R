library(here)
source(here("day1/src/global.R"))


#' 
#' # Rifampin DDI with midazolam
#' 
#' - model file: `day1/model/rifampicin_midazolam.cpp`
#' - look at the model, parameters and compartments
#' - Dosing regimen 1
#'   - rifampin 600 mg po daily x7 `then`
#'   - midazolam 3 mg x1 12 hours after last rifampin dose
#' - Dosing regimen 2
#'   - midazolam 3 mg x1 administered at the same time as 
#'     we administered in regimen 1
#'     
#' - Rifampin goes in `Xgutlumen`, midazolam goes in `mgutlumen`
#'     
#' - Plot rifampin concentraion over time for regimen 1
#' - Plot midazolam concentration over time for regimen 1 and 2
#' 

mod <- mread(here("day1/model/rifampicin_midazolam.cpp"))


rif <- ev(amt = 600, ii = 24, addl = 6, cmt = 1)
mid <- ev(amt = 3, cmt = 2, time = 12)

dose1 <- seq(rif, mid) %>% mutate(reg = "ddi")
dose2 <- filter(dose1, amt==3) %>% mutate(reg = "mono")


dose <- as_data_set(dose2, dose1)


out <- mrgsim(mod, events = dose, end = 240, delta = 0.1, 
              recover = "reg")

plot(out, "Ccentral")

out2 <- filter_sims(out, time >= 168)

plot(out2, "mCcentral")

mod %>% ev(dose) %>% 
  mrgsim(end =240, delta = 0.1) %>% 
  plot("Ccentral,mCcentral")


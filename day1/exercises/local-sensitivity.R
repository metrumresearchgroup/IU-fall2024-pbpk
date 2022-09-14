
# Local Sensitivity Analysis -----------

#library(FME)
requireNamespace("FME")

library(mrgsolve)
library(tidyverse)

mod <- mread("cipro", "model") 

mrgsim(mod)

e <- ev(amt = 400)

mod <- update(mod, end = 12)


func <- function(p) {
  mod %>% 
    param(p) %>% 
    ev(amt = 400) %>%
    obsonly() %>%
    mrgsim(end = 12, delta = 1, rtol = 1E-14, atol = 1E-14)  %>%
    as_tibble()
}

set1 <- c("CRCL", "WT", "klun", "FUP")

set2 <- c("klun",  "kmus", "kadi", "CRCL")

this_set <- set2

parms <- as.list(param(mod))[this_set]

x <- func(parms)

x

out <- FME::sensFun(func = func, parms = parms, sensvar = "CP", map=2, tiny = 1E-6)

plot(out)




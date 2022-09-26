
# Cipro PK: Case-Study -------------------------------------------

#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3555057/pdf/bcp0075-0180.pdf

library(mrgsolve)
library(dplyr)
library(PKPDmisc)
library(purrr)

# Cipro PK: warm-up ------------------------------------

mod <- mread("cipro", here("day1/model"))

mod

param(mod)

init(mod)


# Cipro PK: infusion doses ------------------------------------

#' - Simulate 400 mg as IV infusion over 1 hour
#' - Plot ciprofloxacin concentration over 12 hours

#' 400 mg IV over 1 hour
single <- ev(amt = 400, tinf = 1)

single

out <- mrgsim(mod, events = single, end = 12, delta = 0.1)

plot(out, CP~time)

#' - Simulate 400 mg q12h infusion doses for one day
#' - Plot over 48 hours of dosing

daily <- ev(amt = 400, ii = 12, tinf = 5, addl = 1)

out <- mrgsim(mod, events = daily, end = 48, delta = 0.1, Req = "CP")

plot(out)


# Cipro PK: covariates ----------------------------------

#' - Weight: 40 to 140 kg
#' - Creatinine clearance: 90, 60,45, 30, 15 mL/min
#' - Sex: male / female
#' - Inputs: idata_set, event object
#' - Outputs: ciprofloxacin concentrations & covariate vallues

idata <- expand.idata(SEX = c(0,1), WT = c(40,100), CRCL = c(90,30))

out <- 
  mod %>%
  ev(daily) %>% 
  idata_set(idata) %>% 
  carry_out(SEX,WT,CRCL) %>%
  mrgsim(Req = "CP", delta = 0.1, end = 36)

plot(out, CP ~time | factor(SEX)*factor(WT), scales="same")  

count(as.data.frame(out),ID)

# Cipro PK: AUC ------------------------------

idata <- expand.idata(WT = c(70), CRCL = c(81))

out <- 
  mod %>%
  ev(daily) %>% 
  idata_set(idata) %>% 
  carry_out(SEX,WT,CRCL) %>%
  mrgsim(Req = "CP", delta = 0.1, end = 24)


out %>% 
  as_tibble() %>%
  group_by(WT,CRCL) %>% 
  summarise(AUC = auc_partial(time, CP))
  
out <- 
  mod %>%
  ev(amt = 400, ii = 12, addl = 1, ss = 1) %>% 
  idata_set(idata) %>% 
  carry_out(SEX,WT,CRCL) %>%
  mrgsim(Req = "CP", delta = 0.1, end = 24)


out %>% 
  as_tibble() %>%
  group_by(WT,CRCL) %>% 
  summarise(AUC = auc_partial(time, CP))


# Cipro PK: organ-specific exposure -------------------------------
mod <- mread("cipro_conc", "model", req="")

out <- mrgsim(mod, events = single, end = 120, delta = 0.1, output="df")

long <- tidyr::gather(out, variable, value,CLUN:CKID)

summ <- 
  long %>% 
  group_by(variable) %>% 
  summarise(auc = auc_partial(time,value), Cmax = max(value))

library(ggplot2)
ggplot(summ, aes(x = variable, y = auc)) + geom_col() + theme_bw()


# Cipro PK: simulation with uncertainty ---------------------------

#' - ciprofloxacin 400 mg IV over 2 hours
#' - What is the distribution of AUC after the last dose on day 3 of 
#'   q12h dosing?


mod <- mread("cipro", "model") %>% zero_re() %>% update(end = 72, delta = 0.25)

data <- readRDS("data/cipro_post.RDS") %>% filter(irep <= 100)

e <- ev(amt = 400, tinf = 1.5, ii = 12, until = 72)

out <- mrgsim_ei(mod, e, data, carry_out = "irep", add=0.05)

out %>% 
  filter(time >= 48) %>% 
  group_by(irep) %>% 
  mutate(DV = CP) %>%
  summarise(auc = auc_partial(time,DV)) %>%
  ungroup() %>% 
  summarise(med = median(auc), lo = quantile(auc,0.05), hi = quantile(auc,0.95))

sims <- 
  out %>% 
  filter(time > 0, time <=12) %>%
  mutate(DV = Ckid) %>%
  group_by(time) %>%
  summarise(med = median(DV), lo = quantile(DV,0.05), hi = quantile(DV,0.95))

ggplot(sims, aes(time)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.4, fill = "firebrick") +
  geom_line(aes(y = med), lwd = 1, col = "firebrick") + 
  scale_y_log10() +  theme_bw()

# Cipro PK: AUC:MIC -----------------------------------------------

#' First, make all of the dosing regimens
e1 <- ev(amt = 200, ii = 12, tinf = 1, until = 48, CRCL = 30)
e2 <- ev(amt = 400, ii = 12, tinf = 1, until = 48)
e3 <- ev(amt = 400, ii = 8,  tinf = 1, until = 48)
e4 <- ev(amt = 600, ii = 8,  tinf = 1, until = 48)

lbls <- c("200q12 rf", "400q12", "400q8", "600q8")

#' Test it out
mrgsim_e(mod,e1, delta=0.1, end=48) %>% plot(CP~time)

#' Pull them together into a list
x <- list(e1,e2,e3,e4)

x

#' The posterior
data <- readRDS("data/cipro_post.RDS") %>% filter(irep <= 1000)

#' For regimen (x) simulate from 24 to 48 hours

out <- imap_dfr(x, function(e,i) {
  mod %>% 
    ev(e) %>% 
    idata_set(data) %>%
    update(delta = 0.25, start = 24, end = 48) %>%
    mrgsim(carry.out="irep,ARM",Req="CP") %>% 
    mutate(ID = i*1000 + irep, group=i)
})

#' Add a factor
out <- mutate(out, regimen = factor(group,labels=lbls))


#' mic 0.25
sims <- filter(out, group  %in% c(1,2))
sum25 <- 
  sims %>% 
  group_by(ID,regimen) %>%
  summarise(auc = auc_partial(time,CP)) %>% 
  mutate(aucr = auc/0.5)

sum25


ggplot(sum25, aes(x = aucr)) + 
  geom_histogram(alpha = 0.7,col="grey") + 
  facet_wrap(~regimen) + theme_bw() + 
  geom_vline(xintercept = 125, col = "firebrick", lty =2,lwd=2)

#' mic 1
sims <- filter(out, group %in% c(2,3,4))
sum1 <- 
  sims %>% 
  group_by(ID,regimen) %>%
  summarise(auc = auc_partial(time,CP)) %>% 
  mutate(aucr = auc/1)

sum1

ggplot(sum1, aes(x = aucr)) + 
  geom_histogram(alpha = 0.7,col="grey") + 
  facet_wrap(~regimen) + theme_bw() + 
  geom_vline(xintercept = 125, col = "firebrick", lty =2,lwd=2)


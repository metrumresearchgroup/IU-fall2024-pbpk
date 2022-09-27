# load libraries
library(tidyverse)
library(mrgsolve)

# Chunks 1-2: simplePBPK
# Chunks 3-10: voriPBPK model
# Chunks 11-12: sensitivity analysis
# Chunks 13-14: parameter estimation
# Chunks 15-18: population simulation
# Chunk 19: interactive presentation / shiny
# Chunk 20-21: mAb PBPK


#------------------------------------------------------------------------------#

### Simple PBPK ###

################################################################################
################################## Chunk 1  ####################################
################################################################################

# Populate and compile simplePBPK model
mod <- mread("models/simplePBPK.mod")

################################################################################
################################################################################

################################################################################
################################## Chunk 2 #####################################
################################################################################

## Simulate a 100 mg dose given as an IV bolus
mod %>%
  ev(amt=100, cmt="VEN") %>%
  mrgsim(end=100) %>%
  plot()

################################################################################
################################################################################


#------------------------------------------------------------------------------#

### voriconazole PBPK ###

################################################################################
################################## Chunk 3 #####################################
################################################################################

# Compile voriPBPK model
modA <- mread("models/voriPBPK.mod")

################################################################################
################################################################################

################################################################################
##################################  Chunk 4 ####################################
################################################################################

# Use `calcKp_PT.R` function to calculate voriconzole tissue:plasma partition 
# coefficients according to Poulin and Theil method https://www.ncbi.nlm.nih.gov/pubmed/11782904.

# source partition coefficient calculation script
source("calcKp_PT.R")

#voriconazole physicochemical properties
logP <- 2.56  #lipophilicity
pKa <- 1.76  
fup <- 0.42   #unbound fraction in plasma
type <- 3     #monoprotic base
BP <- 1       #blood:plasma concentration ratio
dat <- read.csv("data/tissue_comp_PT.csv")

#calculate partition coefficients
Kp <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)

#update model parameters partition coefficients
modA <- param(modA, Kp)

################################################################################
################################################################################

################################################################################
##################################  Chunk 5 ####################################
################################################################################

# Simulate the steady state after a 4 mg/kg IV infusion dose of voriconazole given 
# to an adult male with a rate of 4 mg/kg/h twice a day for 7 days. Compare the 
# steady state plasma drug concentration-time profiles from previous simulation to 
# the observed data in `Adult_IV.csv`. (N.B.: observed data were digitized from 
# Zane and Thakker (2014) paper using WebPlotDigitizer https://automeris.io/WebPlotDigitizer/):
  
#load observed IV infusion data
obs <- read.csv("data/Adult_IV.csv")

#set simulation conditions
bw   <- 73
amt  <- 4*bw
rate <- 4*bw
cmt  <- "VEN"
ii   <- 12
addl <- 13
ss   <- 1

#run simulation
sim <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)  

#plot prediction and compare to observed data
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
##################################  Chunk 6 ####################################
################################################################################

# simulate the steady state after a 200 mg PO dose of voriconazole given to an adult male 
# with a rate of twice a day for 7 days. Compare the steady state plasma drug 
# concentration-time profile to the observed data in `Adult_PO.csv`.

obs <- read.csv("data/Adult_PO.csv")

bw   <- 73
amt  <- 200
cmt  <- "GUTLUMEN"
ii   <- 12
addl <- 13
ss   <- 1

sim <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, ss=ss) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)  

gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg PO", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
##################################  Chunk 7 ####################################
################################################################################

# Integrate the model with an additional gut wall enterocyte compartment to account for 
# intestinal clearance, intestinal transit and lumen solubility effects on absorption rate. 
# Note: use about an intestinal clearance that is 30 times lower than hepatic clearance.
# Recompile and re-run the previous step. Any change!!
# https://www.jstage.jst.go.jp/article/yakushi/123/5/123_5_369/_article/-char/en

modA <- mread("models/voriPBPK_ext.mod") %>%
  param(MPPGI = 30.3/30) %>%
  param(Kp)

sim <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, ss=ss) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)  

gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg PO", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################


################################################################################
##################################  Chunk 8 ####################################
################################################################################

# Generate pediatric model

# pediatric (5 yo) male physiology; https://www.ncbi.nlm.nih.gov/pubmed/14506981
pedPhys <- list(WEIGHT = 19,
                Vad = 5.5,
                Vbo = 2.43,
                Vbr = 1.31,
                VguWall = 0.22,
                VguLumen = 0.117,
                Vhe = 0.085,
                Vki = 0.11,
                Vli = 0.467,
                Vlu = 0.125,
                Vmu = 5.6,
                Vsp = 0.05,
                Vbl = 1.5,
                Qad = 0.05*3.4*60,
                Qbo = 0.05*3.4*60,
                Qbr = 0.12*3.4*60,
                Qgu = 0.15*3.4*60, 
                Qhe = 0.04*3.4*60,
                Qki = 0.19*3.4*60,
                Qmu = 0.17*3.4*60,
                Qsp = 0.03*3.4*60,
                Qha = 0.065*3.4*60, 
                Qlu = 3.4*60,
                MPPGL = 26,
                VmaxH = 120.5,
                KmH = 11,
                MPPGI = 0,
                VmaxG = 120.5,
                KmG = 11)

modP <- param(modA, pedPhys)

################################################################################
################################################################################

################################################################################
##################################  Chunk 9 ####################################
################################################################################

# Simulate the steady state after a 4 mg/kg voriconazole IV infusion dosing in a male child 
# subject infused at a rate of 3 mg/kg/h twice a day for seven days. Compare the steady state 
# to the observed data in `Pediatric_IV.csv`.

obs <- read.csv("data/Pediatric_IV.csv")  #load observed data

wt   <- 19  #pediatric body weight
amt  <- 4*wt  
rate <- 3*wt
cmt  <- "VEN"  #intravenous infusion
ii   <- 12
addl <- 13
ss   <- 1

# simulate
sim <- 
  modP %>%
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=1) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Pediatric 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
#################################  Chunk 10 ####################################
################################################################################

# Simulate the steady state after a 4 mg/kg voriconazole PO dosing in a male child subject 
# twice a day for seven days. Compare to obsreved data in `Pediatric_PO.csv`
# Note: Include a similar 30-fold lower intestinal clearance than hepatic clearance.
# https://www.jstage.jst.go.jp/article/yakushi/123/5/123_5_369/_article/-char/en

obs <- read.csv("data/Pediatric_PO.csv")  #load observed data

# adjust intestinal clearance
modP <- modP %>% param(MPPGI = 26 / 30)

# simulation conditions
bw   <- 19
amt  <- 4 * bw
cmt  <- "GUTLUMEN"
ii   <- 12
addl <- 13
ss   <- 1

# simulate
sim <- 
  modP %>%
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, ss=1) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Pediatric 4 mg/kg PO", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Sensitivity Analysis ###

################################################################################
#################################  Chunk 11 ####################################
################################################################################

## graphical ##

# Run graphical sensitivity analysis for the muscle:plasma (`Kpmu`) and 
# lung:plasma (`Kplu`) partition coefficients using adult IV data.

#load observed IV infusion data
obs <- read.csv("data/Adult_IV.csv")

#set simulation conditions
bw   <- 73
amt  <- 4*bw
rate <- 4*bw
cmt  <- "VEN"
ii   <- 12
addl <- 13
ss   <- 1

##' Define an intervention
e <- ev(cmt=cmt, amt=amt, rate=rate, ii= ii, addl=addl, ss=1)

## Sensitivity analysis on Kpmu
idata <- expand.idata(Kpmu = c(3/2, 3, 3*2))

modA %>%
  carry_out(Kpmu) %>% 
  mrgsim_ei(e, idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>% 
  mutate(Kpmu = factor(Kpmu)) %>%
  ggplot(aes(x=time, y=CP, col=Kpmu)) +
  geom_line() +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()

##' Sensitivity analysis on Kplu
idata <- expand.idata(Kplu = c(1/2, 1, 1*2))

modA %>%
  carry_out(Kplu) %>% 
  mrgsim_ei(e, idata, delta = 0.1, recsort=3, obsonly=TRUE, end = 12) %>% 
  mutate(Kplu = factor(Kplu)) %>%
  ggplot(aes(x=time, y=CP, col=Kplu)) +
  geom_line() +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()

################################################################################
################################################################################

################################################################################
#################################  Chunk 12 ####################################
################################################################################

## mrgsim.sa package ##

# The package mrgsim.sa simplifies the process of running local sensitivity analysis 
# https://github.com/kylebaron/mrgsim.sa. The package can do graphical sensitivity as 
# well as obtaining local sensitivity coefficients.
# What is the most ing=fluential parameter?

library(mrgsim.sa)

### sensitivity analysis

#set the output variable of interest
sensvar <- c("CP")  

# graphical sensitivity for each parameter
out_sens <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss) %>% 
  select_par(all_of(names(Kp[-11]))) %>% 
  parseq_fct(.n=3) %>% 
  sens_each(delta = 0.1, recsort=3, obsonly=TRUE, end = 12)

sens_plot(out_sens, "CP") 

# graphical sensitivity for a grid of parameters
out_sens_grid <- 
  modA %>% 
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss) %>% 
  parseq_cv(Kpmu, Kplu, .n=3) %>% 
  sens_grid(delta = 0.1, recsort=3, obsonly=TRUE, end = 12)

sens_plot(out_sens_grid, "CP")

# local sensitivity analysis 
out_lsa <- lsa(modA, var = "CP", par = names(all_of(Kp[-11])), events = ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss), end=12, delta=0.1, eps=1e-4)

lsa_plot(out_lsa, pal = NULL)

# get a summary of sensitivity coefficients
out_lsa_summ <- out_lsa %>% 
  group_by(p_name) %>%
  summarise(mean_sens = mean(sens)) %>%
  ungroup()

ggplot(data=out_lsa_summ, aes(x=reorder(p_name, mean_sens), y=mean_sens)) + 
  geom_col() + 
  labs(x="Parameter", y="Coefficient") +
  coord_flip() +
  geom_hline(yintercept = 0, lty=2) +
  theme_bw()

################################################################################
################################################################################

################################################################################
################################## Chunk 13 ####################################
################################################################################

library(sensobol)
library(future)
library(mrgsim.parallel)

# set parallelization options
nCores <- future::availableCores()
options(future.fork.enable=TRUE, mc.cores = nCores)
plan(multicore, workers = nCores)

# generate parameter sets
N <- 100  # initial sample number
mat <- sobol_matrices(N = N, params = c("Kpmu","Kplu","Kpli","Kpsp"))  # creates sample variates between 0-1
head(mat)

#Transform and groom
params <- unlist(as.list(param(modA))[c("Kpmu","Kplu","Kpli","Kpsp")])
umin <- params / 3
umax <- params * 3
umin
umax

# get samples by applying quantile to the uniform distributed matrix
mat <- as_tibble(mat)
mat2 <- imodify(mat, ~ qunif(.x, umin[[.y]], umax[[.y]]))
mat2 <- mutate(mat2, ID = row_number())
head(mat2)

# run the analysis
## serial
system.time(out <- modA %>% future_mrgsim_ei(e, mat2, nchunk = nCores, end = 1, delta = 1, Request = c("CP"), rtol = 1e-4, output = "df"))

# calculate Cmax
y <- out %>% group_by(ID) %>% summarise(cmax = max(CP)) %>% ungroup()
y <- as.numeric(y$cmax)

# get indices
ind <- sobol_indices(Y = y, N = N, params = names(params), boot = TRUE, R = 1000, first = "jansen")
ind.dummy <- sobol_dummy(Y = y, N = N, params = names(params), boot = TRUE, R = 1000)  # dummy to assess numerical approximation error

# plot indices
plot(ind, dummy = ind.dummy) + ylim(0,1)

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Parameter estimation ###

################################################################################
################################## Chunk 13 ####################################
################################################################################

# Use package `nloptr` (https://cran.r-project.org/web/packages/nloptr/index.html) 
# to optimize for most influential partition coefficient parameter and compare prediction 
# before and after optimization to observed data
# Use package `numDeriv` (https://cran.r-project.org/web/packages/numDeriv/index.html) 
# to generate the 95% CI around the parameter estimates
# For more background on optimization methods, check out Metrum's open science course 
# MI210 https://metrumrg.com/course/mi210-essentials-population-pk-pd-modeling-simulation/

library(nloptr)
library(numDeriv)
library(kableExtra)

### Do some data assembly to create dataset for optimization

#load observed IV infusion data
bw <- 73
obs <- read.csv("data/Adult_IV.csv")

# create an nm-tran dataset
nm_dv <- obs %>%
  select(time, dv=obs) %>%
  mutate(ID = 1,
         amt = 0,
         rate = 0,
         cmt = "VEN",
         evid = 0,
         ii = 0,
         addl = 0,
         ss = 0)

nm_dose <- nm_dv %>%
  slice(1) %>%
  mutate(amt = 4*bw,
         rate = 4*bw,
         ii = 12,
         addl = 13,
         ss = 1,
         dv = NA,
         evid = 1)

nm <- bind_rows(nm_dose, nm_dv) %>% arrange(time)

##set up objective function
OF <- function(pars, dat, pred=F){
  pars <- lapply(pars,exp)  #Get out of log domain for MLE
  pars <- as.list(pars)
  names(pars) <- names(theta)
  
  ## Get a prediction
  out <- modA %>% param(pars) %>% mrgsim_d(dat, carry_out=c("dv"), output="df")
  
  if(pred) return(out)
  
  ##OLS
  #return(sum((out$CP - out$dv)^2))
  
  ##maximum likelihood
  return(-1*sum(dnorm(log(out$dv),
                      mean=log(out$CP),
                      sd=pars$sigma, log=TRUE), na.rm = T))
}

# set initial values
theta <- log(c(Kpmu=2.94, sigma=1))  #initial parameter for MLE; mean and standard deviation

##Fit with nloptr package
##derivative-free optimizers
#fit <- neldermead(theta, OF)  #Nelder-Mead simplex
fit <- newuoa(theta, OF, dat=nm)  #New Unconstrained Optimization with quadratic Approximation

##gradient-basd optimizers
#fit <- tnewton(theta, OF)  #Local optimizer; Nelder-Mead simplex
#fit <- mlsl(theta, OF, lower=log(c(0.1)), upper=log(c(10)))  #global optimizer; multi-level single-linkage; takes very long ~ 10 min

##global optimizers
#fit <- direct(OF, lower=log(c(0.1, 0.1)), upper=log(c(10, 2)))  #DIviding RECTangles algorithm; takes ~ 3 min

fit

p <- as.list(exp(fit$par))  #get the parameters on the linear scale
names(p) <- names(theta)
p

# get standard error and confidence intervals around the estimated parameters
h <- hessian(OF, fit$par, dat=nm)
vc_log <- solve(h)  #variance-covariance matrix
SE_log <- sqrt(diag(vc_log))  #standard error on log scale

# create dataframe with parameter summary
sig <- function(x) signif(x, 3)
paramSumm <- tibble(Parameter = names(theta),
                    Estimate = exp(fit$par),
                    lb = exp(fit$par - (1.96 * SE_log)),
                    ub = exp(fit$par + (1.96 * SE_log))) %>%
  mutate(`Estimate (95% CI)` = paste0(sig(Estimate), " (", sig(lb), ", ", sig(ub), ")")) %>%
  select(Parameter, `Estimate (95% CI)`)

paramSumm %>%
  knitr::kable() %>%
  kable_styling()

# compare initial predictions to those using the optimized parameters
predB4 <- OF(theta, dat=nm, pred=T)
predAfter <- OF(fit$par, dat=nm, pred=T)

gp <- ggplot() +
  geom_point(data=obs, aes(x=time, y=obs)) + 
  geom_errorbar(data=obs, aes(x=time, y=obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data=predB4, aes(x=time, y=CP), lty=2) +
  geom_line(data=predAfter, aes(x=time, y=CP)) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") + 
  theme_bw()
gp


################################################################################
################################################################################

################################################################################
################################## Chunk 14 ####################################
################################################################################

# Simulate a 4 mg/kg voriconazole IV infusion dosing in a male child subject infused 
# at a rate of 3 mg/kg/h twice a day for seven days. Compare the steady state prediction 
# before and after optimization to the observed data in `Pediatric_IV.csv`.

obs <- read.csv("data/Pediatric_IV.csv")  #load observed data

wt <- 19  #adult body weight
amt <- 4*wt  
rate <- 3*wt
cmt <- "VEN"  #intravenous infusion
ii = 12
addl = 13
ss = 1

# simulate
simB4 <- as.data.frame(modP %>%
                         ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=1) %>% 
                         mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

simAfter <- as.data.frame(modP %>%
                            param(Kpmu=p$Kpmu) %>%
                            ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = simB4, aes(x=time, y=CP), lty=2) + 
  geom_line(data = simAfter, aes(x=time, y=CP)) + 
  labs(title="Pediatric 4 mg/kg IV", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Population simulation ###

################################################################################
################################## Chunk 15 ####################################
################################################################################

# Simulate a simple voriconazole dosing scenario for a population of 100 individuals, 
# 50 females and 50 males, with ages ranging between 20-80, weights ranging between 
# 50-100 kg and heights between 1.5 and 1.9 m. Use population saved in file 
# `data/popPars_100.rds`

# load population params
popPars <- readRDS("data/popPars_100.rds")

# add IIV on CL/VmaxH
set.seed(192898)
iVmaxH <- rlnorm(100, meanlog=log(40), sdlog=0.2)

# add to popPars
popPars2 <- lapply(1:length(popPars), function(i){
  pars <- c(popPars[[i]],list(VmaxH = iVmaxH[i]))
  return(pars)
})

# simulate
modA2 <- mread("models/voriPBPK2.mod")
Kpmu <- 0.56
modA2 <- param(modA2, Kpmu=Kpmu)

#set simulation conditions
bw   <- 73
amt  <- 4*bw
rate <- 4*bw
cmt  <- "VEN"
ii   <- 12
addl <- 13
ss   <- 1
delta <- 0.1
end <- 12

# prepare simulation dataset
idata <- popPars2 %>% bind_rows()

e <- ev(amt=amt, cmt=cmt, ii=ii, addl=addl, rate=rate, ss=ss)

data <- e %>% 
  as_tibble %>% 
  bind_cols(idata) %>%
  mutate(amt = 4*BW, rate = 4*BW) %>%
  select(ID, everything())

#run simulation
system.time(sims <- modA2 %>% 
              mrgsim_d(data=data, delta = delta, end = end, obsonly=T, outvars = c("CP"), output="df"))

# get summary stats for population prediction
hi95 <- function(x) quantile(x, probs = c(0.95))
lo05 <- function(x) quantile(x, probs = c(0.05))

sims2 <- sims %>%
  group_by(time) %>%
  mutate(loCP = lo05(CP),
         medCP = median(CP),
         hiCP = hi95(CP)) %>%
  ungroup() %>%
  filter(ID == first(ID))

#plot population predictions 
gp <- ggplot(data = sims2 %>% filter(row_number() != 1), aes(x=time)) + 
  geom_line(aes(y=medCP), col="black") +
  geom_ribbon(aes(ymin=loCP, ymax=hiCP), alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  labs(title="Population simulation", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

################################################################################
################################## Chunk 16 ####################################
################################################################################

# Use `PKNCA` package to calculate summary statistics for the population

library(PKNCA)

# create dose object
dose_obj <-
  PKNCAdose(
    data,
    amt~time|ID
  )

# create conc object
conc_obj <-
  PKNCAconc(
    sims,
    CP~time|ID
  )

# combine
data_obj <- PKNCAdata(conc_obj, dose_obj)

## Calculate the NCA parameters
results_obj <- pk.nca(data_obj)

## Summarize the results
summary(results_obj)


################################################################################
################################################################################

################################################################################
################################## Chunk 17 ####################################
################################################################################

# Compare doses 3 and 4 mg/kg for adults. Which dose keeps voriconazole concentration 
# at steady state above MIC (1 mg/L) for >= 90% of subjects?

library(cowplot)

data_3mg <- e %>% 
  as_tibble %>% 
  bind_cols(idata) %>%
  mutate(amt = 3*BW, rate = 3*BW) %>%
  select(ID, everything())

#run simulation
sims_3mg <- modA2 %>% 
  mrgsim_d(data=data_3mg, delta = delta, end = end, obsonly=T, outvars = c("CP")) %>% 
  filter(row_number() != 1)

sims2_3mg <- sims_3mg %>%
  group_by(time) %>%
  mutate(loCP = lo05(CP),
         medCP = median(CP),
         hiCP = hi95(CP)) %>%
  ungroup() %>%
  filter(ID == first(ID))

#plot population predictions 
gp_3mg <- ggplot(data = sims2_3mg, aes(x=time)) + 
  geom_line(aes(y=medCP), col="black") +
  geom_ribbon(aes(ymin=loCP, ymax=hiCP), alpha = 0.5) +
  scale_y_continuous(trans = "log10", limits=c(0.3,7)) +
  labs(title="Population simulation (3 mg/kg)", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw() + 
  geom_hline(aes(yintercept = 1), lty=2) 

gp_4mg <- gp + 
  labs(title="Population simulation (4 mg/kg)", x="time (h)", y="Plasma concentration (mg/L)") +
  geom_hline(aes(yintercept = 1), lty=2) +
  scale_y_continuous(trans = "log10", limits=c(0.3,7))

gp_34 <- plot_grid(gp_3mg, gp_4mg, ncol=2) 
gp_34

# calculate stats
stats_4mg <- sims %>%
  filter(time != 0) %>%
  group_by(ID) %>%
  mutate(Cmin = min(CP)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(FLG = ifelse(Cmin < 1, 1, 0))
percHi_4mg <- (nrow(stats_4mg[stats_4mg$FLG==0,]) / nrow(stats_4mg)) * 100

stats_3mg <- sims_3mg %>%
  filter(time != 0) %>%
  group_by(ID) %>%
  mutate(Cmin = min(CP)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(FLG = ifelse(Cmin < 1, 1, 0))
percHi_3mg <- (nrow(stats_3mg[stats_3mg$FLG==0,]) / nrow(stats_3mg)) * 100 

stats_df <- tibble(Dose = c("3 mg/kg","4 mg/kg"),
                   `Percent of subjects with SS Cmin > MIC` = c(percHi_3mg, percHi_4mg))

stats_df %>%
  knitr::kable() %>%
  kable_styling()

################################################################################
################################################################################

################################################################################
################################## Chunk 18 ####################################
################################################################################

## parallelization ##

# Use the package `mrgsim.apply` (vignette: https://github.com/kylebaron/mrgsim.parallel) 
# to run the same population simulation over parallel cores. How long does it take 
# compared to the single core simulation?

library(future)
library(mrgsim.parallel)

# set parallelization options
nCores <- future::availableCores()
options(future.fork.enable=TRUE, mc.cores = nCores)
plan(multicore, workers = nCores)

#run simulation
system.time(sims <- modA2 %>% 
              future_mrgsim_d(data=data, nchunk = nCores, delta = delta, end = end, obsonly=T, outvars = c("CP")) %>% 
              filter(row_number() != 1))

# get summary stats for population
hi95 <- function(x) quantile(x, probs = c(0.95))
lo05 <- function(x) quantile(x, probs = c(0.05))

sims2 <- sims %>%
  group_by(time) %>%
  mutate(loCP = lo05(CP),
         medCP = median(CP),
         hiCP = hi95(CP)) %>%
  ungroup() %>%
  filter(ID == first(ID))

#plot population predictions 
gp <- ggplot(data = sims2, aes(x=time)) + 
  geom_line(aes(y=medCP), col="black") +
  geom_ribbon(aes(ymin=loCP, ymax=hiCP), alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  labs(title="Population simulation - parallel", x="time (h)", y="Plasma concentration (mg/L)") +
  theme_bw()
gp

################################################################################
################################################################################

#------------------------------------------------------------------------------#

### Shiny interactive document ###

################################################################################
################################## Chunk 19 ####################################
################################################################################

# Explore the simple shiny app saved as `app.R`

library(shiny)
runApp("app.R")

################################################################################
################################################################################

################################################################################
################################## Chunk 20 ####################################
################################################################################

# Compile the mAb bamlanivimab mAb model and simulate a single 700 mg IV infusion dose infused over 2 hours. 
# What are the drug concentrations in plasma and lung on day 28? How do they compare to bamlanivimab IC90 (0.41481 mcg/mL)?

mod_mab <- mread("models/mAb.mod")

# set up simulation conditions
dose <- 700/150 #700 mg to umol
dur <- 2
rate <- dose/dur
cmt <- 4
end <- 28*24
e <- ev(amt=dose, rate=rate, cmt=cmt)

sim_mab <- mod_mab %>%
  mrgsim_e(e, end=end, outvars=c("Cexg_Plasma", "Cexg_Lung_IS")) %>%
  as_tibble() 

sim_mab2 <- sim_mab %>%
  mutate(time = time/24,
         Plasma = Cexg_Plasma*150,
         Lung = Cexg_Lung_IS*150) %>%
  select(-Cexg_Plasma, -Cexg_Lung_IS) %>%
  gather(tissue, conc, -ID, -time)

# plot
p_mab <- ggplot(data=sim_mab2, aes(x=time, y=conc, col=tissue)) +
  geom_line(lwd=1) +
  geom_hline(yintercept = 0.41481, lty=2) +
  scale_x_continuous(breaks=seq(0, 28, 7)) +
  scale_y_continuous(trans = "log10") +
  labs(x="Time (d)", y=expression("mAb concentration ("*mu*"g/mL)")) +
  theme_bw()
p_mab

################################################################################
################################################################################

################################################################################
################################## Chunk 21 ####################################
################################################################################

# Use the mAb model as a template to add a simple tumor model. 
# Save a new model file, compile, run a quick simulation, and plot mAb tumor 
# concentration-time profile for 28 days.

lines_mod <- read_lines("models/mAb.mod")

lines_muscle <- lines_mod[str_detect(lines_mod, "Muscle")]
lines_tumor <- str_replace_all(lines_muscle, "Muscle", "Tumor")

# tumor-specific modifications
## param
lines_tumor[str_detect(lines_tumor, "^V_Tumor =")] <- "V_Tumor = 0.02 // [L] Total Volume"
lines_tumor[str_detect(lines_tumor, "^V_Tumor_IS =")] <- "V_Tumor_IS = 0.55*0.02  // [L] Interstitial Volume"
lines_tumor[str_detect(lines_tumor, "^V_Tumor_V =")] <- "V_Tumor_V = 0.07*0.02 // [L] Vascular Volume"
lines_tumor[str_detect(lines_tumor, "^PLQ_Tumor =")] <- "PLQ_Tumor = 12.7*0.02 // [L/h] blood flow rate" 
lines_tumor[str_detect(lines_tumor, "^Endothelial_Cell_Frac_Tumor =")] <- "Endothelial_Cell_Frac_Tumor = 0.005"
lines_tumor[str_detect(lines_tumor, "^SV_Tumor =")] <- "SV_Tumor = 0.842"
lines_tumor[str_detect(lines_tumor, "^SIS_Tumor =")] <- "SIS_Tumor = 0.2"

## MAIN
lines_tumor[str_detect(lines_tumor, "^double FcRn_Conc =")] <- ""

## ODE
lines_tumor[str_detect(lines_tumor, "\\(1.0 - SIS_Tumor\\)")] <- "" 
lines_tumor[str_detect(lines_tumor, "\\+ \\(PLQ_Tumor - LF_Tumor\\)\\*Cexg_Tumor_V")] <- ""  # remove the Tumor entry to Aexg; will be patched later

# original model modifications
## PARAM
lines_mod[str_detect(lines_mod, "^PLQ_Lung = 181.913000000")] <- "PLQ_Lung = 181.9130000000000109 + (12.7*0.02) // [L/h]"

## MAIN
lines_mod[str_detect(lines_mod, "^double FcRn_Conc =")] <- paste0(str_remove(lines_mod[str_detect(lines_mod, "^double FcRn_Conc =")], "\\);"), "+V_endosomal_Tumor);")

## ODE
lines_mod[str_detect(lines_mod, "\\- L_LymphNode\\*Cedg_LN\\)/V_LN;")] <- "- L_LymphNode*Cedg_LN + (1.0 - SIS_Tumor)*LF_Tumor*Cedg_Tumor_IS)/V_LN;"
lines_mod[str_detect(lines_mod, "\\- L_LymphNode\\*Cexg_LN\\)/V_LN;")] <- "- L_LymphNode*Cexg_LN + (1.0 - SIS_Tumor)*LF_Tumor*Cexg_Tumor_IS)/V_LN;"
lines_mod[str_detect(lines_mod, "\\+ L_LymphNode\\*Cexg_LN\\);")] <- "+ L_LymphNode*Cexg_LN + (PLQ_Tumor - LF_Tumor)*Cexg_Tumor_V);"

# join lines
lines_mod_new <- 
  # PROB, SET, and PARAM
  lines_mod[1:str_which(lines_mod, "\\[ MAIN \\]") - 1] %>%  
  append(lines_tumor[1:str_which(lines_tumor, "SIS_Tumor = 0.2")]) %>%
  # MAIN
  append(lines_mod[str_which(lines_mod, "\\[ MAIN \\]"):(str_which(lines_mod, "\\[ CMT \\]") - 1)]) %>% 
  append(lines_tumor[str_which(lines_tumor, "double LF_Tumor = PLQ_Tumor\\*0.002"):str_which(lines_tumor, "CFcRn_Tumor_IM_0 = FcRn_Conc\\*1e-4;")]) %>%
  # CMT
  append(lines_mod[str_which(lines_mod, "\\[ CMT \\]"):(str_which(lines_mod, "\\[ ODE \\]") - 1)]) %>%
  append(lines_tumor[str_which(lines_tumor, "^Cedg_Tumor_V$"):str_which(lines_tumor, "^CFcRn_Tumor_IM$")]) %>%
  # ODE
  append(lines_mod[str_which(lines_mod, "\\[ ODE \\]"):(str_which(lines_mod, "\\[ CAPTURE \\]") - 1)]) %>%
  append(lines_tumor[str_which(lines_tumor, "dxdt_Cexg_Tumor_V ="):str_which(lines_tumor, "2.0\\*kdeg_FcRn_Ab\\*Cedg_Tumor_b2IM\\)\\)/V_Tumor_IM;")]) %>%
  # CAPTURE
  append(lines_mod[str_which(lines_mod, "\\[ CAPTURE \\]"):length(lines_mod)])

# save
write_lines(lines_mod_new, "mAb_tumor.mod")

# compile
mod_mab_tumor <- mread("mAb_tumor.mod")

# run simulation

# set up simulation conditions
dose <- 700/150 #700 mg to umol
dur <- 2
rate <- dose/dur
cmt <- 4
end <- 28*24
e <- ev(amt=dose, rate=rate, cmt=cmt)

sim_mab_tumor <- mod_mab_tumor %>%
  mrgsim_e(e, end=end, outvars=c("Cexg_Plasma", "Cexg_Tumor_IS")) %>%
  as_tibble() 
sim_mab_tumor2 <- sim_mab_tumor %>%
  mutate(time = time/24,
         Plasma = Cexg_Plasma*150,
         Tumor = Cexg_Tumor_IS*150) %>%
  select(-Cexg_Plasma, -Cexg_Tumor_IS) %>%
  gather(tissue, conc, -ID, -time)

# plot
p_mab_tumor <- ggplot(data=sim_mab_tumor2 %>% filter(tissue != "Lung"), aes(x=time, y=conc, col=tissue)) +
  geom_line(lwd=1) +
  scale_x_continuous(breaks=seq(0, 28, 7)) +
  scale_y_continuous(trans = "log10") +
  labs(x="Time (d)", y=expression("mAb concentration ("*mu*"g/mL)")) +
  theme_bw()
p_mab_tumor

################################################################################
################################################################################


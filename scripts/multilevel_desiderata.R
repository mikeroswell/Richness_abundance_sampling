
#multilevel model desiderata: fit RE by group (e.g., treatment) s.t. I esitmate
#different RE variances for different groups. Pretty sure lme4 doesn't think
#like that, but this seems easy to do in a language like STAN?
library(ggplot2)
library(dplyr)
library(RTMB)
# imagine 5 tx levels
# 15 species
# n is observed species abundance
# species ID is nested within treatment
ntx <- 5
nsp <- 35
sigma <- seq(1.1, 1.5, 0.1)
mu <- seq(1.2, 2.0, 0.2)
ss = 75
ntst <- rnorm(nsp * ntx, mean = rep(mu, each = nsp), sd = rep(sigma, each = nsp))

dd <- data.frame(
  sp = factor(rep(1:nsp,  ntx*ss))
  , tx = factor(rep(1:ntx, each = ss*nsp))
 ,  n = rpois(ntx*nsp*ss, lambda = exp(rep(ntst, each = ss)))
)

dd %>% 
  ggplot(aes(tx, n, color = sp)) + 
  geom_jitter() + 
  theme_classic() + 
  scale_y_log10()

# want to model expected abundance with a random effect of species ID; where 
# the RE variance varies by tx

sigvar_mod <- lme4::glmer(n ~  (1|tx/sp)
                          , family = poisson
                         , data = dd)

# I don't think that's what happened here. 
summary(sigvar_mod)

# So now I try with RTMB

library(RTMB)
# vignette("RTMB-introduction")
# vignette("RTMB-advanced")
# define parameters
parameters <- list(
  muj = rep(0, ntx) # site intercept
  , sdj = rep(1, ntx) # site sd in intensity parameter
  , b = matrix(rep(0, ntx*nsp), nrow = nsp, ncol = ntx) # species specific intensity parameter (latent)
)

# objective function that takes parameters as input
# data is hard-coded in (for now)
oF <- function(parms){
  getAll(dd, parms, warn = TRUE)
  # set up prediction features of RTMB with the OBS function
  OBS(n) # this is response variable

  ## Random intensity
  b %~% dnorm(mean=muj, sd=sdj)
 
  ## Data
  predLogNum <- b[sp,tx]
  n %~% dpois(exp(predLogNum))
  ## Get predicted weight uncertainties
  ADREPORT(predLogNum)
}
################################
## NOT RUN: Takes a long time ##
################################ 

## process objective function into model
# obj <- MakeADFun(oF, parameters, random =c("b"))
## that took a while. 

## try fitting the model

# opt <- nlminb(obj$par, obj$fn, obj$gr)
# 
# opt

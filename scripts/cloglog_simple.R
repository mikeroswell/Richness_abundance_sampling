# beware count models for richness;
# hazard models show promise. 

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
# library(lme4)

# simulate data
set.seed(101823)
# set # of reps
reps <- 999

# set all spp with same abundance
snn <- 10

# nsp sets richness in most samples
nsp <- 100
nsr <- 200 # double diversity for one smaple

# create three sample sizes... large, half = medium, third = small
lg <- 450
md <-lg/2
sm <- lg/3
xt <- md # medium sample for the xtra diverse comm.


# random sampling
ab_lg <- rmultinom(reps, lg, rep(snn, nsp))
ab_md <- rmultinom(reps, md, rep(snn, nsp))
ab_sm <- rmultinom(reps, sm, rep(snn, nsp))
ab_xt <- rmultinom(reps, md, rep(snn, nsr)) # double richness in "xt"

# format into data frame
dd <- map_dfr(c("sm", "md", "lg", "xt"), function(samp){
     data.frame(abundance = apply(get(paste("ab", samp, sep = "_"))
                                 , 2, sum)
               , species = apply(get(paste("ab", samp, sep = "_"))
                                 , 2, function(x){sum(x>0)})
               , samp = samp
    ) %>% mutate(true_rich = ifelse(samp == "xt", 200, 100)
                 , rperc = species/true_rich)
  })


# look at richness across replicate samples
dd %>% 
  ggplot(aes(samp, species)) +
  geom_violin() +
  scale_y_log10(limits = c(1, NA), n.breaks = 6) +
  theme_classic() +
  labs(x = "sample")

# look at sample size (commonly, "abundance")
dd %>% 
  ggplot(aes(samp, abundance)) +
  geom_violin() +
  scale_y_log10(limits = c(1, NA), n.breaks = 6) +
  theme_classic() +
  labs(x = "sample")

# fit some glms
# proof of concept: count model fine with abundance:
ab_mod <- glm(abundance ~ samp
              , data = dd
              , family = "poisson")
# summary(ab_mod)
exp(coef(ab_mod)) # lg as baseline, 450 individuals, md is 1/2, sm is 1/3

# can't fit the negative binomial as theta ~ Inf here (reasonably!)
# nb_mod <- glm.nb(species ~ samp, data = dd)
# poisson model fits fine; does it make sense?  

pois_mod <- glm(species ~ samp
                , family = "poisson"
                , data = dd)
#  summary(pois_mod)
exp(coef(pois_mod)) # not sure what those numbers mean...


# one issue is that basically all spp. possible detected in large sample,
# but only a fraction in the others. 

# Maybe we do better to think about fraction of spp observed, logistic
logi_mod <- glm(rperc ~ samp 
                        , data = dd
                        , family = "binomial"
                        , weights = true_rich)
  

# summary(logi_mod)
exp(coef(logi_mod))
exp(coef(pois_mod))

# go to hazard models, cloglog

# defining link and inverse
# cloglog <- function(x) log(-log(1-x))
# gompertz <- function(x) 1-exp(-exp(x))

# fit model: similar to logistic but with cloglog link instead of logit
species_cloglog <- glm(rperc ~ samp 
                       , data = dd
                       , family = binomial(link = cloglog)
                       , weights = true_rich
                       )

# summary(species_cloglog)

exp(coef(species_cloglog)) # nearly getting sampling effect right, and noticing
# extra diversity in sample xt!

# that is, half the sampling and _double_ diversity leads to effect size ~1/4


# and, if we consider number of individuals sampled to be exposure, 
# we can include it as an offset. 
# this is nice b/c it's not of direct interest, and is directly measured/known

species_cloglog_off <-  glm(rperc ~ samp + offset(log(abundance))
                            , data = dd
                            , family = binomial(link = cloglog)
                            , weights = true_rich)

# summary(species_cloglog_off) 
exp(coef(species_cloglog_off)) 
# and now we see no effect of sampling at all , but do see the effect of
# diversity in sample xt!



# show breakdown with heterogeneity in species abundances

#################################
# here, lognormal SAD
sim_ab <- MeanRarity::fit_SAD(rich = 100, simpson = 40)$rel_abundances
# here, double the diversity, still lognormal SAD
sim_xt <- MeanRarity::fit_SAD(rich = 200, simpson = 80)$rel_abundances

# remake sample data
ab_lg <- rmultinom(reps, lg, sim_ab)
ab_md <- rmultinom(reps, md, sim_ab)
ab_sm <- rmultinom(reps, sm, sim_ab)
ab_xt <- rmultinom(reps, md, sim_xt)

dd2 <- map_dfr(c("sm", "md", "lg", "xt"), function(samp){
  data.frame(abundance = apply(get(paste("ab", samp, sep = "_"))
                               , 2, sum)
             , species = apply(get(paste("ab", samp, sep = "_"))
                               , 2, function(x){sum(x>0)})
             , samp = samp
  ) %>% mutate(true_rich = ifelse(samp == "xt", 200, 100)
               , rperc = species/true_rich)
})

##################################

# does cloglog work now?
species_cloglog <- glm(rperc ~ samp 
                       , data = dd2
                       , family = binomial(link = cloglog)
                       , weights = true_rich)

#  summary(species_cloglog)

exp(coef(species_cloglog)) # these aren't the right ratios!

# does offset save us?

species_cloglog_off <-  glm(rperc ~ samp   + offset(log(abundance))
                            , data = dd2
                            , family = binomial(link = cloglog)
                            , weights = true_rich)
#  summary(species_cloglog_off) 
exp(coef(species_cloglog_off)) #pushing in the right direction??

# thought question: does offset have to be number of individuals?
# would generalizing to "effort" violate any key assumptions? 

# [NOT RUN] does random effect solve the problem?
species_glmm <-  glmer(rperc ~ samp   + offset(log(abundance)) + (1|sp)
                            , data = dd2
                            , family = binomial(link = cloglog)
                            , weights = true_rich)


# next step wd be generating data that better match assumptions for the the
# random effects using Dirichlet, Logistic-Normal, or other approach to model
# species-level heterogeneity.
#Simplify cloglog story 

library(dplyr)
library(tidyr)
library(ggplot2)

# simulate data
set.seed(101823)
# set # of reps
reps <- 999

# set all spp with same abundance
snn <- 10
nsp <- 100

# create two samples, with one 3x as many individuals as other
lg <- 450
md <-lg/2
sm <- lg/3

ab_lg <- rmultinom(reps, lg, rep(snn, nsp))
ab_md <- rmultinom(reps, md, rep(snn, nsp))
ab_sm <- rmultinom(reps, sm, rep(snn, nsp))

dd <- map_dfr(c("sm", "md", "lg"), function(samp){
     data.frame(abundance = apply(get(paste("ab", samp, sep = "_"))
                                 , 2, sum)
               , species = apply(get(paste("ab", samp, sep = "_"))
                                 , 2, function(x){sum(x>0)})
               , samp = samp
    ) %>% mutate(rperc = species/100)
  })


dd %>% 
  ggplot(aes(samp, species))+
  geom_violin()+
  scale_y_log10(limits = c(1, NA), n.breaks = 6) +
  theme_classic()


dd %>% 
  ggplot(aes(samp, abundance)) +
  geom_violin() +
  scale_y_log10(limits = c(1, NA), n.breaks = 6) +
  theme_classic() 

# fit some glms
# proof of concept with abundance:
ab_mod <- glm(abundance ~ samp
              , data = dd
              , family = "poisson")
summary(ab_mod)
exp(coef(ab_mod)) # lg as baseline, 450 individuals, md is 1/2, sm is 1/3

# can't fit the negative binomial as theta ~ Inf here (reasonably!)
# nb_mod <- glm.nb(species ~ samp, data = dd)
# poisson model fits fine, but might be a bit weird:  

pois_mod <- glm(species ~ samp, family = "poisson", data = dd)
summary(pois_mod)
exp(coef(pois_mod))

# think about fraction of spp observed, logistic
logi_mod <- glm(rperc ~ samp 
                        , data = dd
                        , family = "binomial"
                        , weights = rep.int(100, length(dd[,1])))

summary(logi_mod)
exp(coef(logi_mod))
exp(coef(pois_mod))

# go to hazard models, cloglog
#defining link and inverse
cloglog <- function(x) log(-log(1-x))
gompertz <- function(x) 1-exp(-exp(x))


species_cloglog <- glm(rperc ~ samp 
                       , data = dd
                       , family = binomial(link = cloglog)
                       , weights = rep.int(100, length(dd[,1])))


summary(species_cloglog)

exp(coef(species_cloglog)) # nearly getting sampling effect right!


# and, if we consider number of individuals sampled to be exposure, 
# we can include it as an offset. 

species_cloglog_off <-  glm(rperc ~ samp   + offset(log(abundance))
                            , data = dd
                            , family = binomial(link = cloglog)
                            , weights = rep.int(100, length(dd[,1])))
exp(coef(species_cloglog_off))
summary(species_cloglog_off) # and now we see no effect of sampling at all!


# next step wd be the random effects using Dirichlet, Logistic-Normal, or other approach to model species-level heterogeneity. 
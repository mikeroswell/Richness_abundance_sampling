# Here, I simulate an overall species abundance distribution, and then sample
# from it more intensively in the "treatment" ("tx") than in the control ("ct").
# I look at the upper and lower chunk of the overall abundance distribution,
# defining these chunks as comprised of common and rare spp. respectively.

# I want to model two versions of the treatment effect: the relative difference
# in the number of *species* and of *individuals* of each group, i.e., of the
# "rare" spp and the "common" spp, between "tx" and "ct".

# I feel like (since nothing biological is happening that moderates the
# treatment effect by how rare or common species are), I shouldn't see
# differences in how the rare and common groups of spp. respond to the
# treatment. But when I fit negative binomial models to the counts of species
# (and not to the counts of individuals), it looks like there is a **stronger**
# effect for the **rare** spp than the common ones.

# I think this is expected in some sense... for example, the probability that a
# species is not detected is high for rare spp and low for common spp, so it's
# relatively easy to imagine how a rare species is not detected in the control
# but is in the treatment. It's less likely a common species is absent at all.

# but I don't see why this isn't already accounted for in a negative binomial
# glm.

# do you???


# install.packages(MASS)
library(tidyverse)

set.seed(3333)


# get 400 regional spp
rich <- 400

# create a species abundance distribution. Details aren't that important here;
# key is that some species are very rare and others very common, going with a 
# vaguely realistic distributionn

SAD <- rnbinom(rich*2, size = 0.5, mu = 100) # start with more and drop the 
# things that are practically or actually 0. 
SAD <- (SAD[order(SAD, decreasing = TRUE)]/sum(SAD))[1:rich]

# plot(SAD)

# create a "treatment" effect. Literally just sampling more individuals from the 
# SAD in tx than in control ("ct")

ab_ct <- 100

# iterate
nrep <- 999 # replicates
effect_by_effect <- map_dfr(seq(1.5, 3.5, 0.25), function(tx_eff){
  ab_tx <- ab_ct*tx_eff

  # sample SAD
  tx_obs <- rmultinom(nrep, ab_tx, prob = SAD )
  ct_obs <- rmultinom(nrep, ab_ct, prob = SAD )
  
  # individuals of common, then rare spp in each tx
  common_tx_ab <- apply(tx_obs[1:80,], 2, sum)
  common_ct_ab <- apply(ct_obs[1:80,], 2, sum)
  rare_tx_ab <- apply(tx_obs[220:300,], 2, sum)
  rare_ct_ab <- apply(ct_obs[220:300,], 2, sum)
  
  # number of unique common or rare spp found in teach treatment
  pos <- function(x){sum(x>0)}
  # chao1 estimator
  mychao <- function(x){iNEXT::ChaoRichness(x)$Estimator}
  
  common_tx_r <- apply(tx_obs[1:80,], 2, pos)
  common_ct_r <- apply(ct_obs[1:80,], 2, pos)
  rare_tx_r <- apply(tx_obs[220:300,], 2, pos)
  rare_ct_r <- apply(ct_obs[220:300,], 2, pos)
  
  # proportion of rare or common spp found in each treatment
  common_tx_rp <- apply(tx_obs[1:80,], 2, pos)/80
  common_ct_rp <- apply(ct_obs[1:80,], 2, pos)/80
  rare_tx_rp <- apply(tx_obs[220:300,], 2, pos)/80
  rare_ct_rp <- apply(ct_obs[220:300,], 2, pos)/80
  
  # chao expected proportions
  common_tx_pc <- apply(tx_obs[1:80,], 2, mychao)/80
  common_ct_pc <- apply(ct_obs[1:80,], 2, mychao)/80
  rare_tx_pc <- apply(tx_obs[220:300,], 2, mychao)/80
  rare_ct_pc <- apply(ct_obs[220:300,], 2, mychao)/80
  
  # par(mfrow=c(2,1))
  # hist(common_tx_ab, xlim = c(0, max(common_tx_ab)*1.1))
  # hist(common_ct_ab, xlim = c(0, max(common_tx_ab)*1.1))
  # 
  # 
  # hist(rare_tx_ab, xlim = c(0, max(rare_tx_ab)*1.1))
  # hist(rare_ct_ab, xlim = c(0, max(rare_tx_ab)*1.1))
  # 
  # hist(common_tx_r, xlim = c(0, max(common_tx_r)*1.1))
  # hist(common_ct_r, xlim = c(0, max(common_tx_r)*1.1))
  # 
  # 
  # hist(rare_tx_r, xlim = c(0, max(rare_tx_r)*1.1))
  # hist(rare_ct_r, xlim = c(0, max(rare_tx_r)*1.1))
  
  dd <- bind_rows(list("common_tx_r"= common_tx_r
                       , "common_tx_rp"= common_tx_rp
                       , "common_tx_pc"= common_tx_pc
                       , "common_tx_ab"= common_tx_ab
                       , "common_ct_r" = common_ct_r
                       , "common_ct_rp" = common_ct_rp
                       , "common_ct_pc" = common_ct_pc
                       , "common_ct_ab"= common_ct_ab
                       , "rare_tx_r" = rare_tx_r
                       , "rare_tx_rp" = rare_tx_rp
                       , "rare_tx_pc" = rare_tx_pc
                       , "rare_tx_ab" = rare_tx_ab
                       , "rare_ct_r" = rare_ct_r
                       , "rare_ct_rp" = rare_ct_rp
                       , "rare_ct_pc" = rare_ct_pc
                       , "rare_ct_ab"= rare_ct_ab)) %>% 
    pivot_longer(cols = 1:12 
                 , names_to = c("rarity", "treatment", "response")
                 , names_sep = "_")
  
  
  
  # Why don't I see the same treatment effect (including same effect size) for 
  # rare and common species, and why doesn't the negative binomial model account 
  # for the lower abundance of rare spp?
  # summary(MASS::glm.nb(formula = value~treatment
  #                      , data = dd %>% filter(response == "r" & rarity == "common") ))
  # summary(MASS::glm.nb(formula = value~treatment
  #                      , data = dd %>% filter(response == "r" & rarity == "rare") ))
  # 
  # # do logistic regressions instead
  # summary(glm(formula = value~treatment
  #                      , data = dd %>% filter(response == "rp" & rarity == "common")
  #             , family = "binomial"))
  # summary(glm(formula = value~treatment
  #                      , data = dd %>% filter(response == "rp" & rarity == "rare")
  #             , family = "binomial"))
  
  # now test interaction
  # br <- betareg::betareg(formula = value~treatment*rarity
  #                        , data = dd %>% filter(response == "rp" )
  #                        , link = "logit")
  # 
  # summary(br)
  mymod <- glm(formula = value~ treatment * rarity 
               , data = dd %>% 
                 filter(response == "pc" ) %>% 
                 mutate(value = ifelse(value >1, 0.9999999999999999, value))
                        
               , family = "binomial")
  modsum <- emmeans::contrast(emmeans::emmeans(mymod, ~treatment * rarity)
                              , interaction = "trt.vs.ctrl")
  
  data.frame(emm = modsum, sampling_effect_size = tx_eff)
})

effect_by_effect %>% 
  ggplot(aes(sampling_effect_size, emm.estimate, color = emm.p.value)) + 
  geom_point() +
  theme_classic() +
  labs(x = "treatment effect:\n (sampling intensity in treatment)/\n(sampling intensity in control)"
       , y = "logistic model estimate for \n treatment by rarity_group interaction")

# Tx effect is stronger for rare spp. It's more likely a rare species will be 
# both absent in the control and present in the treatment.




# of course, we expect the same treatment effect when looking only at 
# individuals (and not collapsing them into spp)
# summary(MASS::glm.nb(formula = value~treatment
#                      , data = dd %>% filter(response == "ab" & rarity == "common") ))
# summary(MASS::glm.nb(formula = value~treatment
#                      , data = dd %>% filter(response == "ab" & rarity == "rare") ))

# not much difference in mean Tx effect. SE obviously still higher for rare spp. 

# more like real model 
# summary(MASS::glm.nb(formula = value~treatment*rarity
#                      , data = dd %>% filter(response == "r" )))
# summary(MASS::glm.nb(formula = value~treatment*rarity
#                      , data = dd %>% filter(response == "ab" )))
# 
# 


# # add a sampling effect
# ab_tx <- 3000
# ab_ct <-1000
# tx_true <- rmultinom(nrep, ab_tx, prob = SAD)
# ct_true <- rmultinom(nrep, ab_ct, prob = SAD )
# tx_obs <- apply(tx_true, 2, function(x){rbinom(rich, size = x, prob= 0.2 )})
# ct_obs <- apply(ct_true, 2, function(x){rbinom(rich, size = x, prob= 0.2 )})
# common_tx_ab <- apply(tx_obs[1:80,], 2, sum)
# common_ct_ab <- apply(ct_obs[1:80,], 2, sum)
# rare_tx_ab <- apply(tx_obs[220:300,], 2, sum)
# rare_ct_ab <- apply(ct_obs[220:300,], 2, sum)
# 


# common_tx_r <- apply(tx_obs[1:80,], 2, pos)
# common_ct_r <- apply(ct_obs[1:80,], 2, pos)
# rare_tx_r <- apply(tx_obs[220:300,], 2, pos)
# rare_ct_r <- apply(ct_obs[220:300,], 2, pos)
# 
# par(mfrow=c(2,1))
# hist(common_tx_ab, xlim = c(0, max(common_tx_ab)*1.1))
# hist(common_ct_ab, xlim = c(0, max(common_tx_ab)*1.1))
# 
# 
# hist(rare_tx_ab, xlim = c(0, max(rare_tx_ab)*1.1))
# hist(rare_ct_ab, xlim = c(0, max(rare_tx_ab)*1.1))
# 
# hist(common_tx_r, xlim = c(0, max(common_tx_r)*1.1))
# hist(common_ct_r, xlim = c(0, max(common_tx_r)*1.1))
# 
# 
# hist(rare_tx_r, xlim = c(0, max(rare_tx_r)*1.1))
# hist(rare_ct_r, xlim = c(0, max(rare_tx_r)*1.1))
# 
# # remove abundance effect of treatment. add sampling effect by treatment
# 
# ab_tx <- 3000
# ab_ct <- 3000
# tx_true <- rmultinom(nrep, ab_tx, prob = SAD)
# ct_true <- rmultinom(nrep, ab_ct, prob = SAD )
# tx_obs <- apply(tx_true, 2, function(x){rbinom(rich, size = x, prob= 0.2 )})
# ct_obs <- apply(ct_true, 2, function(x){rbinom(rich, size = x, prob= 0.1 )})
# common_tx_ab <- apply(tx_obs[1:80,], 2, sum)
# common_ct_ab <- apply(ct_obs[1:80,], 2, sum)
# rare_tx_ab <- apply(tx_obs[220:300,], 2, sum)
# rare_ct_ab <- apply(ct_obs[220:300,], 2, sum)
# 
# pos <- function(x){sum(x>0)}
# common_tx_r <- apply(tx_obs[1:80,], 2, pos)
# common_ct_r <- apply(ct_obs[1:80,], 2, pos)
# rare_tx_r <- apply(tx_obs[220:300,], 2, pos)
# rare_ct_r <- apply(ct_obs[220:300,], 2, pos)
# 
# par(mfrow=c(2,1))
# hist(common_tx_ab, xlim = c(0, max(common_tx_ab)*1.1))
# hist(common_ct_ab, xlim = c(0, max(common_tx_ab)*1.1))
# 
# 
# hist(rare_tx_ab, xlim = c(0, max(rare_tx_ab)*1.1))
# hist(rare_ct_ab, xlim = c(0, max(rare_tx_ab)*1.1))
# 
# hist(common_tx_r, xlim = c(0, max(common_tx_r)*1.1))
# hist(common_ct_r, xlim = c(0, max(common_tx_r)*1.1))
# 
# 
# hist(rare_tx_r, xlim = c(0, max(rare_tx_r)*1.1))
# hist(rare_ct_r, xlim = c(0, max(rare_tx_r)*1.1))
# 
# 
# 
# 
# dd <- bind_rows(list("common_tx_r"= common_tx_r, "common_tx_ab"= common_tx_ab, "common_ct_r" = common_ct_r, "common_ct_ab"=common_ct_ab, "rare_tx_r" = rare_tx_r, "rare_tx_ab" = rare_tx_ab, "rare_ct_r" =rare_ct_r, "rare_ct_ab"= rare_ct_ab)) %>% pivot_longer(cols = 1:8 , names_to = c("rarity", "treatment", "response"), names_sep = "_")
# 
# 
# # mods <- map(unique(dd$response), function(resp){
# #   map(unique(dd$rarity), function(rar){
# #     mod_dat <- dd %>% filter(response==resp & rarity == rar)
# #     MASS::glm.nb(formula = value~treatment, data = mod_dat )
# #   })
# # })
# #   
# # map(flatten(mods), summary)
# 
# summary(MASS::glm.nb(formula = value~treatment, data = dd %>% filter(response == "r" & rarity == "common") ))
# summary(MASS::glm.nb(formula = value~treatment, data = dd %>% filter(response == "r" & rarity == "rare") ))
# 
# # Tx effect is stronger for rare spp
# 
# summary(MASS::glm.nb(formula = value~treatment, data = dd %>% filter(response == "ab" & rarity == "common") ))
# summary(MASS::glm.nb(formula = value~treatment, data = dd %>% filter(response == "ab" & rarity == "rare") ))
# 
# # no difference in Tx effect
# 
# # more like real model 
# 
# summary(MASS::glm.nb(formula = value~treatment*rarity, data = dd %>% filter(response == "ab" )))
# 
# summary(MASS::glm.nb(formula = value~treatment*rarity, data = dd %>% filter(response == "r" )))
# 
# 
# 

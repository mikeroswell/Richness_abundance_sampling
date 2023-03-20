# simulate CIG like data to see if my intuition about abundance effects being
# cryptic but driving richness effects makes any sense.

# devtools::install_github("mikeroswell/MeanRarity")
# install.packages("glmmTMB")
# install.packages("emmeans")
library(tidyverse)

set.seed(3)


# get 400 regional spp
rich <- 400
SAD <- MeanRarity::fit_SAD(rich = rich, simpson = 25)

#SAD$rel_abundances
MeanRarity::radplot(SAD$rel_abundances)

# create 16 sites. 
n_sites <- 16
# No average Tx effect besides detection
# sample spp. independently at each site
# variation in total abundance @ site and Tx levels
# make tx effect vary by sp, and make this only a detection effect

dd <- map_dfr(1:n_sites, function(site){
  mu <- sample(1e3:5e3, 1) # big site or small site?
  map_dfr(c("tx", "ctrl"), function(treat){
    # TRUE abundance varies by treatment, but not systematically
    sz <- rnbinom(1, size = 40, mu = mu)
    print(sz)
    site_abs <- c(rmultinom(n= 1, size = sz, prob = SAD$rel_abundances))
    print(paste("tot = ",  sum(site_abs)))
    print(paste("tot_rich = ", sum(site_abs>0)))
    n_obs <- rbinom(n = length(site_abs), size = site_abs, prob = 
                      # Tx effect is a per-species detection difference
                      # as implemented, it is not conserved among sites
                      ifelse(treat == "tx"
                             , sample(c(0.09, 0.02), size = rich, prob = c(7,2), replace = TRUE)
                             , sample(c(0.09, 0.02), size = rich, prob = c(2,7), replace = TRUE)
                      )
                  )
    print(paste("nobs = ", sum(n_obs)))
    print(paste("rich obs = ", sum(n_obs>0)))
    data_frame(site, treatment =treat, sp = 1:rich, site_abs, n_obs)
  })
})

dd <- dd %>% 
  mutate(rarity = ifelse(sp <80, "common" # this is like CIG
                         , ifelse(sp >300, "rare" # sounds like more but I suspect many spp are just 0
                                  , "intermediate")))

summ <- dd %>% group_by(site, treatment, rarity) %>% 
  summarize(ab =  sum(n_obs), rich = sum(n_obs>0)) 


summ %>% ggplot(aes(ab, rich)) + 
  geom_point() +
  theme_classic()

summ %>% 
  pivot_longer(cols = c("ab", "rich"), names_to= "response", values_to = "val") %>% 
  ggplot(aes(treatment, val, color = rarity)) +
  geom_boxplot() +
  facet_wrap(~response) +
  scale_y_log10() + 
  theme_classic()


richmod <- glmmTMB::glmmTMB(val ~ treatment * rarity + (1|site)
                            , family = glmmTMB::nbinom2
                            , data = summ %>% 
                              pivot_longer(cols = c("ab", "rich"), names_to= "response", values_to = "val") %>% 
                              filter(response == "rich" & rarity != "intermediate"))  

richmod <- glmmTMB::glmmTMB(val ~ treatment * rarity + (1|site)
                            , family = 
                            , data = summ %>% 
                              pivot_longer(cols = c("ab", "rich"), names_to= "response", values_to = "val") %>% 
                              filter(response == "rich" & rarity != "intermediate"))  

abmod <- glmmTMB::glmmTMB(val ~ treatment * rarity + (1|site)
                          , family = glmmTMB::nbinom2
                          , data = summ %>% 
                            pivot_longer(cols = c("ab", "rich"), names_to= "response", values_to = "val") %>% 
                            filter(response == "ab" & rarity != "intermediate"))  


emmeans::emmeans(richmod, pairwise~treatment|rarity)
# Richness effect is stronger for rare spp (b/c they have more 0s)

emmeans::emmeans(abmod, pairwise~treatment|rarity)
# abundance effect is the same (B/c it is)
# it has higher SD for rare spp, though (B/c they're rare!)











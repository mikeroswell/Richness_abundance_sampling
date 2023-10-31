# notes for ESA talk
# lose most common species: 55% of individuals and 20% of species, diversity increases

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(RTMB)
library(patchwork)

# rhetorical: lose 56% of individuals and 20% species; Shannon and Simpson increase
MeanRarity::rarity(c(5, rep(1, 4)), l = 0 )
MeanRarity::rarity(c(rep(1, 4)), l = 0 )

#Richness is a relative abundance measure too

# similarly, lose 67% of individuals and 20% of species and expected sample
# richness increases by 20%:
obs1 <- rmultinom(999, 15, c(5, 4, rep(1,3)))
obs2 <- rmultinom(999, 15, c(rep(1,4)))

# par(mfrow= c(2,1))
h <- mean(apply(obs1, 2, function(x){sum(x>1)}))
l <- mean(apply(obs2, 2, function(x){sum(x>1)}))
l/h

# load cloglog solver
source("scripts/cloglog_solver.R")

first_SAD <- cloglog_SAD(rich = 100, sd = 2)
second_SAD <- cloglog_SAD(rich = 100, sd = 1.5)
# plot sads
first_SAD_plot <- MeanRarity::radplot(first_SAD$ab) + ylim(c(0, 0.3))
second_SAD_plot <- MeanRarity::radplot(second_SAD$ab) + ylim(c(0, 0.3))
first_SAD_plot + second_SAD_plot


# load simulated data for cloglog explorations
source("scripts/cloglog_sim_only.R")

dd %>% ggplot(
  aes(x = sample_size, y = species, color = diversity) )+
    geom_violin() +
  theme_classic() +
  scale_y_log10()

dd %>% ggplot(
  aes(x = sample_size, y = rperc, color = diversity) )+
  geom_violin() +
  theme_classic() +
  scale_y_log10()



gn_pois <- glm(species ~ sample_size * diversity
               , family = "poisson"
               , data = dd)

summary(gn_pois)

exp(coef(gn_pois))

gn_nb <- MASS::glm.nb(species ~ sample_size * diversity
               # , family = "poisson"
               , data = dd)

summary(gn_nb)

exp(coef(gn_nb))

gn_logi <- glm(rperc ~ sample_size * diversityb 
                          , family = "binomial"
               , weights = rep(100, 4*reps)
                          , data = dd)

summary(gn_logi)

exp(coef(gn_logi))

gn_clog_flat <- glm(rperc ~ sample_size * diversity
               , family = binomial(link = "cloglog")
               , weights = rep(100, 4*reps)
               , data = dd)

summary(gn_clog_flat)

exp(coef(gn_clog_flat))

gn_clog_off <- glm(rperc ~ log(abundance) + diversity
                   , family = binomial(link = "cloglog")
                   , weights = rep(100, 4*reps)
                   , data = dd)

summary(gn_clog_off)

exp(coef(gn_clog_off))

gn_clog_off <- glm(rperc ~ log(abundance) + diversity * sample_size
                   , family = binomial(link = "cloglog")
                   , weights = rep(100, 4*reps)
                   , data = dd)

gn_clog_off <- glm(pres ~ offset(log(tot_ab)) + diversity *sample_size
                   , family = binomial(link = "cloglog")
                   , data = ddl)

exp(coef(gn_clog_off))

gn_clog_re <- lme4::glmer(pres ~ offset(log(tot_ab)) +  (1|diversity/sp)
                          , family = binomial(link = "cloglog")
                          , data = ddl
                          , start = list(fixef = -log(2))
                          )
exp(lme4::fixef(gn_clog_re))


gn_logi_re <- lme4::glmer(pres ~  diversity + rplct + 
                            sample_size + (1|sp)
                          , family = binomial(link = "logit")
                          , data = ddl)


gn_cloglog_re <- lme4::glmer(pres ~  diversity + sample_size + rplct +  (1|sp)
                          , family = binomial(link = "cloglog")
                          , data = ddl)


exp(lme4::fixef(gn_logi_re))
exp(lme4::fixef(gn_cloglog_re))

# gn_cloglog_re_off <- lme4::glmer(pres ~ +offset(log(tot_ab)) + diversity * sample_size + (1|sp)
#                              , family = binomial(link = cloglog)
#                              , data = ddl)

gn_cloglog_re_off <- lme4::glmer(pres ~ +offset(log(tot_ab)) +
                                   rplct + 
                                   (1|sp/diversity:sample_size)
                                 , family = binomial(link = cloglog)
                                 , data = ddl)



exp(lme4::fixef(gn_cloglog_re_off))


sub <- ddl %>%
  ungroup() %>%  
  filter(rplct %in% c("X1", "X2", "X3", "X4", "X5", "X6", "X7")
               , diversity == "l"
               , sample_size == "s") %>% 
  
  select(-c("diversity", "sample_size", "tot_r"))



dput(sub )

sub <- structure(list(sp = c("1", "1", "1", "1", "1", "2", "2", "2", 
                             "2", "2", "3", "3", "3", "3", "3", "4", "4", "4", "4", "4", "5", 
                             "5", "5", "5", "5"), rplct = c("X1", "X2", "X3", "X4", "X5", 
                                                            "X1", "X2", "X3", "X4", "X5", "X1", "X2", "X3", "X4", "X5", "X1", 
                                                            "X2", "X3", "X4", "X5", "X1", "X2", "X3", "X4", "X5"), abundance = c(131L, 
                                                                                                                                 124L, 128L, 114L, 108L, 62L, 57L, 62L, 65L, 54L, 45L, 45L, 35L, 
                                                                                                                                 38L, 33L, 29L, 31L, 30L, 32L, 30L, 21L, 20L, 25L, 23L, 28L), 
                      pres = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                               1, 1, 1, 1, 1, 1, 1, 1, 1), tot_ab = c(288L, 277L, 280L, 
                                                                      272L, 253L, 288L, 277L, 280L, 272L, 253L, 288L, 277L, 280L, 
                                                                      272L, 253L, 288L, 277L, 280L, 272L, 253L, 288L, 277L, 280L, 
                                                                      272L, 253L)), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, 
                                                                                                                                            -25L))

smaller <- lme4::glmer(pres ~ offset(log(tot_ab)) + (1|sp)
               , family = binomial(link = cloglog)
               , data = sub)

no_off <- lme4::glmer(pres ~ (1|sp)
                       , family = binomial(link = cloglog)
                      , data = sub)

no_re <- glm(pres ~ offset(log(tot_ab)) 
                             , family = binomial(link = cloglog)
                             , data = sub)

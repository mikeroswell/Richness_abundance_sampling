# notes for ESA talk
# lose most common species: 55% of individuals and 20% of species, diversity increases

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

MeanRarity::rarity(c(5, rep(1, 4)), l = 0 )
MeanRarity::rarity(c(rep(1, 4)), l = 0 )

#Richness is a relative abundance measure too

obs1 <- rmultinom(999, 15, c(5, rep(1,4)))
obs2 <- rmultinom(999, 15, c(rep(1,5)))

par(mfrow= c(2,1))
mean(apply(obs1, 2, function(x){sum(x>1)}))
mean(apply(obs2, 2, function(x){sum(x>1)}))

# Create a new SAD solver

cloglog <- function(x) log(-log(1-x))
gompertz <- function(x) 1-exp(-exp(x))

ur_cloglog <- function(x
  , rich = rich
  , sd = sd
  , ...){
    1 - sum(try_cloglog(x, sd, rich))
            }

rich <- 100
sd <- 2
int_uppr <- 0
int_lwr <- -1e5

try_cloglog <- function(x, sd, rich){
  myab <- gompertz(stats::qnorm(seq(from = (1 / rich) / 2
                            , to = 1 - (1 / rich) / 2
                            , by = (1 / rich)
  )
  , mean =x
  , sd = sd))
  myab[order(myab, decreasing = TRUE)]
}

cloglog_SAD <- function(rich, sd, int_lwr = -1e5, ...){
  ur =  tryCatch(
   stats::uniroot(function(x) {
      ur_cloglog(
        rich = rich,
        x = x,
        sd = sd)}
      , lower = int_lwr
      , upper = 0)
   )
    ab = try_cloglog(x = ur$root, sd = sd, rich = rich)
    return(
      list(distr = data.frame(form = "gompertz-normal"
                              , mean = ur$root
                              , sd = sd
                              , rich = rich)
                , ab = ab ))
    }

first_SAD <- cloglog_SAD(rich = 100, sd = 2)
second_SAD <- cloglog_SAD(rich = 100, sd = 1.5)
hist(first_SAD$ab)
hist(second_SAD$ab)
reps <- 999
n_s <- 100 # small sample
n_b <- 300 # big sample
obs_s_l <- rmultinom(size = n_s, n = reps, prob = first_SAD$ab)
obs_b_l <- rmultinom(size = n_b, n = reps, prob = first_SAD$ab)
obs_s_h <- rmultinom(size = n_s, n = reps, prob = second_SAD$ab)
obs_b_h <- rmultinom(size = n_b, n = reps, prob = second_SAD$ab)


dd <- map_dfr(c("s", "b"), function(sample_size){
    map_dfr(c("h", "l"), function(diversity){
      data.frame(abundance = apply(get(paste("obs"
                                             , sample_size, diversity
                                             , sep = "_"))
                                   , 2
                                   , sum)
             , species = apply(get(paste("obs"
                                         , sample_size, diversity
                                         , sep = "_"))
                               , 2
                               , function(x){sum(x>0)})
             , sample_size = sample_size
             , diversity = diversity) %>% 
        mutate(rperc = species/rich)
  })
})

ddl <- map_dfr(c("s", "b"), function(sample_size){
  map_dfr(c("h", "l"), function(diversity){
    
    data.frame(pivot_longer(data.frame(get(paste("obs"
                                                 , sample_size, diversity
                                                 , sep = "_"))) %>% tibble::rownames_to_column(var = "sp")
                            , cols = 2:(reps+1)
                            , names_to = "rplct"
                            , values_to = "abundance"
    ) %>% mutate(pres = as.numeric(abundance >0))
    , sample_size = sample_size
    , diversity = diversity
    ) %>% 
      group_by(rplct
               , sample_size
               , diversity) %>% 
      mutate(tot_ab = sum(abundance)
             , tot_r = sum(pres))
  })
})



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

gn_logi <- glm(rperc ~ sample_size * diversity
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

gn_clog_off <- glm(rperc ~ offset(log(abundance)) + diversity
                   , family = binomial(link = "cloglog")
                   , weights = rep(100, 4*reps)
                   , data = dd)

summary(gn_clog_off)

exp(coef(gn_clog_off))

gn_clog_off <- glm(rperc ~ offset(log(abundance)) + diversity * sample_size
                   , family = binomial(link = "cloglog")
                   , weights = rep(100, 4*reps)
                   , data = dd)

gn_clog_off <- glm(pres ~ offset(log(tot_ab)) + diversity *sample_size
                   , family = binomial(link = "cloglog")
                   , data = ddl)

exp(coef(gn_clog_off))

gn_clog_re <- lme4::glmer(pres ~ offset(log(tot_ab)) + diversity *sample_size + (1|sp)
                          , family = binomial(link = "cloglog")
                          , data = ddl)

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
  filter(rplct %in% c("X1", "X2", "X3", "X4")
               , diversity == "l"
               , sample_size == "s") %>% 
  
  select(-c("diversity", "sample_size"))

str(sub)

dput(sub)

sub <- structure(list(sp = c("1", "1", "1", "1", "2", "2", "2", "2", 
                      "3", "3", "3", "3", "4", "4", "4", "4", "5", "5", "5", "5", "6", 
                      "6", "6", "6", "7", "7", "7", "7", "8", "8", "8", "8", "9", "9", 
                      "9", "9", "10", "10", "10", "10", "11", "11", "11", "11", "12", 
                      "12", "12", "12", "13", "13", "13", "13", "14", "14", "14", "14", 
                      "15", "15", "15", "15", "16", "16", "16", "16", "17", "17", "17", 
                      "17", "18", "18", "18", "18", "19", "19", "19", "19", "20", "20", 
                      "20", "20", "21", "21", "21", "21", "22", "22", "22", "22", "23", 
                      "23", "23", "23", "24", "24", "24", "24", "25", "25", "25", "25", 
                      "26", "26", "26", "26", "27", "27", "27", "27", "28", "28", "28", 
                      "28", "29", "29", "29", "29", "30", "30", "30", "30", "31", "31", 
                      "31", "31", "32", "32", "32", "32", "33", "33", "33", "33", "34", 
                      "34", "34", "34", "35", "35", "35", "35", "36", "36", "36", "36", 
                      "37", "37", "37", "37", "38", "38", "38", "38", "39", "39", "39", 
                      "39", "40", "40", "40", "40", "41", "41", "41", "41", "42", "42", 
                      "42", "42", "43", "43", "43", "43", "44", "44", "44", "44", "45", 
                      "45", "45", "45", "46", "46", "46", "46", "47", "47", "47", "47", 
                      "48", "48", "48", "48", "49", "49", "49", "49", "50", "50", "50", 
                      "50", "51", "51", "51", "51", "52", "52", "52", "52", "53", "53", 
                      "53", "53", "54", "54", "54", "54", "55", "55", "55", "55", "56", 
                      "56", "56", "56", "57", "57", "57", "57", "58", "58", "58", "58", 
                      "59", "59", "59", "59", "60", "60", "60", "60", "61", "61", "61", 
                      "61", "62", "62", "62", "62", "63", "63", "63", "63", "64", "64", 
                      "64", "64", "65", "65", "65", "65", "66", "66", "66", "66", "67", 
                      "67", "67", "67", "68", "68", "68", "68", "69", "69", "69", "69", 
                      "70", "70", "70", "70", "71", "71", "71", "71", "72", "72", "72", 
                      "72", "73", "73", "73", "73", "74", "74", "74", "74", "75", "75", 
                      "75", "75", "76", "76", "76", "76", "77", "77", "77", "77", "78", 
                      "78", "78", "78", "79", "79", "79", "79", "80", "80", "80", "80", 
                      "81", "81", "81", "81", "82", "82", "82", "82", "83", "83", "83", 
                      "83", "84", "84", "84", "84", "85", "85", "85", "85", "86", "86", 
                      "86", "86", "87", "87", "87", "87", "88", "88", "88", "88", "89", 
                      "89", "89", "89", "90", "90", "90", "90", "91", "91", "91", "91", 
                      "92", "92", "92", "92", "93", "93", "93", "93", "94", "94", "94", 
                      "94", "95", "95", "95", "95", "96", "96", "96", "96", "97", "97", 
                      "97", "97", "98", "98", "98", "98", "99", "99", "99", "99", "100", 
                      "100", "100", "100"), rplct = c("X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", 
                                                      "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", 
                                                      "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", 
                                                      "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", 
                                                      "X3", "X4", "X1", "X2", "X3", "X4", "X1", "X2", "X3", "X4"), 
               abundance = c(131L, 124L, 128L, 114L, 62L, 57L, 62L, 65L, 
                             45L, 45L, 35L, 38L, 29L, 31L, 30L, 32L, 21L, 20L, 25L, 23L, 
                             21L, 21L, 25L, 18L, 18L, 17L, 22L, 23L, 20L, 16L, 13L, 15L, 
                             8L, 14L, 9L, 9L, 13L, 15L, 10L, 6L, 9L, 8L, 12L, 12L, 11L, 
                             14L, 11L, 12L, 10L, 4L, 7L, 14L, 4L, 4L, 9L, 6L, 9L, 5L, 
                             5L, 4L, 4L, 8L, 3L, 8L, 4L, 4L, 8L, 7L, 4L, 8L, 4L, 4L, 6L, 
                             2L, 6L, 2L, 3L, 8L, 5L, 4L, 6L, 4L, 5L, 3L, 6L, 5L, 2L, 5L, 
                             3L, 1L, 2L, 3L, 0L, 5L, 3L, 3L, 2L, 6L, 4L, 4L, 5L, 1L, 4L, 
                             6L, 2L, 3L, 3L, 4L, 2L, 3L, 5L, 2L, 2L, 3L, 3L, 2L, 3L, 2L, 
                             4L, 6L, 2L, 2L, 2L, 3L, 0L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 
                             3L, 1L, 4L, 3L, 0L, 2L, 2L, 0L, 2L, 0L, 1L, 1L, 1L, 0L, 1L, 
                             1L, 2L, 0L, 2L, 0L, 1L, 3L, 1L, 2L, 0L, 2L, 3L, 3L, 2L, 0L, 
                             1L, 1L, 2L, 2L, 1L, 1L, 3L, 0L, 0L, 1L, 1L, 2L, 0L, 0L, 0L, 
                             1L, 1L, 2L, 1L, 1L, 1L, 1L, 0L, 3L, 1L, 0L, 0L, 1L, 4L, 0L, 
                             1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 
                             3L, 1L, 1L, 0L, 2L, 3L, 0L, 1L, 0L, 1L, 2L, 1L, 0L, 2L, 1L, 
                             0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 2L, 0L, 0L, 0L, 0L, 1L, 1L, 
                             0L, 0L, 0L, 2L, 1L, 0L, 0L, 1L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 
                             0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 
                             0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 
                             1L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
                             0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
                             0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
                             0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                             0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
                             0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                             0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 
                             0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), pres = c(1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                                                                       1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 
                                                                                       1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 
                                                                                       1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 
                                                                                       1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 
                                                                                       1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 
                                                                                       1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 
                                                                                       1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 
                                                                                       0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 
                                                                                       1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 
                                                                                       0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                       0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                       0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                                       0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
               tot_ab = c(500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 500L, 
                          500L, 500L), tot_r = c(55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 
                                                 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 
                                                 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 
                                                 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 
                                                 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 
                                                 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 
                                                 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 
                                                 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 
                                                 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 61, 55, 58, 56, 
                                                 61)), row.names = c(NA, -400L), class = c("tbl_df", "tbl", 
                                                                                           "data.frame"))

smaller <- lme4::glmer(pres ~ offset(log(tot_ab)) + (1|sp)
               , family = binomial(link = cloglog)
               , data = sub)

no_off <- lme4::glmer(pres ~ (1|sp)
                       , family = binomial(link = cloglog)
                      , data = sub)

no_re <- glm(pres ~ offset(log(tot_ab)) 
                             , family = binomial(link = cloglog)
                             , data = sub)

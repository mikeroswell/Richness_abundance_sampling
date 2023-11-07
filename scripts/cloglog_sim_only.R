# simulate data for cloglog link-ed "GLM" with random variance by group


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


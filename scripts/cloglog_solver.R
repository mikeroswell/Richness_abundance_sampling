# Create a new SAD solver for gompertz-normal

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
  , mean = x
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
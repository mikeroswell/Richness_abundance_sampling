# putzing around with some cloglog stuff 

# from Bolker https://rpubs.com/bbolker/hazodds
## same as ff.cloglog$linkfun:
cloglog <- function(x) log(-log(1-x))
## same as ff.cloglog$linkinv
gompertz <- function(x) 1-exp(-exp(x))
all.equal(cloglog(gompertz(1.375)),1.375)

set.seed(101223)
ndev <- rnorm(500)
some_abundances <- gompertz(ndev) # so these can't be rel.abundances b/c sum is >1

# need a gompertz-normal that is a density. 

# is this translation of some sort?

ndev
gn<- gompertz(ndev)

gns <- gn/sum(gn)

hist(cloglog(gns))

bulk<-sum(gn)

trans1 <- gn +(1-bulk)     
sum(gompertz(trans1))

install.packages("flexsurv")
library(flexsurv)


gdev <- rgompertz(100)
hist(gdev)
hist(cloglog(gdev))
sum(gdev)

point2s<- rep(0.2, 5)
point1s <- rep(0.1, 5)
sum(cloglog(point2s))

sum(cloglog(point1s))

rdev <- rnorm(100)
sum(gompertz(rdev))
log(65.7)
sum(gompertz(rdev -5.0919))
hist(gompertz(rdev -5.0919))
    
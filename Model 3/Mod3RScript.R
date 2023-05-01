library(rstan)
library("gtools")

# Firstly, prepare the data for stan

# B is the reference population data (that we mix our samples from) 
B = matrix(
  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0.3, 0.4, 0.3),nrow=3,ncol=10,byrow=T)+0.05
B=B/rowSums(B)

# These are the old frequencies of the populations
f_old=c(0.7,0.2,0.1)
#------------ALL FOR GENERATING X---------------

set.seed(2) # create simulated values which are reproducible

# M is number of individuals in our sample X
M = 25
# Number of features (should be large, > 1000?)
#L=1000

## Grandparent model: # is L meant to be N
simpop<-function(nindivs,B,gens,f,sigma_fchange=1000){
    twototheG=2^gens
    f_true=rdirichlet(1,f*sigma_fchange)
    gpps=t(rmultinom(nindivs,twototheG,f_true))
    indB=gpps %*% B
    indB=indB/rowSums(indB)
    X=t(apply(indB,1,function(b)rmultinom(1,twototheG,b))) # N, specifying the total number of objects that are put into 
    # K boxes in the typical multinomial experiment
    ## part of list used by stan:
    list(K=dim(B)[1],
         M=nindivs,
         N=dim(B)[2],
         x=X,
         f_old=f,
         B=B,
         twototheG=twototheG,
         inverse_sigma=1,
         ## additional information we want
         sigma_fchange=sigma_fchange,
         indB=indB,
         gpps=gpps,
         f_true=f_true)
}

#create pop data for 2 gen passed and 4 gen passed
set.seed(1)
population_data2g = simpop(25,B,gens=2,f_old,sigma_fchange =1000)
population_data4g = simpop(25,B,gens=4,f_old,sigma_fchange=1000)
population_data7g = simpop(25,B,gens=7,f_old,sigma_fchange=1000)
## our model is ran here
fit2g <- stan(
  file = "Mod3Stan.stan",  # Stan program
  data = population_data2g,  # named list of data
  chains = 2,              # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
  cores = 1,               # number of cores (could use one per chain)
)
fit4g <- stan(
  file = "Mod3Stan.stan",  # Stan program
  data = population_data4g,  # named list of data
  chains = 2,              # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
  cores = 1,               # number of cores (could use one per chain)
)
fit7g <- stan(
  file = "Mod3Stan.stan",  # Stan program
  data = population_data7g,  # named list of data
  chains = 2,              # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
  cores = 1,               # number of cores (could use one per chain)
)

## We correlate the true B estimates with the observed ones
## and see the model performs well
## Extract all the parameters we will check

## Here is the output from the model
pars2g=as.data.frame(fit2g)
pars4g=as.data.frame(fit4g)
pars7g=as.data.frame(fit7g)

#col means of each param
mpars2g=colMeans(pars2g)
mpars4g=colMeans(pars4g)
mpars7g=colMeans(pars7g)

# estimates of indB, indF, sigma_fchange
indBest2g=matrix(mpars2g[grep("indB",names(mpars2g))],nrow=M)
indBest4g=matrix(mpars4g[grep("indB",names(mpars4g))],nrow=M)
indBest7g=matrix(mpars7g[grep("indB",names(mpars7g))],nrow=M)

indFest2g=matrix(mpars2g[grep("indf",names(mpars2g))],nrow=M)
indFest4g=matrix(mpars4g[grep("indf",names(mpars4g))],nrow=M)
indFest7g=matrix(mpars7g[grep("indf",names(mpars7g))],nrow=M)

sigmaest4g=mpars4g[grep("sigma_fchange",names(mpars4g))]
sigmaest2g=mpars2g[grep("sigma_fchange",names(mpars2g))]
sigmaest7g=mpars7g[grep("sigma_fchange",names(mpars7g))]


fest2g=mpars2g[grep("f_true",names(mpars2g))]
fest4g=mpars4g[grep("f_true",names(mpars4g))]
fest7g=mpars7g[grep("f_true",names(mpars7g))]

## Check the f_true estimates for 2g
data.frame(f_old=f_old,
           f_true=as.numeric(population_data2g$f_true),
           f_estimate_2g=fest2g,
           observed_f_2g=colMeans(indFest2g))


## Check the f_true estimates for 4g
data.frame(f_old=f_old,
           f_true=as.numeric(population_data4g$f_true),
           f_estimate_4g=fest4g,
           observed_f_4g=colMeans(indFest4g))

## Check the f_true estimates for 7g
data.frame(f_old=f_old,
           f_true=as.numeric(population_data7g$f_true),
           f_estimate_7g=fest7g,
           observed_f_7g=colMeans(indFest7g))

## Check the correlations
cor(as.numeric(indFest2g),as.numeric(population_data2g$gpps))
cor(as.numeric(indFest4g),as.numeric(population_data4g$gpps))
cor(as.numeric(indFest7g),as.numeric(population_data7g$gpps))

cor(as.numeric(indBest2g),as.numeric(population_data2g$indB))
cor(as.numeric(indBest4g),as.numeric(population_data4g$indB))
cor(as.numeric(indBest7g),as.numeric(population_data7g$indB))

plot(indBest2g,population_data2g$indB)
plot(indBest4g,population_data4g$indB)
plot(indBest7g,population_data7g$indB)

plot(colMeans(population_data2g$indB),colMeans(indBest2g))
plot(colMeans(population_data4g$indB),colMeans(indBest4g))

est_gpps2g = data.frame(f_estimate_1 = indFest2g[,1],
                    f_estimate_2 =indFest2g[,2],
                    f_estimate_3 =indFest2g[,3])
est_gpps7g = data.frame(f_estimate_1 = indFest7g[,1],
                        f_estimate_2 =indFest7g[,2],
                        f_estimate_3 =indFest7g[,3])

est_gpps4g = data.frame(f_estimate_1 = indFest4g[,1],
                        f_estimate_2 =indFest4g[,2],
                        f_estimate_3 =indFest4g[,3])


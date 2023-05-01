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
sigma_fchange=10
#sigma=1

# This function produces our X data (M individuals with N features each), 
# based on the ref pop profiles B
# and the true frequencies of the populations, sampled from multinomial
# distribution.
# rmultinorm simulates data from multinom dist 

# M is number of individuals in our sample X
M = 25
# Number of features (should this be 10 or 1000?)
N=1000

#--------------FUNCTIONS FOR GENERATING X DATA --------------

make_sample_mixtures<-function(f_true,M){
  ## Returns M individuals, each of which contains N features sampled from f_true
  xi=rdirichlet(M,f_true)
  xi
}
make_population_data<-function(A,B,N){ #different N to the donor number (N=1000 samples)
  AB=true_A %*% B # matrix multiplication
  X=sapply(1:nrow(AB),function(i){ ## i is the individual
    rmultinom(1,N,AB[i,]) 
  })
  t(X) # transpose of X
}

#-------------------------------------------------

#------------ALL FOR GENERATING X---------------

set.seed(2) # create simulated values which are reproducible

f_true=rdirichlet(1,f_old*sigma_fchange)[1,] ## cf rdirichlet(100,f_old*sigma)) and look at the column means
true_A=make_sample_mixtures(f_true,M=M) # dimension M by K
# simulate data from new generation with some change in f
X=make_population_data(true_A,B,N) # dimension M by N
# B is dimension K by N

#-----------------------------------------------
# SANITY CHECK
## A check: Do individuals with more from the first component have more from the corresponding feature in B?
plot(X[,1],true_A[,1]) # yes - pos correlation

## In expectation, it should be that X/N = A B
plot(X/N,true_A %*% B) # yes, linear x=y line plotted
print(X)

# Ultimately, we needed X for our model

library(rstan)

population_data = list( # need to give f_old and B
  K = 3,
  M = 25,
  N = 10,
  x = X,
  f_old = c(0.7,0.2,0.1),
  B = B
  
)

# our model is ran here
fit1 <- stan(
  file = "02-modelA.stan",  # Stan program
  data = population_data,  # named list of data
  chains = 4,              # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
  cores = 1,               # number of cores (could use one per chain)
)

# Here is the output from the model, we plot the recovered mixture probabilities
# and see the model performs well 
pars=as.data.frame(fit1) 
print.data.frame(pars)

library(Ternary)
TernaryPlot(alab = expression(bold('f_true'[1])), 
            blab = expression(bold('f_true'[2])),
            clab = expression(bold('f_true'[3])))
TernaryPoints(pars[c("f_true[1]", "f_true[2]","f_true[3]")] ,
              col='blue', pch=20, cex=0.5)
TernaryPoints(f_true, col='red', pch=4, cex=0.5)
TernaryPoints(f_old, col='black', pch=4, cex=0.5)
library(ggplot2)
stan_hist(fit1, pars = c("sigma_fchange")) +
  geom_vline(xintercept = 10, linetype='dotted', color ='blue') + 
    ggtitle("Posterior distribution of sigma_fchange")

# Firstly, prepare the data for stan

B = matrix(
  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0.3, 0.4, 0.3),nrow=3,ncol=10,byrow=T)+0.05
B=B/rowSums(B)

# This is the actual frequencies of the populations
ftrue=c(0.7,0.2,0.1)


# This function produces a sample child, based on the parent pop profiles 
# and the true frequencies of the populations, sampled from multinomial
# distribution.
make_population_data<-function(B,ftrue,N=100){
  xi=rmultinom(1,N,ftrue)[,1]
  x=rowSums(do.call("cbind",sapply(1:length(xi),function(i){
    rmultinom(xi[i],1,B[i,])
  })))
  list(K=dim(B)[1],N=dim(B)[2],B=B,x=x)
}

# Here we assign the output of the above function to the pop data we
# will be feeding into our stan model.
population_data=make_population_data(B,ftrue,N=1000)

# Next, we need to call stan function to draw posterior 
# samples

#my_file <- file.path("C:", "Users", "Joach", "Desktop", "my_file.csv")
#"C:\Users\Team Knowhow\Documents\YEAR 4\Project\STAN models\01-script.R"
# use forward slashes in file path otherwise leads to an error
library(rstan)
#library(Ternary)

# our model is ran here
fit1 <- stan(
  file = "Mod1Stan.stan",  # Stan program
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
TernaryPlot(alab = expression(bold('f'[1])), 
            blab = expression(bold('f'[2])),
            clab = expression(bold('f'[3])))
TernaryPoints(pars[c("f[1]", "f[2]","f[3]")] ,
              col='black', pch=20, cex=0.5)
TernaryPoints(ftrue, col='red', pch=4, cex=0.5)



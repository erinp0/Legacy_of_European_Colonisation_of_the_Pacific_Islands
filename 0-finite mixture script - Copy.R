library(dplyr)
library(ggplot2)

# Simulating some observations

N = 1000 # number of observations

mean = c(1,5,9)
sd = c(1,3,2)

lambda = c(0.2,0.5,0.3)

#sample 1000 from our categorical variable Z
z <- sample(1:3, size = N, prob = lambda, replace = T)

#add noise
epsilon = rnorm(N)

#simulate observations
y = mean[z] + sd[z]*epsilon

data_frame(y, z = as.factor(z)) %>% 
  ggplot(aes(x = y, fill = z)) +
  geom_density(alpha = 0.3) +
  ggtitle("Three data generating processes")


pop_data = list(N= N, y = y, n_groups = 3)
library(rstan)
fit1 <- stan(
  file = "0-modelstan.stan",  # Stan program
  data = pop_data,  # named list of data
  chains = 4,              # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
  cores = 1,               # number of cores (could use one per chain)
)

 
pars=as.data.frame(fit1)
print.data.frame(pars)
library(Ternary)
TernaryPlot(alab = expression(bold('lambda'[1])), 
            blab = expression(bold('lambda'[2])),
            clab = expression(bold('lambda'[3])))
TernaryPoints(pars[c("lambda[1]", "lambda[2]","lambda[3]")] ,
              col='black', pch=20, cex=0.5)
TernaryPoints(lambda, col='red', pch=4, cex=0.5)

#plot(fit1) + 
 # stat_summary(fun.y='mean', geom="point", shape=18,
  #             size=3, color="red")
#ggtitle("")

stan_hist(fit1, pars = c("lambda")) +
  geom_vline(xintercept = 10, linetype='dotted', color ='blue') +
  xlim(0,1)+
  ggtitle("Posterior distribution of lambda")


stan_hist(fit1, pars = c("sd[1]", "sd[2]","sd[3]")) +
  #xlim(-2,15)
  ggtitle("Posterior distributions of the standard deviations")

stan_hist(fit1, pars = c("mean[1]", "mean[2]","mean[3]"),binwidth=0.25) +
  xlim(-2,12) +
  ggtitle("Posterior distributions of the mean")


#saveRDS(fit1, file = "C:/Users/Team Knowhow/Documents/YEAR 4/Project/STAN models/Model 0")

#model_old = readRDS("C:/Users/Team Knowhow/Documents/YEAR 4/Project/STAN models/Model 0")
#ls()

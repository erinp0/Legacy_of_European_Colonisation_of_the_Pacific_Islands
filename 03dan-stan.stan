//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
// This Stan program defines a model, for estimating the
// change in frequencies of populations from
// generation to generation. E.g. a population may become
// more or less prevalent as time goes on. 

// k is number of reference populations, N is number of features. f_old holds
// frequencies of the ref pops in the past generation. (f_new holds the frequencies
// in the current generation.) There are M individuals in the old generation.

data {
  int<lower = 1> K; 
  int<lower = 1> N;
  int<lower = 1> M;
  matrix[K,N] B;
  int<lower = 0> x[M,N];
  real<lower=0> twototheG;
  simplex[K] f_old;
  real<lower=0> inverse_sigma; // 1/mean sigma_fchange value
}

// The parameters accepted by the model. Our model
// accepts one parameters 'f_new' which is a vector representing
// the mixture probabilities of each surrogate pop. in the current generation
parameters {
  simplex[K] f_true;  // the true f that the samples now come from
  real<lower=0> sigma_fchange; // how much f changed by
  simplex[K] indf[M]; // the individual f's that the samples have
}
transformed parameters {
  vector<lower = 0>[N] indB[M];
  for(m in 1:M){ // indB = F %*% B
    for(n in 1:N){
      indB[m][n]=0;
      for(k in 1:K){
        indB[m][n] += indf[m][k]*B[k,n];
      }
    }
  }
}

model {
  // priors
  sigma_fchange~exponential(inverse_sigma); // argument is rate parameter
  f_true~dirichlet(f_old*sigma_fchange); //change in frequency
  for(m in 1:M){
    indf[m]~dirichlet(f_true*twototheG); // dirichlet approximation to the gpp distribution (which is multinomial)
  }
  // likelihood    
  for(m in 1:M){
    x[m,] ~ multinomial(indB[m]); // multinomial likelihood
  }
}


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
  simplex[K] f_old;
  matrix[K,N] B;
  int<lower = 0> x[M,N]; 
}

// The parameters accepted by the model. Our model
// accepts one parameters 'f_new' which is a vector representing
// the mixture probabilities of each surrogate pop. in the current generation
parameters {
  simplex[K] f_true; 
  //real<lower = 0>sigma; 
  real<lower = 0>sigma_fchange;
}

transformed parameters {
  vector[N] alpha[M];

  for(m in 1:M){
    for(i in 1:N){
      alpha[m][i]=0;
      for(j in 1:K){
        alpha[m][i] = alpha[m][i] + f_true[j]*B[j,i];

      }
    }
  }
}

// The model to be estimated. We model the output
// 'f_new' (new mixture frequencies) to be dirichlet distributed 
// with parameter (1,1,1), with length p (=3 in this case).
// x is dirichlet distributed with parameter alpha
model {
  // priors
  //sigma~exponential(1);
  
  f_true~dirichlet(f_old*sigma_fchange); 
  
  for(m in 1:M){
    x[m,]~multinomial(alpha[m]); 
  }
}


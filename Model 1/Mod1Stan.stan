//
// int<lower=0> N means N can only take values 0,1,2...


//This Stan program defines a simple model, for recovering
// the mixture probabilities (f_i) from 1 generation of
// ancestry.

// The input data is a matrix 'B' with P rows and k 
// columns, where p is the number of reference populations
// and k is the number of features. In this simple example
// we will consider 3 reference populations, each with 10
// features. Ref pop. are dirichlet distributed with 
// parameter alpha (1,1,...,1) (length N=10).

// X is the observed data (children), in this case we only
// have 1 child and this child has 10 features. Obs data is 
// multinomial distributed with parameter (sum of f_k.B_ki, from k=1 to 3)
// and is vector length N.

data {
  int<lower = 1> K;
  int<lower = 1> N;
  matrix[K,N] B;
  int<lower = 0> x[N]; //changed from vector to simplex
}

// The parameters accepted by the model. Our model
// accepts one parameters 'f' which is a simplex representing
// the mixture probabilities of each ref pop.
parameters {
  simplex[K] f;
}

transformed parameters {
  vector[N] alpha;
  for(i in 1:N){
    alpha[i]=0;
    for(j in 1:K){
      alpha[i] = alpha[i] + f[j]*B[j,i];
    }
  }
}

// The model to be estimated. We model the output
// 'f' (mixture probability) to be dirichlet distributed 
// with parameter (1,1,1), with length p (=3 in this case).
// x is dirichlet distributed with parameter alpha
model {
  f~dirichlet(rep_vector(1, K));
  x~multinomial(alpha);
}


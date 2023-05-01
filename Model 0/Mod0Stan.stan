// saved as 0-modelstan.stan
data {
  int N;
  vector[N] y;
  int n_groups; // number of populations
}
parameters {
  ordered[n_groups] mean; // ordered prevents label switching whilst model is running
  vector<lower = 0>[n_groups] sd;
  simplex[n_groups] lambda;
}
model {
  vector[n_groups] contributions;
  // priors
  mean ~ normal(0, 10);
  sd ~ cauchy(0, 2);
  lambda ~ dirichlet(rep_vector(2.0, n_groups));
  
  
  // likelihood
  for(i in 1:N) {
    for(k in 1:n_groups) {
      contributions[k] = log(lambda[k]) + normal_lpdf(y[i] | mean[k], sd[k]);
    }
    target += log_sum_exp(contributions);
  }
}

data {
  int<lower=0> N;
  int n[N]; //larvae count
  vector[N] r; //radial distance
  vector[N] theta; //angle
  vector[N] leaf; //leaf count
  int indices[N]; //treatment group indices
  //int year[N];
  int<lower=0> K; //treatments * years
  vector[K] rel_density;
  matrix<lower = -50, upper = 50>[K,2] drift;
  real c_order;
  real standard_N; //minimum number of larvae released (maximum density possible);
  int<lower=0> M; //treatments
  int tmt[K];
}

parameters {
  //real<lower = -log(standard_N), upper = log(standard_N)> log_a_mean;
  //real<lower = 0, upper = 10> log_a_sigma;
  vector<lower = -log(standard_N), upper = log(standard_N)>[M] log_a_intercept;
  real<upper = log(100)> log_b_coef;
  real<upper = log(10)> log_resid_sigma;
}



transformed parameters {
  vector[N] expect;
  vector[N] shifted_r;
  real b_coef = exp(log_b_coef);
  vector[M] a_intercept = exp(log_a_intercept);
  vector[N] residual;
  real<lower = 0, upper = 10> resid_sigma = exp(log_resid_sigma);
  
  for (i in 1:N){
    shifted_r[i] = sqrt((r[i]*sin(theta[i]) - drift[indices[i], 2])^2 + (r[i]*cos(theta[i]) - drift[indices[i], 1])^2);
    expect[i] = rel_density[indices[i]] * a_intercept[tmt[indices[i]]] * exp(-b_coef * shifted_r[i]^c_order) / (2 * pi());
    residual[i] = (n[i] / leaf[i]) - expect[i];
  }
}

model {
  //target += uniform_lpdf(exp(log_a_mean) | 0, standard_N);
  //log_a_sigma ~ uniform(0,10);
  //log_a_intercept ~ normal(log_a_mean, log_a_sigma);
  target += uniform_lpdf(exp(log_a_intercept) | 0, standard_N);
  
  log_b_coef ~ normal(0, 4);
  target += uniform_lpdf(resid_sigma | 0,10);
  //for (i in 1:N) target += poisson_lpmf(n[i] | expect[i] * leaf[i]);
  for (i in 1:N) target += normal_lpdf(residual[i] | 0, resid_sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = normal_lpdf(residual[i] | 0, resid_sigma);
}

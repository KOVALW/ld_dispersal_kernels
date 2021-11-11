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
  
  int tmt[K];
  int<lower=0> qs; //1 + number of qs (either 1 +1 or 2 +1)
  
  vector[qs+2] q_sds;
  vector[qs] q_sca;
  vector[qs] q_ctr;
  matrix[qs, qs] q_rots; //rotations from PCA
}


parameters {
  vector[qs] q_pcp;
  real<upper= log(10)> log_resid_sigma;
}



transformed parameters {
  
  real b_coef;
  vector[qs-1] a_intercept;
  real<lower = 0, upper = 10> resid_sigma = exp(log_resid_sigma);
  
  vector[N] expect;
  vector[N] shifted_r;
  vector[N] residual;
  
  vector[qs] used_qs;
  for (j in 1:qs) {
    row_vector[qs] curr = q_rots[j];
    used_qs[j] = q_ctr[j] + dot_product(curr, q_pcp) * q_sca[j]; //this transfroms from pca to useable parms
  }
  
  b_coef = exp(used_qs[1]);
  for (i in 2:qs) a_intercept[i-1] = exp(used_qs[i]);

  for (i in 1:N){
    shifted_r[i] = sqrt((r[i]*sin(theta[i]) - drift[indices[i], 2])^2 + (r[i]*cos(theta[i]) - drift[indices[i], 1])^2);
    expect[i] = rel_density[indices[i]] * a_intercept[tmt[indices[i]]] * exp(-b_coef * shifted_r[i]^c_order) / (2 * pi());
    residual[i] = (n[i] / leaf[i]) - expect[i];
  }
}

model {
  //target += uniform_lpdf(a_intercept | 0, 1e6);
  
  //target += normal_lpdf(log(b_coef) | 0, 4);
  //target += uniform_lpdf(resid_sigma | 0,10);
  for (i in 1:qs){
  q_pcp[i] ~ normal(0, q_sds[i]);
  }
  
  log_resid_sigma ~ normal(q_sds[qs+1], q_sds[qs+2]);
  
  for (i in 1:N) target += normal_lpdf(residual[i] | 0, resid_sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = normal_lpdf(residual[i] | 0, resid_sigma);
}

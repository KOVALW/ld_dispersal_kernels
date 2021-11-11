
functions {
  real distort(real x, real y, real beta, real gamma){
    real term1 = (2*gamma^2 + 1)*((x/beta)^2 + 2*gamma*x/beta);
    real term2 = (gamma^2 + 1)*(y^2 + 2*gamma^2);
    real term3 = 2*gamma*(x/beta + gamma);
    real rootterm = (gamma^2 + 1)*((x/beta + gamma)^2 + y^2 + 1);
    return(sqrt(term1+term2-term3*sqrt(rootterm)));
  }
  
  real distort_distance(real x, real y, real beta, real gamma, real psi){
    real x_adj = x*cos(psi) + y*sin(psi);
    real y_adj = x*sin(psi) - y*cos(psi);
    return(distort(x_adj, y_adj, beta, gamma));
  }
}
data {
  int<lower=0> N; // number data observations
  vector[N] r; //radial distance
  int n[N]; //observed larvae
  vector[N] leaf; //leaf count
  vector[N] theta; // angle
  int indices[N]; //treatment group indices 
  
  int<lower=0> K; //number of treatment groups * timepoints
//  int tmt[K];
  int qpoint[K];
  int windpoint[K];
  vector[K] rel_density;
  real c_order;
  
  int<lower=0> winds; // number of wind prior sets ( either 1 or 2)
  int<lower=0> qs; //1 + number of qs (either 1 +1 or 2 +1)
  
  matrix[winds, 3] wind_sds;
  matrix[winds, 3] wind_sca;
  matrix[winds, 3] wind_ctr;
  vector[qs+2] q_sds;
  vector[qs] q_sca;
  vector[qs] q_ctr;
  
  matrix[3, 3] wind_rots[winds]; //rotations from PCA
  matrix[qs, qs] q_rots; //rotations from PCA

} 

transformed data {
  vector[N] x;
  vector[N] y;
  for (i in 1:N){
    x[i] = r[i]*cos(theta[i]);
    y[i] = r[i]*sin(theta[i]);
  }
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
//  vector[K] releasedE5; //number released in each treatment group, in hundreds of thousands
  matrix[winds, 3] wind_pcp;
  vector[qs] q_pcp;
  real<upper= log(10)> log_resid_sigma;
  
//  matrix[K,2] shift; //shifted center of dispersal
}

transformed parameters {
  vector[qs] used_qs;
  matrix[winds, 3] used_winds;
  
  
  vector[winds] logbeta;
  vector[winds] loggamma;
  vector[winds] beta;
  vector[winds] gamma;
  vector[winds] psi;
  
  
  real b_coef;
  vector[qs-1] a_intercept;
  real<lower = 0, upper = 10> resid_sigma;
  
  vector[N] expect;
  vector[N] shifted_r;
  vector[N] residual;
  
  
  for (j in 1:qs) {
    row_vector[qs] curr = q_rots[j];
    used_qs[j] = q_ctr[j] + dot_product(curr, q_pcp) * q_sca[j]; //this transfroms from pca to useable parms
  }
  
  for (w in 1:winds) {
    for (j in 1:3) {
    row_vector[3] curr = wind_rots[w, j];
    used_winds[w,j] = wind_ctr[w, j] + dot_product(curr, wind_pcp[w]) * wind_sca[w, j];
  }
  
  resid_sigma = exp(log_resid_sigma);
  
  logbeta[w] = used_winds[w,1];
  loggamma[w] = used_winds[w,2];
  psi[w] = used_winds[w,3];
  }
  
  gamma = exp(loggamma);
  beta = exp(logbeta);

  b_coef = exp(used_qs[1]);
  for (i in 2:qs) a_intercept[i-1] = exp(used_qs[i]);
  
  for (i in 1:N){
    shifted_r[i] = distort_distance(x[i], y[i], 
      beta[windpoint[indices[i]]], 
      gamma[windpoint[indices[i]]], 
      psi[windpoint[indices[i]]]);
    expect[i] = rel_density[indices[i]] * a_intercept[qpoint[indices[i]]] * exp(-b_coef * shifted_r[i]^c_order)/(2*pi()*beta[windpoint[indices[i]]]*sqrt(gamma[windpoint[indices[i]]]^2+1));
    residual[i] = (n[i] / leaf[i]) - expect[i];
//print("r ", r[i], ", shifted ", shifted_r[i], ", expect ", expect[i], ", resid ", residual[i]);
  }
}

model {
  
  for (j in 1:winds) {
    wind_pcp[j] ~ normal(0, wind_sds[j]);
  }
  
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

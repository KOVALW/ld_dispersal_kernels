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
  int<lower=0> N;
  vector[N] u; //wind velocity in x direction
  vector[N] v; //wind velocity in y direction
  int release[N];
  
  int<lower=0> P; //number of parameters from PCA
  int<lower=0> winds; // number of wind prior sets ( either 1 or 2)

  vector[P] wind_sca;
  vector[P] wind_ctr;
  
  matrix[P, P] wind_rots; //rotations from PCA
  vector[2] mag_penalty;
  //vector[P] sds; //std devs from PCA
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  matrix[winds, P] wind_pcp;
}

transformed parameters {
  vector[N] dist;
  matrix[winds, P] used_winds;
//  real gamma_adj = gamma;
//  vector[P] used_par;
  vector[winds] logbeta;
  vector[winds] loggamma;
  vector[winds] beta;
  vector[winds] gamma;
  vector<lower = -pi()*1/6, upper = pi()*2/3>[winds] psi;
  vector[winds] penalty_term;
  
  for (w in 1:winds) {
    for (j in 1:P) {
    row_vector[P] curr = wind_rots[j];
    used_winds[w,j] = wind_ctr[j] + dot_product(curr, wind_pcp[w]) * wind_sca[j];
  }
  
  logbeta[w] = used_winds[w,1];
  loggamma[w] = used_winds[w,2];
  psi[w] = used_winds[w,3];
  }
  
  beta = exp(logbeta);
  gamma = exp(loggamma);
  
  for (i in 1:N){ 
    dist[i] = log(distort_distance(u[i], v[i], beta[release[i]], gamma[release[i]], psi[release[i]]));
}

  for (i in 1:winds){
  penalty_term[i] = sqrt(logbeta[i]^2 + loggamma[i]^2);
  }
}

model {
  for (i in 1:N){
    target += normal_lpdf(dist[i] | 0, 1);
  }

//  target += uniform_lpdf(logbeta | -5, 5);
//  target += uniform_lpdf(loggamma | -5, 5);
  for (i in 1:winds) target += von_mises_lpdf(psi[i] | pi()/4, 4);
  for (i in 1:winds) target += gamma_lpdf(penalty_term[i] | mag_penalty[1], mag_penalty[2]);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] =  normal_lpdf(dist[i] | 0, 1);
}

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
  
  int<lower=0> P; //number of parameters from PCA
  matrix[P, P] rot; //rotations from PCA
  vector[P] sca; //scales from PCA
  vector[P] ctr; //centers from PCA
  vector[2] mag_penalty;
  //vector[P] sds; //std devs from PCA
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[P] pc_par;
}

transformed parameters {
  vector[N] dist;
//  real gamma_adj = gamma;
  vector[P] used_par;
    real beta;
    real gamma;
    real logbeta;
    real loggamma;
    real penalty_term;
    real<lower = -pi()*1/6, upper = pi()*2/3> psi;
  
    for (j in 1:P) {
    row_vector[P] curr = rot[j];
    used_par[j] = ctr[j] + dot_product(curr, pc_par) * sca[j]; //this transfroms from pca to useable parms
  }
  
  logbeta = used_par[1];
  loggamma = used_par[2];
  psi = used_par[3];
  
  beta = exp(logbeta);
  gamma = exp(loggamma);
  
  for (i in 1:N){ 
    dist[i] = log(distort_distance(u[i], v[i], beta, gamma, psi));
}
  penalty_term = sqrt(logbeta^2 + loggamma^2);
}

model {
  for (i in 1:N){
    target += normal_lpdf(dist[i] | 0, 1);
  }

//  target += uniform_lpdf(logbeta | -5, 5);
//  target += uniform_lpdf(loggamma | -5, 5);
  target += von_mises_lpdf(psi | pi()/4, 4);
  target += gamma_lpdf(penalty_term | mag_penalty[1], mag_penalty[2]);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] =  normal_lpdf(dist[i] | 0, 1);
}

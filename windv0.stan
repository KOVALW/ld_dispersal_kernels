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
  //int<lower=0> K; //number of treatment groups
  //int year[N];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real logbeta; //elliptical contour coherency
  real loggamma; //elliptical drift
  real<lower = -pi()*1/6, upper = pi()*2/3> psi; //rotation angle
}

transformed parameters {
  vector[N] dist;
//  real gamma_adj = gamma;
  real beta = exp(logbeta);
  real gamma = exp(loggamma);
  for (i in 1:N){ 
    dist[i] = log(distort_distance(u[i], v[i], beta, gamma, psi));
}
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (i in 1:N){
    target += normal_lpdf(dist[i] | 0, 1);
  }
  logbeta ~ uniform(-5, 5);
  loggamma ~ uniform(-5, 5);

  psi ~ von_mises(pi()/4,4);

}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] =  normal_lpdf(dist[i] | 0, 1);
}

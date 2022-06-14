
data{
  int<lower=0> N;
  int<lower=0> deaths[N];
  real mob[N];
  int T;
  int M;
  vector[M] mobility_counterfactual;
}

transformed data{
  //pre-calculate values of the gamma pdf
  int deaths0 = 0; // assume 0 deaths on day 0
  real omega[N];
  real h[N];
  real omega0; 
  real h0;
  for (n in 1:N){
    // convert mean and sd to alpha and beta
    h[n] = exp(gamma_lpdf(n | (18.8)^2 / (8.46)^2, 18.8 / (8.46)^2)); 
    omega[n] = gamma_cdf(n + 0.5, (6.48)^2 / (3.83)^2, 6.48 / (3.83)^2) - gamma_cdf(n - 0.5, (6.48)^2 / (3.83)^2, 6.48 / (3.83)^2); 
  }
  h0 = exp(gamma_lpdf(0 | (18.8)^2 / (8.46)^2, 18.8 / (8.46)^2)); 
  omega0 = gamma_cdf(0.5, (6.48)^2 / (3.83)^2, 6.48 / (3.83)^2) - gamma_cdf(0, (6.48)^2 / (3.83)^2, 6.48 / (3.83)^2);
}

parameters{
  real R01;
  real R02;
  real beta1;
  real beta2;
  real<lower=0> delta;
}

transformed parameters{
  real R[N];
  real RD[N];

  for (n in 1:N){
    if (n <= T){
      R[n] = exp(log(R01) - beta1 * (1-mob[n]));
    }
    else {
      R[n] = exp(log(R02) - beta2 * (1-mob[n]));
    }
  }
  for (n in 1:N){
    RD[n] = R01 * h[n]; // s = 0
    for (s in 1:(n-1)){ // s = 1, ..., n-1
      RD[n] += R[s] * h[n-s];
    }
    RD[n] += R[n] * h0; // s = n
  }
}

model{
  real sumD_omega;
  
  // Priors: 
  R01 ~ uniform(0,5);
  R02 ~ uniform(0,5);
  beta1 ~ uniform(-100,100);
  beta2 ~ uniform(-100,100);
  delta ~ exponential(1);
  
  // Model:
  deaths[1] ~ neg_binomial_2(RD[1] * (deaths0 * omega[1] + deaths[1] * omega0), delta);
  for (n in 2:N){
    sumD_omega = deaths0 * omega[n]; // s = 0
    for (s in 1:(n-1)){ // s = 1, ..., n-1
      sumD_omega += deaths[s] * omega[n-s];
    }
    sumD_omega += deaths[n] * omega0; // s = n
    deaths[n] ~ neg_binomial_2(RD[n] * sumD_omega, delta);
  }
}

generated quantities {
  vector[M] R1_counterfactual;
  vector[M] R2_counterfactual;
  real R1_threshold;
  real R2_threshold;
  
  R1_counterfactual = exp(log(R01) - beta1 * (1-mobility_counterfactual));
  R2_counterfactual = exp(log(R02) - beta2 * (1-mobility_counterfactual));
  
  R1_threshold = log(R01) / beta1;
  R2_threshold = log(R02) / beta2;
}

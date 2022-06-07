
data{
  int<lower=0> N;
  int<lower=0> deaths[N];
  real mob[N];
}

transformed data{
  //pre-calculate values of the gamma pdf
  real omega[N];
  real h[N];
  real omega0;
  for (n in 1:N){
    h[n] = exp(gamma_lpdf(n | 18.8, 8.46)); //this parameterisation is wrong!
    omega[n] = exp(gamma_lpdf(n | 6.48, 3.83)); //this parameterisation is wrong!
  }
  omega0 = exp(gamma_lpdf(0 | 6.48, 3.83));
}

parameters{
  real R01;
  real R02;
  real beta1;
  real beta2;
  real<lower=0> delta;
  real<lower=2, upper=N-1> T;
}

transformed parameters{
  real R[N];
  real RD[N];

  for (n in 1:N){
    if (n < T){
      R[n] = exp(log(R01) - beta1 * (1-mob[n]));
    }
    else {
      R[n] = exp(log(R02) - beta2 * (1-mob[n]));
    }
  }
  for (n in 1:N){
    RD[n] = R01 * h[n];
    for (s in 1:(n-1)){
      RD[n] = RD[n] + R[s] * h[n-s];
    }
  }
}

model{
  real sumD_omega;
  
  //Priors: 
  R01 ~ uniform(0,5);
  R02 ~ uniform(0,5);
  beta1 ~ uniform(-100,100);
  beta2 ~ uniform(-100,100);
  delta ~ exponential(1);
  
  //Model:
  //need to include deaths[1]!
  for (n in 2:N){
    //assume D[0] = 0
    sumD_omega = deaths[n] * omega0;
    for (s in 1:(n-1)){
      sumD_omega += deaths[s] * omega[n-s];
    }
    deaths[n] ~ neg_binomial_2(RD[n] * sumD_omega, delta);
  }
}

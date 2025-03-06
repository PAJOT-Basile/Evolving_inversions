data {
  int<lower=0> N;     // number of data points
  vector[N] d_x;      // distance along transect
  vector[N] gen;     // genotype
  real upc;           // upper limit for centre
  real loc;           // lower limit for centre
  real uplw;          // upper limit for width
}

parameters {
  real<lower=loc, upper=upc> centre;            // centre, must be within transect
  real<upper=uplw> log_width;                   // log_width 
  real<lower=-10,upper=10> left;                // left AF (logit)
  real<lower=-10,upper=10> right;               // right AF (logit)
  real<lower=-1, upper=1> f;                    // Fis at cline centre
}

transformed parameters {
  real w = exp(log_width);                               // width back on natural scale
  real lp = exp(left)/(1+exp(left));                     // end frequencies as proportions
  real rp = exp(right)/(1+exp(right));
  vector[N] p_x;
  vector[N] f_x;
  vector[N] z_x;
  vector[N] z1;
  vector[N] z2;
  vector[N] z3;
  vector[N] L;
  for (i in 1:N) {
    p_x[i] = 1/(1 + exp(0 - 4*((d_x[i] - centre)/w)));   // underlying frequency cline, proportion wave
    f_x[i] = f - f*2*fabs(p_x[i] - 0.5);                  // linear decline in Fis from f at centre to 0 at each end
    if(p_x[i]/(p_x[i]-1) > (p_x[i]-1)/p_x[i]){
      if(f_x[i] < p_x[i]/(p_x[i]-1)){f_x[i] = p_x[i]/(p_x[i]-1);}}
    if(p_x[i]/(p_x[i]-1) < (p_x[i]-1)/p_x[i]){
      if(f_x[i] < (p_x[i]-1)/p_x[i]){f_x[i] = (p_x[i]-1)/p_x[i];}}  // these conditionals keep p_x in 0,1
    z_x[i] = lp + ((rp - lp)*p_x[i]);              // AF cline
    z1[i] = z_x[i]^2 + f_x[i]*z_x[i]*(1-z_x[i]);
    z2[i] = 2*z_x[i]*(1-z_x[i])*(1 - f_x[i]);
    z3[i] = (1-z_x[i])^2 + f_x[i]*z_x[i]*(1-z_x[i]);
    if(gen[i] == 0){L[i] = z3[i];}         // likelihood for each genotype
    if(gen[i] == 1){L[i] = z2[i];}
    if(gen[i] == 2){L[i] = z1[i];}
    }
}

model {
  target += sum(log(L));  // log-likelihood of the data - Stan will maximise this quantity
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = L[n];
  }
}




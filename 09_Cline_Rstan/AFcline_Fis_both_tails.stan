data {
  int<lower=0> N;     // number of data points
  vector[N] d_x;      // distance along transect
  vector[N] gen;     // genotype
  real upc;
  real loc;
  real uplw;
  real upd;         // max distance of tails from cline centre
}

parameters {
  real<lower=loc, upper=upc> centre;            // centre, must be within transect
  real<upper=uplw> log_width;                   // log_width 
  real<lower=-1, upper=1> f;                    // Fis at cline centre
  real<lower=0,upper=upd> dl;                      // distance of left tail from cline centre
  real<lower=0,upper=1> tl;                         // ratio of left tail slope to centre slope
  real<lower=0,upper=upd> dr;                      // distance of right tail from cline centre
  real<lower=0,upper=1> tr;                         // ratio of right tail slope to centre slope
}

transformed parameters {
  real w = exp(log_width);                               // width back on natural scale
  real zl = 4*dl/w;
  real Al = 1/(1+exp(zl));
  real kl = tl*(1-Al);
  real zr = 4*dr/w;
  real Ar = 1/(1+exp(zr));
  real kr = tr*(1-Ar);
  real dp = (1/(1+exp(0-4*((centre+dr)/w))))-(1/(1+exp(0-4*((centre-dl)/w))));
  real pl = Al*exp(kl*(zl+4*((centre-dl)/w))) - Al*exp(kl*(zl+4*(((centre-dl)-1)/w)));
  real pr = (1-Ar*exp(kr*(zr-4*(((centre+dr)+1)/w)))) - (1-Ar*exp(kr*(zr-4*((centre+dr)/w))));
  real Bl = dp/pl;
  real Br = dp/pr;
  vector[N] p_x;
  vector[N] u_x;
  vector[N] f_x;
  vector[N] z_x;
  vector[N] z1;
  vector[N] z2;
  vector[N] z3;
  vector[N] L;
  for (i in 1:N) {
    u_x[i] = 0 - 4*((d_x[i] - centre)/w);
    p_x[i] = 1/(1 + exp(u_x[i]));   // underlying frequency cline, proportion wave, centre
    if(d_x[i] < (centre-dl)){p_x[i] = Al*exp(kl*(zl-u_x[i]));}
    if(d_x[i] > (centre+dr)){p_x[i] = 1-Ar*exp(kr*(zr+u_x[i]));}
    f_x[i] = f - f*2*fabs(p_x[i] - 0.5);                  // linear decline in Fis from f at centre to 0 at each end
    if(p_x[i]/(p_x[i]-1) > (p_x[i]-1)/p_x[i]){
      if(f_x[i] < p_x[i]/(p_x[i]-1)){f_x[i] = p_x[i]/(p_x[i]-1);}}
    if(p_x[i]/(p_x[i]-1) < (p_x[i]-1)/p_x[i]){
      if(f_x[i] < (p_x[i]-1)/p_x[i]){f_x[i] = (p_x[i]-1)/p_x[i];}}  // these conditionals keep p_x in 0,1
    z_x[i] = p_x[i];              // AF cline, no end frequencies when adding tails
    z1[i] = z_x[i]^2 + f_x[i]*z_x[i]*(1-z_x[i]);
    z2[i] = 2*z_x[i]*(1-z_x[i])*(1 - f_x[i]);
    z3[i] = (1-z_x[i])^2 + f_x[i]*z_x[i]*(1-z_x[i]);
    if(gen[i] == 0){L[i] = z3[i];}         // likelihood for each genotype
    if(gen[i] == 1){L[i] = z2[i];}
    if(gen[i] == 2){L[i] = z1[i];}
    }
}

model {
  target += sum(log(L));
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = L[n];
  }
}




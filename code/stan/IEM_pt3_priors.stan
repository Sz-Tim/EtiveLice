// Project: E(tive)Lice
// Tim Szewczyk
// tim.szewczyk@sams.ac.uk
// Integrated ensemble model


data {
  int<lower=0> nDays;
  int<lower=0> nFarms;
  int<lower=0> nSims;
  int<lower=0> nStages;
  int<lower=1> nAttachCov;
  int<lower=1> nSurvCov;
  int<lower=1> nSamples;
  array[nFarms] int<lower=0> nPens;
  array[nStages, 2, nDays, nFarms] int y;
  matrix[nDays, nFarms] N_attach;
  array[nStages, 2] matrix[nDays, nFarms] mu_true;
  array[nFarms] matrix<lower=0>[nDays, nSims] IP_mx;
  array[nFarms] matrix[nDays, nAttachCov] attach_env_mx;
  array[nFarms] matrix[nDays, nSurvCov] sal_mx;
  matrix[nDays, nFarms] temp_mx;
  matrix[nDays, nFarms] nFish_mx;
  array[nFarms, nSamples] int sampledDays;
  matrix[nDays, nFarms] nFishSampled_mx;
  // priors: [mean, sd]
  array[nAttachCov, 2] real prior_attach_beta;
  array[nSurvCov, nStages, 2] real prior_surv_beta;
  array[nStages-1, 2] real prior_thresh_GDD_F;
  array[nStages-2, 2] real prior_thresh_GDD_M;
  array[2, nStages-1, 2] real prior_pMolt_F;
  array[2, nStages-2, 2] real prior_pMolt_M;
  array[2] real prior_lifespan;
  array[nStages, 2] real prior_logit_detect_p;
}

transformed data {

}

parameters {
  vector[nStages] logit_detect_p; // p(detect) by stage
  real<lower=0> nb_prec; // neg_binom precision
}

transformed parameters {
  vector[nStages] detect_p;
  real lprior = 0;  // prior contributions to the log posterior

  // re-scale and de-center parameters
  for(stage in 1:nStages) {
    detect_p[stage] = inv_logit(fma(logit_detect_p[stage],
                                    prior_logit_detect_p[stage,2],
                                    prior_logit_detect_p[stage,1]));
  }

  // priors
  lprior += std_normal_lpdf(logit_detect_p);
  lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
}

model {
  target += lprior;
}


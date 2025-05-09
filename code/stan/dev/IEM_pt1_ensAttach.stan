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
  array[nFarms] matrix<lower=0>[nDays, nSims] IP_mx;
  array[nFarms] matrix[nDays, nAttachCov] attach_env_mx;
  array[nFarms] matrix[nDays, nSurvCov] sal_mx;
  matrix[nDays, nFarms] temp_mx;
  matrix[nDays, nFarms] nFish_mx;
  array[nFarms, nSamples] int sampledDays;
  matrix[nDays, nFarms] nFishSampled_mx;
  array[nDays, nFarms] int y_attach;
  // priors: [mean, sd]
  array[nAttachCov, 2] real prior_attach_beta;
  array[nSurvCov, nStages, 2] real prior_surv_beta;
  array[nStages-1, 2] real prior_thresh_GDD_F;
  array[nStages-1, 2] real prior_thresh_GDD_M;
  array[2, nStages-1, 2] real prior_pMolt_F;
  array[2, nStages-1, 2] real prior_pMolt_M;
  array[2] real prior_lifespan;
  array[nStages, 2] real prior_logit_detect_p;
}

parameters {
  vector[nSims] ensWts_p_uc;  // unconstrained mixture proportions
  vector[nAttachCov] attach_beta_z; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  real<lower=0> nb_prec; // neg_binom 1/sqrt(precision)
  real<lower=0> IP_bg_m3; // background infection pressure per m3
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  simplex[nSims] ensWts_p = softmax(ensWts_p_uc);  // mixture proportions
  real<lower=0> IP_bg = IP_bg_m3 * pi() * square(30) * 20;
  vector[nAttachCov] attach_beta; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nDays, nFarms] ensIP;
  matrix[nDays, nFarms] pr_attach;
  matrix[nDays, nFarms] N_attach;
  matrix[nDays, nFarms] ybar_attach;

  // re-scale and de-center parameters
  for(i in 1:nAttachCov) {
    attach_beta[i] = fma(attach_beta_z[i],
                         prior_attach_beta[i,2],
                         prior_attach_beta[i,1]);
  }

  profile("ensIP") {
    // calculate ensIP, pr_attach, N_attach, ybar_attach
    for(farm in 1:nFarms) {
      ensIP[,farm] = IP_mx[farm] * ensWts_p + IP_bg * nPens[farm];
      pr_attach[,farm] = inv_logit(attach_env_mx[farm] * attach_beta);
    }
    N_attach = ensIP .* pr_attach * 0.5;
    for(farm in 1:nFarms) {
      for(day in 1:nDays) {
        if(nFish_mx[day, farm] == 0) {
          N_attach[day, farm] = 0;
          ybar_attach[day, farm] = 0.000001;
        } else {
          ybar_attach[day, farm] = N_attach[day, farm] ./ nFish_mx[day, farm];
        }
      }
    }
  }
  // priors
  lprior += std_normal_lpdf(ensWts_p_uc);
  lprior += std_normal_lpdf(attach_beta_z);
  lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
  lprior += normal_lpdf(IP_bg_m3 | 0, 0.5) - normal_lccdf(0 | 0, 0.5);
}

model {
  profile("likelihood") {
    for(farm in 1:nFarms) {
      target += neg_binomial_2_lpmf(y_attach[, farm] | ybar_attach[, farm], nb_prec);
    }
  }
  target += lprior;
}

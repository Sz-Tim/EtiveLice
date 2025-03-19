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
  real<lower=0> IP_volume;
  real<lower=0> y_bar_minimum;
  array[nFarms] int<lower=0> nPens;
  array[nStages, 2, nDays, nFarms] int y;
  // matrix[nDays, nFarms] N_attach;
  // array[nStages, 2] matrix[nDays, nFarms] mu_true;
  array[nFarms] matrix<lower=0>[nDays, nSims] IP_mx;
  array[nFarms] matrix[nDays, nAttachCov] attach_env_mx;
  array[nFarms] matrix[nDays, nSurvCov] sal_mx;
  matrix[nDays, nFarms] temp_z_mx;
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
  array[2] int nStagesSex = {nStages, nStages-1};
  array[nFarms] matrix[nDays, 2] temp_X;
  matrix[nStages, nStages] zero_matrix_stages_stages = rep_matrix(0, nStages, nStages);
  matrix[nDays, nDays] zero_matrix_days_days = rep_matrix(0, nDays, nDays);
  array[nDays, nStages, 2] int days_array = rep_array(nDays, nDays, nStages, 2);
  array[2, nFarms, nDays] matrix[nStages, nStages] trans_mx_init;
  array[2, nFarms] matrix[nStages, nDays] N_init;
  matrix[nFarms, nDays] t_nFish_mx = nFish_mx';
  matrix[nFarms, nDays] t_nFishSampled_mx = nFishSampled_mx';
  // array[nStages, 2] matrix[nFarms, nDays] mu_true2;
  array[2, nStages, nFarms, nDays] int y2; // reshaped y

  for(farm in 1:nFarms) {
    temp_X[farm, ,1] = ones_vector(nDays);
    temp_X[farm, ,2] = temp_z_mx[, farm];
  }
  for(sex in 1:2) {
    for(farm in 1:nFarms) {
      N_init[sex, farm] = rep_matrix(0, nStages, nDays);
    }
  }
  for(sex in 1:2) {
    for(farm in 1:nFarms) {
      for(day in 1:nDays) {
        trans_mx_init[sex, farm, day] = zero_matrix_stages_stages;
      }
    }
  }
  for(sex in 1:2) {
    for(stage in 1:nStages) {
      for(farm in 1:nFarms) {
        // mu_true2[stage, sex, farm] = mu_true[stage, sex, , farm]';
        y2[sex, stage, farm] = y[stage, sex, , farm];
      }
    }
  }
}

parameters {
  real<lower=0> IP_bg_m3; // background infection pressure per m3
  vector[nSims] ensWts_p_uc;  // unconstrained mixture proportions
  vector[nAttachCov] attach_beta_z; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta_z; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  matrix[2, nStages-1] pMoltF_beta_z; // [Int, temp][Ch-Pr, Pr-Ad, Ad-Gr]
  matrix[2, nStages-2] pMoltM_beta_z; // [Int, temp][Ch-Pr, Pr-Ad]
  vector[nStages] logit_detect_p; // p(detect) by stage
  real<lower=0> nb_prec; // neg_binom precision
}

transformed parameters {
  // parameters
  real lprior = 0;  // prior contributions to the log posterior
  real<lower=0> IP_bg = IP_bg_m3 * IP_volume;
  simplex[nSims] ensWts_p = softmax(ensWts_p_uc);  // mixture proportions
  vector[nAttachCov] attach_beta; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  matrix[2, nStages-1] pMoltF_beta; // [Int, temp][Ch-Pr, Pr-Ad, Ad-Gr]
  matrix[2, nStages-2] pMoltM_beta; // [Int, temp][Ch-Pr, Pr-Ad]
  vector[nStages] detect_p;
  // intermediate quantities
  matrix[nDays, nFarms] ensIP;
  matrix[nDays, nFarms] pr_attach;
  matrix[nDays, nFarms] N_attach;
  array[nFarms] matrix[nDays, nStages] stage_Surv; // survival rates
  array[2, nFarms] matrix[nDays, nStages-1] pMolt; // transition probabilities
  array[2, nFarms, nDays] matrix[nStages, nStages] trans_mx; // transition matrix

  array[2, nFarms] matrix[nStages, nDays] N; // total daily N
  array[2, nFarms] matrix[nStages, nDays] mu; // total daily mu
  array[2, nFarms] matrix[nStages, nDays] y_bar; // expected mean accounting for detection

profile("param_rescale") {
  // re-scale and de-center parameters
  for(i in 1:nAttachCov) {
    attach_beta[i] = fma(attach_beta_z[i],
                         prior_attach_beta[i,2],
                         prior_attach_beta[i,1]);
  }
  for(i in 1:nSurvCov) {
    for(stage in 1:nStages) {
      surv_beta[i,stage] = fma(surv_beta_z[i,stage],
                               prior_surv_beta[i,stage,2],
                               prior_surv_beta[i,stage,1]);
    }
  }
  for(i in 1:2) {
    for(stage in 1:(nStages-1)) {
      pMoltF_beta[i,stage] = fma(pMoltF_beta_z[i,stage],
                                 prior_pMolt_F[i,stage,2],
                                 prior_pMolt_F[i,stage,1]);
    }
  }
  for(i in 1:2) {
    for(stage in 1:(nStages-2)) {
      pMoltM_beta[i,stage] = fma(pMoltM_beta_z[i,stage],
                                 prior_pMolt_M[i,stage,2],
                                 prior_pMolt_M[i,stage,1]);
    }
  }
  for(stage in 1:nStages) {
    detect_p[stage] = inv_logit(fma(logit_detect_p[stage],
                                    prior_logit_detect_p[stage,2],
                                    prior_logit_detect_p[stage,1]));
  }
}
profile("attachment") {
  // calculate ensIP, pr_attach, N_attach, stage_logSurv
  for(farm in 1:nFarms) {
    ensIP[,farm] = IP_mx[farm] * ensWts_p + IP_bg * nPens[farm];
    pr_attach[,farm] = inv_logit(attach_env_mx[farm] * attach_beta);
    stage_Surv[farm] = inv_logit(sal_mx[farm] * surv_beta);
    for(stage in 1:(nStages-2)) {
      pMolt[1, farm, , stage] = inv_logit(temp_X[farm] * pMoltF_beta[, stage]);
      pMolt[2, farm, , stage] = inv_logit(temp_X[farm] * pMoltM_beta[, stage]);
    }
    pMolt[1, farm, , nStages-1] = inv_logit(temp_X[farm] * pMoltF_beta[, nStages-1]);
    pMolt[2, farm, , nStages-1] = zeros_vector(nDays);
  }
  N_attach = ensIP .* pr_attach * 0.5;
}
profile("trans_mx") {
  trans_mx = trans_mx_init;
  for(sex in 1:2) {
    for(farm in 1:nFarms) {
      for(day in 1:nDays) {
        for(stage in 1:(nStages-1)) {
          trans_mx[sex, farm, day, stage, stage] = (1-pMolt[sex, farm, day, stage]) * stage_Surv[farm, day, stage];
          trans_mx[sex, farm, day, stage+1, stage] = pMolt[sex, farm, day, stage] * stage_Surv[farm, day, stage];
        }
        trans_mx[sex, farm, day, nStages, nStages] = stage_Surv[farm, day, nStages] * (sex-1);
      }
    }
  }
}

  N = N_init;

profile("N") {
  for(sex in 1:2) {
    for(farm in 1:nFarms) {
      N[sex, farm, 1, 1] = N_attach[1, farm];
      for(day in 1:(nDays-1)) {
        N[sex, farm, , day+1] = trans_mx[sex, farm, day] * N[sex, farm, , day];
        N[sex, farm, 1, day+1] += N_attach[day+1, farm];
      }
    }
  }
}
profile("observations") {
  for(sex in 1:2) {
    for(farm in 1:nFarms) {
      for(stage in 1:nStagesSex[sex]) {
        mu[sex, farm, stage] = N[sex, farm, stage] ./ t_nFish_mx[farm, ];
        y_bar[sex, farm, stage] = mu[sex, farm, stage] .* t_nFishSampled_mx[farm, ] * detect_p[stage] + y_bar_minimum;
      }
    }
  }
}
profile("priors") {
  // priors
  lprior += normal_lpdf(IP_bg_m3 | 0, 0.5) - normal_lccdf(0 | 0, 0.5);
  lprior += std_normal_lpdf(ensWts_p_uc);
  lprior += std_normal_lpdf(attach_beta_z);
  for(i in 1:nSurvCov) {
    lprior += std_normal_lpdf(surv_beta_z[i,]);
  }
  for(i in 1:2) {
    lprior += std_normal_lpdf(pMoltF_beta_z[i,]);
    lprior += std_normal_lpdf(pMoltM_beta_z[i,]);
  }
  lprior += std_normal_lpdf(logit_detect_p);
  lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
}
}

model {

profile("likelihood") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      for(farm in 1:nFarms) {
        target += neg_binomial_2_lpmf(y2[sex, stage, farm, sampledDays[farm]] |
                                      y_bar[sex, farm, stage, sampledDays[farm]],
                                      nb_prec);
      }
    }
  }
}
  target += lprior;
}


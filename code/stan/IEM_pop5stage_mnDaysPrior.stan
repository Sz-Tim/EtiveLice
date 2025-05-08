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
  array[nFarms] matrix<lower=0>[nDays, nSims] IP_mx;
  array[nFarms] matrix[nDays, nAttachCov] attach_env_mx;
  array[nFarms] matrix[nDays, nSurvCov] sal_mx;
  matrix[nDays, nFarms] temp_z_mx;
  matrix[nDays, nFarms] nFish_mx;
  array[nSamples, 2] int sample_i; // rows: sample; cols: farm, day; sorted by farm, day
  array[nFarms, 2] int sample_ii; // start/end indexes for each farm in sample_i
  matrix[nDays, nFarms] nFishSampled_mx;
  int<lower=0,upper=1> sample_prior_only;
  // priors: [mean, sd]
  array[nAttachCov, 2] real prior_attach_beta;
  array[nSurvCov, nStages, 2] real prior_surv_beta;
  array[2, nStages-1, 2] real prior_mnDaysStage_F;
  array[nStages-1, 2] real prior_logit_detect_p;
}

transformed data {
  array[2] int nStagesSex = {nStages, nStages-1};
  array[nFarms] matrix[nDays, 2] temp_X;
  matrix[nStages+2, nStages+2] zero_matrix_stages_stages = rep_matrix(0, nStages+2, nStages+2);
  array[1, nFarms, nDays] matrix[nStages+2, nStages+2] trans_mx_init;
  array[1, nFarms] matrix[nStages+2, nDays] N_init;
  matrix[nFarms, nDays] t_nFishSampled_mx = nFishSampled_mx';
  array[1, nStages, nSamples] int y2;
  matrix[nDays, nFarms] nFishNoZero_mx = nFish_mx;
  matrix[nDays, nFarms] fishPresent;
  int sex = 1;

  for(farm in 1:nFarms) {
    temp_X[farm, ,1] = ones_vector(nDays);
    temp_X[farm, ,2] = temp_z_mx[, farm];
  }
  for(farm in 1:nFarms) {
    N_init[sex, farm] = rep_matrix(0, nStages+2, nDays);
  }
  for(farm in 1:nFarms) {
    for(day in 1:nDays) {
      trans_mx_init[sex, farm, day] = zero_matrix_stages_stages;
      fishPresent[day, farm] = nFish_mx[day, farm] > 0;
      if(nFish_mx[day, farm] == 0) {
        nFishNoZero_mx[day, farm] = 1;
      }
    }
  }
  for(sample in 1:nSamples) {
    for(stage in 1:nStages) {
      y2[sex, stage, sample] = y[stage, sex, sample_i[sample, 2], sample_i[sample, 1]];
    }
  }
}

parameters {
  real<lower=0> IP_bg_m3; // background infection pressure per m3
  vector[nSims] ensWts_p_uc;  // unconstrained mixture proportions
  vector[nAttachCov] attach_beta_z; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta_z; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  matrix[2, nStages-1] mnDaysStage_beta_z; // [Int, temp][Ch-Pr, Pr-Ad, Ad-Gr]
  vector[nStages] logit_detect_p; // p(detect) by stage
  real<lower=0> nb_prec; // neg_binom precision
}

transformed parameters {
  // parameters
  real lprior = 0;  // prior contributions to the log posterior
  real<lower=0> IP_bg = IP_bg_m3 * IP_volume;
  simplex[nSims] ensWts_p = softmax(ensWts_p_uc);  // mixture proportions
  vector[nAttachCov] attach_beta; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta; // p(surv) [Int, sal][Ch, Pr, Ad]
  matrix[2, nStages-1] pMoltF_beta; // [Int, temp][Ch-Pr, Pr-Ad]
  vector[nStages] detect_p;
  // intermediate quantities
  matrix[nDays, nFarms] ensIP;
  matrix[nDays, nFarms] pr_attach;
  matrix[nDays, nFarms] N_attach;
  array[nFarms] matrix[nDays, nStages] stage_Surv; // survival rates
  array[1, nFarms] matrix[nDays, nStages-1] pMolt; // transition probabilities
  array[1, nFarms, nDays] matrix[nStages+2, nStages+2] trans_mx; // transition matrix
  array[1, nFarms] matrix[nStages+2, nDays] mu; // latent daily mean lice per fish
  array[1, nStages] row_vector[nSamples] y_bar;

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
      pMoltF_beta[i,stage] = fma(mnDaysStage_beta_z[i,stage],
                                   prior_mnDaysStage_F[i,stage,2],
                                   prior_mnDaysStage_F[i,stage,1]/2);
    }
  }
  for(stage in 1:(nStages-1)) {
    detect_p[stage] = inv_logit(fma(logit_detect_p[stage],
                                    prior_logit_detect_p[stage,2],
                                    prior_logit_detect_p[stage,1]));
  }
  detect_p[nStages] = 1;
}
profile("attachment") {
  // calculate ensIP, pr_attach, N_attach, stage_logSurv
  for(farm in 1:nFarms) {
    ensIP[,farm] = IP_mx[farm] * ensWts_p + IP_bg * nPens[farm];
    pr_attach[,farm] = inv_logit(attach_env_mx[farm] * attach_beta);
    stage_Surv[farm] = inv_logit(sal_mx[farm] * surv_beta);
    for(stage in 1:(nStages-1)) {
      pMolt[1, farm, , stage] = 1 / (temp_X[farm] * pMoltF_beta[, stage]);
    }
    for(stage in 1:nStages) {
      stage_Surv[farm, , stage] = stage_Surv[farm, , stage] .* fishPresent[, farm];
    }
  }
  N_attach = ensIP .* pr_attach .* fishPresent * 0.5 ./ nFishNoZero_mx;
}
profile("trans_mx") {
  trans_mx = trans_mx_init;
  for(farm in 1:nFarms) {
    for(day in 1:nDays) {
      // Chalimus 1-2
      trans_mx[sex, farm, day, 1, 1] = (1-pMolt[sex, farm, day, 1]) * stage_Surv[farm, day, 1];
      trans_mx[sex, farm, day, 2, 1] = pMolt[sex, farm, day, 1] * stage_Surv[farm, day, 1];
      trans_mx[sex, farm, day, 2, 2] = (1-pMolt[sex, farm, day, 1]) * stage_Surv[farm, day, 1];
      trans_mx[sex, farm, day, 3, 2] = pMolt[sex, farm, day, 1] * stage_Surv[farm, day, 1];
      // Pre-adult 1-2
      trans_mx[sex, farm, day, 3, 3] = (1-pMolt[sex, farm, day, 2]) * stage_Surv[farm, day, 2];
      trans_mx[sex, farm, day, 4, 3] = pMolt[sex, farm, day, 2] * stage_Surv[farm, day, 2];
      trans_mx[sex, farm, day, 4, 4] = (1-pMolt[sex, farm, day, 2]) * stage_Surv[farm, day, 2];
      trans_mx[sex, farm, day, 5, 4] = pMolt[sex, farm, day, 2] * stage_Surv[farm, day, 2];
      // Adult
      trans_mx[sex, farm, day, 5, 5] = stage_Surv[farm, day, 3];
    }
  }
}

  mu = N_init;

profile("N") {
  for(farm in 1:nFarms) {
    mu[sex, farm, 1, 1] = N_attach[1, farm];
    for(day in 1:(nDays-1)) {
      mu[sex, farm, , day+1] = trans_mx[sex, farm, day] * mu[sex, farm, , day];
      mu[sex, farm, 1, day+1] += N_attach[day+1, farm];
    }
  }
}
profile("observations") {
  for(farm in 1:nFarms) {
    array[2] int farmInd = sample_ii[farm, ];
    // Chalimus: (mu[Ch1] + mu[Ch2]) * nFish * pDet[Ch]
    y_bar[sex, 1, farmInd[1]:farmInd[2]] =
      ones_row_vector(2) * mu[sex, farm, 1:2, sample_i[farmInd[1]:farmInd[2], 2]] .*
      t_nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
      detect_p[1] + y_bar_minimum;
    // Pre-adult: (mu[PA1] + mu[PA2]) * nFish * pDet[PA]
    y_bar[sex, 2, farmInd[1]:farmInd[2]] =
      ones_row_vector(2) * mu[sex, farm, 3:4, sample_i[farmInd[1]:farmInd[2], 2]] .*
      t_nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
      detect_p[2] + y_bar_minimum;
    // Adult: mu[Ad] * nFish * pDet[Ad]
    y_bar[sex, 3, farmInd[1]:farmInd[2]] =
      mu[sex, farm, 5, sample_i[farmInd[1]:farmInd[2], 2]] .*
      t_nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
      detect_p[3] + y_bar_minimum;
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
    lprior += std_normal_lpdf(mnDaysStage_beta_z[i,]);
  }
  lprior += std_normal_lpdf(logit_detect_p);
  lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
}
}

model {

profile("likelihood") {
  for(stage in 1:nStagesSex[sex]) {
    target += neg_binomial_2_lpmf(y2[sex, stage] | y_bar[sex, stage], nb_prec);
  }
}
  target += lprior;
}

generated quantities {
  // array[1, nStages, nFarms, nSamples] int y_pred;
  array[1, nStages, nSamples] int y_pred;
  for(stage in 1:nStagesSex[sex]) {
      y_pred[sex, stage] = neg_binomial_2_rng(y_bar[sex, stage], nb_prec);
  }
}


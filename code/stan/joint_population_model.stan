// Project: E(tive)Lice
// Tim Szewczyk
// tim.szewczyk@sams.ac.uk
// Integrated ensemble model

data {
  int<lower=0> nDays; // number of days from start to end
  int<lower=0> nHours; // number of hours from start to end
  int<lower=0> nFarms; // number of farms
  int<lower=0> nSims; // number of IP simulations for ensembling
  int<lower=0> nStages; // number of on-fish stages: 3 = Ch, PA, Ad
  int<lower=1> nAttachCov; // number of covariates for p(Attach)
  int<lower=1> nSurvCov; // number of covariates for p(Surv) incl. intercept
  int<lower=1> nSamples; // total number of sampling bouts across all farms
  real<lower=0> IP_volume; // per-pen volume to use for scaling IP_m3 to N_copepodids
  real<lower=0> y_bar_minimum; // minimum possible expected lice counts (0 not allowed)
  array[nFarms] int<lower=0> nPens; // number of pens per farm
  array[nDays,24] int<lower=1,upper=nHours> day_hour; // hour numbers corresponding to each day
  array[nStages, nDays, nFarms] int y_F; // FEMALE lice counts for each day 1:nDays
  array[nFarms] matrix<lower=0>[nHours, nSims] IP_mx; // IP from particle tracking simulations
  array[nFarms] matrix[nHours, nAttachCov] attach_env_mx; // covariates for p(Attach)
  array[nFarms] matrix[nDays, nSurvCov] surv_env_mx; // covariates for p(Surv) incl. intercept
  matrix[nDays, nFarms] temp_z_mx; // z-transformed temperature
  matrix<lower=0,upper=1>[nDays, nFarms] treatDays; // indicator for treatment application
  matrix[nDays, nFarms] nFish_mx; // number of fish on each farm each day
  array[nSamples, 2] int sample_i; // rows: sample; cols: farm, day; sorted by farm, day
  array[nFarms, 2] int sample_ii; // start/end indexes for each farm in sample_i
  matrix[nDays, nFarms] nFishSampled_mx; // number of fish sampled on each farm each day
  int<lower=0,upper=1> sample_prior_only; // 0: include likelihood, 1: prior only
  // priors: [mean, sd]
  array[nAttachCov, 2] real prior_attach_beta; // p(attach) slopes
  array[nSurvCov, nStages, 2] real prior_surv_beta; // p(surv) intercept & slopes
  array[2, nStages-1, 2] real prior_mnDaysStage_F; // mean days per stage (Ch, PA)
  array[nStages-1, 2] real prior_logit_detect_p; // detection probability
}

transformed data {
  array[nFarms] matrix[nDays, 2] temp_X;
  matrix[nStages+2, nStages+2] zero_matrix_stages_stages = rep_matrix(0, nStages+2, nStages+2);
  array[nFarms, nDays] matrix[nStages+2, nStages+2] trans_mx_init;
  array[nFarms] matrix[nStages+2, nDays] N_init;
  matrix[nFarms, nDays] t_nFishSampled_mx = nFishSampled_mx';
  array[nStages, nSamples] int y2_F;
  matrix[nDays, nFarms] nFishNoZero_mx = nFish_mx;
  matrix[nDays, nFarms] fishPresent;
  matrix[nDays, nFarms] fishScale;

  for(farm in 1:nFarms) {
    temp_X[farm, ,1] = ones_vector(nDays);
    temp_X[farm, ,2] = temp_z_mx[, farm];
  }
  for(farm in 1:nFarms) {
    N_init[farm] = rep_matrix(0, nStages+2, nDays);
  }
  for(farm in 1:nFarms) {
    for(day in 1:nDays) {
      trans_mx_init[farm, day] = zero_matrix_stages_stages;
      fishPresent[day, farm] = nFish_mx[day, farm] > 0;
      if(nFish_mx[day, farm] == 0) {
        nFishNoZero_mx[day, farm] = 1;
      }
    }
  }
  fishScale = fishPresent ./ nFishNoZero_mx;
  for(sample in 1:nSamples) {
    for(stage in 1:nStages) {
      y2_F[stage, sample] = y_F[stage, sample_i[sample, 2], sample_i[sample, 1]];
    }
  }
}

parameters {
  real<lower=0> IP_bg_m3; // background infection pressure per m3
  vector[nSims] ensWts_p_uc;  // unconstrained mixture proportions
  vector[nAttachCov] attach_beta_z; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta_z; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  real<lower=0,upper=1> treatEfficacy; // mortality induced by treatment
  matrix[2, nStages-1] mnDaysStage_beta_z; // [Int, temp][Ch-Pr, Pr-Ad, Ad-Gr]
  vector[nStages] logit_detect_p; // p(detect) by stage
  // real<lower=1,upper=4> IP_scale; // scaling factor for IP -> N_attach
  real<lower=0> IP_halfSat_m3; // half-saturation constant for attachment rate (cop/m3)
  real<lower=0> nb_prec; // neg_binom precision
}

transformed parameters {
  // parameters
  real lprior = 0;  // prior contributions to the log posterior
  real<lower=0> IP_bg = IP_bg_m3 * IP_volume; // background N_copepodids per pen
  real<lower=0> IP_halfSat = IP_halfSat_m3 * IP_volume; // half-saturation constant for attachment rate per pen
  simplex[nSims] ensWts_p = softmax(ensWts_p_uc);  // mixture proportions
  vector[nAttachCov] attach_beta; // p(attach) [RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta; // p(surv) [Int, sal][Ch, Pr, Ad]
  matrix[2, nStages-1] mnDaysStage_beta; // [Int, temp][Ch-Pr, Pr-Ad]
  vector[nStages] detect_p;
  // intermediate quantities
  matrix[nDays, nFarms] ensIP;
  // matrix[nDays, nFarms] pr_attach;
  matrix[nHours, nFarms] pr_attachSaturated;
  matrix[nHours, nFarms] N_attach;
  array[nFarms] matrix[nDays, nStages] stage_Surv; // survival rates
  array[nFarms] matrix[nDays, nStages-1] pMolt; // transition probabilities
  array[nFarms, nDays] matrix[nStages+2, nStages+2] trans_mx; // transition matrix
  array[nFarms] matrix[nStages+2, nDays] mu; // latent daily mean lice per fish
  array[nStages] row_vector[nSamples] y_bar;

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
  for(i in 1:(nStages-1)) {
    for(stage in 1:(nStages-1)) {
      mnDaysStage_beta[i,stage] = fma(mnDaysStage_beta_z[i,stage],
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

  // calculate ensIP, pr_attach, N_attach, stage_logSurv
  for(farm in 1:nFarms) {
    ensIP[,farm] = IP_mx[farm] * ensWts_p + IP_bg * nPens[farm];
    // pr_attach[,farm] = inv_logit(attach_env_mx[farm] * attach_beta);
    stage_Surv[farm] = inv_logit(surv_env_mx[farm] * surv_beta);
    for(stage in 1:(nStages-1)) {
      pMolt[farm, , stage] = 1 / (temp_X[farm] * mnDaysStage_beta[, stage]);
    }
    for(stage in 1:nStages) {
      stage_Surv[farm, , stage] = stage_Surv[farm, , stage] .*
                                    fishPresent[, farm] .*
                                    (1 - (treatDays[, farm] * treatEfficacy));
    }
  }
  // N_attach = (ensIP.^(1/IP_scale) .* pr_attach).^IP_scale .* fishPresent * 0.5 ./ nFishNoZero_mx;
  for(farm in 1:nFarms) {
    pr_attachSaturated[,farm] = (inv_logit(attach_env_mx[farm] * attach_beta) * IP_halfSat * nPens[farm])
    ./ (ensIP[,farm] + IP_halfSat * nPens[farm]);
  }
  N_attach = 0.5 * ensIP .* pr_attachSaturated;

  trans_mx = trans_mx_init;
  for(farm in 1:nFarms) {
    for(day in 1:nDays) {
      // Chalimus 1-2
      trans_mx[farm, day, 1, 1] = (1-pMolt[farm, day, 1]) * stage_Surv[farm, day, 1];
      trans_mx[farm, day, 2, 1] = pMolt[farm, day, 1] * stage_Surv[farm, day, 1];
      trans_mx[farm, day, 2, 2] = (1-pMolt[farm, day, 1]) * stage_Surv[farm, day, 1];
      trans_mx[farm, day, 3, 2] = pMolt[farm, day, 1] * stage_Surv[farm, day, 1];
      // Pre-adult 1-2
      trans_mx[farm, day, 3, 3] = (1-pMolt[farm, day, 2]) * stage_Surv[farm, day, 2];
      trans_mx[farm, day, 4, 3] = pMolt[farm, day, 2] * stage_Surv[farm, day, 2];
      trans_mx[farm, day, 4, 4] = (1-pMolt[farm, day, 2]) * stage_Surv[farm, day, 2];
      trans_mx[farm, day, 5, 4] = pMolt[farm, day, 2] * stage_Surv[farm, day, 2];
      // Adult
      trans_mx[farm, day, 5, 5] = stage_Surv[farm, day, 3];
    }
  }

  mu = N_init;
  // Population projection: Newly attached Ch1 + existing population transitions
  for(farm in 1:nFarms) {
    mu[farm, 1, 1] = sum(N_attach[day_hour[1], farm]) * fishScale[1,farm];
    for(day in 1:(nDays-1)) {
      mu[farm, , day+1] = trans_mx[farm, day] * mu[farm, , day];
      mu[farm, 1, day+1] += sum(N_attach[day_hour[day+1], farm]) * fishScale[day+1,farm];
    }
  }
  for(farm in 1:nFarms) {
    array[2] int farmInd = sample_ii[farm, ];
    // Chalimus: (mu[Ch1] + mu[Ch2]) * nFishSampled * pDet[Ch]
    y_bar[1, farmInd[1]:farmInd[2]] =
      ones_row_vector(2) * mu[farm, 1:2, sample_i[farmInd[1]:farmInd[2], 2]] .*
      t_nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
      detect_p[1] + y_bar_minimum;
    // Pre-adult: (mu[PA1] + mu[PA2]) * nFishSampled * pDet[PA]
    y_bar[2, farmInd[1]:farmInd[2]] =
      ones_row_vector(2) * mu[farm, 3:4, sample_i[farmInd[1]:farmInd[2], 2]] .*
      t_nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
      detect_p[2] + y_bar_minimum;
    // Adult: mu[Ad] * nFishSampled * pDet[Ad]
    y_bar[3, farmInd[1]:farmInd[2]] =
      mu[farm, 5, sample_i[farmInd[1]:farmInd[2], 2]] .*
      t_nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
      detect_p[3] + y_bar_minimum;
  }
  // priors
  {
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
    lprior += student_t_lpdf(IP_halfSat_m3 | 3, 5, 10) - student_t_lccdf(0 | 3, 5, 10);
  }
}

model {
  if(sample_prior_only==0) {
    for(stage in 1:nStages) {
      target += neg_binomial_2_lpmf(y2_F[stage] | y_bar[stage], nb_prec);
    }
  }
  target += lprior;
}

generated quantities {
  array[nStages, nSamples] int y_pred;
  for(stage in 1:nStages) {
      y_pred[stage] = neg_binomial_2_rng(y_bar[stage], nb_prec);
  }
}


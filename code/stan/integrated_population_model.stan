// Project: E(tive)Lice
// Tim Szewczyk
// tim.szewczyk@sams.ac.uk
// Integrated ensemble model

functions {
  array[,] matrix make_trans_mx(array[,] matrix trans_mx_init,
                                int nFarms,
                                int nDays,
                                array[] int stg_grp_ii,
                                array[] matrix stage_Surv,
                                array[] matrix pMolt) {
    array[nFarms, nDays] matrix[5, 5] trans_mx = trans_mx_init;
    for(farm in 1:nFarms) {
      for(day in 1:nDays) {
        // Chalimus 1-2
        trans_mx[farm, day, 1, 1] = (1-pMolt[farm, day, stg_grp_ii[1]]) * stage_Surv[farm, day, stg_grp_ii[1]];
        trans_mx[farm, day, 2, 1] = pMolt[farm, day, stg_grp_ii[1]] * stage_Surv[farm, day, stg_grp_ii[1]];
        trans_mx[farm, day, 2, 2] = (1-pMolt[farm, day, stg_grp_ii[2]]) * stage_Surv[farm, day, stg_grp_ii[2]];
        trans_mx[farm, day, 3, 2] = pMolt[farm, day, stg_grp_ii[2]] * stage_Surv[farm, day, stg_grp_ii[2]];
        // Pre-adult 1-2
        trans_mx[farm, day, 3, 3] = (1-pMolt[farm, day, stg_grp_ii[3]]) * stage_Surv[farm, day, stg_grp_ii[3]];
        trans_mx[farm, day, 4, 3] = pMolt[farm, day, stg_grp_ii[3]] * stage_Surv[farm, day, stg_grp_ii[3]];
        trans_mx[farm, day, 4, 4] = (1-pMolt[farm, day, stg_grp_ii[4]]) * stage_Surv[farm, day, stg_grp_ii[4]];
        trans_mx[farm, day, 5, 4] = pMolt[farm, day, stg_grp_ii[4]] * stage_Surv[farm, day, stg_grp_ii[4]];
        // Adult
        trans_mx[farm, day, 5, 5] = stage_Surv[farm, day, stg_grp_ii[5]];
      }
    }
    return trans_mx;
  }

  array[] row_vector calc_y_bar(array[] matrix mu,
                                int nFarms,
                                int nStageGroups,
                                int nSamples,
                                array[,] int sample_i,
                                array[,] int sample_ii,
                                matrix nFishSampled_mx,
                                vector detect_p,
                                real y_bar_minimum) {
    array[nStageGroups] row_vector[nSamples] y_bar;
    for(farm in 1:nFarms) {
      array[2] int farmInd = sample_ii[farm, ];
      if(farmInd[1] > 0 && farmInd[2] > 0) {
        // Chalimus: (mu[Ch1] + mu[Ch2]) * nFishSampled * pDet[Ch]
        y_bar[1, farmInd[1]:farmInd[2]] =
          ones_row_vector(2) * mu[farm, 1:2, sample_i[farmInd[1]:farmInd[2], 2]] .*
          nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
          detect_p[1] + y_bar_minimum;
        // Pre-adult: (mu[PA1] + mu[PA2]) * nFishSampled * pDet[PA]
        y_bar[2, farmInd[1]:farmInd[2]] =
          ones_row_vector(2) * mu[farm, 3:4, sample_i[farmInd[1]:farmInd[2], 2]] .*
          nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
          detect_p[2] + y_bar_minimum;
        // Adult: mu[Ad] * nFishSampled * pDet[Ad]
        y_bar[3, farmInd[1]:farmInd[2]] =
          mu[farm, 5, sample_i[farmInd[1]:farmInd[2], 2]] .*
          nFishSampled_mx[farm, sample_i[farmInd[1]:farmInd[2], 2]] *
          detect_p[3] + y_bar_minimum;
      }
    }
    return(y_bar);
  }

}

data {
  // control switches
  int<lower=0,upper=1> sample_prior_only; // 0: include likelihood, 1: prior only
  int<lower=0,upper=1> GQ_ypred; // 0: don't, 1: generate ypred for each sample
  int<lower=0,upper=1> GQ_new; // 0: don't, 1: simulate new dates
  // info
  int<lower=0> nDays; // number of days from start to end
  int<lower=0> nHours; // number of hours from start to end
  int<lower=0> nFarms; // number of farms
  int<lower=0> nSims; // number of IP simulations for ensembling
  int<lower=0> nStages; // number of on-fish stages: 5 = Ch1, Ch2, PA1, PA2, Ad
  int<lower=0> nStageGroups; // number of on-fish stages: 3 = Ch, PA, Ad
  int<lower=1> nAttachCov; // number of covariates for p(Attach)
  int<lower=1> nSurvCov; // number of covariates for p(Surv) incl. intercept
  int<lower=1> nSamples; // total number of sampling bouts across all farms
  int<lower=1> nTrtMethods; // number of treatment methods (bath, wellboat, feed, mechanical)
  int<lower=1> nTrtTypes; // number of treatment types (method:substance)
  vector<lower=0>[nFarms] IP_volume; // per-farm volume to use for scaling IP_m3 to N_copepodids
  real<lower=0> y_bar_minimum; // minimum possible expected lice counts (0 not allowed)
  array[nDays,24] int<lower=1,upper=nHours> day_hour; // hour numbers corresponding to each day
  array[nStages] int<lower=1,upper=nStageGroups> stg_grp_ii; // index for stages into stage groups (1,1,2,2,3)
  // data
  array[nStageGroups, nDays, nFarms] int y_F; // FEMALE lice counts for each day 1:nDays
  array[nFarms] matrix<lower=0>[nHours, nSims] IP_mx; // IP from particle tracking simulations
  array[nFarms] matrix[nHours, nAttachCov] attach_env_mx; // covariates for p(Attach)
  array[nFarms] matrix[nDays, nSurvCov] surv_env_mx; // covariates for p(Surv) incl. intercept
  matrix[nDays, nFarms] temp_z_mx; // z-transformed temperature
  // array[nDays, nFarms] int<lower=0,upper=nTrtTypes> treatApplied;
  array[nFarms] matrix<lower=0,upper=1>[nDays, nTrtTypes] trtApplied; // indicator for treatment application
  matrix[nDays, nFarms] nFish_mx; // number of fish on each farm each day
  array[nSamples, 2] int sample_i; // rows: sample; cols: farm, day; sorted by farm, day
  array[nFarms, 2] int sample_ii; // start/end indexes for each farm in sample_i
  matrix[nFarms, nDays] nFishSampled_mx; // number of fish sampled on each farm each day
  // priors: [mean, sd] unless otherwise stated
  vector[2] prior_IP_bg_m3; // background IP
  array[nAttachCov, 2] real prior_attach_beta; // p(attach) slopes
  array[nSurvCov, nStageGroups, 2] real prior_surv_beta; // p(surv) intercept & slopes
  vector<lower=0>[3] prior_surv_int_farm_sd; // p(surv) int sd: student_t(nu, mu, sd)
  vector[2] prior_logit_trtEff_global; // treatment efficacy: normal(mu, sd)
  vector<lower=0>[3] prior_trtEff_sd_types; // treatment efficacy: sd among types within methods
  array[2, nStageGroups-1, 2] real prior_mnDaysStage_F; // mean days per stage (Ch, PA)
  array[nStageGroups-1, 2] real prior_logit_detect_p; // detection probability
  vector[2] prior_inv_sqrt_nb_prec; // negative binomial precision: normal(mu, sd)
  // Generated Quantities: For predicting new observations starting nDays+1
  int<lower=0> nDays_GQ; // number of days from start to end
  int<lower=0> nHours_GQ; // number of hours from start to end
  int<lower=1> nSamples_GQ; // total number of sampling bouts across all farms
  array[nDays_GQ,24] int<lower=1,upper=nHours_GQ> day_hour_GQ; // hour numbers corresponding to each day
  array[nFarms] matrix<lower=0>[nHours_GQ, nSims] IP_mx_GQ; // IP from particle tracking simulations
  array[nFarms] matrix[nHours_GQ, nAttachCov] attach_env_mx_GQ; // covariates for p(Attach)
  array[nFarms] matrix[nDays_GQ, nSurvCov] surv_env_mx_GQ; // covariates for p(Surv) incl. intercept
  matrix[nDays_GQ, nFarms] temp_z_mx_GQ; // z-transformed temperature
  matrix[nDays_GQ, nFarms] nFish_mx_GQ; // number of fish on each farm each day
  array[nFarms] matrix<lower=0,upper=1>[nDays_GQ, nTrtTypes] trtApplied_GQ; // indicator for treatment application
  array[nSamples_GQ, 2] int sample_i_GQ; // rows: sample; cols: farm, day; sorted by farm, day
  array[nFarms, 2] int sample_ii_GQ; // start/end indexes for each farm in sample_i
  matrix[nFarms, nDays_GQ] nFishSampled_mx_GQ; // number of fish sampled on each farm each day
}

transformed data {
  array[nStageGroups, nSamples] int y_F_sparse; // sparse: observations only
  array[nFarms] matrix[nDays, 2] temp_X; // reshape and add intercept 1's
  matrix[nStages, nStages] zero_matrix_stages_stages = rep_matrix(0, nStages, nStages);
  array[nFarms, nDays] matrix[nStages, nStages] trans_mx_init; // initial 0's
  array[nFarms] matrix[nStages, nDays] N_init; // initial 0's
  matrix[nDays, nFarms] nFishNoZero_mx = nFish_mx; // nFish with 0's replaced by 1's
  matrix[nDays, nFarms] fishPresent; // indicator; 1=yes, 0=no
  matrix[nDays, nFarms] fishScale; // number of attached copepodids --> mean per fish
  // Generated Quantities: For predicting new observations starting nDays+1
  array[nFarms] matrix[nDays_GQ, 2] temp_X_GQ;
  array[nFarms, nDays_GQ] matrix[nStages, nStages] trans_mx_init_GQ;
  array[nFarms] matrix[nStages, nDays_GQ] N_init_GQ;
  matrix[nDays_GQ, nFarms] nFishNoZero_mx_GQ = nFish_mx_GQ;
  matrix[nDays_GQ, nFarms] fishPresent_GQ;
  matrix[nDays_GQ, nFarms] fishScale_GQ;

  for(farm in 1:nFarms) {
    temp_X[farm, ,1] = ones_vector(nDays);
    temp_X[farm, ,2] = temp_z_mx[, farm];
  }
  for(farm in 1:nFarms) {
    N_init[farm] = rep_matrix(0, nStages, nDays);
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
    for(grp in 1:nStageGroups) {
      y_F_sparse[grp, sample] = y_F[grp, sample_i[sample, 2], sample_i[sample, 1]];
    }
  }

  if(GQ_new==1) {
    for(farm in 1:nFarms) {
      temp_X_GQ[farm, ,1] = ones_vector(nDays_GQ);
      temp_X_GQ[farm, ,2] = temp_z_mx_GQ[, farm];
    }
    for(farm in 1:nFarms) {
      N_init_GQ[farm] = rep_matrix(0, nStages, nDays_GQ);
    }
    for(farm in 1:nFarms) {
      for(day in 1:nDays_GQ) {
        trans_mx_init_GQ[farm, day] = zero_matrix_stages_stages;
        fishPresent_GQ[day, farm] = nFish_mx_GQ[day, farm] > 0;
        if(nFish_mx_GQ[day, farm] == 0) {
          nFishNoZero_mx_GQ[day, farm] = 1;
        }
      }
    }
    fishScale_GQ = fishPresent_GQ ./ nFishNoZero_mx_GQ;
  }
}

parameters {
  // NB: Prior mean and sd are applied using fma(param_z, prior_sd, prior_mean)
  // with param_z ~ Norm(0, 1), solely for more efficient execution in Stan.
  // This is equivalent to param ~ N(prior_mean, prior_sd)
  real<lower=0> IP_bg_m3; // background infection pressure per m3
  vector[nSims] ensWts_p_uc;  // unconstrained mixture proportions
  vector[nAttachCov-1] attach_beta_z; // p(attach) [RW_logit, sal_z, temp_z, uv_z]
  real<upper=0> attach_betaUV2_z; // p(attach) uv_z^2 coefficient: constrain to concave down
  matrix[nSurvCov, nStageGroups] surv_beta_z; // p(surv) [Int, sal][Ch1-2, PA1-2, Ad]
  row_vector<lower=0>[nStageGroups] surv_int_farm_sd; // p(surv) sd for farm-level Int; logit-scale
  matrix[nFarms, nStageGroups] surv_int_farm_z; // p(surv) [farm][Ch1-2, PA1-2, Ad]
  real logit_trtEff_global; // global treatment efficacy mean
  vector[nTrtTypes] logit_trtEff_type; // treatment efficacies for each type
  real<lower=0> trtEff_sd_types; // treatment efficacy sd among types
  matrix[2, nStageGroups-1] mnDaysStage_beta_z; // [Int, temp][Ch-PA, PA-Ad]
  vector[nStageGroups] logit_detect_p; // p(detect) by stage group
  real<lower=0> inv_sqrt_nb_prec; // negative binomial precision (1/sqrt(prec))
}

transformed parameters {
  // parameters
  real lprior = 0;  // prior contributions to the log posterior
  vector<lower=0>[nFarms] IP_bg = IP_bg_m3 * IP_volume; // background N_copepodids per pen
  simplex[nSims] ensWts_p = softmax(ensWts_p_uc);  // mixture proportions
  vector[nAttachCov] attach_beta; // p(attach) [RW_logit, sal_z, temp_z, uv_z, uv_z^2]
  vector<lower=0,upper=1>[nTrtTypes] trtEff_type = inv_logit(logit_trtEff_type); // treatment efficacy (0-1)
  vector[nTrtTypes] log1m_trtEff_type = log1m(trtEff_type); // log(1 - trtEff_type)
  matrix[nSurvCov, nStageGroups] surv_beta; // p(surv) [Int, sal][Ch1-2, PA1-2, Ad]
  array[nFarms] matrix[nSurvCov, nStageGroups] surv_beta_farm; // p(surv) [farmInt, sal][Ch1-2, PA1-2, Ad]
  matrix[2, nStageGroups-1] mnDaysStage_beta; // [Int, temp][Ch-PA, PA-Ad]
  vector[nStageGroups] detect_p;
  // intermediate quantities
  matrix<lower=0>[nHours, nFarms] ensIP; // ensemble IP
  matrix<lower=0,upper=1>[nHours, nFarms] pr_attach; // attachment probability
  matrix<lower=0>[nHours, nFarms] N_attach; // number of copepodids that attach
  array[nFarms] matrix<lower=0,upper=1>[nDays, nStageGroups] stage_Surv; // survival rates
  array[nFarms] matrix<lower=0,upper=1>[nDays, nStageGroups-1] pMolt; // transition probabilities
  array[nFarms, nDays] matrix<lower=0,upper=1>[nStages, nStages] trans_mx; // transition matrix
  array[nFarms] matrix<lower=0>[nStages, nDays] mu; // latent daily mean lice per fish
  array[nStageGroups] row_vector<lower=0>[nSamples] y_bar; // expected observed (mu*prDet)
  real<lower=0> nb_prec = 1 / inv_sqrt_nb_prec^2; // neg_binom precision

  // re-scale and de-center parameters
  for(i in 1:(nAttachCov-1)) {
    attach_beta[i] = fma(attach_beta_z[i],
                         prior_attach_beta[i,2],
                         prior_attach_beta[i,1]);
  }
  attach_beta[nAttachCov] = fma(attach_betaUV2_z,
                                prior_attach_beta[nAttachCov,2],
                                prior_attach_beta[nAttachCov,1]);
  for(i in 1:nSurvCov) {
    for(grp in 1:nStageGroups) {
      surv_beta[i,grp] = fma(surv_beta_z[i,grp],
                               prior_surv_beta[i,grp,2],
                               prior_surv_beta[i,grp,1]);
    }
  }
  for(farm in 1:nFarms) {
      surv_beta_farm[farm,1,] = surv_beta[1,] + surv_int_farm_z[farm,] .* surv_int_farm_sd;
      surv_beta_farm[farm,2,] = surv_beta[2,];
  }
  for(i in 1:2) {
    for(grp in 1:(nStageGroups-1)) {
      mnDaysStage_beta[i,grp] = fma(mnDaysStage_beta_z[i,grp],
                                      prior_mnDaysStage_F[i,grp,2],
                                      prior_mnDaysStage_F[i,grp,1]);
    }
  }
  for(grp in 1:(nStageGroups-1)) {
    detect_p[grp] = inv_logit(fma(logit_detect_p[grp],
                                    prior_logit_detect_p[grp,2],
                                    prior_logit_detect_p[grp,1]));
  }
  detect_p[nStageGroups] = 1; // all adults are detected

  // calculate ensIP, pr_attach, N_attach, stage_Surv
  for(farm in 1:nFarms) {
    ensIP[,farm] = IP_mx[farm] * ensWts_p + IP_bg[farm];
    pr_attach[,farm] = inv_logit(attach_env_mx[farm] * attach_beta);
    stage_Surv[farm] = inv_logit(surv_env_mx[farm] * surv_beta_farm[farm]);
    for(grp in 1:(nStageGroups-1)) {
      pMolt[farm, , grp] = 1 / (temp_X[farm] * mnDaysStage_beta[, grp]);
    }
    for(grp in 1:nStageGroups) {
      stage_Surv[farm, , grp] = stage_Surv[farm, , grp] .* fishPresent[, farm] .* exp(trtApplied[farm] * log1m_trtEff_type);
    }
  }
  N_attach = 0.5 * ensIP .* pr_attach;
  trans_mx = make_trans_mx(trans_mx_init, nFarms, nDays, stg_grp_ii, stage_Surv, pMolt);

  mu = N_init;
  // Population projection: Newly attached Ch1 + existing population transitions
  for(farm in 1:nFarms) {
    mu[farm, 1, 1] = sum(N_attach[day_hour[1], farm]) * fishScale[1,farm];
    for(day in 1:(nDays-1)) {
      mu[farm, , day+1] = trans_mx[farm, day] * mu[farm, , day];
      mu[farm, 1, day+1] += sum(N_attach[day_hour[day+1], farm]) * fishScale[day+1,farm];
    }
  }
  y_bar = calc_y_bar(mu, nFarms, nStageGroups, nSamples, sample_i, sample_ii, nFishSampled_mx, detect_p, y_bar_minimum);
  // priors
  {
    lprior += normal_lpdf(IP_bg_m3 | prior_IP_bg_m3[1], prior_IP_bg_m3[2]) -
      normal_lccdf(0 | prior_IP_bg_m3[1], prior_IP_bg_m3[2]);
    lprior += std_normal_lpdf(ensWts_p_uc);
    lprior += std_normal_lpdf(attach_beta_z);
    lprior += std_normal_lpdf(attach_betaUV2_z) - std_normal_lcdf(0);
    lprior += normal_lpdf(logit_trtEff_global | prior_logit_trtEff_global[1], prior_logit_trtEff_global[2]);
    lprior += normal_lpdf(logit_trtEff_type | logit_trtEff_global, trtEff_sd_types);
    lprior += student_t_lpdf(trtEff_sd_types | prior_trtEff_sd_types[1], prior_trtEff_sd_types[2], prior_trtEff_sd_types[3]) -
      student_t_lccdf(0 | prior_trtEff_sd_types[1], prior_trtEff_sd_types[2], prior_trtEff_sd_types[3]);
    for(i in 1:nSurvCov) {
        lprior += std_normal_lpdf(surv_beta_z[i,]);
    }
    for(grp in 1:nStageGroups) {
      lprior += std_normal_lpdf(surv_int_farm_z[,grp]);
    }
    for(i in 1:2) {
      lprior += std_normal_lpdf(mnDaysStage_beta_z[i,]);
    }
    lprior += std_normal_lpdf(logit_detect_p);
    lprior += normal_lpdf(inv_sqrt_nb_prec | prior_inv_sqrt_nb_prec[1], prior_inv_sqrt_nb_prec[2]) -
      normal_lccdf(0 | prior_inv_sqrt_nb_prec[1], prior_inv_sqrt_nb_prec[2]);
    lprior += student_t_lpdf(surv_int_farm_sd | prior_surv_int_farm_sd[1], prior_surv_int_farm_sd[2], prior_surv_int_farm_sd[3]) -
      student_t_lccdf(0 | prior_surv_int_farm_sd[1], prior_surv_int_farm_sd[2], prior_surv_int_farm_sd[3]);
  }
}

model {
  if(sample_prior_only==0) {
    for(grp in 1:nStageGroups) {
      target += neg_binomial_2_lpmf(y_F_sparse[grp] | y_bar[grp], nb_prec);
    }
  }
  target += lprior;
}

generated quantities {
  array[nStages, nSamples] int y_pred;
  matrix[nHours_GQ, nFarms] ensIP_GQ;
  matrix[nHours_GQ, nFarms] pr_attach_GQ;
  matrix[nHours_GQ, nFarms] N_attach_GQ;
  array[nFarms] matrix[nDays_GQ, nStageGroups] stage_Surv_GQ; // survival rates
  array[nFarms] matrix[nDays_GQ, nStageGroups-1] pMolt_GQ; // transition probabilities
  array[nFarms, nDays_GQ] matrix[nStages, nStages] trans_mx_GQ; // transition matrix
  array[nFarms] matrix[nStages, nDays_GQ] mu_GQ; // latent daily mean lice per fish
  array[nStageGroups] row_vector[nSamples_GQ] y_bar_GQ; // expected observed mean lice per fish

  // predicted y for each observation in training data
  if(GQ_ypred==1) {
    for(grp in 1:nStageGroups) {
        y_pred[grp] = neg_binomial_2_rng(y_bar[grp], nb_prec);
    }
  }
  // predictions for NEW observations
  if(GQ_new==1) {
    // calculate ensIP_GQ, pr_attach_GQ, N_attach_GQ, stage_Surv_GQ
    for(farm in 1:nFarms) {
      ensIP_GQ[,farm] = IP_mx_GQ[farm] * ensWts_p + IP_bg[farm];
      pr_attach_GQ[,farm] = inv_logit(attach_env_mx_GQ[farm] * attach_beta);
      stage_Surv_GQ[farm] = inv_logit(surv_env_mx_GQ[farm] * surv_beta_farm[farm]);
      for(grp in 1:(nStageGroups-1)) {
        pMolt_GQ[farm, , grp] = 1 / (temp_X_GQ[farm] * mnDaysStage_beta[, grp]);
      }
      for(grp in 1:nStageGroups) {
        stage_Surv_GQ[farm, , grp] = stage_Surv_GQ[farm, , grp] .* fishPresent_GQ[, farm] .* exp(trtApplied_GQ[farm] * log1m_trtEff_type);
      }
    }
    N_attach_GQ = 0.5 * ensIP_GQ .* pr_attach_GQ;
    trans_mx_GQ = make_trans_mx(trans_mx_init_GQ, nFarms, nDays_GQ, stg_grp_ii, stage_Surv_GQ, pMolt_GQ);
    // Population projection: Newly attached Ch1 + existing population transitions
    mu_GQ = N_init_GQ;
    for(farm in 1:nFarms) {
      // day 1 = final training day + 1
      mu_GQ[farm, , 1] = trans_mx[farm, nDays] * mu[farm, , nDays];
      mu_GQ[farm, 1, 1] += sum(N_attach_GQ[day_hour_GQ[1], farm]) * fishScale_GQ[1,farm];
      for(day in 1:(nDays_GQ-1)) {
        mu_GQ[farm, , day+1] = trans_mx_GQ[farm, day] * mu_GQ[farm, , day];
        mu_GQ[farm, 1, day+1] += sum(N_attach_GQ[day_hour_GQ[day+1], farm]) * fishScale_GQ[day+1,farm];
      }
    }
    y_bar_GQ = calc_y_bar(mu_GQ, nFarms, nStageGroups, nSamples_GQ, sample_i_GQ, sample_ii_GQ, nFishSampled_mx_GQ, detect_p, y_bar_minimum);
  }
}


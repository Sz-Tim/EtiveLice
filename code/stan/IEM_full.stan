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
  array[2] int nStagesSex = {nStages, nStages-1};
  array[nFarms, nDays] vector[nDays] cohort_GDD;
  matrix[nDays, nDays] zero_matrix_days_days = rep_matrix(0, nDays, nDays);
  array[nDays, nStages, 2] int days_array = rep_array(nDays, nDays, nStages, 2);
  array[nFarms, 2] matrix[nDays, nDays] cohort_zero_array;
  array[nFarms, 2, nDays, nStages, 2] int cohort_molts_init;
  array[2, nStages] matrix[nDays, nFarms] N_init;
  // array[nStages, 2] matrix[nFarms, nDays] mu_true2;
  array[2, nStages, nFarms, nDays] int y2; // reshaped y


  for(sex in 1:2) {
    for(stage in 1:nStages) {
      N_init[sex, stage] = rep_matrix(0, nDays, nFarms);
      for(farm in 1:nFarms) {
        // mu_true2[stage, sex, farm] = mu_true[stage, sex, , farm]';
        y2[sex, stage, farm] = y[stage, sex, , farm];
      }
    }
  }
  for(farm in 1:nFarms) {
    for(sex in 1:2) {
      cohort_zero_array[farm, sex] = zero_matrix_days_days;
      cohort_molts_init[farm, sex, , , ] = days_array;
      cohort_molts_init[farm, sex, , 1, 1] = linspaced_int_array(nDays, 1, nDays); // chalimus start index
    }
    for(cohort in 1:nDays) {
      cohort_GDD[farm, cohort] = zeros_vector(nDays);
      cohort_GDD[farm, cohort, cohort:nDays] = cumulative_sum(temp_mx[cohort:nDays, farm]);
    }
  }
}

parameters {
  real<lower=0> IP_bg_m3; // background infection pressure per m3
  vector[nSims] ensWts_p_uc;  // unconstrained mixture proportions
  vector[nAttachCov] attach_beta_z; // p(attach) [Int, RW_logit, sal_z, uv, uv^2]
  matrix[nSurvCov, nStages] surv_beta_z; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  real lifespan_z;
  positive_ordered[3] thresh_GDD_F_z;
  positive_ordered[2] thresh_GDD_M_z;
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
  real lifespan;
  matrix[nStages-1, 2] thresh_GDD; // GDD molting thresholds [Pr, Ad, Gr][F, M]
  vector[nStages] detect_p;
  // intermediate quantities
  matrix[nDays, nFarms] ensIP;
  matrix[nDays, nFarms] pr_attach;
  matrix[nDays, nFarms] N_attach;
  array[nFarms] matrix[nDays, nStages] stage_logSurv;

profile("param_rescale") {
  // re-scale and de-center parameters: implies param ~ Norm(prior[2], prior[1])
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
  lifespan = fma(lifespan_z,
                 prior_lifespan[2],
                 prior_lifespan[1]);
  for(stage in 1:(nStages-1)) {
    thresh_GDD[stage,1] = fma(thresh_GDD_F_z[stage],
                              prior_thresh_GDD_F[stage,2],
                              prior_thresh_GDD_F[stage,1]);
  }
  for(stage in 1:(nStages-2)) {
    thresh_GDD[stage,2] = fma(thresh_GDD_M_z[stage],
                              prior_thresh_GDD_M[stage,2],
                              prior_thresh_GDD_M[stage,1]);
  }
  thresh_GDD[3,2] = lifespan + 10; // M to gravid, not possible
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
    stage_logSurv[farm] = log_inv_logit(sal_mx[farm] * surv_beta);
  }
  N_attach = ensIP .* pr_attach * 0.5;
}
profile("priors") {
  // priors
  lprior += normal_lpdf(IP_bg_m3 | 0, 0.5) - normal_lccdf(0 | 0, 0.5);
  lprior += std_normal_lpdf(ensWts_p_uc);
  lprior += std_normal_lpdf(attach_beta_z);
  for(i in 1:nSurvCov) {
    lprior += std_normal_lpdf(surv_beta_z[i,]);
  }
  lprior += std_normal_lpdf(lifespan_z);
  lprior += std_normal_lpdf(thresh_GDD_F_z);
  lprior += std_normal_lpdf(thresh_GDD_M_z);
  lprior += std_normal_lpdf(logit_detect_p);
  lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
}
}

model {
  array[nFarms, 2] matrix[nDays, nDays] cohort_N; // cohort daily N
  array[nFarms, 2] matrix[nDays, nDays] cohort_logSurv; // cohort daily s
  array[nFarms, 2, nDays, nStages, 2] int cohort_molts; // day indexes for stages
  array[2, nStages] matrix[nDays, nFarms] N; // total daily N
  array[2, nStages] matrix[nDays, nFarms] mu; // total daily mu
  array[2, nStages] matrix[nDays, nFarms] y_bar; // expected mean accounting for detection

  cohort_N = cohort_zero_array;
  cohort_logSurv = cohort_zero_array;
  cohort_molts = cohort_molts_init;
  N = N_init;

profile("cohort_stage") {
  for(farm in 1:nFarms) {
    for(sex in 1:2) {
      for(cohort in 1:nDays) {
        int current_stage = 1;
        for(day in (cohort+1):nDays) {
          if(current_stage <= nStages) { // is it still alive?
            if(cohort_GDD[farm, cohort, day] > lifespan) { // has it reached its lifespan?
              cohort_molts[farm, sex, cohort, current_stage:nStages, 2] = rep_array(day, nStages-current_stage+1);
              current_stage = nStages + 1;
            }
            if(current_stage < nStages) { // is there another life stage?
              if(cohort_GDD[farm, cohort, day] > thresh_GDD[current_stage, sex]) { // can it molt?
                  cohort_molts[farm, sex, cohort, current_stage, 2] = day;
                current_stage += 1;
                cohort_molts[farm, sex, cohort, current_stage, 1] = min(day+1, nDays);
              }
            }
          }
        }
      }
    }
  }
}
profile("cohort_surv") {
  for(farm in 1:nFarms) {
    for(sex in 1:2) {
      for(cohort in 1:nDays) {
        for(stage in 1:nStagesSex[sex]) {
          cohort_logSurv[farm, sex, cohort_molts[farm, sex, cohort, stage, 1]:cohort_molts[farm, sex, cohort, stage, 2], cohort] =
            stage_logSurv[farm, cohort_molts[farm, sex, cohort, stage, 1]:cohort_molts[farm, sex, cohort, stage, 2], stage];
        }
      }
    }
  }
}
// map_rect this since each is returning a (row) vector?
profile("cohort_N") {
  for(farm in 1:nFarms) {
    for(sex in 1:2) {
      // initial cohort N: 0 + N_attach for diagonal = cohorts' day 1
      cohort_N[farm, sex] = add_diag(cohort_N[farm, sex], N_attach[, farm]);
      for(cohort in 1:nDays) {
          // N[start:end] = N[start] * cumulprod(surv[1:end])
          cohort_N[farm, sex, cohort:cohort_molts[farm, sex, cohort, nStages, 2], cohort] =
            cohort_N[farm, sex, cohort, cohort] *
            exp(cumulative_sum(cohort_logSurv[farm, sex, cohort:cohort_molts[farm, sex, cohort, nStages, 2], cohort]));
      }
    }
  }
}
profile("N") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      for(farm in 1:nFarms) {
        for(cohort in 1:nDays) {
          N[sex, stage, cohort_molts[farm, sex, cohort, stage, 1]:cohort_molts[farm, sex, cohort, stage, 2], farm] +=
            cohort_N[farm, sex, cohort_molts[farm, sex, cohort, stage, 1]:cohort_molts[farm, sex, cohort, stage, 2], cohort];
        }
      }
    }
  }
}
profile("observations") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      mu[sex, stage] = N[sex, stage] ./ nFish_mx;
      y_bar[sex, stage] = mu[sex, stage] .* nFishSampled_mx * detect_p[stage] + y_bar_minimum;
    }
  }
}
profile("likelihood") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      for(farm in 1:nFarms) {
        target += neg_binomial_2_lpmf(y2[sex, stage, farm, sampledDays[farm]] |
                                      y_bar[sex, stage, sampledDays[farm], farm],
                                      nb_prec);
      }
    }
  }
}
  target += lprior;
}


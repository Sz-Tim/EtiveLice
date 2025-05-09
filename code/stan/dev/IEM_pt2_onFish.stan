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
  array[2] int nStagesSex = {nStages, nStages-1};
  array[nFarms, nDays] vector[nDays] cohort_GDD;
  matrix[nDays, nDays] zero_matrix_days_days = rep_matrix(0, nDays, nDays);
  array[nDays, nStages, 2] int days_array = rep_array(nDays, nDays, nStages, 2);
  array[nFarms, 2] matrix[nDays, nDays] cohort_zero_array;
  array[nFarms, 2, nDays, nStages, 2] int cohort_molts_init;
  array[2, nStages] matrix[nDays, nFarms] N_init;

  for(sex in 1:2) {
    for(stage in 1:nStages) {
      N_init[sex, stage] = rep_matrix(0, nDays, nFarms);
    }
  }
  for(farm in 1:nFarms) {
    for(sex in 1:2) {
      cohort_zero_array[farm, sex] = zero_matrix_days_days;
      cohort_molts_init[farm, sex, , , ] = days_array;
      cohort_molts_init[farm, sex, , 1, 1] = linspaced_int_array(nDays, 1, nDays); // chalimus start index
    }
  }
  for(farm in 1:nFarms) {
    for(cohort in 1:nDays) {
      cohort_GDD[farm, cohort] = zeros_vector(nDays);
      cohort_GDD[farm, cohort, cohort:nDays] = cumulative_sum(temp_mx[cohort:nDays, farm]);
    }
  }
}

parameters {
  matrix[nSurvCov, nStages] surv_beta_z; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  real lifespan_z;
  positive_ordered[3] thresh_GDD_F_z;
  positive_ordered[2] thresh_GDD_M_z;
  real<lower=0> nb_prec; // neg_binom precision
}

transformed parameters {
  matrix[nSurvCov, nStages] surv_beta; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  real lifespan;
  matrix[nStages-1, 2] thresh_GDD; // GDD molting thresholds [Pr, Ad, Gr][F, M]
  real lprior = 0;  // prior contributions to the log posterior
  array[nFarms] matrix[nDays, nStages] stage_logSurv;

  profile("trans_params") {
    // re-scale and de-center parameters
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

    for(farm in 1:nFarms) {
      stage_logSurv[farm] = log_inv_logit(sal_mx[farm] * surv_beta);
    }

    // priors
    for(i in 1:nSurvCov) {
      lprior += std_normal_lpdf(surv_beta_z[i,]);
    }
    lprior += std_normal_lpdf(thresh_GDD_F_z);
    lprior += std_normal_lpdf(thresh_GDD_M_z);
    lprior += std_normal_lpdf(lifespan_z);
    lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
  }
}

model {
  array[nFarms, 2] matrix[nDays, nDays] cohort_N; // cohort daily N
  array[nFarms, 2] matrix[nDays, nDays] cohort_logSurv; // cohort daily s
  array[nFarms, 2, nDays, nStages, 2] int cohort_molts; // day indexes for stages
  array[2, nStages] matrix[nDays, nFarms] N; // total daily N
  array[2, nStages] matrix[nDays, nFarms] mu; // total daily mu

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
          if(current_stage <= nStagesSex[sex]) { // is it still alive?
            if(cohort_GDD[farm, cohort, day] > lifespan) { // has it reached its lifespan?
              cohort_molts[farm, sex, cohort, current_stage:nStages, 2] = rep_array(day, nStages-current_stage+1);
              current_stage = nStages + 1;
            }
            if(current_stage < nStagesSex[sex]) { // is there another life stage?
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
profile("mu") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      mu[sex, stage] = N[sex, stage] ./ nFish_mx;
    }
  }
}
profile("likelihood") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      for(farm in 1:nFarms) {
        target += normal_lpdf(mu_true[stage, sex, sampledDays[farm], farm] |
                              mu[sex, stage, sampledDays[farm], farm],
                              nb_prec);
      }
    }
  }
}

  target += lprior;
}


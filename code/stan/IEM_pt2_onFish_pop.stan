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
  array[nDays, nStages, 2] int days_array = rep_array(nDays, nDays, nStages, 2);
  array[2, nFarms, nDays] matrix[nStages, nStages] trans_mx_init;
  array[2, nFarms] matrix[nStages, nDays] N_init;
  matrix[nFarms, nDays] t_nFish_mx = nFish_mx';

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
}

parameters {
  matrix[nSurvCov, nStages] surv_beta_z; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  matrix[2, nStages-1] pMoltF_beta_z; // [Int, temp][Ch-Pr, Pr-Ad, Ad-Gr]
  matrix[2, nStages-2] pMoltM_beta_z; // [Int, temp][Ch-Pr, Pr-Ad]
  real<lower=0> nb_prec; // neg_binom precision
}

transformed parameters {
  matrix[nSurvCov, nStages] surv_beta; // p(surv) [Int, sal][Ch, Pr, Ad, Gr]
  matrix[2, nStages-1] pMoltF_beta; // [Int, temp][Ch-Pr, Pr-Ad, Ad-Gr]
  matrix[2, nStages-2] pMoltM_beta; // [Int, temp][Ch-Pr, Pr-Ad]
  real lprior = 0;  // prior contributions to the log posterior
  array[nFarms] matrix[nDays, nStages] stage_Surv; // survival rates
  array[2, nFarms] matrix[nDays, nStages-1] pMolt; // transition probabilities
  array[2, nFarms, nDays] matrix[nStages, nStages] trans_mx; // transition matrix
  array[2, nFarms] matrix[nStages, nDays] N; // total daily N
  array[2, nFarms] matrix[nStages, nDays] mu; // total daily mu

  profile("trans_params") {
    // re-scale and de-center parameters
    for(i in 1:nSurvCov) {
    for(stage in 1:nStages) {
      surv_beta[i,stage] = fma(surv_beta_z[i,stage],
                               prior_surv_beta[i,stage,2],
                               prior_surv_beta[i,stage,1]);
    }
  }
    pMoltF_beta[1,1] = fma(pMoltF_beta_z[1,1], 0.5, -2.6); // 1/(15 days)
    pMoltF_beta[1,2] = fma(pMoltF_beta_z[1,2], 0.5, -2.9); // 1/(20 days)
    pMoltF_beta[1,3] = fma(pMoltF_beta_z[1,3], 0.5, -2.9); // 1/(20 days)
    pMoltM_beta[1,1] = fma(pMoltM_beta_z[1,1], 0.5, -2.6); // 1/(15 days)
    pMoltM_beta[1,2] = fma(pMoltM_beta_z[1,2], 0.5, -2.9); // 1/(20 days)
    pMoltF_beta[2,] = fma(pMoltF_beta_z[2,], 0.5, 0.5);
    pMoltM_beta[2,] = fma(pMoltM_beta_z[2,], 0.5, 0.5);

    for(farm in 1:nFarms) {
      stage_Surv[farm] = inv_logit(sal_mx[farm] * surv_beta);
      for(stage in 1:(nStages-2)) {
        pMolt[1, farm, , stage] = inv_logit(temp_X[farm] * pMoltF_beta[, stage]);
        pMolt[2, farm, , stage] = inv_logit(temp_X[farm] * pMoltM_beta[, stage]);
      }
      pMolt[1, farm, , nStages-1] = inv_logit(temp_X[farm] * pMoltF_beta[, nStages-1]);
      pMolt[2, farm, , nStages-1] = zeros_vector(nDays);
    }
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

profile("N") {
  N = N_init;
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
profile("mu") {
  for(sex in 1:2) {
    for(farm in 1:nFarms) {
      for(stage in 1:nStagesSex[sex]) {
        mu[sex, farm, stage] = N[sex, farm, stage] ./ t_nFish_mx[farm, ];
      }
    }
  }
}

profile("priors") {
    // priors
    for(i in 1:nSurvCov) {
      lprior += std_normal_lpdf(surv_beta_z[i,]);
    }
    for(i in 1:2) {
      lprior += std_normal_lpdf(pMoltF_beta_z[i,]);
      lprior += std_normal_lpdf(pMoltM_beta_z[i,]);
    }
    lprior += normal_lpdf(nb_prec | 0, 2) - normal_lccdf(0 | 0, 2);
  }
}

model {
profile("likelihood") {
  for(sex in 1:2) {
    for(stage in 1:nStagesSex[sex]) {
      for(farm in 1:nFarms) {
        target += normal_lpdf(mu_true[stage, sex, sampledDays[farm], farm] |
                              mu[sex, farm, stage, sampledDays[farm]],
                              nb_prec);
      }
    }
  }
}
  target += lprior;
}

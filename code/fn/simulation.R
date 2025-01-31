# Simulation functions


simulate_farm_pops <- function(params, info, influx_df, farm_env, out_dir) {
  library(tidyverse)

  cat(format(now(), "%F %T"), "Initializing", out_dir, "\n")
  dir.create(out_dir, recursive=T)

  #---- transformed data
  # IP_mx[day, farm, sim]
  IP_mx <- array(influx_df$influx_100m, dim=c(info$nDays, info$nFarms, info$nSims))
  # attach_env_mx[day, farm, (1, RW, sal, uv, uv^2)]
  attach_env_mx <- array(1, dim=c(info$nDays, info$nFarms, 5))
  attach_env_mx[,,2] <- farm_env$RW_logit
  attach_env_mx[,,3] <- farm_env$salinity_z
  attach_env_mx[,,4] <- farm_env$uv
  attach_env_mx[,,5] <- farm_env$uv_sq
  # sal_mx[day, farm, (1, sal)]
  sal_mx <- array(1, dim=c(info$nDays, info$nFarms, 2))
  sal_mx[,,2] <- farm_env$salinity
  # temp_mx[day, farm, (1, temp)]
  temp_mx <- array(1, dim=c(info$nDays, info$nFarms, 2))
  temp_mx[,,2] <- farm_env$temperature
  # nFish_mx[day, farm]
  nFish_mx <- matrix(farm_env$nFish_est, nrow=info$nDays)
  # sampledDays[day, farm]
  sampledDays_mx <- matrix(as.numeric(farm_env$sampled), nrow=info$nDays)
  # nFishSampled[day, farm]
  nFishSampled_mx <- matrix(as.numeric(farm_env$nFishSampled), nrow=info$nDays)


  #---- Data structures
  #### Ensemble infection pressure: Number of copepodids in radius
  # ensIP[day, farm]
  ensIP <- array(0, dim=c(info$nDays, info$nFarms))

  #### Attachment: Probability of attachment per copepodid
  # pr_attach[day, farm]
  pr_attach <- array(0, dim=c(info$nDays, info$nFarms))

  #### For each daily cohort of newly attached copepodids by farm, track daily
  #### abundance, stage, and accumulated GDD for each farm, separating M/F
  # N[cohort, day, farm, M/F]
  cohort_N <- array(0, dim=c(info$nDays, info$nDays, info$nFarms, 2))
  # stage[cohort, day, farm, M/F]
  cohort_stage <- cohort_N
  # GDD[cohort, day, farm, M/F]
  cohort_GDD <- array(0, dim=c(info$nDays, info$nDays, info$nFarms))

  #### Daily survival rate on each farm for each stage
  # surv_rate[day, farm, stage]
  stage_survRate <- array(0, dim=c(info$nDays, info$nFarms, info$nStages+2))

  #### Total true and observed N by farm, stage, and day
  # N[day, farm, stage, M/F]  --- total abundance
  N <- array(0, dim=c(info$nDays, info$nFarms, info$nStages, 2))
  # mu[day, farm, stage, M/F]  --- latent mean per fish
  mu <- array(0, dim=c(info$nDays, info$nFarms, info$nStages, 2))
  # y[day, farm, stage, M/F]  --- sampled mean per fish
  y <- array(0, dim=c(info$nDays, info$nFarms, info$nStages, 2))

  cat(format(now(), "%F %T"), "Starting simulation \n")

  #---- Ensemble IP
  for(day in 1:info$nDays) {
    for(farm in 1:info$nFarms) {
      ensIP[day, farm] <- IP_mx[day, farm, ] %*% params$ensWts_p + params$IP_bg[farm]
    }
  }

  #---- Attachment
  for(day in 1:info$nDays) {
    for(farm in 1:info$nFarms) {
      pr_attach[day, farm] <- plogis(attach_env_mx[day, farm, ] %*%
                                       params$attach_beta)
    }
  }

  #---- Survival rates
  for(day in 1:info$nDays) {
    for(stage in 1:info$nStages) {
      stage_survRate[day, , stage] <- plogis(sal_mx[day,,] %*%
                                               params$surv_beta[stage,])
    }
  }

  #---- Cohort progression
  for(cohort in 1:info$nDays) {
    for(day in cohort:info$nDays) {
      if(cohort == day) {
        # initialize as attached copepodids, equal M/F distribution
        cohort_stage[cohort, day, , ] <- 1
        cohort_N[cohort, day, , ] <- ((((ensIP)^0.5 * pr_attach)[day,])^2)/2
        # assume GDD applies first day (attachment at 00:00:00)
        cohort_GDD[cohort, day, ] <- temp_mx[day, , 2]
      } else {
        # apply survival rate based on stage at start of day
        for(farm in 1:info$nFarms) {
          cohort_N[cohort, day, farm, ] <- cohort_N[cohort, day-1, farm, ] *
            stage_survRate[day-1, farm, cohort_stage[cohort, day-1, farm, ]]
        }
        # determine stage at END of day based on GDD accumulation
        cohort_stage[cohort, day, , 1] <- cohort_stage[cohort, day-1, , 1] +
          (cohort_GDD[cohort, day-1, ] > params$stageGDD_M[cohort_stage[cohort, day-1, , 1]])
        cohort_stage[cohort, day, , 2] <- cohort_stage[cohort, day-1, , 2] +
          (cohort_GDD[cohort, day-1, ] > params$stageGDD_F[cohort_stage[cohort, day-1, , 2]])
        # accumulate GDD
        cohort_GDD[cohort, day,] <- cohort_GDD[cohort, day-1, ] + temp_mx[day, , 2]
      }
    }
  }

  #---- Sum cohort abundances by stage for each day within each farm
  for(day in 1:info$nDays) {
    for(stage in 1:info$nStages) {
      for(sex in 1:2) {
        N[day, , stage, sex] <- apply(cohort_N[, day, , sex] *
                                        (cohort_stage[, day, , sex] == stage),
                                      2, sum)
      }
    }
  }

  #--- Calculate latent mean lice per fish for each day, farm, stage, and sex
  for(day in 1:info$nDays) {
    for(stage in 1:info$nStages) {
      for(sex in 1:2) {
        mu[day, , stage, sex] <- N[day, , stage, sex]/nFish_mx[day, ]
      }
    }
  }
  mu[is.infinite(mu)] <- 0  # when nFish == 0

  #--- Sample fish and calculate mean
  for(day in 1:info$nDays) {
    for(farm in 1:info$nFarms) {
      if(sampledDays_mx[day, farm] == 1) {
        for(stage in 1:info$nStages) {
          for(sex in 1:2) {
            y[day, farm, stage, sex] <- rnbinom(
              n=nFishSampled_mx[day, farm],
              mu=mu[day, farm, stage, sex] * params$detect_p[stage],
              size=params$nb_prec) |>
              mean()
          }
        }
      }
    }
  }

  out_df <- farm_env |>
    select(date, sepaSite, sampled) |>
    mutate(mu_chal=c(apply(mu[,,1,], 1:2, sum)),
           mu_prea=c(apply(mu[,,2,], 1:2, sum)),
           mu_af_ng=c(mu[,,3,2]),
           mu_af_gr=c(mu[,,4,2]),
           mu_af=c(apply(mu[,,3:4,2], 1:2, sum)),
           mu_am=c(apply(mu[,,3:4,1], 1:2, sum)),
           y_chal=c(apply(y[,,1,], 1:2, sum)),
           y_prea=c(apply(y[,,2,], 1:2, sum)),
           y_af_ng=c(y[,,3,2]),
           y_af_gr=c(y[,,4,2]),
           y_af=c(apply(y[,,3:4,2], 1:2, sum)),
           y_am=c(apply(y[,,3:4,1], 1:2, sum)))

  cat(format(now(), "%F %T"), "Writing simulation to", out_dir, "\n")

  cat(format(now(), "%F %T"), "  Storing inputs \n")
  saveRDS(params, glue("{out_dir}/params.rds"))
  saveRDS(IP_mx, glue("{out_dir}/IP_mx.rds"))
  saveRDS(attach_env_mx, glue("{out_dir}/attach_env_mx.rds"))
  saveRDS(IP_mx, glue("{out_dir}/IP_mx.rds"))
  saveRDS(sal_mx, glue("{out_dir}/sal_mx.rds"))
  saveRDS(temp_mx, glue("{out_dir}/temp_mx.rds"))
  saveRDS(nFish_mx, glue("{out_dir}/nFish_mx.rds"))

  cat(format(now(), "%F %T"), "  Storing calculated variables  \n")
  saveRDS(ensIP, glue("{out_dir}/ensIP.rds"))
  saveRDS(pr_attach, glue("{out_dir}/pr_attach.rds"))
  saveRDS(stage_survRate, glue("{out_dir}/stage_survRate.rds"))

  cat(format(now(), "%F %T"), "  Storing cohort structures  \n")
  saveRDS(cohort_N, glue("{out_dir}/cohort_N.rds"))
  saveRDS(cohort_stage, glue("{out_dir}/cohort_stage.rds"))
  saveRDS(cohort_GDD, glue("{out_dir}/cohort_GDD.rds"))

  cat(format(now(), "%F %T"), "  Storing lice structures  \n")
  saveRDS(N, glue("{out_dir}/N.rds"))
  saveRDS(mu, glue("{out_dir}/mu.rds"))
  saveRDS(y, glue("{out_dir}/y.rds"))
  write_csv(out_df, glue("{out_dir}/out_df.csv"))

  return(out_df)

}

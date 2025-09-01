# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Simulation functions


simulate_farm_pops_mn_lpf <- function(params, info, influx_df, farm_env, out_dir) {
  # Attached copepodids are immediately translated to 'mean per fish'
  library(tidyverse)

  cat(format(now(), "%F %T"), "Initializing", out_dir, "\n")
  dir.create(out_dir, recursive=T)
  cat(format(now(), "%F %T"), "  Creating and storing inputs \n")
  saveRDS(info, glue("{out_dir}/info.rds"))
  saveRDS(params, glue("{out_dir}/params.rds"))


  #---- transformed data
  # IP_mx[farm, day, sim]
  IP_mx <- make_IP_mx(influx_df, info, out_dir)

  # attach_env_mx[farm, day, (1, RW, sal, uv, uv^2)]
  attach_env_mx <- make_attach_env_mx(farm_env, info, params, out_dir)

  # sal_mx[day, farm, (1, sal)]
  sal_mx <- make_sal_mx(farm_env, info, params, out_dir)

  # temp_mx[day, farm, (1, temp)]
  temp_mx <- make_temp_mx(farm_env, info, out_dir)
  temp_z_mx <- make_temp_z_mx(farm_env, info, out_dir)

  # nFish_mx[day, farm]
  nFish_mx <- make_nFish_mx(farm_env, info, out_dir)

  # sampledDays[sample, c(farm, day)] -- database-like structure for Stan
  sampledDays <- make_sampledDays(farm_env, out_dir)

  # nFishSampled[day, farm]
  nFishSampled_mx <- make_nFishSampled_mx(farm_env, info, nFish_mx, out_dir)

  # treatDays_mx[day, farm] -- initialize as 0s
  treatDays_mx <- matrix(0, nrow=info$nDays, ncol=info$nFarms, out_dir)


  #---- Data structures
  #### Ensemble infection pressure: Number of copepodids in radius
  # ensIP[day, farm]
  ensIP <- array(0, dim=c(info$nDays, info$nFarms))

  #### Attachment: Probability of attachment per copepodid
  # pr_attach[day, farm]
  pr_attach <- array(0, dim=c(info$nDays, info$nFarms))

  #### Number of attached copepodids
  # N_attach[day, farm]
  N_attach <- array(0, dim=c(info$nDays, info$nFarms))
  y_attach <- N_attach

  #### For each daily cohort of newly attached copepodids by farm, track daily
  #### abundance, stage, and accumulated GDD for each farm, separating M/F
  # N[cohort, day, farm, M/F]
  cohort_N <- array(0, dim=c(info$nDays, info$nDays, info$nFarms, 2))
  # stage[cohort, day, farm, M/F]
  cohort_stage <- cohort_N
  # GDD[cohort, day, farm]
  cohort_GDD <- array(0, dim=c(info$nDays, info$nDays, info$nFarms))

  #### Daily survival rate on each farm for each stage
  # surv_rate[farm, day, stage]
  stage_survRate <- array(0, dim=c(info$nFarms, info$nDays, info$nStages))

  #### Total true and observed N by farm, stage, and day
  # mu[stage, M/F, day, farm]  --- latent mean per fish
  mu <- array(0, dim=c(info$nStages, 2, info$nDays, info$nFarms))
  # y_bar[stage, M/F, day, farm]  --- expected mean observed lice per fish
  y_bar <- array(0, dim=c(info$nStages, 2, info$nDays, info$nFarms))
  # y[stage, M/F, day, farm]  --- sampled mean per fish
  y <- array(0, dim=c(info$nStages, 2, info$nDays, info$nFarms))

  cat(format(now(), "%F %T"), "Starting simulation \n")

  #---- Ensemble IP, attachement, survival rates
  for(farm in 1:info$nFarms) {
    ensIP[,farm] <- IP_mx[farm,,] %*% params$ensWts_p + params$IP_bg*info$nPens[farm]
    if(length(params$attach_beta) > 1) {
      pr_attach[,farm] <- plogis(attach_env_mx[farm,,] %*% params$attach_beta)
    } else {
      pr_attach[,farm] <- plogis(attach_env_mx[farm,,] * params$attach_beta)
    }
    stage_survRate[farm,,] <- plogis(sal_mx[farm,,] %*% params$surv_beta)
  }
  # assume even sex ratio
  N_attach <- (ensIP^(1/params$IP_scale) * pr_attach)^params$IP_scale / nFish_mx / 2
  N_attach[nFish_mx==0] <- 0
  for(day in 1:info$nDays) {
    for(farm in 1:info$nFarms) {
      y_attach[day, farm] <- rnbinom(1, mu=N_attach[day, farm], size=params$nb_prec)
    }
  }

  #---- Calculate GDD by cohort through simulation period
  for(cohort in 1:info$nDays) {
    for(day in cohort:info$nDays) {
      if(day == 1) {
        cohort_GDD[cohort, day, ] <-  temp_mx[day, ]
      } else {
        cohort_GDD[cohort, day, ] <- cohort_GDD[cohort, day-1, ] + temp_mx[day, ]
      }
    }
  }

  #---- Stage calculation
  for(cohort in 1:info$nDays) {
    for(farm in 1:info$nFarms) {
      for(sex in 1:2) {
        cohort_stage[cohort, 1:cohort, farm, sex] <- 0
        cohort_stage[cohort, cohort:info$nDays, farm, sex] <- 1
        for(day in cohort:info$nDays) {
          if(nFish_mx[day, farm] == 0) {
            cohort_stage[cohort, day:info$nDays, farm, sex] <- 0
          }
          if(cohort_stage[cohort, day, farm, sex] > 0) {
            if(cohort_stage[cohort, day, farm, sex] < info$nStages) {
              if(cohort_GDD[cohort, day, farm] > params$thresh_GDD[cohort_stage[cohort, day, farm, sex],sex]) {
                cohort_stage[cohort, day:info$nDays, farm, sex] <- cohort_stage[cohort, day, farm, sex] + 1
              }
            }
            if(cohort_GDD[cohort, day, farm] > params$lifespan) {
              cohort_stage[cohort, day:info$nDays, farm, sex] <- 0
            }
          }
        }
      }
    }
  }

  #---- Daily abundance progression
  for(day in 1:info$nDays) {
    for(cohort in 1:day) {
      if(cohort == day) {
        # initialize cohort
        cohort_N[cohort, day, , ] <- N_attach[day,]
      } else {
        # update cohorts
        for(farm in 1:info$nFarms) {
          # update abundance
          if(cohort_stage[cohort, day-1, farm, sex] == 0) {
            cohort_N[cohort, day, farm, ] <- 0
          } else {
            cohort_N[cohort, day, farm, ] <- cohort_N[cohort, day-1, farm, ] *
              stage_survRate[farm, day-1, cohort_stage[cohort, day-1, farm,]]
          }
        }
      }
      #---- Sum cohort abundances by stage for each day within each farm
      for(farm in 1:info$nFarms) {
        for(sex in 1:2) {
          mu[cohort_stage[cohort, day, farm, sex], sex, day, farm] <-
            mu[cohort_stage[cohort, day, farm, sex], sex, day, farm] +
            cohort_N[cohort, day, farm, sex]
        }
      }
    }
    #--- Calculate latent mean lice per fish for each day, farm, stage, and sex
    for(stage in 1:info$nStages) {
      for(sex in 1:2) {
        y_bar[stage, sex, day, ] <- mu[stage, sex, day, ] * nFishSampled_mx[day,] * params$detect_p[stage]
      }
    }
    #--- Sample fish and calculate mean
    for(farm in 1:info$nFarms) {
      if(day %in% sampledDays[sampledDays[,1]==farm,2]) {
        for(sex in 1:2) {
          y[, sex, day, farm] <- rnbinom(info$nStages,
                                       mu=y_bar[, sex, day, farm],
                                       size=params$nb_prec)
        }
        #--- Treat if AF > threshold
        if((y[info$nStages, 1, day, farm]/nFishSampled_mx[day, farm] > params$treat_thresh) &
            day < info$nDays)  {
          treatDays_mx[day+1, farm] = 1
          stage_survRate[farm, day,] <- stage_survRate[farm, day,] * (1 - params$treat_efficacy)
        }
      }
    }
  }

  out_df <- farm_env |>
    select(date, sepaSite, pen, sampled, nFishSampled) |>
    arrange(day, sepaSite, pen) |>
    mutate(treat=c(treatDays_mx)) |>
    mutate(mu_chal=c(apply(mu[1,,,], 2:3, sum)),
           mu_prea=c(apply(mu[2,,,], 2:3, sum)),
           mu_af=c(mu[3,1,,]),
           mu_am=c(mu[3,2,,]),
           ybar_chal=c(apply(y_bar[1,,,], 2:3, sum)),
           ybar_prea=c(apply(y_bar[2,,,], 2:3, sum)),
           ybar_af=c(y_bar[3,1,,]),
           ybar_am=c(y_bar[3,2,,]),
           y_chal=c(apply(y[1,,,], 2:3, sum)),
           y_prea=c(apply(y[2,,,], 2:3, sum)),
           y_af=c(y[3,1,,]),
           y_am=c(y[3,2,,])) |>
    mutate(across(starts_with("y"), ~.x/nFishSampled))


  cat(format(now(), "%F %T"), "Writing simulation to", out_dir, "\n")

  cat(format(now(), "%F %T"), "  Storing calculated variables  \n")
  saveRDS(treatDays_mx, glue("{out_dir}/treatDays_mx.rds"))
  saveRDS(ensIP, glue("{out_dir}/ensIP.rds"))
  saveRDS(pr_attach, glue("{out_dir}/pr_attach.rds"))
  saveRDS(stage_survRate, glue("{out_dir}/stage_survRate.rds"))

  cat(format(now(), "%F %T"), "  Storing cohort structures  \n")
  saveRDS(cohort_N, glue("{out_dir}/cohort_N.rds"))
  saveRDS(cohort_stage, glue("{out_dir}/cohort_stage.rds"))
  saveRDS(cohort_GDD, glue("{out_dir}/cohort_GDD.rds"))

  cat(format(now(), "%F %T"), "  Storing lice structures  \n")
  saveRDS(N_attach, glue("{out_dir}/N_attach.rds"))
  saveRDS(y_attach, glue("{out_dir}/y_attach.rds"))
  saveRDS(mu, glue("{out_dir}/mu.rds"))
  saveRDS(y_bar, glue("{out_dir}/y_bar.rds"))
  saveRDS(y, glue("{out_dir}/y.rds"))
  write_csv(out_df, glue("{out_dir}/out_df.csv"))

  return(out_df)

}




make_stan_data <- function(dat_dir, source="sim", GQ_ypred=TRUE, GQ_start=NULL, priors_only=FALSE, prior_ls=NULL) {
  library(tidyverse)
  library(glue)

  info <- readRDS(glue("{dat_dir}info.rds"))
  params <- readRDS(glue("{dat_dir}params.rds"))
  if(is.null(GQ_start)) {
    dates <- 1:info$nDays
    hours <- 1:info$nHours
    dates_GQ <- 1
    hours_GQ <- 1
  } else {
    info$nDays <- as.numeric(ymd(GQ_start) - info$dateRange[1])
    info$nDays_GQ <- as.numeric(info$dateRange[2] - ymd(GQ_start))
    info$nHours <- info$nDays*24
    info$nHours_GQ <- info$nDays_GQ*24
    dates <- 1:info$nDays
    hours <- 1:info$nHours
    dates_GQ <- info$nDays + (1:info$nDays_GQ)
    hours_GQ <- info$nHours + (1:info$nHours_GQ)
  }

  stan_dat <- list(
    GQ_ypred=as.numeric(GQ_ypred),
    GQ_new=as.numeric(!is.null(GQ_start)),
    nDays=info$nDays,
    nHours=info$nHours,
    nFarms=info$nFarms,
    nSims=info$nSims,
    nStages=info$nStages,
    nStageGroups=info$nStageGroups,
    stg_grp_ii=info$stg_grp_ii,
    nPens=info$nPens,
    nTrtMethods=info$nTrtMethods,
    nTrtTypes=info$nTrtTypes,
    trt_meth_ii=info$trt_meth_ii,
    nAttachCov=length(params$attach_beta),
    nSurvCov=nrow(params$surv_beta),
    IP_volume=info$IP_penVolume,
    y_bar_minimum=1e-5,
    day_hour=readRDS(glue("{dat_dir}day_hour.rds"))[dates,],
    # farm data, treatments, sampling info
    nFish_mx=readRDS(glue("{dat_dir}nFish_mx.rds"))[dates,],
    trtApplied=readRDS(glue("{dat_dir}trtApplied_mx.rds"))[,dates,],
    sample_i=readRDS(glue("{dat_dir}sampledDays.rds")) |> as_tibble() |> filter(day %in% dates) |> as.matrix(),
    nFishSampled_mx=t(readRDS(glue("{dat_dir}nFishSampled_mx.rds"))[dates,]),
    # IP from biotracker
    IP_mx=readRDS(glue("{dat_dir}IP_mx.rds"))[,hours,],
    # farm environment
    attach_env_mx=readRDS(glue("{dat_dir}attach_env_mx.rds"))[,hours,],
    surv_env_mx=readRDS(glue("{dat_dir}sal_mx.rds"))[,dates,],
    temp_mx=readRDS(glue("{dat_dir}temp_mx.rds"))[dates,],
    temp_z_mx=readRDS(glue("{dat_dir}temp_z_mx.rds"))[dates,],
    ydayh_mx=readRDS(glue("{dat_dir}ydayh_mx.rds"))[hours,],
    # priors
    sample_prior_only=as.numeric(priors_only),
    # attach_beta: [c(RW, Sal, Temp, UV, UV^2), c(mu, sigma)]; normal (logit scale)
    prior_attach_beta=cbind(c(1, rep(0.25, length(params$attach_beta)-2), 0),
                            c(rep(0.25, length(params$attach_beta)-1), 0.25)),
    # treatment: trtEff[global] ~ N(mu, sigma) (logit scale)
    prior_logit_trtEff_global=c(0, 1),
    # treatment: trtEff[method] ~ N(trtEff[global], trtEff_sd_methods)  (logit scale)  student_t(nu, mu, sd)
    prior_trtEff_sd_methods=c(3, 0, 0.75),
    # treatment: trtEff[type] ~ N(trtEff[method], trtEff_sd_types)  (logit scale)  student_t(nu, mu, sd)
    prior_trtEff_sd_types=c(3, 0, 0.75),
    # surv_beta: [c(Int, Temp), c(Ch1/Ch2, PA1/PA2, Ad), c(mu, sigma)]; normal (logit scale)
    prior_surv_beta=array(c(rep(c(4, rep(0.2, nrow(params$surv_beta)-1)), info$nStageGroups),
                            rep(c(1, rep(0.1, nrow(params$surv_beta)-1)), info$nStageGroups)),
                          dim=c(nrow(params$surv_beta), info$nStageGroups, 2),
                          dimnames=list(c("int", "temp"),
                                        c("Ch1,2", "PA1,2", "Ad"),
                                        c("mu", "sd"))),
    # surv beta: sd for farm-level intercepts: student_t(nu, mu, sd)
    prior_surv_int_farm_sd=c(3, 0, 0.75),
    # mnDaysStage: [c(Int, Temp), c(Ch1/Ch2, PA1/PA2), c(mu, sd)]; normal
    prior_mnDaysStage_F=array(c(params$mnDaysStageCh[,1]/2, params$mnDaysStagePA[,1]/2,
                                params$mnDaysStageCh[,2], params$mnDaysStagePA[,2]),
                              dim=c(2, info$nStageGroups-1, 2),
                              dimnames=list(c("int", "temp"),
                                            c("Ch1,2", "PA1,2"),
                                            c("mu", "sd"))),
    # logit_detect_p: [c(Ch, Pr), c(mu, sd)]
    prior_logit_detect_p=cbind(c(-1, 1),
                               c(0.5, 0.5)),
    # IP_bg_m3: c(mu, sigma); normal, T(0, )
    prior_IP_bg_m3=c(0.05, 0.05),
    # inv_sqrt_nb_prec: df; normal(mu, sd); nb_prec = 1/inv_sqrt_nb_prec^2
    prior_inv_sqrt_nb_prec=c(0, 1)
  )
  # reformat sample info for Stan
  stan_dat$nSamples <- nrow(stan_dat$sample_i)
  stan_dat$sample_ii <- make_sample_ii(stan_dat$sample_i, info$nFarms)
  # update any priors specified in prior_ls argument
  if(!is.null(prior_ls)) {
    for(i in 1:length(prior_ls)) {
      stan_dat[[names(prior_ls)[i]]] <- prior_ls[[i]]
    }
  }
  if(source=="sim") {
    # add male Ch/PA, re-divide assuming a 50:50 ratio
    stan_dat$y <- readRDS(glue("{dat_dir}y_obs.rds"))[,,dates,]
    stan_dat$y_F <- stan_dat$y[,1,,]
    stan_dat$y_F[1,,] <- round((stan_dat$y_F[1,,] + stan_dat$y[1,2,,])/2)
    stan_dat$y_F[2,,] <- round((stan_dat$y_F[2,,] + stan_dat$y[2,2,,])/2)
  } else {
    stan_dat$y <- readRDS(glue("{dat_dir}y_obs.rds"))[,dates,]
    stan_dat$y_F <- stan_dat$y
    stan_dat$y_F[1,,] <- round(stan_dat$y_F[1,,]/2)
    stan_dat$y_F[2,,] <- round(stan_dat$y_F[2,,]/2)
  }
  if(is.null(GQ_start)) {
    stan_dat <- c(
      stan_dat,
      list(nDays_GQ=1,
           nHours_GQ=1,
           nSamples_GQ=1,
           day_hour_GQ=matrix(1, nrow=1, ncol=24),
           IP_mx_GQ=array(0, dim=c(stan_dat$nFarms, 1, stan_dat$nSims)),
           attach_env_mx_GQ=array(0, dim=c(stan_dat$nFarms, 1, stan_dat$nAttachCov)),
           surv_env_mx_GQ=array(0, dim=c(stan_dat$nFarms, 1, stan_dat$nSurvCov)),
           temp_z_mx_GQ=matrix(0, nrow=1, ncol=stan_dat$nFarms),
           ydayh_mx_GQ=matrix(0, nrow=1, ncol=3),
           nFish_mx_GQ=matrix(0, nrow=1, ncol=stan_dat$nFarms),
           trtApplied_GQ=array(0, dim=c(stan_dat$nFarms, 1, stan_dat$nTrtTypes)),
           sample_i_GQ=matrix(0, nrow=1, ncol=2),
           sample_ii_GQ=matrix(0, nrow=stan_dat$nFarms, ncol=2),
           nFishSampled_mx_GQ=matrix(0, nrow=stan_dat$nFarms, ncol=1),
           y_F_GQ=array(0, dim=c(stan_dat$nStageGroups, stan_dat$nDays_GQ, stan_dat$nFarms))
      ))
  } else {
    stan_dat <- c(
      stan_dat,
      list(nDays_GQ=info$nDays_GQ,
           nHours_GQ=info$nHours_GQ,
           day_hour_GQ=readRDS(glue("{dat_dir}day_hour.rds"))[dates_GQ,] - info$nHours,
           IP_mx_GQ=readRDS(glue("{dat_dir}IP_mx.rds"))[,hours_GQ,],
           attach_env_mx_GQ=readRDS(glue("{dat_dir}attach_env_mx.rds"))[,hours_GQ,],
           surv_env_mx_GQ=readRDS(glue("{dat_dir}sal_mx.rds"))[,dates_GQ,],
           temp_z_mx_GQ=readRDS(glue("{dat_dir}temp_z_mx.rds"))[dates_GQ,],
           ydayh_mx_GQ=readRDS(glue("{dat_dir}ydayh_mx.rds"))[hours_GQ,],
           nFish_mx_GQ=readRDS(glue("{dat_dir}nFish_mx.rds"))[dates_GQ,],
           trtApplied_GQ=readRDS(glue("{dat_dir}trtApplied_mx.rds"))[,dates_GQ,],
           sample_i_GQ=readRDS(glue("{dat_dir}sampledDays.rds")) |>
             as_tibble() |>
             filter(day %in% dates_GQ) |>
             mutate(day=day - info$nDays) |>
             as.matrix(),
           nFishSampled_mx_GQ=t(readRDS(glue("{dat_dir}nFishSampled_mx.rds"))[dates_GQ,])
      ))
    stan_dat$nSamples_GQ <- nrow(stan_dat$sample_i_GQ)
    stan_dat$sample_ii_GQ <- make_sample_ii(stan_dat$sample_i_GQ, info$nFarms)
    if(source=="sim") {
      # add male Ch/PA, re-divide assuming a 50:50 ratio
      stan_dat$y_GQ <- readRDS(glue("{dat_dir}y_obs.rds"))[,,dates_GQ,]
      stan_dat$y_F_GQ <- stan_dat$y_GQ[,1,,]
      stan_dat$y_F_GQ[1,,] <- round((stan_dat$y_F_GQ[1,,] + stan_dat$y_GQ[1,2,,])/2)
      stan_dat$y_F_GQ[2,,] <- round((stan_dat$y_F_GQ[2,,] + stan_dat$y_GQ[2,2,,])/2)
    } else {
      stan_dat$y_GQ <- readRDS(glue("{dat_dir}y_obs.rds"))[,dates_GQ,]
      stan_dat$y_F_GQ <- stan_dat$y_GQ
      stan_dat$y_F_GQ[1,,] <- round(stan_dat$y_F_GQ[1,,]/2)
      stan_dat$y_F_GQ[2,,] <- round(stan_dat$y_F_GQ[2,,]/2)
    }
  }

  return(list(dat=stan_dat, params=params))
}




make_sample_ii <- function(sample_i, nFarms) {
  sample_i |>
    as_tibble() |>
    mutate(index=row_number()) |>
    group_by(sepaSite) |>
    summarise(start=min(index),
              end=max(index)) |>
    full_join(tibble(sepaSite=1:nFarms), by=join_by(sepaSite)) |>
    replace_na(list(start=0, end=0)) |>
    arrange(sepaSite) |>
    select(-sepaSite) |>
    as.matrix()
}



make_IP_mx <- function(influx_df, info, out_dir=NULL) {
  IP_mx <- arrange(influx_df, sim, day, sepaSite, pen)$influx |>
    array(dim=c(info$nFarms, info$nHours, info$nSims))
  if(!is.null(out_dir)) {
    saveRDS(IP_mx, glue("{out_dir}/IP_mx.rds"))
  }
  return(IP_mx)
}



make_attach_env_mx <- function(farm_env, info, params, out_dir=NULL) {
  if("pen" %notin% names(farm_env)) {
    farm_env$pen <- "a"
  }
  farm_env <- farm_env |> arrange(time, sepaSite, pen)
  attach_env_mx <- array(1, dim=c(info$nFarms, info$nHours, 5))
  attach_env_mx[,,1] <- farm_env$RW_logit
  attach_env_mx[,,2] <- farm_env$salinity_z
  attach_env_mx[,,3] <- farm_env$temperature_z
  attach_env_mx[,,4] <- farm_env$uv_z
  attach_env_mx[,,5] <- farm_env$uv_z_sq
  attach_env_mx <- attach_env_mx[,,1:length(params$attach_beta), drop=F]
  if(!is.null(out_dir)) {
    saveRDS(attach_env_mx, glue("{out_dir}/attach_env_mx.rds"))
  }
  return(attach_env_mx)
}

make_sal_mx <- function(farm_env, info, params, out_dir=NULL) {
  farm_env <- farm_env |> arrange(date, sepaSite)
  sal_mx <- array(1, dim=c(info$nFarms, info$nDays, 2))
  sal_mx[,,2] <- farm_env$salinity_m30
  sal_mx <- sal_mx[,,1:nrow(params$surv_beta), drop=F]
  if(!is.null(out_dir)) {
    saveRDS(sal_mx, glue("{out_dir}/sal_mx.rds"))
  }
  return(sal_mx)
}

make_trtApplied_mx <- function(trt_df, info, out_dir=NULL) {
  trtApp_mx <- array(0, dim=c(info$nFarms, info$nDays, info$nTrtTypes))
  for(trt in 1:info$nTrtTypes) {
    trtApp_mx[,,trt] <- trt_df[[paste0("t_", trt)]]
  }
  if(!is.null(out_dir)) {
    saveRDS(trtApp_mx, glue("{out_dir}/trtApp_mx.rds"))
  }
  return(trtApp_mx)
}

make_temp_mx <- function(farm_env, info, out_dir=NULL) {
  farm_env <- farm_env |> arrange(date, sepaSite)
  temp_mx <- matrix(farm_env$temperature, nrow=info$nDays, byrow=T)
  if(!is.null(out_dir)) {
    saveRDS(temp_mx, glue("{out_dir}/temp_mx.rds"))
  }
  return(temp_mx)
}


make_temp_z_mx <- function(farm_env, info, out_dir=NULL) {
  farm_env <- farm_env |> arrange(date, sepaSite)
  temp_z_mx <- matrix(farm_env$temperature_z, nrow=info$nDays, byrow=T)
  if(!is.null(out_dir)) {
    saveRDS(temp_z_mx, glue("{out_dir}/temp_z_mx.rds"))
  }
  return(temp_z_mx)
}


make_nFish_mx <- function(farm_env, info, out_dir=NULL) {
  farm_env <- farm_env |> arrange(date, sepaSite)
  nFish_mx <- matrix(farm_env$nFish_est, nrow=info$nDays, byrow=T)
  if(!is.null(out_dir)) {
    saveRDS(nFish_mx, glue("{out_dir}/nFish_mx.rds"))
  }
  return(nFish_mx)
}

make_sampledDays <- function(farm_env, out_dir=NULL) {
  if("pen" %notin% names(farm_env)) {
    farm_env$pen <- "a"
  }
  sampledDays <- farm_env |>
    arrange(sepaSite, pen, day) |>
    filter(sampled) |>
    filter(day > 28) |> # 4 week burn-in
    select(sepaSite, day) |>
    mutate(sepaSite=as.numeric(sepaSite)) |>
    as.matrix()
  if(!is.null(out_dir)) {
    saveRDS(sampledDays, glue("{out_dir}/sampledDays.rds"))
  }
  return(sampledDays)
}

make_nFishSampled_mx <- function(farm_env, info, nFish_mx, out_dir=NULL) {
  farm_env <- farm_env |> arrange(date, sepaSite)
  nFishSampled_mx <- matrix(as.numeric(farm_env$nFishSampled), nrow=info$nDays, byrow=T) *
    (nFish_mx > 0)
  if(!is.null(out_dir)) {
    saveRDS(nFishSampled_mx, glue("{out_dir}/nFishSampled_mx.rds"))
  }
  return(nFishSampled_mx)
}


make_ydayh_mx <- function(farm_env=NULL, out_dir=NULL) {
  n_h <- 366*24
  if(is.null(farm_env)) {
    hour_vec <- 1:n_h
  } else {
    hour_vec <- 1:max(farm_env$elapsedHours)
  }
  ydayh_mx <- cbind(1,
                   cos(hour_vec/n_h*2*pi),
                   sin(hour_vec/n_h*2*pi))
  if(!is.null(out_dir)) {
    saveRDS(ydayh_mx, glue("{out_dir}/ydayh_mx.rds"))
  }
  return(ydayh_mx)
}

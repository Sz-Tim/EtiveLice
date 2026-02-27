make_stan_data <- function(dat_dir, source="sim", GQ_ypred=TRUE, GQ_new=FALSE, priors_only=FALSE, prior_ls=NULL) {
  library(tidyverse)
  library(glue)

  info <- readRDS(glue("{dat_dir}info.rds"))
  params <- readRDS(glue("{dat_dir}params.rds"))
  dates <- 1:info$nDays
  hours <- 1:info$nHours

  stan_dat <- list(
    GQ_ypred = as.numeric(GQ_ypred),
    GQ_new = as.numeric(GQ_new),
    nDays = info$nDays,
    nHours = info$nHours,
    nFarms = info$nFarms,
    nSims = info$nSims,
    nStages = info$nStages,
    nPens = info$nPens,
    nAttachCov = length(params$attach_beta),
    nSurvCov = nrow(params$surv_beta),
    IP_volume = info$IP_penVolume,
    y_bar_minimum = 1e-5,
    day_hour = readRDS(glue("{dat_dir}day_hour.rds"))[dates,],
    # farm counts, treatments, sampling info
    y = readRDS(glue("{dat_dir}y.rds"))[,,dates,],
    sample_i = readRDS(glue("{dat_dir}sampledDays.rds")),
    nFishSampled_mx = readRDS(glue("{dat_dir}nFishSampled_mx.rds"))[dates,],
    nFish_mx = readRDS(glue("{dat_dir}nFish_mx.rds"))[dates,],
    treatDays = readRDS(glue("{dat_dir}treatDays_mx.rds"))[dates,],
    # IP from biotracker
    IP_mx = readRDS(glue("{dat_dir}IP_mx.rds"))[,hours,],
    # farm environment
    attach_env_mx = readRDS(glue("{dat_dir}attach_env_mx.rds"))[,hours,],
    surv_env_mx = readRDS(glue("{dat_dir}sal_mx.rds"))[,dates,],
    temp_mx = readRDS(glue("{dat_dir}temp_mx.rds"))[dates,],
    temp_z_mx = readRDS(glue("{dat_dir}temp_z_mx.rds"))[dates,],
    # priors
    sample_prior_only = as.numeric(priors_only),
    # attach_beta: [c(RW, Sal, Temp, UV, UV^2), c(mu, sigma)]; normal (logit scale)
    # prior_attach_beta = cbind(c(1, 0.25, 0.25, 0, 0),
    #                           c(0.5, 0.5, 0.5, 0.5, 0.5)),
    prior_attach_beta = cbind(c(1, rep(0.25, length(params$attach_beta)-2), 0),
                              c(rep(0.25, length(params$attach_beta)-1), 0.25)),
    # surv_beta: [c(Int, Temp), c(Ch, Pr, Ad), c(mu, sigma)]; normal (logit scale)
    prior_surv_beta = array(c(rep(c(4, rep(0.2, nrow(params$surv_beta)-1)), info$nStages),
                              rep(c(1, rep(0.1, nrow(params$surv_beta)-1)), info$nStages)),
                            dim=c(nrow(params$surv_beta), info$nStages, 2),
                            dimnames=list(c("int", "temp"),
                                          c("Ch", "Pr", "Ad"),
                                          c("mu", "sd"))),
    # surv beta: sd for farm-level intercepts: student_t(nu, mu, sd)
    prior_surv_int_farm_sd = c(3, 0, 0.75),
    # mnDaysStage: [c(Int, Temp), c(Ch-Pr, Pr-Ad), c(mu, sd)]; normal
    prior_mnDaysStage_F = array(c(params$mnDaysStageCh[,1], params$mnDaysStagePA[,1],
                                  params$mnDaysStageCh[,2], params$mnDaysStagePA[,2]),
                                dim=c(2, info$nStages-1, 2),
                                dimnames=list(c("int", "temp"),
                                              c("Pr", "Ad"),
                                              c("mu", "sd"))),
    # logit_detect_p: [c(Ch, Pr), c(mu, sd)]
    prior_logit_detect_p = cbind(c(-1, 1),
                                 c(0.5, 0.5)),
    # IP_bg_m3: c(mu, sigma); normal, T(0, )
    prior_IP_bg_m3 = c(0.05, 0.05),
    # inv_sqrt_nb_prec: df; normal(mu, sd); nb_prec = 1/inv_sqrt_nb_prec^2
    prior_inv_sqrt_nb_prec = c(0, 1),
    # treatEfficacy: beta(alpha, beta)
    prior_treatEfficacy = c(2, 2),
    # IP_halfStat_m3: c(nu, mu, sigma); student_t, T(0, )
    prior_IP_halfSat_m3 = c(3, 20, 10)
  )
  # reformat sample info for Stan
  stan_dat$nSamples <- nrow(stan_dat$sample_i)
  stan_dat$sample_ii <- stan_dat$sample_i |>
    as_tibble() |>
    mutate(index=row_number()) |>
    group_by(sepaSite) |>
    summarise(start=min(index),
              end=max(index)) |>
    select(-sepaSite) |>
    as.matrix()
  # update any priors specified in prior_ls argument
  if(!is.null(prior_ls)) {
    for(i in 1:length(prior_ls)) {
      stan_dat[[names(prior_ls)[i]]] <- prior_ls[[i]]
    }
  }
  if(source=="sim") {
    # add male Ch/PA, re-divide assuming a 50:50 ratio
    stan_dat$y_F <- stan_dat$y[,1,,]
    stan_dat$y_F[1,,] <- round((stan_dat$y_F[1,,] + stan_dat$y[1,2,,])/2)
    stan_dat$y_F[2,,] <- round((stan_dat$y_F[2,,] + stan_dat$y[2,2,,])/2)
  }

  return(list(dat=stan_dat, params=params))
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
  sal_mx <- array(1, dim=c(info$nFarms, info$nDays, 2))
  sal_mx[,,2] <- farm_env$salinity_m30
  sal_mx <- sal_mx[,,1:nrow(params$surv_beta), drop=F]
  if(!is.null(out_dir)) {
    saveRDS(sal_mx, glue("{out_dir}/sal_mx.rds"))
  }
  return(sal_mx)
}

make_temp_mx <- function(farm_env, info, out_dir=NULL) {
  temp_mx <- matrix(farm_env$temperature, nrow=info$nDays, byrow=T)
  if(!is.null(out_dir)) {
    saveRDS(temp_mx, glue("{out_dir}/temp_mx.rds"))
  }
  return(temp_mx)
}


make_temp_z_mx <- function(farm_env, info, out_dir=NULL) {
  temp_z_mx <- matrix(farm_env$temperature_z, nrow=info$nDays, byrow=T)
  if(!is.null(out_dir)) {
    saveRDS(temp_z_mx, glue("{out_dir}/temp_z_mx.rds"))
  }
  return(temp_z_mx)
}


make_nFish_mx <- function(farm_env, info, out_dir=NULL) {
  nFish_mx <- matrix(farm_env$nFish_est, nrow=info$nDays, byrow=T)
  if(!is.null(out_dir)) {
    saveRDS(nFish_mx, glue("{out_dir}/nFish_mx.rds"))
  }
  return(nFish_mx)
}

make_sampledDays <- function(farm_env, out_dir=NULL) {
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
  nFishSampled_mx <- matrix(as.numeric(farm_env$nFishSampled), nrow=info$nDays, byrow=T) *
    (nFish_mx > 0)
  if(!is.null(out_dir)) {
    saveRDS(nFishSampled_mx, glue("{out_dir}/nFishSampled_mx.rds"))
  }
  return(nFishSampled_mx)
}



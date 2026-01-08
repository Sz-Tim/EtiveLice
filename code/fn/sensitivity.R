


sample_parameter_distributions <- function(n_sim=30, out_dir,
                                           mode="random",
                                           egg_post=NULL, mort_post=NULL, sink_post=NULL,
                                           swim_Sigma=diag(1, nrow=4, ncol=4),
                                           light_Sigma=matrix(c(1, 0.4, 0.4, 1), nrow=2),
                                           salMax_Sigma=matrix(c(1, 0.4, 0.4, 1), nrow=2),
                                           salSpan_Sigma=matrix(c(1, 0.4, 0.4, 1), nrow=2),
                                           alt_ls=NULL
                                           ) {
  swim_mx <- MASS::mvrnorm(n_sim, c(0, 0, 0, 0), swim_Sigma)
  light_mx <- MASS::mvrnorm(n_sim, c(0,0), light_Sigma)
  salMax_mx <- MASS::mvrnorm(n_sim, c(0,0), salMax_Sigma)
  salSpan_mx <- MASS::mvrnorm(n_sim, c(0,0), salSpan_Sigma)

  # Define minimum and maximum for each parameter
  bounds <- list(
    variableDh=c("true", "false"),
    variableDhV=c("true", "false"),
    D_h=c(1e-4, 1e1),
    D_hVert=c(1e-6, 1e0),
    mortSal_fn=c("constant", "logistic"),
    eggTemp_fn=c("constant", "logistic"),
    lightN=c(0.05, 0.5),
    lightC=c(2e-6, 2e-4),
    coldPrefN=c("true", "false"),
    swimUpN=c(0.05, 5)*1e-3, # mm/s -> m/s
    swimDownN=c(0.05, 2.5)*1e-3, # mm/s -> m/s
    swimUpC=c(0.1, 10)*1e-3, # mm/s -> m/s
    swimDownC=c(0.1, 5)*1e-3, # mm/s -> m/s
    passiveSinkSal=c(T, F),
    salThreshMaxN=c(20, 32),
    salThreshSpanN=c(0, 15),
    salThreshMaxC=c(20, 32),
    salThreshSpanC=c(0, 15),
    viableDD=c(30, 55),
    maxDepth=c(10, 300),
    connectRadius=c(15, 300)
  )
  # Replace any defaults
  if(!is.null(alt_ls)) {
    for(i in seq_along(alt_ls)) {
      bounds[[names(alt_ls)[i]]] <- alt_ls[[i]]
    }
  }

  if(mode=="lhs") {
    library(lhs)
    LHS <- randomLHS(n_sim, length(bounds))
    # Several parameters are sampled on a transformed scale. This is to counteract
    # the over-representation of larger values (and lack of resolution among low
    # values) when the plausible range spans orders of magnitude.
    sim.i <- tibble(
      variableDh=bounds$variableDh[qinteger(LHS[,1], 1, 2)],
      variableDhV=bounds$variableDhV[qinteger(LHS[,2], 1, 2)],
      # Diffusion coefficients: sample on a log scale
      D_h=log(bounds$D_h) |>
        qunif_minmax(LHS[,3]) |>
        exp(),
      D_hVert=log(bounds$D_hVert) |>
        qunif_minmax(LHS[,4]) |>
        exp(),
      # Mortality function
      mortSal_fn=bounds$mortSal_fn[qinteger(LHS[,5], 1, length(bounds$mortSal_fn))],
      # Egg production function
      eggTemp_fn=bounds$eggTemp_fn[qinteger(LHS[,6], 1, length(bounds$eggTemp_fn))],
      # Light responses: sample on a sqrt scale
      lightThreshNauplius=sqrt(bounds$lightN) |>
        qunif_minmax(LHS[,7]) |>
        pow(2),
      lightThreshCopepodid=sqrt(bounds$lightC) |>
        qunif_minmax(LHS[,8]) |>
        pow(2),
      # Nauplius only swim up if surface is colder
      swimColdNauplius=bounds$coldPrefN[qinteger(LHS[,9], 1, length(bounds$coldPrefN))],
      # Swim speeds: sample on a sqrt scale
      swimUpSpeedNaupliusMean=sqrt(bounds$swimUpN) |>
        qunif_minmax(LHS[,10]) |>
        pow(2) |>
        multiply(-1),
      swimDownSpeedNaupliusMean=sqrt(bounds$swimDownN) |>
        qunif_minmax(LHS[,11]) |>
        pow(2),
      swimUpSpeedCopepodidMean=sqrt(bounds$swimUpC) |>
        qunif_minmax(LHS[,12]) |>
        pow(2) |>
        multiply(-1),
      swimDownSpeedCopepodidMean=sqrt(bounds$swimDownC) |>
        qunif_minmax(LHS[,13]) |>
        pow(2),
      # Sink rate: salinity-dependent or defined by swimDownSpeed ?
      passiveSinkRateSal=bounds$passiveSinkSal[qinteger(LHS[,14], 1, 2)],
      # Salinity thresholds: define max, then define psu span of 0-100% sinking
      salinityThreshNaupliusMax=bounds$salThreshMaxN |>
        qunif_minmax(LHS[,15]),
      salinityThreshNaupliusMin=salinityThreshNaupliusMax -
        (bounds$salThreshSpanN |>
           qunif_minmax(LHS[,16])),
      salinityThreshCopepodidMax=bounds$salThreshMaxC |>
        qunif_minmax(LHS[,17]),
      salinityThreshCopepodidMin=salinityThreshCopepodidMax -
        (bounds$salThreshSpanC |>
           qunif_minmax(LHS[,18])),
      # Degree days for transition to copepodid
      viableDegreeDays=bounds$viableDD |>
        qunif_minmax(LHS[,19]),
      # Maximum preferred depth: sample on a log scale
      maxDepth=log(bounds$maxDepth) |>
        qunif_minmax(LHS[,20]) |>
        exp(),
      # Connectivity radius around pens
      connectivityThresh=sqrt(bounds$connectRadius) |>
        qunif_minmax(LHS[,21]) |>
        pow(2)
    )  |>
      rowwise() |>
      mutate(eggTemp_b=sample(egg_post[[eggTemp_fn]], 1),
             mortSal_b=sample(mort_post[[mortSal_fn]], 1)) |>
      ungroup()
  }

  if(mode=="random") {
    # Several parameters are sampled on a transformed scale. This is to counteract
    # the over-representation of larger values (and lack of resolution among low
    # values) when the plausible range spans orders of magnitude.
    sim.i <- tibble(
      variableDh=bounds$variableDh |>
        sample(n_sim, replace=T),
      variableDhV=bounds$variableDhV |>
        sample(n_sim, replace=T),
      # Diffusion coefficients: sample on a log scale
      D_h=log(bounds$D_h) |>
        runif_minmax(n_sim) |>
        exp(),
      D_hVert=log(bounds$D_hVert) |>
        runif_minmax(n_sim) |>
        exp(),
      # Mortality function
      mortSal_fn=bounds$mortSal_fn |>
        sample(n_sim, replace=T),
      # Egg production function
      eggTemp_fn=bounds$eggTemp_fn |>
        sample(n_sim, replace=T),
      # Light responses: sample on a sqrt scale
      lightThreshNauplius=sqrt(bounds$lightN) |>
        qunif_minmax(pnorm(light_mx[,1])) |>
        pow(2),
      lightThreshCopepodid=sqrt(bounds$lightC) |>
        qunif_minmax(pnorm(light_mx[,2])) |>
        pow(2),
      # Nauplius only swim up if surface is colder
      swimColdNauplius=bounds$coldPrefN |>
        sample(n_sim, replace=T),
      # Swim speeds: sample on a sqrt scale
      swimUpSpeedNaupliusMean=sqrt(bounds$swimUpN) |>
        qunif_minmax(pnorm(swim_mx[,1])) |>
        pow(2) |>
        multiply(-1),
      swimDownSpeedNaupliusMean=sqrt(bounds$swimDownN) |>
        qunif_minmax(pnorm(swim_mx[,2])) |>
        pow(2),
      swimUpSpeedCopepodidMean=sqrt(bounds$swimUpC) |>
        qunif_minmax(pnorm(swim_mx[,3])) |>
        pow(2) |>
        multiply(-1),
      swimDownSpeedCopepodidMean=sqrt(bounds$swimDownC) |>
        qunif_minmax(pnorm(swim_mx[,4])) |>
        pow(2),
      # Sink rate: salinity-dependent or defined by swimDownSpeed ?
      passiveSinkRateSal=bounds$passiveSinkSal |>
        sample(n_sim, replace=T),
      # Salinity thresholds: define max, then define psu span of 0-100% sinking
      salinityThreshNaupliusMax=bounds$salThreshMaxN |>
        qunif_minmax(pnorm(salMax_mx[,1])),
      salinityThreshNaupliusMin=salinityThreshNaupliusMax -
        (bounds$salThreshSpanN |>
           qunif_minmax(pnorm(salSpan_mx[,1]))),
      salinityThreshCopepodidMax=bounds$salThreshMaxC |>
        qunif_minmax(pnorm(salMax_mx[,2])),
      salinityThreshCopepodidMin=salinityThreshCopepodidMax -
        (bounds$salThreshSpanC |>
           qunif_minmax(pnorm(salSpan_mx[,2]))),
      # Degree days for transition to copepodid
      viableDegreeDays=bounds$viableDD |>
        runif_minmax(n_sim),
      # Maximum preferred depth: sample on a log scale
      maxDepth=log(bounds$maxDepth) |>
        runif_minmax(n_sim) |>
        exp(),
      # Connectivity radius around pens
      connectivityThresh=sqrt(bounds$connectRadius) |>
        runif_minmax(n_sim) |>
        pow(2)
    )  |>
      rowwise() |>
      mutate(eggTemp_b=sample(egg_post[[eggTemp_fn]], 1),
             mortSal_b=sample(mort_post[[mortSal_fn]], 1)) |>
      ungroup()
  }
  # Add posterior samples from sink rate regression
  if(is.null(sink_post)) {
    sim.i <- sim.i |>
      mutate(passiveSinkInt=swimDownSpeedNaupliusMean,
             passiveSinkSlope=0)
  } else {
    sim.i <- sim.i |>
      mutate(passiveSinkInt=sample(sink_post$Intercept, n_sim, replace=T)*passiveSinkRateSal +
               swimDownSpeedNaupliusMean*(!passiveSinkRateSal),
             passiveSinkSlope=sample(sink_post$slope, n_sim, replace=T)*passiveSinkRateSal)
  }
  sim.i <- sim.i |>
    mutate(across(where(is.numeric), ~signif(.x, 5))) |>
    mutate(i=str_pad(row_number(), nchar(n_sim), "left", "0"),
           outDir=glue("{out_dir}/sim_{i}/"))
  return(sim.i)
}




calc_GSA_connectivity <- function(f, farms, sim, site_areas, cores=4) {
  library(furrr)
  if(get_os()=="windows") {
    plan(multisession, workers=cores)
  } else {
    plan(multicore, workers=cores)
  }
  c_i <- future_map_dfr(f,
                        ~load_connectivity(.x,
                                           source_names=farms,
                                           dest_names=farms,
                                           liceScale=1)) |>
    mutate(sim=sim)
  plan(sequential)
  c_i_daily <- calc_daily_fluxes(c_i, site_areas)
  c_summary <- c_i_daily |> calc_sensitivity_outcomes(sim)
  c_farm <- c_i_daily |> group_by(sepaSite) |> calc_sensitivity_outcomes(sim)

  return(list(og=c_i, daily=c_i_daily, summary=c_summary, farm=c_farm))
}





calc_daily_fluxes <- function(c_long, site_areas) {
  list(
    c_long |>
      calc_influx(destination, value, sim, date) |>
      rename(sepaSite=destination),
    c_long |>
      calc_self_infection(source, destination, value, sim, date) |>
      rename(sepaSite=source, selfflux=self),
    c_long |>
      calc_outflux(source, value, dest_areas=NULL, sim, date) |>
      rename(sepaSite=source)
  ) |>
    reduce(full_join) |>
    complete(sepaSite, sim, date,
             fill=list(influx=0, N_influx=0,
                       selfflux=0,
                       outflux=0, N_outflux=0)) |>
    left_join(site_areas, by=join_by(sepaSite)) |>
    mutate(influx_m2=influx/area,
           selfflux_m2=selfflux/area,
           outflux_m2=outflux/area) |>
    select(-area)
}




calc_sensitivity_outcomes <- function(c_daily, sim) {
  flux_ls <- list(
    c_daily |>
      summarise(across(contains("flux"), list(mn=mean, md=median))),
    c_daily |>
      filter(between(date, ymd("2023-03-01"), ymd("2023-05-31"))) |>
      summarise(across(contains("flux"), list(MAM_mn=mean, MAM_md=median))),
    c_daily |>
      filter(between(date, ymd("2023-06-01"), ymd("2023-08-31"))) |>
      summarise(across(contains("flux"), list(JJA_mn=mean, JJA_md=median))),
    c_daily |>
      filter(between(date, ymd("2023-09-01"), ymd("2023-11-30"))) |>
      summarise(across(contains("flux"), list(SON_mn=mean, SON_md=median)))
  )
  if(is_grouped_df(c_daily)) {
    flux_ls |>
      reduce(full_join) |>
      mutate(sim=sim)
  } else {
    flux_ls |>
      reduce(bind_cols) |>
      mutate(sim=sim)
  }
}



make_sensitivity_sum_df <- function(sum_ls, sim.i) {
  sum_ls |>
    reduce(bind_rows) |>
    inner_join(sim.i |> select(-ends_with("_b"), -passiveSinkInt, -passiveSinkSlope, -outDir),
               by=join_by(sim==i)) |>
    mutate(mortSal_fn=as.numeric(mortSal_fn=="logistic"),
           eggTemp_fn=as.numeric(eggTemp_fn=="logistic"),
           passiveSinkRateSal=as.numeric(passiveSinkRateSal),
           across(starts_with("swim"), abs))
}



run_sensitivity_ML <- function(outcome_names, predictor_names,
                               sum_df, sim.i,
                               method="rf") {
  df_ls <- map(outcome_names,
      ~sum_df |>
        rename(outcome=.x) |>
        mutate(outcome=outcome^0.25) |>
        select(outcome, any_of(predictor_names)) |>
        drop_na())
  if(method=="rf") {
    rf_ls <- map(df_ls, ~randomForest(outcome ~ ., data=.x, importance=T))
  }
  exp_ls <- list(x_valid=df_ls |> map(~.x |> select(-outcome)),
                 y_valid=df_ls |> map(~.x$outcome))
  return(list(rf=rf_ls, exp=exp_ls))
}




summarise_importance <- function(ML_ls, outcome_names) {
  imap_dfr(ML_ls,
           ~importance(.x) |>
             as_tibble(rownames="param") |>
             mutate(outcome=outcome_names[.y])) |>
    mutate(flux=case_when(grepl("influx", outcome) ~ "influx",
                          grepl("outflux", outcome) ~ "outflux",
                          grepl("selfflux", outcome) ~ "self"),
           mn_md=if_else(grepl("_mn", outcome), "mean", "median"),
           type=case_when(grepl("^N_", outcome) ~ "n_farms",
                          grepl("_m2", outcome) ~ "IP_m2",
                          .default="IP total"),
           season=case_when(grepl("MAM", outcome) ~ "MAM",
                            grepl("JJA", outcome) ~ "JJA",
                            grepl("SON", outcome) ~ "SON",
                            .default="all"))
}



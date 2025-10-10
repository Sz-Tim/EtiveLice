




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
    inner_join(sim.i |> select(1:20), by=join_by(sim==i)) |>
    mutate(mortSal_fn=as.numeric(mortSal_fn=="logistic"),
           eggTemp_fn=as.numeric(eggTemp_fn=="logistic"),
           passiveSinkRateSal=as.numeric(passiveSinkRateSal),
           across(starts_with("swim"), abs))
}



run_sensitivity_ML <- function(outcome_names, sum_df, sim.i, method="rf") {
  df_ls <- map(outcome_names,
      ~sum_df |>
        rename(outcome=.x) |>
        mutate(outcome=outcome^0.25) |>
        select(outcome, any_of(names(sim.i)[1:19])) |>
        drop_na())
  if(method=="rf") {
    map(df_ls, ~randomForest(outcome ~ ., data=.x, importance=T))
  }
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



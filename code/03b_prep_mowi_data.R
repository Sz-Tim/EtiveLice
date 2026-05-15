# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Prep for Stan runs with Mowi data



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR@dev")



# farm data ---------------------------------------------------------------

mowi_df_ext <- read_csv("data/aquaculture/mowi_cleaned.csv") |>
  complete(sepaSite, date, fill=list(nFish=0, biomass=0)) |>
  arrange(date, sepaSite) |>
  rename(nFish_est=nFish) |>
  mutate(day=as.numeric(factor(date)),
         sepaSite=factor(sepaSite),
         pen=sepaSite,
         nFishSampled=if_else(is.na(nFishSampled), 0, nFishSampled),
         sampled=nFishSampled > 0)
farm_i <- read_csv("data/farm_sites_widerLinnhe_2022-2025.csv") |>
  filter(sepaSite %in% names(sepa_key)) |>
  inner_join(read_csv(glue("data/sim/inputs/farm_influx_vols.csv"), show_col_types=F) |>
               filter(sim=="01") |>
               mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1]) |>
               group_by(sepaSite) |>
               summarise(vol=sum(vol)) |>
               ungroup(),
             by="sepaSite") |>
  arrange(sepaSite)
trt_df <- read_csv("data/aquaculture/lice_data_compiled.csv") |>
  rename(date=weekBeginning) |>
  inner_join(read_csv("data/aquaculture/mowi_cleaned.csv") |> select(sepaSite, date)) |>
  filter(!is.na(mitigation)) |>
  filter(mitigation %in% c("Bath", "Bio Reduction", "In Feed", "Physical")) |>
  select(sepaSite, date, mitigation) |>
  right_join(mowi_df_ext |> select(sepaSite, date),
             by=join_by(sepaSite, date)) |>
  mutate(trt=as.numeric(!is.na(mitigation))) |>
  arrange(sepaSite, date)

info <- list(nDays=as.numeric(diff(range(mowi_df_ext$date)))+1,
             nHours=(as.numeric(diff(range(mowi_df_ext$date)))+1)*24,
             nFarms=n_distinct(mowi_df_ext$sepaSite),
             nSims=10,
             nStages=5, # cII-V, Ad (Piasecki 2023)
             nStageGroups=3, # ch, mobile, adult
             stg_grp_ii=c(1, 1, 2, 2, 3), # index for matching stages to stage groups
             nPens=1,
             IP_penVolume=farm_i$vol,
             dateRange=range(mowi_df_ext$date))
day_hour <- matrix(1:info$nHours, ncol=24, byrow=T)
nFish_mx <- make_nFish_mx(mowi_df_ext, info)
sampledDays <- make_sampledDays(mowi_df_ext |> mutate(pen=sepaSite))
nFishSampled_mx <- make_nFishSampled_mx(mowi_df_ext, info, nFish_mx)
treatDays_mx <- matrix(trt_df$trt, nrow=info$nDays, ncol=info$nFarms)
mowi_y_obs <- mowi_df_ext |>
  mutate(across(all_of(c("Ch", "PA", "AF")), ~if_else(is.na(.x), 0, .x))) |>
  select(sepaSite, date, Ch, PA, AF) |>
  pivot_longer(all_of(c("Ch", "PA", "AF")), names_to="stage", values_to="N") |>
  mutate(stage=factor(stage, levels=c("Ch", "PA", "AF"))) |>
  arrange(sepaSite, date, stage)
y_obs <- array(c(mowi_y_obs$N), dim=c(info$nStageGroups, info$nDays, info$nFarms))

saveRDS(info, "data/aquaculture/mowi_stan/info.rds")
saveRDS(params, "data/aquaculture/mowi_stan/params.rds")
saveRDS(day_hour, "data/aquaculture/mowi_stan/day_hour.rds")
saveRDS(nFish_mx, "data/aquaculture/mowi_stan/nFish_mx.rds")
saveRDS(sampledDays, "data/aquaculture/mowi_stan/sampledDays.rds")
saveRDS(nFishSampled_mx, "data/aquaculture/mowi_stan/nFishSampled_mx.rds")
saveRDS(treatDays_mx, "data/aquaculture/mowi_stan/treatDays.rds")
saveRDS(y_obs, "data/aquaculture/mowi_stan/y_obs.rds")



# environmental data ------------------------------------------------------

farm_env <- read_csv()
farm_env_daily <- read_csv()
attach_env_mx <- make_attach_env_mx(farm_env, info, params)
sal_mx <- make_sal_mx(farm_env_daily, info, params)
temp_mx <- make_temp_mx(farm_env_daily, info)
temp_z_mx <- make_temp_z_mx(farm_env_daily, info)

saveRDS(attach_env_mx, "data/aquaculture/mowi_stan/attach_env_mx.rds")
saveRDS(sal_mx, "data/aquaculture/mowi_stan/sal_mx.rds")
saveRDS(temp_mx, "data/aquaculture/mowi_stan/temp_mx.rds")
saveRDS(temp_z_mx, "data/aquaculture/mowi_stan/temp_z_mx.rds")



# infestation pressure ----------------------------------------------------

influx_df <- read_csv()
make_IP_mx(influx_df, info)
saveRDS(IP_mx, "data/aquaculture/mowi_stan/IP_mx.rds")



# params ------------------------------------------------------------------

# Dimensions only, for setting priors in make_stan_data()
params <- list(attach_beta=rep(0, ncol(attach_env_mx)),
               surv_beta=matrix(0, ncol=info$nStageGroups, nrow=2))
# For calculating priors: mean stage duration
mnDaysStage_df <- expand_grid(Ch=seq(125, 175, length.out=100),
                              PA=seq(325, 375, length.out=100)) |>
  mutate(PA=PA-Ch) |>
  pivot_longer(everything(), names_to="stage", values_to="GDD") |>
  mutate(temperature=list(seq_range(farm_env$temperature, length.out=100))) |>
  unnest(temperature) |>
  mutate(days=GDD/temperature,
         temp_z=(temperature-farm_env_avg$temperature[1])/farm_env_avg$temperature[2])
# split to Ch1, Ch2 from Ch1+2: assume equal duration
params$mnDaysStageCh <- summary(
  lm(days ~ temp_z, data=mnDaysStage_df |> filter(stage=="Ch")))$coefficients[,1:2] *
  rep(c(1, sqrt(nrow(mnDaysStage_df)/2)), each=2)
# split to PA1, PA2 from PA1+2: assume equal duration
params$mnDaysStagePA <- summary(
  lm(days ~ temp_z, data=mnDaysStage_df |> filter(stage=="PA")))$coefficients[,1:2] *
  rep(c(1, sqrt(nrow(mnDaysStage_df)/2)), each=2)

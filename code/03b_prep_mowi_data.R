# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Prep for Stan runs with Mowi data



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR@dev")
dir("code/fn", ".R", full.names=T) |> walk(source)

dat_dir <- "data/aquaculture"
dat_stan_dir <- "data/aquaculture/mowi_stan"
inputs_dir <- "data/biotracker/2022_2025"

sepa_key <- c("APT1"="Etive 4",
              "ARDG1"="Ardgour",
              "CAG1"="Camas Glas",
              "CALL1"="Leven",
              "FFMC84"="Etive 3",
              "GORS1"="Gorsten",
              "INV1"="Invasion Bay",
              "KING1"="Kingairloch",
              "MCLN1"="Macleans Nose",
              "SAR1"="Etive 6")



# farm data ---------------------------------------------------------------

mowi_df_ext <- read_csv(glue("{dat_dir}/mowi_cleaned.csv")) |>
  complete(sepaSite, date, fill=list(nFish=0, biomass=0)) |>
  arrange(date, sepaSite) |>
  rename(nFish_est=nFish) |>
  mutate(day=as.numeric(factor(date)),
         pen=sepaSite,
         nFishSampled=if_else(is.na(nFishSampled), 0, nFishSampled),
         sampled=nFishSampled > 0) |>
  left_join(read_csv(glue("{inputs_dir}/farm_dat.csv"), show_col_types=F) |>
              summarise(.by=c(sepaSite, maximumBiomassAllowedTonnes))) |>
  left_join(read_csv(glue("{inputs_dir}/farm_influx_vols.csv"), show_col_types=F) |>
              mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1]) |>
              summarise(vol=sum(vol),
                        .by=sepaSite)) |>
  mutate(RW=sqrt(biomass/(20*vol)),
         RW_logit=brms::logit_scaled(RW, lb=-1e-5),
         sepaSite=factor(sepaSite))
farm_i <- read_csv("data/farm_sites_widerLinnhe_2022-2025.csv") |>
  filter(sepaSite %in% names(sepa_key)) |>
  inner_join(read_csv(glue("{inputs_dir}/farm_influx_vols.csv"), show_col_types=F) |>
               mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1]) |>
               group_by(sepaSite) |>
               summarise(vol=sum(vol)) |>
               ungroup(),
             by="sepaSite") |>
  arrange(sepaSite)
trt_df <- read_csv(glue("{dat_dir}/mowi_trt_cleaned.csv")) |>
  select(sepaSite, date, TypeNum) |>
  mutate(TypeNum=paste0("t_", TypeNum),
         TypeApplied=1) |>
  inner_join(read_csv(glue("{dat_dir}/mowi_cleaned.csv")) |> select(sepaSite, date)) |>
  arrange(TypeNum, sepaSite, date) |>
  pivot_wider(names_from=TypeNum, values_from=TypeApplied, values_fill=0) |>
  right_join(mowi_df_ext |> select(sepaSite, date),
             by=join_by(sepaSite, date)) |>
  mutate(across(starts_with("t_"), ~replace_na(.x, 0))) |>
  arrange(sepaSite, date)
trt_meth_ii <- read_csv(glue("{dat_dir}/mowi_trt_cleaned.csv")) |>
  summarise(.by=c(MethodNum, TypeNum, Method, Type)) |>
  arrange(TypeNum)

info <- list(nDays=as.numeric(diff(range(mowi_df_ext$date)))+1,
             nHours=(as.numeric(diff(range(mowi_df_ext$date)))+1)*24,
             nFarms=n_distinct(mowi_df_ext$sepaSite),
             nSims=8,
             nStages=5, # cII-V, Ad (Piasecki 2023)
             nStageGroups=3, # ch, mobile, adult
             stg_grp_ii=c(1, 1, 2, 2, 3), # index for matching stages to stage groups
             nPens=1,
             nTrtMethods=max(trt_meth_ii$MethodNum),
             nTrtTypes=max(trt_meth_ii$TypeNum),
             trt_meth_ii=trt_meth_ii$MethodNum,
             IP_penVolume=farm_i$vol,
             dateRange=range(mowi_df_ext$date))
# Dimensions only, for setting priors and covariate structures
params <- list(attach_beta=rep(0, 5),
               surv_beta=matrix(0, ncol=3, nrow=2))

day_hour <- matrix(1:info$nHours, ncol=24, byrow=T)
nFish_mx <- make_nFish_mx(mowi_df_ext, info)
sampledDays <- make_sampledDays(mowi_df_ext |> mutate(pen=sepaSite))
nFishSampled_mx <- make_nFishSampled_mx(mowi_df_ext, info, nFish_mx)
trtApplied_mx <- make_trtApplied_mx(trt_df, info)
mowi_y_obs <- mowi_df_ext |>
  mutate(across(all_of(c("Ch", "PA", "AF")), ~if_else(is.na(.x), 0, .x))) |>
  select(sepaSite, date, Ch, PA, AF) |>
  pivot_longer(all_of(c("Ch", "PA", "AF")), names_to="stage", values_to="N") |>
  mutate(stage=factor(stage, levels=c("Ch", "PA", "AF"))) |>
  arrange(sepaSite, date, stage)
y_obs <- array(c(mowi_y_obs$N), dim=c(info$nStageGroups, info$nDays, info$nFarms))

saveRDS(info, glue("{dat_stan_dir}/info.rds"))
saveRDS(day_hour, glue("{dat_stan_dir}/day_hour.rds"))
saveRDS(nFish_mx, glue("{dat_stan_dir}/nFish_mx.rds"))
saveRDS(sampledDays, glue("{dat_stan_dir}/sampledDays.rds"))
saveRDS(nFishSampled_mx, glue("{dat_stan_dir}/nFishSampled_mx.rds"))
saveRDS(trtApplied_mx, glue("{dat_stan_dir}/trtApplied_mx.rds"))
saveRDS(y_obs, glue("{dat_stan_dir}/y_obs.rds"))



# environmental data ------------------------------------------------------

farm_env <- read_csv(glue("{inputs_dir}/farm_env_hourly.csv"), show_col_types=F) |>
  filter(sepaSite %in% names(sepa_key)) |>
  mutate(sepaSite=factor(sepaSite, levels=names(sepa_key)),
         date=date(time),
         hour=hour(time)) |>
  mutate(day=as.numeric(as.factor(date))) |>
  inner_join(
    mowi_df_ext |>
      select(sepaSite, date, RW_logit),
    by=join_by(sepaSite, date)) |>
  arrange(sepaSite, time) |>
  # calculate farm-level averages
  mutate(u=u*100, # cm/s
         v=v*100, # cm/s
         w=w*100, # cm/s
         uv=uv*100, # cm/s
         uv_sq=uv^2,
         salinity_m30=salinity - 30) |> # recenter so intercept = high salinity
  mutate(across(c(temperature, u, v, w, uv, salinity), ~c(scale(.x)), .names="{.col}_z")) |>
  mutate(uv_z_sq=uv_z^2)
# for back-transforming z-scores
farm_env_avg <- farm_env |>
  reframe(across(where(is.numeric), ~c(mn=mean(.x), sd=sd(.x)))) |>
  select(-ends_with("_z")) |>
  mutate(metric=c("mean", "sd"))
# daily values
farm_env_daily <- farm_env |>
  select(-time, -hour, -elapsedHours) |>
  summarise(across(where(is.numeric), mean),
            .by=c(sepaSite, day, date))
attach_env_mx <- make_attach_env_mx(farm_env, info, params)
sal_mx <- make_sal_mx(farm_env_daily, info, params)
temp_mx <- make_temp_mx(farm_env_daily, info)
temp_z_mx <- make_temp_z_mx(farm_env_daily, info)
ydayh_mx <- make_ydayh_mx(farm_env)

saveRDS(farm_env, glue("{dat_stan_dir}/farm_env.rds"))
saveRDS(farm_env_avg, glue("{dat_stan_dir}/farm_env_avg.rds"))
saveRDS(farm_env_daily, glue("{dat_stan_dir}/farm_env_daily.rds"))
saveRDS(attach_env_mx, glue("{dat_stan_dir}/attach_env_mx.rds"))
saveRDS(sal_mx, glue("{dat_stan_dir}/sal_mx.rds"))
saveRDS(temp_mx, glue("{dat_stan_dir}/temp_mx.rds"))
saveRDS(temp_z_mx, glue("{dat_stan_dir}/temp_z_mx.rds"))
saveRDS(ydayh_mx, glue("{dat_stan_dir}/ydayh_mx.rds"))



# infestation pressure ----------------------------------------------------

influx_df <- read_csv(glue("{inputs_dir}/influx_hourly.csv"), show_col_types=F) |>
  filter(sepaSite %in% names(sepa_key)) |>
  full_join(tibble(time=seq_range(range(farm_env$time), by=3600),
                   sepaSite=names(sepa_key)[1], pen=paste0(names(sepa_key)[1], "_01"), sim="01"),
            by=join_by(sepaSite, sim, pen, time)) |>
  filter(between(time, min(farm_env$time), max(farm_env$time))) |>
  complete(time, pen, sim, fill=list(influx=0, influx_m3=0)) |>
  mutate(date=date(time),
         day=as.numeric(as.factor(date)),
         hour=hour(time),
         sepaSite=str_split_fixed(pen, "_", 2)[,1],
         sepaSite=factor(sepaSite, levels=names(sepa_key))) |>
  arrange(sim, sepaSite, pen, time) |>
  summarise(pen=first(pen),
            influx=sum(influx),
            .by=c(time, date, day, hour, sepaSite, sim))
IP_mx <- make_IP_mx(influx_df, info)
saveRDS(IP_mx, glue("{dat_stan_dir}/IP_mx.rds"))



# params ------------------------------------------------------------------

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

saveRDS(params, glue("{dat_stan_dir}/params.rds"))

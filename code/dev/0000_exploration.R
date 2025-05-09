library(tidyverse)
library(glue)
library(scico)
library(biotrackR)
theme_set(theme_bw())



# load data ---------------------------------------------------------------

etive_sites <- c("FFMC84", "FFMC32", "APT1", "SAR1", "FFMC27")
ens_dir <- "../sealice_ensembling/"
site_i <- read_csv(glue("{ens_dir}/data/farm_sites_100m_areas.csv")) |>
  rename(area=area_m2) |>
  filter(sepaSite %in% etive_sites) |>
  mutate(sepaSite=factor(sepaSite, levels=etive_sites)) |>
  arrange(sepaSite) |>
  mutate(depth=c(52, 69, 31, 42, 28),
         vol=area*depth)
site_env <- read_csv(glue("{ens_dir}/out/siteConditions_2019-2023.csv")) |>
  rename(sepaSite=site) |>
  filter(date >= "2023-01-01" & sepaSite %in% etive_sites) |>
  mutate(sepaSite=factor(sepaSite, levels=etive_sites),
         across(where(is.numeric), ~c(scale(.x)), .names="{.col}_z"))
site_env_avg <- site_env |>
  reframe(across(where(is.numeric), ~c(mn=mean(.x), sd=sd(.x)))) |>
  select(-ends_with("_z"))
copInflux_df <- readRDS(glue("{ens_dir}/out/sim_2019-2023/processed/connectivity_day.rds")) |>
  filter(date >= "2023-01-01" & sepaSite %in% etive_sites) |>
  select(sepaSite, sim, date, influx) |>
  mutate(sepaSite=factor(sepaSite, levels=etive_sites)) |>
  left_join(site_i |> select(sepaSite, area, vol)) |>
  mutate(influx=influx*24, # hourly to daily
         influx_m2=influx/area,
         influx_m3=influx/vol)
farm_dat <- read_csv(glue("{ens_dir}/data/lice_biomass_2017-01-01_2023-12-31.csv")) |>
  filter(date >= "2023-01-01") |>
  inner_join(site_i |> select(sepaSite, vol)) |>
  mutate(RW=sqrt(actualBiomassOnSiteTonnes*1000/(20*vol)),
         RW_logit=pmax(boot::logit(RW), -1e3))
# TODO: TO BE REPLACED WITH ACTUAL NUMBERS!
nFish_dat <- farm_dat |>
  filter(actualBiomassOnSiteTonnes > 0) |>
  mutate(nFish=240*actualBiomassOnSiteTonnes) |>
  group_by(sepaSite) |>
  summarise(nFish=median(nFish)) |>
  full_join(farm_dat |> select(sepaSite, date)) |>
  group_by(sepaSite) |>
  mutate(nFish_est=if_else(row_number()==1, nFish, lag(nFish)*(0.999^row_number())))


# site overviews ----------------------------------------------------------

site_env |>
  select(-ends_with("_z")) |>
  pivot_longer(cols=c("u", "v", "w", "uv", "salinity", "temperature"),
               names_to="variable", values_to="value") |>
  mutate(variable=factor(variable,
                         levels=c("salinity", "temperature", "u", "v", "w", "uv"))) |>
  ggplot(aes(date, value, colour=sepaSite)) +
  geom_line() +
  scale_colour_scico_d("Site", palette="devon", end=0.8, direction=-1) +
  facet_grid(variable~., scales="free_y") +
  labs(x="Date", y="Daily mean, 0-30m (WeStCOMS2)") +
  theme(panel.grid.major.y=element_line(colour="grey70", linewidth=0.2),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.1))

copInflux_df |>
  complete(date, sepaSite, sim, fill=list(influx_m2=0)) |>
  ggplot(aes(date, influx_m3, colour=sim)) +
  geom_line() +
  scale_colour_viridis_d(option="turbo") +
  scale_y_continuous("Mean daily copepodids / hour / m3 (150m radius cylinder)",
                     transform="pseudo_log",
                     breaks=c(0, 1, 10, 20, 40, 60)) +
  facet_grid(sepaSite~.)

copInflux_df |>
  complete(date, sepaSite, sim, fill=list(influx_m2=0)) |>
  ggplot(aes(date, influx_m3, colour=sim)) +
  geom_line() +
  scale_colour_viridis_d(option="turbo") +
  scale_y_continuous("Mean daily copepodids / hour / m3 (150m radius cylinder)") +
  facet_grid(sepaSite~.)

farm_dat |>
  ggplot(aes(date, actualBiomassOnSiteTonnes, colour=sepaSite)) +
  geom_line() +
  scale_colour_scico_d("Site", palette="devon", end=0.8, direction=-1)

farm_dat |>
  ggplot(aes(date, RW, colour=sepaSite)) +
  geom_line() +
  scale_colour_scico_d("Site", palette="devon", end=0.8, direction=-1)


# simulation --------------------------------------------------------------

# . fn (temp) -------------------------------------------------------------

make_composition <- function(x) {
  x/sum(x)
}



# . setup -----------------------------------------------------------------

# On-fish stages: Chalimus, pre-adult, adult, gravid

# Detection rates and attachment rates will not be fully identifiable. We could
# either put strong priors on one of them, or fix one (which amounts to the
# same thing -- an extremely strong spike prior). For now, I'll fix the
# gravid detection rate at 1.

# Kragesteen separates males from females because the development time seems
# to be different. Males are really only relevant in that they appear in counts
# of chalimus and pre-adults (maybe non-gravid adults? Depends on Mowi's data).

set.seed(1000)

# metadata
info <- list(nDays=n_distinct(copInflux_df$date),
             nSites=n_distinct(copInflux_df$sepaSite),
             nSims=n_distinct(copInflux_df$sim),
             nStages=4)

# plausible parameters
params <- list(
  ensWts_p=rbeta(20, 2, 5) |> make_composition(),
  attach_beta=c(-0.5, 1, 0.18, 0.3, -0.01), # intercept, RW_logit, salinity, uv, uv^2
  surv_p=c(0.97, 0.998, 0.997, 0.942, 0, 0),
  surv_beta=cbind(c(-2.5, -1, -1, -1.5),
                  c(0.2, 0.18, 0.18, 0.16)),
  stageDurCumulGDD_F=c(150, 350, 550, 1500, Inf),
  stageDurCumulGDD_M=c(130, 315, 1500, 1500, Inf),
  detect_p=c(0.1, 0.7, 0.99, 1)
)

# PT variant copepodid infection pressure (daily copepodids/m3)
IP_df <- copInflux_df |>
  complete(date, sepaSite, sim, fill=list(influx=0, influx_m3=0)) |>
  mutate(day=as.numeric(as.factor(date))) |>
  select(day, sepaSite, sim, influx_m3, influx) |>
  arrange(sim, sepaSite, day)
# array version for faster calculations: IP[day, site, sim]
IP_mx <- array(IP_df$influx, dim=c(info$nDays, info$nSites, info$nSims))

# site conditions
siteEnv_df <- site_env |>
  mutate(day=as.numeric(as.factor(date))) |>
  inner_join(farm_dat |> select(sepaSite, date, RW_logit)) |>
  inner_join(nFish_dat |> select(sepaSite, date, nFish_est)) |>
  arrange(sepaSite, day) |>
  mutate(uv=uv*100,
         uv_sq=uv^2)
attach_env_mx <- array(1, dim=c(info$nDays, info$nSites, 5))
attach_env_mx[,,2] <- siteEnv_df$RW_logit
attach_env_mx[,,3] <- siteEnv_df$salinity_z
attach_env_mx[,,4] <- siteEnv_df$uv
attach_env_mx[,,5] <- siteEnv_df$uv_sq
sal_mx <- matrix(siteEnv_df$salinity, nrow=info$nDays)
temp_mx <- matrix(siteEnv_df$temperature, nrow=info$nDays)
nFish_mx <- matrix(siteEnv_df$nFish_est, nrow=info$nDays)

# Data structures
# N[cohort, day, site, M/F]
cohort_N <- array(0, dim=c(info$nDays, info$nDays, info$nSites, 2))
# stage[cohort, day, site, M/F]
cohort_stage <- cohort_N
# accumulated GDD[cohort, day, site, M/F]
cohort_GDD <- array(0, dim=c(info$nDays, info$nDays, info$nSites))
# surv_rate[day, site, stage]
stage_survRate <- array(0, dim=c(info$nDays, info$nSites, info$nStages+2))
# N[day, site, stage]
mu <- array(0, dim=c(info$nDays, info$nSites, info$nStages, 2))
y <- array(0, dim=c(info$nDays, info$nSites, info$nStages, 2))

# . ensemble IP -----------------------------------------------------------

ensIP <- apply(IP_mx, 1:2,
               function(x, ensWts_p) {x %*% ensWts_p},
               ensWts_p=params$ensWts_p)
matplot(ensIP, lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))


# . attachment ------------------------------------------------------------

attach_mx <- apply(attach_env_mx, 1:2,
                   function(x, attach_beta) {boot::inv.logit(x %*% attach_beta)},
                   attach_beta=params$attach_beta)

matplot(attach_mx, lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))
matplot(ensIP * attach_mx, lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))



# . cohort progression ----------------------------------------------------


# How does IP relate to mean N per fish?
for(cohort in 1:info$nDays) {
  for(day in 1:info$nDays) {
    # calculate survival rates
    for(stage in 1:info$nStages) {
      stage_survRate[day, , stage] <- boot::inv.logit(
        params$surv_beta[stage,1] + params$surv_beta[stage,2]*sal_mx[day,]
      )
    }
    if(cohort == day) {
      # initialize as attached copepodids, equal M/F distribution
      cohort_stage[cohort, day, , ] <- 1
      cohort_N[cohort, day, , ] <- ((((ensIP)^0.5 * attach_mx)[day,])^2)/2
      # assume GDD applies first day (attachment at 00:00:00)
      cohort_GDD[cohort, day, ] <- temp_mx[day, ]
    }

    if(day > cohort) {
      # apply survival rate based on stage at start of day
      # cohort_N[cohort, day, , ] <- cohort_N[cohort, day-1, , ] *
      #   params$surv_p[cohort_stage[cohort, day-1, , ]]
      for(site in 1:info$nSites) {
        cohort_N[cohort, day, site, ] <- cohort_N[cohort, day-1, site, ] *
          stage_survRate[day-1, site, cohort_stage[cohort, day-1, site, ]]
      }
      # determine stage at END of day based on GDD accumulation
      cohort_stage[cohort, day, , 1] <- cohort_stage[cohort, day-1, , 1] +
        (cohort_GDD[cohort, day-1, ] > params$stageDurCumulGDD_M[cohort_stage[cohort, day-1, , 1]])
      cohort_stage[cohort, day, , 2] <- cohort_stage[cohort, day-1, , 2] +
        (cohort_GDD[cohort, day-1, ] > params$stageDurCumulGDD_F[cohort_stage[cohort, day-1, , 2]])
      # accumulate GDD
      cohort_GDD[cohort, day,] <- cohort_GDD[cohort, day-1, ] + temp_mx[day, ]
    }
  }
}

matplot(cohort_GDD[1, , ], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))
matplot(cohort_GDD[100, , ], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))

matplot(cohort_N[1, , , 1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))
matplot(cohort_N[1, , , 2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))

matplot(cohort_N[200, , , 1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))
matplot(cohort_N[200, , , 2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))

matplot(cohort_stage[1, , , 1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))
matplot(cohort_stage[1, , , 2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))

matplot(cohort_stage[200, , , 1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))
matplot(cohort_stage[200, , , 2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1))

matplot(stage_survRate[, , 1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 1))
matplot(stage_survRate[, , 2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 1))
matplot(stage_survRate[, , 3], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 1))
matplot(stage_survRate[, , 4], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 1))


# . sum stages ------------------------------------------------------------

for(day in 1:info$nDays) {
  for(stage in 1:info$nStages) {
    for(sex in 1:2) {
      mu[day, , stage, sex] <- apply(cohort_N[, day, , sex] * (cohort_stage[, day, , sex] == stage),
                                     2, sum)
    }
  }
}

matplot(mu[, ,1,1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))
matplot(mu[, ,2,1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))
matplot(mu[, ,3,1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))
matplot(mu[, ,4,1], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))

matplot(mu[, ,1,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))
matplot(mu[, ,2,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))
matplot(mu[, ,3,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))
matplot(mu[, ,4,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, max(c(mu))))



# . detection -------------------------------------------------------------

for(day in 1:info$nDays) {
  for(stage in 1:info$nStages) {
    for(site in 1:info$nSites) {
      for(sex in 1:2) {
        y[day, site, stage, sex] <- rnbinom(20,
                                            mu=mu[day, site, stage, sex]/nFish_mx[day, site] *
                                              params$detect_p[stage],
                                            size=5) |> mean()
      }
    }
  }
}
y_wk <- y[seq(1, info$nDays, by=7), , ,]

y_chal <- y[, , 1, 1] + y[, , 1, 2]
y_prea <- y[, , 2, 1] + y[, , 2, 2]
y_adfe <- y[, , 3, 2] + y[, , 4, 2]
matplot(y_chal, lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))
matplot(y_prea, lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))
matplot(y_adfe, lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 15))
matplot(y[, ,3,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))
matplot(y[, ,4,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))


matplot(y_wk[, ,1,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))
matplot(y_wk[, ,2,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))
matplot(y_wk[, ,3,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))
matplot(y_wk[, ,4,2], lty=1, type="l", col=scico(5, palette="devon", end=0.8, direction=-1), ylim=c(0, 10))







# Hamre -------------------------------------------------------------------

tab4 <- read_csv("data/Hamre2019_Table4.csv") |>
  filter(PctLice=="all") |>
  select(Sex, Stage, starts_with("T")) |>
  pivot_longer(starts_with("T"), names_to="temperature", values_to="dpi") |>
  mutate(temperature=as.numeric(str_remove(temperature, "T")),
         Stage=factor(Stage, levels=c("Larvae", "Ch1", "Ch2", "Pa1", "Pa2", "Adult", "Egg_first", "Egg_refract")),
         GDD=dpi*temperature)
tab4 |>
  filter(temperature > 7) |>
  filter(!Stage %in% c("Larvae", "Egg_refract")) |>
  ggplot(aes(temperature, dpi, colour=Stage)) +
  geom_line() +
  scale_colour_viridis_d(option="turbo", begin=0.1) +
  facet_grid(Sex~.)

tab4 |>
  filter(temperature > 7) |>
  filter(!Stage %in% c("Larvae", "Egg_refract")) |>
  ggplot(aes(Stage, GDD, colour=temperature, group=temperature)) +
  scale_colour_viridis_c(option="plasma", end=0.9) +
  geom_line() +
  facet_grid(Sex~.)


rM <- function(temp, beta) {
  beta[2]*temp^2 + beta[3]*temp + beta[4]
}

MnM <- function(temp, beta, dpi) {
  beta[1] + rM(temp, beta) * dpi
}

beta_M <- c(-0.354753, 0.000677, 0.010294, 0.005729)
beta_F <- c(-0.152008, 0.000485, 0.008667, 0.003750)
plot(3:22, rM(3:22, beta_M), ylim=c(0, 0.6))
points(3:22, rM(3:22, beta_F), pch=19)

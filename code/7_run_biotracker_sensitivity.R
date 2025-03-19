# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Run biotracker simulations


# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(doFuture)
library(sf)
dirf("code/fn", ".R") |> walk(source)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
set.seed(101)


# define parameters -------------------------------------------------------

cores_per_sim <- 10
parallel_sims <- 3
start_date <- "2023-01-01"
end_date <- "2023-12-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hydro1="/home/sa04ts/hydro/etive28/Archive",
             hydro2="/home/sa04ts/hydro/WeStCOMS2/Archive",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v1-0-0.jar",
             out=glue("{getwd()}/out/sensitivity")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hydro1="E:/hydro/etive28/Archive",
               hydro2="E:/hydro/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker.jar",
               out=glue("D:/EtiveLice/out/sensitivity"))
)

post_dir <- glue("{dirs$proj}/data/lit/fit")
egg_post <- list(constant=tibble(b=rnorm(4000, 28.2, 2)),
                 logistic=read_csv(glue("{post_dir}/egg_temp_logistic_post.csv"), show_col_types=F),
                 linear=read_csv(glue("{post_dir}/egg_temp_linear_post.csv"), show_col_types=F)) |>
  map(~.x |> mutate(across(everything(), ~signif(.x, 3))) |>
        unite("all", everything(), sep=",") %>%
        .$all)
mort_post <- list(constant=tibble(b=plogis(rnorm(4000, qlogis(0.01), 0.25))),
                 logistic=read_csv(glue("{post_dir}/mort_sal_logistic_post.csv"), show_col_types=F)) |>
  map(~.x |> mutate(across(everything(), ~signif(.x, 3))) |>
        unite("all", everything(), sep=",") %>%
        .$all)
n_sim <- 999
swim_mx <- MASS::mvrnorm(n_sim, c(0,0), matrix(c(1, 0.6, 0.6, 1), nrow=2))
light_mx <- MASS::mvrnorm(n_sim, c(0,0), matrix(c(1, 0.6, 0.6, 1), nrow=2))
salMin_mx <- MASS::mvrnorm(n_sim, c(0,0), matrix(c(1, 0.6, 0.6, 1), nrow=2))
salMax_mx <- MASS::mvrnorm(n_sim, c(0,0), matrix(c(1, 0.6, 0.6, 1), nrow=2))
sim.i <- tibble(D_h=runif(n_sim, 0.05, 0.5),
                D_hVert=runif(n_sim, 0.0005, 0.005),
                mortSal_fn=sample(c("constant", "logistic"), n_sim, replace=T),
                eggTemp_fn=sample(c("constant", "logistic"), n_sim, replace=T),
                lightThreshCopepodid=qunif(pnorm(light_mx[,1]), (2e-6)^0.5, (2e-4)^0.5)^2,
                lightThreshNauplius=qunif(pnorm(light_mx[,2]), (0.05)^0.5, (0.5)^0.5)^2,
                swimUpSpeedMean=-(qunif(pnorm(swim_mx[,1]), (1e-4)^0.5, (2e-2)^0.5))^2,
                swimDownSpeedMean=(qunif(pnorm(swim_mx[,2]), (1e-4)^0.5, (2e-2)^0.5))^2,
                passiveSinkRateSal=sample(c(T, F), n_sim, replace=T),
                salinityThreshCopepodidMin=qunif(pnorm(salMin_mx[,1]), 20, 28),
                salinityThreshNaupliusMin=qunif(pnorm(salMin_mx[,2]), 20, 28),
                salinityThreshCopepodidMax=pmin(salinityThreshCopepodidMin + qunif(pnorm(salMax_mx[,1]), 0.1, 6), 32),
                salinityThreshNaupliusMax=pmin(salinityThreshNaupliusMin + qunif(pnorm(salMax_mx[,2]), 0.1, 6), 32),
                viableDegreeDays=runif(n_sim, 35, 45),
                connectivityThresh=runif(n_sim, 30, 250)) |>
  mutate(across(where(is.numeric), ~signif(.x, 5))) |>
  mutate(i=str_pad(row_number(), 3, "left", "0"),
         outDir=glue("{dirs$out}/sim_{i}/")) |>
  rowwise() |>
  mutate(eggTemp_b=sample(egg_post[[eggTemp_fn]], 1),
         mortSal_b=sample(mort_post[[mortSal_fn]], 1)) |>
  ungroup()
write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))
sim_seq <- 1:nrow(sim.i)
sim_seq <- sim_seq[-(1:188)]


# set properties ----------------------------------------------------------

walk(sim_seq, ~dir.create(sim.i$outDir[.x], showWarnings=F))

walk(sim_seq,
     ~set_biotracker_properties(
       # run settings
       properties_file_path=glue("{dirs$out}/sim_{sim.i$i[.x]}.properties"),
       parallelThreads=cores_per_sim,
       parallelThreadsHD=6,
       start_ymd=as.numeric(str_remove_all(start_date, "-")),
       numberOfDays=nDays,
       nparts=10,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh1=glue("{dirs$mesh}/etive28_mesh.nc"),
       mesh1Domain="etive28",
       datadir=glue("{dirs$hydro1}/"),
       mesh2=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       mesh2Domain="westcoms2",
       datadir2=glue("{dirs$hydro2}/"),
       # sites
       sitefile=glue("D:/EtiveLice/data/farm_sites_2023.csv"),
       sitefileEnd=glue("D:/EtiveLice/data/farm_sites_2023.csv"),
       siteDensityPath=glue("D:/EtiveLice/data/lice_daily_2023-01-01_2023-12-31_FARMS.csv"),
       # dynamics
       D_hVert=sim.i$D_hVert[.x],
       D_h=sim.i$D_h[.x],
       stepsPerStep=30,
       # biology
       fixDepth="false",
       startDepth=10,
       eggTemp_fn=sim.i$eggTemp_fn[.x],
       eggTemp_b=sim.i$eggTemp_b[.x],
       mortSal_fn=sim.i$mortSal_fn[.x],
       mortSal_b=sim.i$mortSal_b[.x],
       salinityThreshCopepodidMin=sim.i$salinityThreshCopepodidMin[.x],
       salinityThreshCopepodidMax=sim.i$salinityThreshCopepodidMax[.x],
       salinityThreshNaupliusMin=sim.i$salinityThreshNaupliusMin[.x],
       salinityThreshNaupliusMax=sim.i$salinityThreshNaupliusMax[.x],
       lightThreshCopepodid=sim.i$lightThreshCopepodid[.x],
       lightThreshNauplius=sim.i$lightThreshNauplius[.x],
       swimUpSpeedCopepodidMean=sim.i$swimUpSpeedMean[.x],
       swimUpSpeedCopepodidStd=abs(sim.i$swimUpSpeedMean[.x]/5),
       swimDownSpeedCopepodidMean=sim.i$swimDownSpeedMean[.x],
       swimDownSpeedCopepodidStd=sim.i$swimDownSpeedMean[.x]/5,
       swimUpSpeedNaupliusMean=sim.i$swimUpSpeedMean[.x]/2,
       swimUpSpeedNaupliusStd=abs(sim.i$swimUpSpeedMean[.x]/10),
       swimDownSpeedNaupliusMean=sim.i$swimDownSpeedMean[.x]/2,
       swimDownSpeedNaupliusStd=sim.i$swimDownSpeedMean[.x]/10,
       passiveSinkingIntercept=if_else(sim.i$passiveSinkRateSal[.x], 0.001527, sim.i$swimDownSpeedMean[.x]),
       passiveSinkingSlope=if_else(sim.i$passiveSinkRateSal[.x], -1.68e-5, 0),
       viableDegreeDays=sim.i$viableDegreeDays[.x],
       connectivityThresh=sim.i$connectivityThresh[.x],
       # recording
       verboseSetUp="true",
       recordConnectivity="true",
       connectivityInterval=24,
       connectDepth1_min=0,
       connectDepth1_max=5,
       connectDepth2_min=5,
       connectDepth2_max=20,
       recordVertDistr="false",
       recordElemActivity="false"))


# run simulations ---------------------------------------------------------

if(os=="linux") {
  plan(multicore, workers=parallel_sims)
} else {
  plan(multisession, workers=parallel_sims)
}
sim_sets <- split(sim_seq, rep(1:parallel_sims, length(sim_seq)/parallel_sims))
foreach(j=1:parallel_sims, .options.future=list(globals=structure(TRUE, add="sim.i"))) %dofuture% {
  for(i in sim_sets[[j]]) {
    setwd(dirs$proj)
    biotrackR::run_biotracker(
      jdk_path=dirs$jdk,
      jar_path=dirs$jar,
      f_properties=glue::glue("{dirs$out}/sim_{sim.i$i[i]}.properties"),
      sim_dir=glue::glue("{sim.i$outDir[i]}")
    )
  }
}
plan(sequential)


# process output ----------------------------------------------------------


sim.i <- read_csv(glue("{dirs$out}/sim_i.csv"))
farms_linnhe <- read_csv(glue("{dirs$proj}/data/farm_sites_2023.csv"))
farms_etive <- c("FFMC84", "FFMC32", "APT1", "SAR1", "FFMC27")

# Influx
mesh_fp <- st_read(glue("{dirs$mesh}/WeStCOMS2_meshFootprint.gpkg"))
sim_dirs <- dirf(dirs$out, "^sim_[0-9][0-9][0-9]$")
c_0_5_summary <- c_5_20_summary <- c_0_20_summary <- vector("list", length(sim_dirs))
c_0_5_etive <- c_5_20_etive <- c_0_20_etive <- vector("list", length(sim_dirs))

# for(i in 1:length(sim_dirs)) {
for(i in 1:188) {
  f_0_5 <- dirrf(sim_dirs[i], "connectivity_0.0-5.0.*csv")
  f_5_20 <- dirrf(sim_dirs[i], "connectivity_5.0-20.0.*csv")
  sim <- str_sub(sim_dirs[i], -3, -1)
  site_areas <- farms_linnhe |>
    st_as_sf(coords=c("easting", "northing"), crs=27700) |>
    st_buffer(dist=sim.i$connectivityThresh[i]) |>
    st_intersection(mesh_fp) %>%
    mutate(area=as.numeric(st_area(.))) |>
    st_drop_geometry()
  if(length(f_0_5) > 0) {
    c_0_5_i <- f_0_5 |>
      map_dfr(~load_connectivity(.x,
                                 source_names=farms_linnhe$sepaSite,
                                 dest_names=farms_linnhe$sepaSite,
                                 liceScale=1) |>
                mutate(sim=sim))
    c_0_5_i_daily <- calc_daily_fluxes(c_0_5_i, site_areas)
    c_0_5_summary[[i]] <- c_0_5_i_daily |> calc_sensitivity_outcomes(sim)
    c_0_5_etive[[i]] <- c_0_5_i_daily |>
      filter(sepaSite %in% farms_etive) |>
      group_by(sepaSite) |>
      calc_sensitivity_outcomes(sim)
  }
  if(length(f_5_20) > 0) {
    c_5_20_i <- f_5_20 |>
      map_dfr(~load_connectivity(.x,
                                 source_names=farms_linnhe$sepaSite,
                                 dest_names=farms_linnhe$sepaSite,
                                 liceScale=1) |>
                mutate(sim=sim))
    c_5_20_i_daily <- calc_daily_fluxes(c_5_20_i, site_areas)
    c_5_20_summary[[i]] <- c_5_20_i_daily |> calc_sensitivity_outcomes(sim)
    c_5_20_etive[[i]] <- c_5_20_i_daily |>
      filter(sepaSite %in% farms_etive) |>
      group_by(sepaSite) |>
      calc_sensitivity_outcomes(sim)
  }
  if(length(f_0_5) > 0 & length(f_5_20) > 0) {
    c_0_20_i <- bind_rows(c_0_5_i, c_5_20_i) |>
      group_by(source, destination, sim, date) |>
      summarise(value=sum(value))
    c_0_20_i_daily <- calc_daily_fluxes(c_0_20_i, site_areas)
    c_0_20_summary[[i]] <- c_0_20_i_daily |> calc_sensitivity_outcomes(sim)
    c_0_20_etive[[i]] <- c_0_20_i_daily |>
      filter(sepaSite %in% farms_etive) |>
      group_by(sepaSite) |>
      calc_sensitivity_outcomes(sim)
  }
}


c_0_5_sum_df <- make_sensitivity_sum_df(c_0_5_summary, sim.i)
saveRDS(c_0_5_sum_df, "out/sensitivity/c_0_5_sum_df.rds")

c_5_20_sum_df <- make_sensitivity_sum_df(c_5_20_summary, sim.i)
saveRDS(c_5_20_sum_df, "out/sensitivity/c_5_20_sum_df.rds")

c_0_20_sum_df <- make_sensitivity_sum_df(c_0_5_summary, sim.i)
saveRDS(c_0_20_sum_df, "out/sensitivity/c_0_20_sum_df.rds")

c_0_5_sum_etive_df <- make_sensitivity_sum_df(c_0_5_etive, sim.i)
saveRDS(c_0_5_sum_etive_df, "out/sensitivity/c_0_5_sum_etive_df.rds")

c_5_20_sum_etive_df <- make_sensitivity_sum_df(c_5_20_etive, sim.i)
saveRDS(c_5_20_sum_etive_df, "out/sensitivity/c_5_20_sum_etive_df.rds")

c_0_20_sum_etive_df <- make_sensitivity_sum_df(c_0_20_etive, sim.i)
saveRDS(c_0_20_sum_etive_df, "out/sensitivity/c_0_20_sum_etive_df.rds")


c_0_5_sum_df <- readRDS("out/sensitivity/c_0_5_sum_df.rds")
c_5_20_sum_df <- readRDS("out/sensitivity/c_5_20_sum_df.rds")
c_0_20_sum_df <- readRDS("out/sensitivity/c_0_20_sum_df.rds")

outcome_names <- grep("^N_|_m2", names(c_0_5_sum_df), value=T)

c_0_5_sum_df |>
  pivot_longer(any_of(names(sim.i)[1:15])) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")
c_5_20_sum_df |>
  pivot_longer(any_of(names(sim.i)[1:15])) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")
c_0_20_sum_df |>
  pivot_longer(any_of(names(sim.i)[1:15])) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")


c_0_5_sum_etive_df |>
  mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(names(sim.i)[1:15])) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_5_20_sum_etive_df |>
  mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(names(sim.i)[1:15])) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_0_20_sum_etive_df |>
  mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(names(sim.i)[1:15])) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")

c_0_20_sum_df |>
  ggplot(aes(salinityThreshCopepodidMin, salinityThreshCopepodidMax, colour=influx_m2_MAM_md^0.25)) +
  geom_point() +
  scale_colour_viridis_c(option="turbo", begin=0.1)

c_0_5_sum_df |>
  select(any_of(c("sim", names(sim.i)[1:15], outcome_names))) |>
  pivot_longer(any_of(names(sim.i)[1:15]), names_to="param", values_to="param_val") |>
  pivot_longer(any_of(outcome_names), names_to="outcome", values_to="outcome_val") |>
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
                          .default="all")) |>
  filter(type=="IP_m2", flux=="influx", season != "all") |>
  ggplot(aes(param_val, outcome_val^0.25, colour=season)) +
  geom_point(shape=1, alpha=0.5) +
  geom_line(stat="smooth", method="loess", se=F) +
  scale_colour_viridis_d(end=0.8) +
  ggtitle("0-5m") +
  facet_wrap(~param, scales="free_x")

c_5_20_sum_df |>
  select(any_of(c("sim", names(sim.i)[1:15], outcome_names))) |>
  pivot_longer(any_of(names(sim.i)[1:15]), names_to="param", values_to="param_val") |>
  pivot_longer(any_of(outcome_names), names_to="outcome", values_to="outcome_val") |>
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
                          .default="all")) |>
  filter(type=="IP_m2", flux=="influx", season != "all") |>
  ggplot(aes(param_val, outcome_val^0.25)) +
  geom_point(shape=1, alpha=0.25) +
  geom_line(stat="smooth", method="loess", se=F, colour="blue", linewidth=1) +
  scale_colour_viridis_d(end=0.8) +
  ggtitle("5-20m") +
  facet_wrap(~param, scales="free_x")

c_0_20_sum_df |>
  select(any_of(c("sim", names(sim.i)[1:15], outcome_names))) |>
  pivot_longer(any_of(names(sim.i)[1:15]), names_to="param", values_to="param_val") |>
  pivot_longer(any_of(outcome_names), names_to="outcome", values_to="outcome_val") |>
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
                          .default="all")) |>
  filter(type=="IP_m2", flux=="influx", season != "all") |>
  ggplot(aes(param_val, (outcome_val/20)^0.25)) +
  geom_point(shape=1, alpha=0.25) +
  geom_line(stat="smooth", method="loess", se=F, colour="blue", linewidth=1) +
  scale_colour_viridis_d(end=0.8) +
  labs(x="Parameter value",
       y="Mean IP among all farms: (cop / m3 / d)^0.25",
       title="0-20m") +
  facet_wrap(~param, scales="free_x", nrow=3)
ggsave("admin/project_meetings/figs_temp/sensitivity_0_20m_scatter.png",
       width=12, height=8, dpi=300)



library(randomForest)

rf_0_5_ls <- run_sensitivity_ML(outcome_names, c_0_5_sum_df, sim.i)
importance_0_5_df <- summarise_importance(rf_0_5_ls, outcome_names)

rf_5_20_ls <- run_sensitivity_ML(outcome_names, c_5_20_sum_df, sim.i)
importance_5_20_df <- summarise_importance(rf_5_20_ls, outcome_names)

rf_0_20_ls <- run_sensitivity_ML(outcome_names, c_0_20_sum_df, sim.i)
importance_0_20_df <- summarise_importance(rf_0_20_ls, outcome_names)

importance_df <- bind_rows(
  importance_0_5_df |> mutate(depth="0-5m"),
  importance_5_20_df |> mutate(depth="5-20m"),
  importance_0_20_df |> mutate(depth="0-20m"),
)


param_df <- tibble(param=c("D_h",
                           "D_hVert",
                           "passiveSinkRateSal",
                           "viableDegreeDays",
                           "salinityThreshNaupliusMax",
                           "salinityThreshNaupliusMin",
                           "salinityThreshCopepodidMax",
                           "salinityThreshCopepodidMin",
                           "connectivityThresh",
                           "swimUpSpeedMean",
                           "swimDownSpeedMean",
                           "lightThreshNauplius",
                           "lightThreshCopepodid",
                           "mortSal_fn",
                           "eggTemp_fn"),
                   pretty=c("Diffusion (horizontal)",
                            "Diffusion (vertical)",
                            "Passive sink rate",
                            "Nauplius to Copepodid GDD",
                            "Upper salinity thresh: Nauplius",
                            "Lower salinity thresh: Nauplius",
                            "Upper salinity thresh: Copepodid",
                            "Lower salinity thresh: Copepodid",
                            "Connectivity radius",
                            "Swim speed (up)",
                            "Swim speed (down)",
                            "Light threshold: Nauplius",
                            "Light threshold: Copepodid",
                            "Mortality fn",
                            "Egg production fn"))
importance_avg <- importance_df |>
  group_by(param) |>
  summarise(mn=mean(`%IncMSE`, na.rm=T)) |>
  arrange(mn) |>
  left_join(param_df)

importance_df |>
  filter(flux=="influx") |>
  filter(type=="IP_m2") |>
  left_join(param_df) |>
  mutate(pretty=factor(pretty, levels=importance_avg$pretty)) |>
  ggplot(aes(pretty, `%IncMSE`, fill=depth)) +
  geom_boxplot() +
  scale_fill_manual(values=c("grey80", "#DEE318FF", "#35608DFF")) +
  scale_x_discrete(labels=label_wrap_gen(23)) +
  labs(x="Parameter", y="Relative importance") +
  theme(panel.grid.major.x=element_line(colour="grey90"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="inside",
        legend.position.inside=c(0.1, 0.85))
ggsave("admin/project_meetings/figs_temp/sensitivity_depth.png", width=8, height=5.5, dpi=300)

importance_df |>
  filter(flux=="influx") |>
  filter(type=="IP_m2") |>
  filter(season != "all") |>
  left_join(param_df) |>
  mutate(pretty=factor(pretty, levels=importance_avg$pretty),
         season=factor(season,
                       levels=c("MAM", "JJA", "SON"),
                       labels=c("Spring: MAM", "Summer: JJA", "Autumn: SON"))) |>
  ggplot(aes(pretty, `%IncMSE`, fill=season)) +
  geom_boxplot() +
  scale_fill_viridis_d(option="turbo", begin=0.3, end=0.9) +
  # scale_fill_manual(values=c("white", "#DEE318FF", "#35608DFF")) +
  scale_x_discrete(labels=label_wrap_gen(23)) +
  labs(x="Parameter", y="Relative importance") +
  facet_grid(depth~.) +
  theme(panel.grid.major.x=element_line(colour="grey90"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

importance_df |>
  mutate(param=factor(param, levels=importance_avg$param)) |>
  ggplot(aes(`%IncMSE`, param, fill=depth)) +
  geom_boxplot() +
  facet_grid(.~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))

importance_df |>
  mutate(param=factor(param, levels=importance_avg$param)) |>
  ggplot(aes(`%IncMSE`, param, fill=season)) +
  geom_boxplot() +
  facet_grid(depth~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))




importance_0_5_df |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_0_5_df |>
  filter(flux == "influx") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))






importance_5_20_df |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_5_20_df |>
  filter(type != "IP total") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_5_20_df |>
  filter(type != "IP total") |>
  filter(flux == "influx") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))

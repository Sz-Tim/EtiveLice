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
set.seed(1003)


# define parameters -------------------------------------------------------

cores_per_sim <- 26
parallel_sims <- 6
start_date <- "2023-01-01"
end_date <- "2023-12-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hf0="/home/sa04ts/hydro/WeStCOMS2/Archive",
             hf2="/home/sa04ts/hydro/WeStCOMS2/SWAN/Archive_daily/",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v2-2-0.jar",
             out=glue("{getwd()}/out/sensitivity")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hf0="E:/hydro/WeStCOMS2/Archive",
               hf2="E:/hydro/WeStCOMS2/SWAN/Archive_daily/",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v2-2-0.jar",
               out=glue("D:/EtiveLice/out/sensitivity"))
)

n_sim <- 2000
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
sink_post <- read_csv(glue("{post_dir}/sink_sal_linear_post.csv"), show_col_types=F) |>
  mutate(across(everything(), ~.x/1e3)) #mm/s to m/s
sim.i <- sample_parameter_distributions(n_sim, dirs$out, "lhs", egg_post, mort_post, sink_post)
write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))
sim_seq <- 1:nrow(sim.i)
sim_seq <- sim_seq


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
       mesh0=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       hfDir0=glue("{dirs$hf0}/"),
       hfFilePrefix0="westcoms2",
       hfDir2="",
       hfFilePrefix2="",
       hfDirPrefix2="",
       hfDirSuffix2="",
       # sites
       sitefile=glue("{dirs$proj}/data/farm_sites_GSA_2023-2024.csv"),
       sitefileEnd=glue("{dirs$proj}/data/farm_sites_GSA_2023-2024.csv"),
       siteDensityPath=glue("{dirs$proj}/data/lice_daily_2023-01-01_2024-12-31_GSA.csv"),
       # dynamics
       variableDh=sim.i$variableDh[.x],
       variableDhV=sim.i$variableDhV[.x],
       D_h=sim.i$D_h[.x],
       D_hVert=sim.i$D_hVert[.x],
       stepsPerStep=30,
       stokesDrift="false",
       # biology
       fixDepth="false",
       startDepth=3,
       maxDepth=sim.i$maxDepth[.x],
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
       swimColdNauplius=sim.i$swimColdNauplius[.x],
       swimUpSpeedCopepodidMean=sim.i$swimUpSpeedCopepodidMean[.x],
       swimUpSpeedCopepodidStd=abs(sim.i$swimUpSpeedCopepodidMean[.x]/5),
       swimDownSpeedCopepodidMean=sim.i$swimDownSpeedCopepodidMean[.x],
       swimDownSpeedCopepodidStd=sim.i$swimDownSpeedCopepodidMean[.x]/5,
       swimUpSpeedNaupliusMean=sim.i$swimUpSpeedNaupliusMean[.x],
       swimUpSpeedNaupliusStd=abs(sim.i$swimUpSpeedNaupliusMean[.x]/5),
       swimDownSpeedNaupliusMean=sim.i$swimDownSpeedNaupliusMean[.x],
       swimDownSpeedNaupliusStd=sim.i$swimDownSpeedNaupliusMean[.x]/5,
       passiveSinkingIntercept=sim.i$passiveSinkInt[.x],
       passiveSinkingSlope=sim.i$passiveSinkSlope[.x],
       viableDegreeDays=sim.i$viableDegreeDays[.x],
       connectivityThresh=sim.i$connectivityThresh[.x],
       # recording
       verboseSetUp="true",
       recordConnectivity="true",
       connectImmature="true",
       connectivityInterval=24,
       connectDepth1_min=0,
       connectDepth1_max=5,
       connectDepth2_min=5,
       connectDepth2_max=15,
       connectDepth3_min=15,
       connectDepth3_max=30,
       recordVertDistr="false",
       recordPsteps="true",
       recordImmature="true",
       pstepsInterval=730,
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
farms_GSA <- read_csv(glue("{dirs$proj}/data/farm_sites_GSA_2023-2024.csv"))
farms_etive <- c("FFMC84", "FFMC32", "APT1", "SAR1", "FFMC27")

# Influx
mesh_fp <- st_read(glue("{dirs$mesh}/WeStCOMS2_meshFootprint.gpkg"))
sim_dirs <- dirf(dirs$out, "^sim_[0-9][0-9][0-9][0-9]$")
c_0_5_summary <- c_5_15_summary <- c_15_30_summary <- c_0_30_summary <- vector("list", length(sim_dirs))
c_0_5_farm <- c_5_15_farm <- c_15_30_farm <- c_0_30_farm <- vector("list", length(sim_dirs))

cores <- 30
for(i in 1:length(sim_dirs)) {
  f_0_5 <- dirrf(sim_dirs[i], "connectivity_0.0-5.0.*csv")
  f_5_15 <- dirrf(sim_dirs[i], "connectivity_5.0-15.0.*csv")
  f_15_30 <- dirrf(sim_dirs[i], "connectivity_15.0-30.0.*csv")
  sim <- str_sub(sim_dirs[i], -4, -1)
  site_areas <- farms_GSA |>
    st_as_sf(coords=c("easting", "northing"), crs=27700) |>
    st_buffer(dist=sim.i$connectivityThresh[i]) |>
    st_intersection(mesh_fp) %>%
    mutate(area=as.numeric(st_area(.))) |>
    st_drop_geometry()
  if(length(f_0_5) > 0) {
    c_0_5 <- calc_GSA_connectivity(f_0_5, farms_GSA$sepaSite, sim, site_areas, cores)
    c_0_5_i <- c_0_5$og
    c_0_5_i_daily <- c_0_5$daily
    c_0_5_summary[[i]] <- c_0_5$summary
    c_0_5_farm[[i]] <- c_0_5$farm
  } else {
    c_0_5_i <- NULL
  }
  if(length(f_5_15) > 0) {
    c_5_15 <- calc_GSA_connectivity(f_5_15, farms_GSA$sepaSite, sim, site_areas, cores)
    c_5_15_i <- c_5_15$og
    c_5_15_i_daily <- c_5_15$daily
    c_5_15_summary[[i]] <- c_5_15$summary
    c_5_15_farm[[i]] <- c_5_15$farm
  } else {
    c_5_15_i <- NULL
  }
  if(length(f_15_30) > 0) {
    c_15_30 <- calc_GSA_connectivity(f_15_30, farms_GSA$sepaSite, sim, site_areas, cores)
    c_15_30_i <- c_15_30$og
    c_15_30_i_daily <- c_15_30$daily
    c_15_30_summary[[i]] <- c_15_30$summary
    c_15_30_farm[[i]] <- c_15_30$farm
  } else {
    c_15_30_i <- NULL
  }
  if(length(f_0_5) > 0 | length(f_5_15) > 0 | length(f_15_30) > 0) {
    c_0_30_i <- bind_rows(c_0_5_i, c_5_15_i, c_15_30_i) |>
      group_by(source, destination, sim, date) |>
      summarise(value=sum(value))
    c_0_30_i_daily <- calc_daily_fluxes(c_0_30_i, site_areas)
    c_0_30_summary[[i]] <- c_0_30_i_daily |> calc_sensitivity_outcomes(sim)
    c_0_30_farm[[i]] <- c_0_30_i_daily |>
      group_by(sepaSite) |>
      calc_sensitivity_outcomes(sim)
  }
  cat("Finished", i, "\n")
}


c_0_5_sum_df <- make_sensitivity_sum_df(c_0_5_summary, sim.i)
saveRDS(c_0_5_sum_df, "out/sensitivity/processed/c_0_5_sum_df.rds")

c_5_15_sum_df <- make_sensitivity_sum_df(c_5_15_summary, sim.i)
saveRDS(c_5_15_sum_df, "out/sensitivity/processed/c_5_15_sum_df.rds")

c_15_30_sum_df <- make_sensitivity_sum_df(c_15_30_summary, sim.i)
saveRDS(c_15_30_sum_df, "out/sensitivity/processed/c_15_30_sum_df.rds")

c_0_30_sum_df <- make_sensitivity_sum_df(c_0_30_summary, sim.i)
saveRDS(c_0_30_sum_df, "out/sensitivity/processed/c_0_30_sum_df.rds")

c_0_5_sum_farm_df <- make_sensitivity_sum_df(c_0_5_farm, sim.i)
saveRDS(c_0_5_sum_farm_df, "out/sensitivity/processed/c_0_5_sum_farm_df.rds")

c_5_15_sum_farm_df <- make_sensitivity_sum_df(c_5_15_farm, sim.i)
saveRDS(c_5_15_sum_farm_df, "out/sensitivity/processed/c_5_15_sum_farm_df.rds")

c_15_30_sum_farm_df <- make_sensitivity_sum_df(c_15_30_farm, sim.i)
saveRDS(c_15_30_sum_farm_df, "out/sensitivity/processed/c_15_30_sum_farm_df.rds")

c_0_30_sum_farm_df <- make_sensitivity_sum_df(c_0_30_farm, sim.i)
saveRDS(c_0_30_sum_farm_df, "out/sensitivity/processed/c_0_30_sum_farm_df.rds")




# plots -------------------------------------------------------------------

logCols <- c("D_h", "D_hVert", "maxDepth")
rtCols <- c("lightThreshCopepodid", "lightThreshNauplius",
            "swimDownSpeedCopepodidMean", "swimDownSpeedNaupliusMean",
            "swimUpSpeedCopepodidMean", "swimUpSpeedNaupliusMean")

c_0_5_sum_df <- readRDS("out/sensitivity/processed/c_0_5_sum_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))
c_5_15_sum_df <- readRDS("out/sensitivity/processed/c_5_15_sum_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))
c_15_30_sum_df <- readRDS("out/sensitivity/processed/c_15_30_sum_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))
c_0_30_sum_df <- readRDS("out/sensitivity/processed/c_0_30_sum_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))

c_0_5_sum_farm_df <- readRDS("out/sensitivity/processed/c_0_5_sum_farm_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))
c_5_15_sum_farm_df <- readRDS("out/sensitivity/processed/c_5_15_sum_farm_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))
c_15_30_sum_farm_df <- readRDS("out/sensitivity/processed/c_15_30_sum_farm_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))
c_0_30_sum_farm_df <- readRDS("out/sensitivity/processed/c_0_30_sum_farm_df.rds") |>
  mutate(across(all_of(logCols), log),
         across(all_of(rtCols), sqrt)) |>
  mutate(D_h=if_else(variableDh, NA_real_, D_h),
         D_hVert=if_else(variableDhV, NA_real_, D_hVert))


outcome_names <- grep("^N_|_m2", names(c_0_5_sum_df), value=T)
param_names <- names(sim.i)[1:21]


c_0_5_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")
c_5_15_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")
c_15_30_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")
c_0_30_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x")

# How does relative importance vary among farms?
# Make maps
# What location-level characteristics might explain the variation?
c_0_5_sum_farm_df |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_5_15_sum_farm_df |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_5_15_sum_farm_df |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_15_30_sum_farm_df |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_0_30_sum_farm_df |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")

c_0_30_sum_df |>
  ggplot(aes(salinityThreshCopepodidMin, salinityThreshCopepodidMax, colour=influx_m2_MAM_md^0.25)) +
  geom_point() +
  scale_colour_viridis_c(option="turbo", begin=0.1)

c_0_5_sum_df |>
  select(any_of(c("sim", param_names, outcome_names))) |>
  pivot_longer(any_of(param_names), names_to="param", values_to="param_val") |>
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
  # geom_line(stat="smooth", method="loess", se=F) +
  scale_colour_viridis_d(end=0.8) +
  ggtitle("0-5m") +
  facet_wrap(~param, scales="free_x")

c_5_15_sum_df |>
  select(any_of(c("sim", param_names, outcome_names))) |>
  pivot_longer(any_of(param_names), names_to="param", values_to="param_val") |>
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
  # geom_line(stat="smooth", method="loess", se=F, colour="blue", linewidth=1) +
  scale_colour_viridis_d(end=0.8) +
  ggtitle("5-15m") +
  facet_wrap(~param, scales="free_x")

c_15_30_sum_df |>
  select(any_of(c("sim", param_names, outcome_names))) |>
  pivot_longer(any_of(param_names), names_to="param", values_to="param_val") |>
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
  # geom_line(stat="smooth", method="loess", se=F, colour="blue", linewidth=1) +
  scale_colour_viridis_d(end=0.8) +
  ggtitle("15-30m") +
  facet_wrap(~param, scales="free_x")

c_0_30_sum_df |>
  select(any_of(c("sim", param_names, outcome_names))) |>
  pivot_longer(any_of(param_names), names_to="param", values_to="param_val") |>
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
  ggplot(aes(param_val, (outcome_val/15)^0.25)) +
  geom_point(shape=1, alpha=0.25) +
  # geom_line(stat="smooth", method="loess", se=F, colour="blue", linewidth=1) +
  scale_colour_viridis_d(end=0.8) +
  labs(x="Parameter value",
       y="Mean IP among all farms: (cop / m3 / d)^0.25",
       title="0-30m") +
  facet_wrap(~param, scales="free_x", nrow=4)
ggsave("admin/project_meetings/figs_temp/sensitivity_0_30m_scatter.png",
       width=12, height=8, dpi=300)



library(randomForest)

rf_0_5_ls <- run_sensitivity_ML(outcome_names, c_0_5_sum_df, sim.i)
importance_0_5_df <- summarise_importance(rf_0_5_ls, outcome_names)

rf_5_15_ls <- run_sensitivity_ML(outcome_names, c_5_15_sum_df, sim.i)
importance_5_15_df <- summarise_importance(rf_5_15_ls, outcome_names)

rf_15_30_ls <- run_sensitivity_ML(outcome_names, c_15_30_sum_df, sim.i)
importance_15_30_df <- summarise_importance(rf_15_30_ls, outcome_names)

rf_0_30_ls <- run_sensitivity_ML(outcome_names, c_0_30_sum_df, sim.i)
importance_0_30_df <- summarise_importance(rf_0_30_ls, outcome_names)

importance_df <- bind_rows(
  importance_0_5_df |> mutate(depth="0-5m"),
  importance_5_15_df |> mutate(depth="5-15m"),
  importance_15_30_df |> mutate(depth="15-30m"),
  importance_0_30_df |> mutate(depth="0-30m")
) |>
  mutate(depth=factor(depth, levels=c("0-5m", "5-15m", "15-30m", "0-30m")))


param_df <- tibble(param=c("D_h",
                           "D_hVert",
                           "passiveSinkRateSal",
                           "viableDegreeDays",
                           "salinityThreshNaupliusMax",
                           "salinityThreshNaupliusMin",
                           "salinityThreshCopepodidMax",
                           "salinityThreshCopepodidMin",
                           "connectivityThresh",
                           "swimUpSpeedNaupliusMean",
                           "swimDownSpeedNaupliusMean",
                           "swimUpSpeedCopepodidMean",
                           "swimDownSpeedCopepodidMean",
                           "lightThreshNauplius",
                           "lightThreshCopepodid",
                           "mortSal_fn",
                           "eggTemp_fn",
                           "maxDepth",
                           "stokesDrift"),
                   pretty=c("Diffusion (horizontal)",
                            "Diffusion (vertical)",
                            "Passive sink rate",
                            "Nauplius to Copepodid GDD",
                            "Upper salinity thresh: Nauplius",
                            "Lower salinity thresh: Nauplius",
                            "Upper salinity thresh: Copepodid",
                            "Lower salinity thresh: Copepodid",
                            "Connectivity radius",
                            "Vert speed Nauplius (up)",
                            "Vert speed Nauplius (down)",
                            "Vert speed Copepodid (up)",
                            "Vert speed Copepodid (down)",
                            "Light threshold: Nauplius",
                            "Light threshold: Copepodid",
                            "Mortality fn",
                            "Egg production fn",
                            "Max depth preferred",
                            "Stokes drift"))
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
  scale_fill_manual(values=c("#FDE725FF", "#21908CFF", "#440154FF", "grey80")) +
  scale_x_discrete(labels=label_wrap_gen(23)) +
  labs(x="Parameter", y="Relative importance for influx") +
  theme(panel.grid.major.x=element_line(colour="grey90"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="inside",
        legend.position.inside=c(0.1, 0.85))
ggsave("admin/project_meetings/figs_temp/sensitivity_depth.png", width=8, height=5.5, dpi=300)


importance_df |>
  filter(flux=="influx") |>
  filter(type=="IP_m2") |>
  left_join(param_df) |>
  mutate(pretty=factor(pretty, levels=importance_avg$pretty)) |>
  ggplot(aes(pretty, `%IncMSE`, colour=depth, linetype=mn_md, shape=mn_md)) +
  geom_point() +
  geom_line(aes(group=paste(depth, mn_md))) +
  scale_shape_manual(values=c(1, 19)) +
  scale_linetype_manual(values=c(2, 1)) +
  scale_colour_manual(values=c("#FDE725FF", "#21908CFF", "#440154FF", "grey")) +
  scale_x_discrete(labels=label_wrap_gen(23)) +
  labs(x="Parameter", y="Relative importance for influx") +
  facet_grid(season~.) +
  theme(panel.grid.major.x=element_line(colour="grey90"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="inside",
        legend.position.inside=c(0.1, 0.85),
        legend.box="horizontal")



importance_df |>
  filter(flux=="influx") |>
  left_join(param_df) |>
  mutate(pretty=factor(pretty, levels=importance_avg$pretty)) |>
  ggplot(aes(pretty, `%IncMSE`, fill=depth)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#FDE725FF", "#21908CFF", "#440154FF", "grey80")) +
  scale_x_discrete(labels=label_wrap_gen(23)) +
  labs(x="Parameter", y="Relative importance") +
  facet_grid(type~.) +
  theme(panel.grid.major.x=element_line(colour="grey90"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="inside",
        legend.position.inside=c(0.1, 0.85))





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
  scale_fill_manual(values=c("#FDE725FF", "#21908CFF", "#440154FF", "grey80")) +
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


importance_5_15_df |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_5_15_df |>
  filter(type != "IP total") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_5_15_df |>
  filter(type != "IP total") |>
  filter(flux == "influx") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))


importance_15_30_df |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_15_30_df |>
  filter(type != "IP total") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))
importance_15_30_df |>
  filter(type != "IP total") |>
  filter(flux == "influx") |>
  ggplot(aes(`%IncMSE`, param)) +
  geom_boxplot() +
  geom_point(alpha=0.5) +
  facet_grid(flux~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))

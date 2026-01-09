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
               out=glue("{getwd()}/out/sensitivity"))
)

n_sim <- 2000
post_dir <- glue("{dirs$proj}/data/lit/fit")
egg_post <- list(constant=read_csv(glue("{post_dir}/egg_temp_constant_post.csv"), show_col_type=F),
                 logistic=read_csv(glue("{post_dir}/egg_temp_logistic_post.csv"), show_col_types=F),
                 linear=read_csv(glue("{post_dir}/egg_temp_linear_post.csv"), show_col_types=F)) |>
  map(~.x |> mutate(across(everything(), ~signif(.x, 3))) |>
        unite("all", everything(), sep=",") %>%
        .$all)
mort_post <- list(constant=read_csv(glue("{post_dir}/mort_sal_constant_post.csv"), show_col_types=F),
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
       parallelThreadsHD=8,
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
       siteDensityPath=glue("{dirs$proj}/data/lice_daily_2023-01-01_2024-12-31_GSA_05lpf_maxB.csv"),
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
       pstepsInterval=365,
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
# sim_dirs <- dirf(dirs$out, "^sim_[0-9][0-9][0-9][0-9]$")
sim_dirs <- dirrf(dirs$out, "pstepsImmature_20231231") |> dirname() |> dirname()
c_0_5_summary <- c_5_15_summary <- c_15_30_summary <- c_0_30_summary <- vector("list", length(sim_dirs))
c_0_5_farm <- c_5_15_farm <- c_15_30_farm <- c_0_30_farm <- vector("list", length(sim_dirs))

cores <- 20
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

farm_meta <- read_csv("data/farm_GSA_metadata.csv")
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
  facet_wrap(~name, scales="free_x", ncol=7)
c_5_15_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x", ncol=7)
c_15_30_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x", ncol=7)
c_0_30_sum_df |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  facet_wrap(~name, scales="free_x", ncol=7)

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
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_MAM_md^0.25, colour=sepaSite)) +
  geom_point(shape=1, alpha=0.5) +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=-1) +
  facet_wrap(~name, scales="free_x")
c_0_30_sum_farm_df |>
  inner_join(farm_meta, by="sepaSite") |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(value, influx_m2_SON_md^0.25, colour=logProximity)) +
  geom_point(shape=1, alpha=0.5) +
  scale_colour_viridis_c(option="turbo", begin=0.05) +
  facet_wrap(~name, scales="free_x", ncol=7) +
  theme(legend.position="bottom")

c_0_30_sum_farm_df |>
  inner_join(farm_meta, by="sepaSite") |>
  mutate(across(any_of(param_names), ~c(scale(.x)))) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(logProximity, influx_m2_SON_md^0.25, colour=value)) +
  geom_point(shape=1, alpha=0.5) +
  scale_colour_viridis_c(option="turbo", begin=0.05) +
  facet_wrap(~name, scales="free_x", ncol=7) +
  theme(legend.position="bottom")


c_0_30_sum_farm_df |>
  inner_join(farm_meta, by="sepaSite") |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(logProximity, value, colour=influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scale_colour_viridis_c(option="turbo", begin=0.05) +
  facet_wrap(~name, scales="free_y", ncol=7)
c_0_30_sum_farm_df |>
  inner_join(farm_meta, by="sepaSite") |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(logProximity, value, colour=influx_m2_SON_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scale_colour_viridis_c(option="turbo", begin=0.05) +
  facet_wrap(~name, scales="free_y", ncol=7)

c_0_30_sum_farm_df |>
  inner_join(farm_meta, by="sepaSite") |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(fetch, value, colour=influx_m2_MAM_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scale_colour_viridis_c(option="turbo", begin=0.05) +
  facet_wrap(~name, scales="free_y", ncol=7)
c_0_30_sum_farm_df |>
  inner_join(farm_meta, by="sepaSite") |>
  # mutate(sepaSite=factor(sepaSite, levels=farms_etive)) |>
  pivot_longer(any_of(param_names)) |>
  ggplot(aes(fetch, value, colour=influx_m2_SON_md^0.25)) +
  geom_point(shape=1, alpha=0.5) +
  # geom_smooth(method="loess", se=F) +
  scale_colour_viridis_c(option="turbo", begin=0.05) +
  facet_wrap(~name, scales="free_y", ncol=7)

farm_bbox <- list(xmin=125000, xmax=225500, ymin=690000, ymax=785000)
linnhe_fp <- mesh_fp |> st_crop(unlist(farm_bbox))
pA <- c_0_30_sum_farm_df |>
  group_by(sepaSite) |>
  summarise(influx_MAM_m2=mean(influx_m2_MAM_md^0.25)^4) |>
  left_join(farms_GSA) |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_point(aes(easting, northing, colour=influx_MAM_m2)) +
  scale_colour_viridis_c(option="turbo", begin=0.05)
pB <- c_0_30_sum_farm_df |>
  group_by(sepaSite) |>
  summarise(sd_influx_MAM_m2=sd(influx_m2_MAM_md^0.25)^4) |>
  left_join(farms_GSA) |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_point(aes(easting, northing, colour=sd_influx_MAM_m2)) +
  scale_colour_viridis_c(option="viridis", end=0.95)
pC <- c_0_30_sum_farm_df |>
  group_by(sepaSite) |>
  summarise(influx_SON_m2=mean(influx_m2_SON_md^0.25, na.rm=T)^4) |>
  left_join(farms_GSA) |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_point(aes(easting, northing, colour=influx_SON_m2)) +
  scale_colour_viridis_c(option="turbo", begin=0.05)
pD <- c_0_30_sum_farm_df |>
  group_by(sepaSite) |>
  summarise(sd_influx_SON_m2=sd(influx_m2_SON_md^0.25, na.rm=T)^4) |>
  left_join(farms_GSA) |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_point(aes(easting, northing, colour=sd_influx_SON_m2)) +
  scale_colour_viridis_c(option="viridis", end=0.95)
cowplot::plot_grid(pA, pB, pC, pD, nrow=2, align="hv", axis="tblr")



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



# random forest -----------------------------------------------------------

library(randomForest)

rf_0_5_ls <- run_sensitivity_ML(outcome_names, param_names,
                                c_0_5_sum_df |>
                                  mutate(D_h=if_else(variableDh, max(D_h)+5, D_h),
                                         D_hVert=if_else(variableDhV, max(D_hVert)+5, D_hVert)),
                                sim.i)
importance_0_5_df <- summarise_importance(rf_0_5_ls, outcome_names)

rf_5_15_ls <- run_sensitivity_ML(outcome_names, param_names,
                                 c_5_15_sum_df |>
                                   mutate(D_h=if_else(variableDh, max(D_h)+5, D_h),
                                          D_hVert=if_else(variableDhV, max(D_hVert)+5, D_hVert)),
                                 sim.i)
importance_5_15_df <- summarise_importance(rf_5_15_ls, outcome_names)

rf_15_30_ls <- run_sensitivity_ML(outcome_names, param_names,
                                  c_15_30_sum_df |>
                                    mutate(D_h=if_else(variableDh, max(D_h)+5, D_h),
                                           D_hVert=if_else(variableDhV, max(D_hVert)+5, D_hVert)),
                                  sim.i)
importance_15_30_df <- summarise_importance(rf_15_30_ls, outcome_names)

rf_0_30_ls <- run_sensitivity_ML(outcome_names, param_names,
                                 c_0_30_sum_df |>
                                   mutate(D_h=if_else(variableDh, max(D_h)+5, D_h),
                                          D_hVert=if_else(variableDhV, max(D_hVert)+5, D_hVert)),
                                 sim.i)
importance_0_30_df <- summarise_importance(rf_0_30_ls, outcome_names)

importance_df <- bind_rows(
  importance_0_5_df |> mutate(depth="0-5m"),
  importance_5_15_df |> mutate(depth="5-15m"),
  importance_15_30_df |> mutate(depth="15-30m"),
  importance_0_30_df |> mutate(depth="0-30m")
) |>
  mutate(depth=factor(depth, levels=c("0-5m", "5-15m", "15-30m", "0-30m")))

param_df <- list("D_h"="Diffusion: horizontal",
                   "D_hVert"="Diffusion: vertical",
                   "variableDh"="WeStCOMS horizontal diffusion",
                   "variableDhV"="WeStCOMS vertical diffusion",
                   "eggTemp_fn"="Gravid egg production",
                   "mortSal_fn"="Larval mortality rate",
                   "swimUpSpeedNaupliusMean"="Upward swimming speed: N",
                   "swimDownSpeedNaupliusMean"="Downward sinking speed: N",
                   "swimUpSpeedCopepodidMean"="Upward swimming speed: C",
                   "swimDownSpeedCopepodidMean"="Downward sinking speed: C",
                   "passiveSinkRateSal"="Salinity-dependent sink rate",
                   "salinityThreshNaupliusMin"="Salinity: Lower threshold: N",
                   "salinityThreshNaupliusMax"="Salinity: Upper threshold: N",
                   "salinityThreshCopepodidMin"="Salinity: Lower threshold: C",
                   "salinityThreshCopepodidMax"="Salinity: Upper threshold: C",
                   "salRangeN"="Salinity range: N",
                   "salRangeC"="Salinity range: C",
                   "lightThreshNauplius"="Light threshold: N",
                   "lightThreshCopepodid"="Light threshold: C",
                   "swimColdNauplius"="Cold preference: N",
                   "viableDegreeDays"="Development: N",
                   "maxDepth"="Maximum preferred depth",
                   "connectivityThresh"="IP radius") |>
  as_tibble() |>
  pivot_longer(everything(), names_to="param", values_to="param_pretty") |>
  mutate(param_order=factor(param,
                          levels=c("D_h", "D_hVert", "variableDh", "variableDhV",
                                   "maxDepth", "connectivityThresh", "passiveSinkRateSal",
                                   "eggTemp_fn", "mortSal_fn", "viableDegreeDays",
                                   "swimColdNauplius",
                                   "salinityThreshNaupliusMin", "salinityThreshNaupliusMax",
                                   "salinityThreshCopepodidMin", "salinityThreshCopepodidMax",
                                   "salRangeN", "salRangeC",
                                   "swimDownSpeedNaupliusMean", "swimUpSpeedNaupliusMean",
                                   "swimDownSpeedCopepodidMean", "swimUpSpeedCopepodidMean",
                                   "lightThreshNauplius", "lightThreshCopepodid"))) |>
  arrange(param_order)

importance_avg <- importance_df |>
  group_by(param) |>
  summarise(mn=mean(`%IncMSE`, na.rm=T)) |>
  arrange(mn) |>
  left_join(param_df)

importance_df |>
  filter(flux=="influx") |>
  filter(type=="IP_m2") |>
  left_join(param_df) |>
  mutate(pretty=factor(param_pretty, levels=importance_avg$param_pretty)) |>
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
  mutate(pretty=factor(param_pretty, levels=importance_avg$param_pretty)) |>
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
  mutate(pretty=factor(param_pretty, levels=importance_avg$param_pretty)) |>
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
  mutate(pretty=factor(param_pretty, levels=importance_avg$param_pretty),
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
  mutate(param=factor(param_pretty, levels=importance_avg$param_pretty)) |>
  ggplot(aes(`%IncMSE`, param, fill=depth)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#FDE725FF", "#21908CFF", "#440154FF", "grey80")) +
  facet_grid(.~type) +
  theme(panel.grid.major.y=element_line(colour="grey90"))

importance_df |>
  mutate(param=factor(param_pretty, levels=importance_avg$param_pretty)) |>
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




library(DALEX)
library(ingredients)
explainers <- map(1:length(rf_0_30_ls$rf),
                  ~explain(rf_0_30_ls$rf[[.x]],
                           data=rf_0_30_ls$exp$x_valid[[.x]],
                           y=rf_0_30_ls$exp$y_valid[[.x]]))
LDs <- map(explainers[grep("influx_m2", outcome_names)],
            ~conditional_dependence(.x))
map_dfr(seq_along(LDs),
        ~LDs[[.x]] |>
          mutate(id=.x,
                 outcome=outcome_names[grep("influx_m2", outcome_names)][.x])) |>
  mutate(season=case_when(grepl("MAM", outcome) ~ "MAM",
                          grepl("JJA", outcome) ~ "JJA",
                          grepl("SON", outcome) ~ "SON",
                          .default="all")) |>
  mutate(x=case_when(`_vname_` %in% logCols ~ exp(`_x_`),
                     `_vname_` %in% rtCols ~ (`_x_`)^2,
                     .default=`_x_`)) |>
  ggplot(aes(`_x_`, `_yhat_`, group=id, colour=season)) +
  geom_line(alpha=0.5) +
  scale_colour_manual(values=c("black", viridis::turbo(3, begin=0.3, end=0.9))) +
  facet_wrap(~`_vname_`, scales="free_x", ncol=7)

ALEs <- map(explainers[grep("influx_m2", outcome_names)],
            ~accumulated_dependence(.x))
map_dfr(seq_along(ALEs),
        ~ALEs[[.x]] |>
          mutate(id=.x,
                 outcome=outcome_names[grep("influx_m2", outcome_names)][.x])) |>
  mutate(season=case_when(grepl("MAM", outcome) ~ "MAM",
                          grepl("JJA", outcome) ~ "JJA",
                          grepl("SON", outcome) ~ "SON",
                          .default="all")) |>
  mutate(x=case_when(`_vname_` %in% logCols ~ exp(`_x_`),
                     `_vname_` %in% rtCols ~ (`_x_`)^2,
                     .default=`_x_`)) |>
  ggplot(aes(`_x_`, `_yhat_`, group=id, colour=season)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_line(alpha=0.5) +
  scale_colour_manual(values=c("black", viridis::turbo(3, begin=0.3, end=0.9))) +
  facet_wrap(~`_vname_`, scales="free_x", ncol=7)

FIs <- map(explainers[grep("influx_m2", outcome_names)],
           ~feature_importance(.x))
map_dfr(seq_along(FIs),
        ~FIs[[.x]] |>
          as_tibble() |>
          group_by(permutation) |>
          mutate(rel_loss=(dropout_loss - min(dropout_loss))/(max(dropout_loss) - min(dropout_loss))) |>
          group_by(variable) |>
          summarise(mean_dropout_loss=mean(rel_loss)) |>
          ungroup() |>
          mutate(id=.x,
                 outcome=outcome_names[grep("influx_m2", outcome_names)][.x])) |>
  mutate(season=case_when(grepl("MAM", outcome) ~ "MAM",
                          grepl("JJA", outcome) ~ "JJA",
                          grepl("SON", outcome) ~ "SON",
                          .default="all"),
         metric=if_else(grepl("_mn", outcome), "mean", "median")) |>
  filter(! variable %in% c("_baseline_", "_full_model_", "variableDh", "variableDhV")) |>
  group_by(variable) |>
  mutate(grand_mean=mean(mean_dropout_loss)) |>
  ungroup() |>
  arrange(grand_mean) |>
  mutate(variable=factor(variable, levels=unique(variable))) |>
  ggplot(aes(mean_dropout_loss, variable, colour=season)) +
  # geom_bar(stat="identity", position="dodge", colour="grey30") +
  geom_path(aes(group=season)) +
  geom_point() +
  scale_colour_manual(values=c("grey", viridis::turbo(3, begin=0.3, end=0.9))) +
  facet_grid(.~metric) +
  theme(panel.grid.major.y=element_line(colour="grey90", linewidth=0.15))

# maps --------------------------------------------------------------------

farms_GSA <- read_csv(glue("{dirs$proj}/data/farm_sites_GSA_2023-2024.csv"))
sim.i <- read_csv(glue("{dirs$out}/sim_i.csv"))
mesh_fp <- st_read(glue("{dirs$mesh}/WeStCOMS2_meshFootprint.gpkg"))
mesh_sf <- st_read(glue("{dirs$mesh}/WeStCOMS2_mesh.gpkg"))
farm_bbox <- list(xmin=125000, xmax=225500, ymin=690000, ymax=785000)
linnhe_fp <- mesh_fp |> st_crop(unlist(farm_bbox))
linnhe_sf <- mesh_sf |> st_crop(unlist(farm_bbox))
sim_dirs <- dirf(dirs$out, "^sim_[0-9][0-9][0-9][0-9]$")

ip_ls <- vector("list", length(sim_dirs))

library(furrr)
if(get_os()=="windows") {
  plan(multisession, workers=12)
} else {
  plan(multicore, workers=12)
}
# TODO: This is obviously too unwieldy for a large number of simulations.
# Read by time step and store temporary files
# Also worth using data.table
for(i in (1:length(sim_dirs))[25:50]) {
  sim <- str_sub(sim_dirs[i], -4, -1)
  fN <- dirrf(sim_dirs[i], "pstepsImmature")
  fC <- dirrf(sim_dirs[i], "pstepsMature")
  ip_ls[[i]] <- full_join(
    map_dfr(fN, ~load_psteps(.x, liceScale=1/730) |>
              rename_with(~"ipN_h", starts_with("t_")) |>
              mutate(date=ymd(str_split_fixed(basename(.x), "_", 3)[,2]))),
    map_dfr(fC, ~load_psteps(.x, liceScale=1/730) |>
              rename_with(~"ipC_h", starts_with("t_")) |>
              mutate(date=ymd(str_split_fixed(basename(.x), "_", 3)[,2]))),
    by=join_by(i, date)) |>
    mutate(sim=sim) |>
    filter(i %in% linnhe_sf$i)
}
plan(sequential)

ip_df <- data.table::rbindlist(ip_ls) |> as_tibble()

# ip_sf <- inner_join(linnhe_sf |> select(i, area, depth, geom),
#            ip_df,
#            by="i") |>
#   mutate(vol=(area*pmin(depth, 30)),
#          ipN_m3=ipN_h/vol,
#          ipC_m3=ipC_h/vol)
# ip_sf |>
#   ggplot() +
#   geom_sf(data=linnhe_fp) +
#   geom_sf(aes(fill=ipC_m3), colour=NA) +
#   geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
#   scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5, 10)) +
#   facet_grid(date~sim) +
#   # facet_wrap(~sim) +
#   theme(legend.position="bottom",
#         legend.key.width=unit(3, "cm"),
#         legend.key.height=unit(0.2, "cm"))

ip_sum_df <- inner_join(
  linnhe_sf |> select(i, geom),
  inner_join(linnhe_sf |> as_tibble() |> select(i, area, depth),
             ip_df |> complete(i, sim, date, fill=list(ipN_h=0, ipC_h=0)),
             by="i") |>
    mutate(vol=(area*pmin(depth, 30)),
           ipN_m3=(ipN_h/vol)^0.25,
           ipC_m3=(ipC_h/vol)^0.25) |>
    summarise(mnN=mean(ipN_m3),
              sdN=sd(ipN_m3),
              mnC=mean(ipC_m3),
              sdC=sd(ipC_m3),
              .by=c("i", "date"))
) |>
  mutate(month=as.numeric(as.factor(date)))
pA <- ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=mnN^4), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5, 10)) +
  # scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
pB <- ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=(sdN)), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
pC <- ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=sdN/mnN), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
pD <- ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=mnC^4), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5, 10)) +
  # scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
pE <- ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=(sdC)), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
pF <- ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=sdC/mnC), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
ggpubr::ggarrange(pA, pD, ncol=2, nrow=1, common.legend=TRUE, legend="right")
ggpubr::ggarrange(pB, pE, ncol=2, nrow=1, common.legend=TRUE, legend="right")
ggpubr::ggarrange(pC, pF, ncol=2, nrow=1, common.legend=TRUE, legend="right")

ip_dates <- unique(ip_sum_df$date)
lims_sdN <- range(ip_sum_df$sdN)
lims_sdC <- range(ip_sum_df$sdC)
lims_logsdN <- c(-5, max(log(filter(ip_sum_df, sdN > 0)$sdN)))
lims_logsdC <- c(-5, max(log(filter(ip_sum_df, sdC > 0)$sdC)))
lims_mnN <- range(ip_sum_df$mnN)
lims_mnC <- range(ip_sum_df$mnC)
mnBreaks <- c(0, 0.01, 0.1, 1)
for(i in seq_along(ip_dates)) {
  date_i <- ip_dates[i]
  date_i_c <- format(date_i, "%Y-%m-%d")
  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=mnN^4), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_b("Mean IP (Nauplii/m3/h) in top 30m",
                         option="turbo",
                         breaks=c(0.001, 0.01, 0.1, 1, 5),
                         labels=c(0.001, 0.01, 0.1, 1, 5),
                         limits=c(0, 5)) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/meanIP_N_{date_i_c}.png"), p, width=5, height=5.5)
  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=mnC^4), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_b("Mean IP (Copepodids/m3/h) in top 30m",
                         option="turbo",
                         breaks=c(0.001, 0.01, 0.1, 1, 5),
                         labels=c(0.001, 0.01, 0.1, 1, 5),
                         limits=c(0, 5)) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/meanIP_C_{date_i_c}.png"), p, width=5, height=5.5)

  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=mnN), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_c("Mean IP (Nauplii/m3/h) in top 30m",
                         option="turbo",
                         breaks=mnBreaks^0.25,
                         labels=mnBreaks,
                         limits=lims_mnN) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/contmeanIP_N_{date_i_c}.png"), p, width=5, height=5.5)
  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=mnC), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_c("Mean IP (Copepodids/m3/h) in top 30m",
                         option="turbo",
                         breaks=mnBreaks^0.25,
                         labels=mnBreaks,
                         limits=lims_mnC) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/contmeanIP_C_{date_i_c}.png"), p, width=5, height=5.5)


  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=sdN), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_c("sd(nauplii IP) in top 30m", option="turbo", limits=lims_sdN) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/sdIP_N_{date_i_c}.png"), p,
         width=5, height=5.5)
  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=sdC), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_c("sd(copepodid IP) in top 30m", option="turbo", limits=lims_sdC) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/sdIP_C_{date_i_c}.png"), p, width=5, height=5.5)

  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=pmax(log(sdN), -5)), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_c("log sd(nauplii IP) in top 30m", option="turbo",
                         limits=lims_logsdN,
                         breaks=seq(-5,max(lims_logsdN), length.out=4),
                         labels=round(exp(seq(-5,max(lims_logsdN), length.out=4)), 3)) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/logsdIP_N_{date_i_c}.png"), p, width=5, height=5.5)
  p <- ip_sum_df |>
    filter(date==date_i) |>
    ggplot() +
    geom_sf(data=linnhe_fp) +
    geom_sf(aes(fill=pmax(log(sdC), -5)), colour=NA) +
    geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
    scale_fill_viridis_c("log sd(copepodid IP) in top 30m", option="turbo",
                         limits=lims_logsdC,
                         breaks=seq(-5,max(lims_logsdC), length.out=4),
                         labels=round(exp(seq(-5,max(lims_logsdC), length.out=4)), 3)) +
    ggtitle(date_i_c) +
    theme(axis.title=element_blank(),
          legend.position="bottom",
          legend.title.position="top",
          legend.key.width=unit(1.5, "cm"),
          legend.key.height=unit(0.2, "cm"))
  ggsave(glue("figs/logsdIP_C_{date_i_c}.png"), p, width=5, height=5.5)
}


library(av)
sets <- c(outer(c("contmeanIP_", "meanIP_", "logsdIP_", "sdIP_"),
              c("N_", "C_"),
              paste0))
for(i in sets) {
  dirf("figs/temp/", glue("^{i}.*png")) |>
    av_encode_video(glue("figs/anim_{i}.mp4"),
                    framerate=2, vfilter="scale=-2:'min(720,ih)'")
}





# vert dist ---------------------------------------------------------------

# I tested this out, but I think it is more worthwhile to run *more* simulations
# rather than track vertical distributions. Analysis of vertical dynamics can
# use bins at farms rather than the full map.

farms_GSA <- read_csv(glue("{dirs$proj}/data/farm_sites_GSA_2023-2024.csv"))
sim.i <- read_csv(glue("{dirs$out}/sim_i.csv"))
mesh_fp <- st_read(glue("{dirs$mesh}/WeStCOMS2_meshFootprint.gpkg"))
mesh_sf <- st_read(glue("{dirs$mesh}/WeStCOMS2_mesh.gpkg"))
farm_bbox <- list(xmin=125000, xmax=225500, ymin=690000, ymax=785000)
linnhe_fp <- mesh_fp |> st_crop(unlist(farm_bbox))
linnhe_sf <- mesh_sf |> st_crop(unlist(farm_bbox))
sim_dirs <- dirf(dirs$out, "^sim_[0-9][0-9][0-9][0-9]$")

ip_ls <- zSum_ls <- zBin_ls <- vector("list", length(sim_dirs))

for(i in 1:length(sim_dirs)) {
  sim <- str_sub(sim_dirs[i], -4, -1)
  f <- dirrf(sim_dirs[i], "vertDistrMature")
  z_i <- map_dfr(f, ~read_csv(.x) |>
                   mutate(sim=sim,
                          date=ymd(str_split_fixed(basename(.x), "_", 3)[,2]))
  ) |>
    filter(!is.na(value)) # TODO: WHY ARE THESE NA???
  ip_ls[[i]] <- z_i |>
    summarise(lice=sum(value)/730, .by=c("i", "sim", "date"))
  zSum_ls[[i]] <- z_i |>
    summarise(lice=sum(value)/730, .by=c("z", "sim", "date")) |>
    mutate(prop=lice/sum(lice), .by=c("sim", "date"))
  zBin_ls[[i]] <- z_i |>
    mutate(depBin=case_when(z < 5 ~ "0-5",
                            z >= 5 & z < 15 ~ "5-15",
                            z >=15 & z < 30 ~ "15-30",
                            z >= 30 ~ "30+") |>
             factor(levels=c("0-5", "5-15", "15-30", "30+"))) |>
    summarise(lice=sum(value)/730, .by=c("i", "depBin", "sim", "date"))
}

zSum_ls |>
  reduce(bind_rows) |>
  arrange(sim, date, z) |>
  ggplot(aes(lice, z)) +
  geom_vline(xintercept=0, colour="grey80") +
  geom_point(alpha=0.5, shape=1) +
  geom_path(aes(group=sim), alpha=0.5) +
  scale_y_reverse() +
  facet_wrap(~date)
zSum_ls |>
  reduce(bind_rows) |>
  arrange(sim, date, z) |>
  ggplot(aes(prop, z)) +
  geom_vline(xintercept=0, colour="grey80") +
  geom_point(alpha=0.5, shape=1) +
  geom_path(aes(group=sim), alpha=0.5) +
  scale_y_reverse() +
  facet_wrap(~date)
zSum_ls |>
  reduce(bind_rows) |>
  arrange(sim, date, z) |>
  mutate(month=month(date)) |>
  ggplot(aes(prop, z)) +
  geom_vline(xintercept=0, colour="grey80") +
  geom_point(alpha=0.5, shape=1) +
  geom_path(aes(group=date, colour=month), alpha=0.5) +
  scale_y_reverse() +
  facet_wrap(~sim)
zSum_ls |>
  reduce(bind_rows) |>
  arrange(sim, date, z) |>
  mutate(month=as.numeric(as.factor(date))) |>
  ggplot(aes(month, z, fill=prop)) +
  geom_raster() +
  scale_fill_viridis_c(option="turbo") +
  scale_y_reverse() +
  facet_wrap(~sim)

inner_join(linnhe_sf |> select(i, area, depth, geom),
           ip_ls |> reduce(bind_rows),
           by="i") |>
  mutate(IP_m2=lice/area,
         IP_m3=lice/(area*pmin(depth, 30))) |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=IP_m3), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5, 10)) +
  facet_grid(date~sim) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))

ip_sum_df <- inner_join(
  linnhe_sf |> select(i, geom),
  inner_join(linnhe_sf |> as_tibble() |> select(i, area, depth),
             ip_ls |> reduce(bind_rows) |> complete(i, sim, date, fill=list(lice=0)),
             by="i") |>
    mutate(IP_m3=lice/(area*pmin(depth, 30))) |>
    summarise(mn=mean(IP_m3), sd=sd(IP_m3), .by=c("i", "date"))
)
ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=mn), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5, 10)) +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=log10(sd)), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))
ip_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=sd/mn), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_wrap(~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))


inner_join(linnhe_sf |> select(i, area, depth, geom),
           zBin_ls |> reduce(bind_rows) |> filter(date=="2023-04-02"),
           by="i") |>
  mutate(depRng=case_when(depBin=="0-5" ~ 5,
                          depBin=="5-15" ~ 10,
                          depBin=="15-30" ~ 30,
                          depBin=="30+" ~ 30-pmin(depth, 30))) |>
  filter(depBin != "30+") |>
  mutate(IP_m2=lice/area,
         IP_m3=lice/(area*depRng)) |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=IP_m3^0.25), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5)) +
  facet_grid(depBin~sim) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))


zBin_sum_df <- inner_join(
  linnhe_sf |> select(i, geom),
  inner_join(linnhe_sf |> as_tibble() |> select(i, area, depth),
             zBin_ls |> reduce(bind_rows) |>
               filter(depBin != "30+") |>
               droplevels() |>
               complete(i, sim, date, depBin, fill=list(lice=0)),
             by="i") |>
    mutate(depRng=case_when(depBin=="0-5" ~ 5,
                            depBin=="5-15" ~ 10,
                            depBin=="15-30" ~ 30)) |>
    mutate(IP_m3=lice/(area*depRng)) |>
    summarise(mn=mean(IP_m3), sd=sd(IP_m3), .by=c("i", "date", "depBin"))
)

zBin_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=mn), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_b(option="turbo", breaks=c(0.001, 0.01, 0.1, 1, 5)) +
  facet_grid(depBin~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))

zBin_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=log10(sd)), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_grid(depBin~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))

zBin_sum_df |>
  ggplot() +
  geom_sf(data=linnhe_fp) +
  geom_sf(aes(fill=sd/mn), colour=NA) +
  geom_point(data=farms_GSA, aes(easting, northing), colour="red", shape=1) +
  scale_fill_viridis_c(option="turbo") +
  facet_grid(depBin~date) +
  theme(legend.position="bottom",
        legend.key.width=unit(3, "cm"),
        legend.key.height=unit(0.2, "cm"))



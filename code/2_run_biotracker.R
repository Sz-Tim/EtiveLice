# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Run biotracker simulations


# setup
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(doFuture)
library(sf)
dirf("code/fn", ".R") |> walk(source)
theme_set(theme_bw() + theme(panel.grid=element_blank()))

set.seed(1001)


# define parameters -------------------------------------------------------

cores_per_sim <- 15
parallel_sims <- 2
start_date <- "2023-01-01"
end_date <- "2024-12-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hf0="/home/sa04ts/hydro/etive28/Archive",
             hf1="/home/sa04ts/hydro/WeStCOMS2/Archive",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v2-2-1.jar",
             out=glue("{getwd()}/out/biotracker")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hf0="E:/hydro/etive28/Archive",
               hf1="E:/hydro/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v2-2-1.jar",
               out=glue("E:/EtiveLice/out/biotracker"))
)

n_sim <- 20
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
sim.i <- sample_parameter_distributions(n_sim, dirs$out, "lhs", egg_post, mort_post, sink_post) |>
  mutate(connectivityThresh=30,
         variableDh="false",
         variableDhV="false")
write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))
sim_seq <- 1:nrow(sim.i)


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
       releaseInterval=1,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh0=glue("{dirs$mesh}/etive28_mesh.nc"),
       meshType0="FVCOM",
       hfDir0=glue("{dirs$hf0}/"),
       hfDirPrefix0="netcdf_",
       hfFilePrefix0="etive28",
       mesh1=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       meshType1="FVCOM",
       hfDir1=glue("{dirs$hf1}/"),
       hfDirPrefix1="netcdf_",
       hfFilePrefix1="westcoms2",
       # sites
       sitefile=glue("E:/EtiveLice/data/pen_sites_linnhe_2023-2024.csv"),
       sitefileEnd=glue("E:/EtiveLice/data/pen_sites_linnhe_2023-2024.csv"),
       siteDensityPath=glue("E:/EtiveLice/data/lice_daily_2023-01-01_2024-12-31_05lpf_maxB.csv"),
       # sitefile="D:/EtiveLice/data/pen_sites_linnhe_2023-2024.csv",
       # sitefileEnd="D:/EtiveLice/data/pen_sites_etive_2023-2024.csv",
       # siteDensityPath="D:/EtiveLice/data/lice_daily_2023-01-01_2024-12-31.csv",
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
       connectivityThresh=30,
       # recording
       verboseSetUp="true",
       recordConnectivity="true",
       connectivityInterval=1,
       connectDepth1_min=0,
       connectDepth1_max=20,
       recordPsteps= "false",
       recordVertDistr="false",
       recordSiteEnv="true",
       siteEnvInterval=1,
       siteEnvMaxDepth=20,
       recordElemActivity="false"))


# run simulations ---------------------------------------------------------

plan(multisession, workers=parallel_sims)
sim_sets <- split(sim_seq, rep(1:parallel_sims, length(sim_seq)/parallel_sims))
foreach(j=1:parallel_sims, .options.future=list(globals=structure(TRUE, add="sim.i"))) %dofuture% {
  for(i in sim_sets[[j]]) {
    setwd(dirs$proj)
    cat("Starting", sim.i$i[i], "\n")
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
pens_linnhe <- read_csv("data/pen_sites_linnhe_2023-2024.csv")
pens_etive <- read_csv("data/pen_sites_etive_2023-2024.csv")
etive_farms <- read_csv("data/pen_sites_etive_2023-2024.csv") |>
  mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1])
linnhe_farms <- read_csv("data/pen_sites_linnhe_2023-2024.csv") |>
  mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1])

farm_dat <- read_csv("data/lice_biomass_2017-01-01_2024-12-31.csv") |>
  filter(between(date, ymd(start_date), ymd(end_date))) |>
  select(sepaSite, date, weeklyAverageAf, actualBiomassOnSiteTonnes, maximumBiomassAllowedTonnes)
write_csv(farm_dat, "data/sim/inputs/farm_dat.csv")

# Site conditions
# TODO: Add depth from siteEnv_Average.csv
site_env_df <- dirf(glue("{dirs$out}/sim_01"), "siteEnv_2") |>
  map_dfr(~read_csv(.x, show_col_types=F) |>
            mutate(elapsedHours=str_sub(str_split_fixed(basename(.x), "_", 3)[,3], 1, -5))) |>
  rename(pen=site) |>
  mutate(pen=paste0(str_split_fixed(pen, "_", 2)[,1],
                    "_",
                    str_pad(str_split_fixed(pen, "_", 2)[,2], 2, "left", "0")),
         time=ymd(start_date) + dhours(as.numeric(elapsedHours))) |>
  inner_join(linnhe_farms |> select(sepaSite, pen)) #|>
  # select(sepaSite, pen, date, depth, u, v, w, uv, salinity, temperature)
write_csv(site_env_df, "data/sim/inputs/farm_env.csv")


site_env_df |>
  select(sepaSite, time, pen, uv, salinity, temperature, light) |>
  pivot_longer(4:7) |>
  group_by(name) |>
  mutate(valRel=(value-min(value))/(max(value)-min(value)),
         value=c(scale(value))) |>
  ungroup() |>
  mutate(pen=factor(pen,
                    levels=c(t(outer(paste0(c("FFMC27", "SAR1", "APT1", "FFMC32", "FFMC84"), "_"),
                                     str_pad(1:20, 2, "left", "0"),
                                     paste0)))),
         sepaSite=factor(sepaSite, levels=c("FFMC27", "SAR1", "APT1", "FFMC32", "FFMC84"))) |>
  ggplot(aes(time, pen, fill=valRel)) +
  geom_raster() +
  facet_grid(name~., axes="all_x", axis.labels="margins") +
  scale_fill_viridis_c(option="turbo") +
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day") +
  theme(axis.text.y=element_text(size=6))

site_env_df |>
  select(sepaSite, time, pen, uv, salinity, temperature, light) |>
  pivot_longer(4:7) |>
  group_by(name, sepaSite, time) |>
  summarise(value=mean(value)) |>
  group_by(name) |>
  mutate(valRel=(value-min(value))/(max(value)-min(value)),
         value=c(scale(value))) |>
  ungroup() |>
  mutate(sepaSite=factor(sepaSite, levels=c("FFMC27", "SAR1", "APT1", "FFMC32", "FFMC84"))) |>
  ggplot(aes(time, sepaSite, fill=valRel)) +
  geom_raster() +
  facet_grid(name~., axes="all_x", axis.labels="margins") +
  scale_fill_viridis_c(option="turbo") +
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day")

site_env_df |>
  select(sepaSite, time, pen, uv, salinity, temperature, light) |>
  pivot_longer(4:7) |>
  group_by(name) |>
  ungroup() |>
  mutate(pen=factor(pen,
                    levels=c(t(outer(paste0(c("FFMC27", "SAR1", "APT1", "FFMC32", "FFMC84"), "_"),
                                     str_pad(1:20, 2, "left", "0"),
                                     paste0)))),
         sepaSite=factor(sepaSite, levels=c("FFMC27", "SAR1", "APT1", "FFMC32", "FFMC84"))) |>
  ggplot(aes(time, value, group=pen, colour=sepaSite)) +
  geom_line() +
  scico::scale_colour_scico_d(palette="devon", end=0.8, direction=1, guide="none") +
  facet_grid(name~., axes="all_x", axis.labels="margins", scales="free_y") +
  scale_x_datetime(date_breaks="2 days", date_labels="%b-%d", date_minor_breaks="1 day")

# Influx
mesh_fp <- st_read(glue("{dirs$mesh}/etive28_mesh_footprint.gpkg"))
c_long <- dirrf(dirs$out, "connectivity.*csv") |>
  map_dfr(~load_connectivity(.x,
                             source_names=pens_linnhe$pen,
                             dest_names=pens_linnhe$pen,
                             liceScale=1) |>
            mutate(sim=str_sub(str_split_fixed(.x, "sim_", 2)[,2], 1, 2),
                   elapsedHours=str_sub(str_split_fixed(basename(.x), "_", 4)[,4], 1, -5))) |>
  mutate(time=ymd(start_date) + dhours(as.numeric(elapsedHours)))
site_areas <- sim.i |>
  mutate(connectivityThresh=30) |>
  select(i, connectivityThresh) |>
  mutate(site_df=list(pens_linnhe)) |>
  unnest(site_df) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700) %>%
  st_buffer(dist=.$connectivityThresh) |>
  st_intersection(mesh_fp) %>%
  mutate(area=as.numeric(st_area(.))) |>
  st_drop_geometry() |>
  rename(sim=i) |>
  select(-connectivityThresh)
write_csv(site_areas, "data/sim/inputs/farm_influx_areas.csv")

# mean hourly IP for each day
c_daily <- list(
  c_long |>
    calc_influx(destination, value, sim, date) |>
    rename(pen=destination) |>
    inner_join(site_areas, by=join_by(pen, sim)) |>
    mutate(influx_m2=influx/area) |>
    select(-area),
  c_long |>
    mutate(destination=factor(destination, levels=levels(source))) |>
    calc_self_infection(source, destination, value, sim, date) |> rename(pen=source) |>
    inner_join(site_areas, by=join_by(pen, sim)) |>
    mutate(self_m2=self/area) |>
    select(-area)
) |>
  reduce(full_join) |>
  complete(pen, sim, date,
           fill=list(influx=0, influx_m2=0, N_influx=0,
                     self=0, self_m2=0)) |>
  mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1])

write_csv(c_daily, "data/sim/inputs/influx.csv")



c_long |>
  group_by(destination, time, sim) |>
  summarise(value=sum(value)) |>
  ungroup() |>
  mutate(sepaSite=str_split_fixed(destination, "_", 2)[,1]) |>
  filter(sepaSite=="SAR1") |>
  inner_join(site_env_df, by=join_by(sepaSite==sepaSite, time==time, destination==pen)) |>
  ggplot(aes(time)) +
  geom_rug(aes(colour=uv), linewidth=1.5, length=unit(1, "npc")) +
  geom_point(aes(y=value, group=destination), shape=1) +
  geom_line(aes(y=value, group=destination), linewidth=1) +
  scale_colour_viridis_c(option="turbo") +
  facet_grid(destination~.)

c_long |>
  group_by(destination, time, sim) |>
  summarise(value=sum(value)) |>
  ungroup() |>
  mutate(sepaSite=str_split_fixed(destination, "_", 2)[,1]) |>
  group_by(sepaSite, time, sim) |>
  summarise(lice=mean(value)) |>
  ungroup() |>
  inner_join(site_env_df |>
               group_by(sepaSite, time) |>
               summarise(across(where(is.numeric), mean)) |>
               ungroup(),
             by=join_by(sepaSite==sepaSite, time==time)) |>
  ggplot(aes(time)) +
  geom_rug(aes(colour=uv), linewidth=1.5, length=unit(1, "npc")) +
  geom_point(aes(y=lice), shape=1) +
  geom_line(aes(y=lice), linewidth=1) +
  scale_colour_viridis_c(option="turbo") +
  facet_wrap(~sepaSite, scales="free_y")


c_long |>
  group_by(destination, time, sim) |>
  summarise(value=sum(value)) |>
  ungroup() |>
  mutate(sepaSite=str_split_fixed(destination, "_", 2)[,1]) |>
  group_by(sepaSite, time, sim) |>
  summarise(lice=mean(value)) |>
  ungroup() |>
  inner_join(site_env_df |>
               group_by(sepaSite, time) |>
               summarise(across(where(is.numeric), mean)) |>
               ungroup(),
             by=join_by(sepaSite==sepaSite, time==time)) |>
  ggplot(aes(salinity, lice)) +
  geom_point(shape=1, alpha=0.5) +
  facet_wrap(~sepaSite, scales="free_y")

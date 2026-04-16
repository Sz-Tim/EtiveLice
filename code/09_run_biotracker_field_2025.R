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

cores_per_sim <- 15
parallel_sims <- 2
start_date <- "2025-03-01"
end_date <- "2025-10-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hf0="/home/sa04ts/hydro/etive28/Archive",
             hf1="/home/sa04ts/hydro/WeStCOMS2/Archive",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v2-2-2.jar",
             dat="/home/sa04ts/EtiveLice/data",
             out=glue("{getwd()}/out/biotracker")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hf0="E:/hydro/etive28/Archive",
               hf1="E:/hydro/WeStCOMS3/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v2-2-2.jar",
               dat="E:/EtiveLice/data",
               out="E:/EtiveLice/out/biotracker/field_2025")
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
sim.i <- sample_parameter_distributions(
  n_sim, dirs$out, "lhs", egg_post, mort_post, sink_post,
  alt_ls=list(
    D_h=c(1e-4, 1e0),
    D_hVert=c(1e-6, 1e-2),
    lightN=c(0.05, 0.5),
    lightC=c(2e-6, 2e-4),
    swimUpN=c(0.05, 1)*1e-3, # mm/s -> m/s
    swimDownN=c(0.05, 1)*1e-3, # mm/s -> m/s
    swimUpC=c(0.5, 10)*1e-3, # mm/s -> m/s
    swimDownC=c(0.1, 5)*1e-3, # mm/s -> m/s
    passiveSinkSal=c(T, F),
    salThreshMaxN=c(20, 32),
    salThreshSpanN=c(0, 15),
    salThreshMaxC=c(20, 32),
    salThreshSpanC=c(0, 15),
    viableDD=c(30, 55),
    maxDepth=c(40, 300)
  )) |>
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
       nparts=20,
       releaseInterval=1,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh0=glue("{dirs$mesh}/etive28_mesh.nc"),
       meshType0="FVCOM",
       hfDir0=glue("{dirs$hf0}/"),
       hfDirPrefix0="netcdf_",
       hfFilePrefix0="etive28",
       mesh1=glue("{dirs$mesh}/WeStCOMS3_mesh.nc"),
       meshType1="FVCOM",
       hfDir1=glue("{dirs$hf1}/"),
       hfDirPrefix1="netcdf_",
       hfFilePrefix1="westcoms3",
       # sites
       sitefile=glue("{dirs$dat}/pen_sites_linnhe_2025.csv"),
       sitefileEnd=glue("{dirs$dat}/SAMS_etive_sites_OS.csv"),
       siteDensityPath=glue("{dirs$dat}/lice_daily_2025-03-01_2025-12-31.csv"),
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
       connectImmature="true",
       connectivityInterval=1,
       connectDepth1_min=0,
       connectDepth1_max=2,
       connectDepth2_min=11,
       connectDepth2_max=13,
       connectDepth3_min=19,
       connectDepth3_max=21,
       recordPsteps= "true",
       pstepsInterval=1,
       recordImmature="true",
       recordVertDistr="false",
       recordSiteEnv="false",
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
pens_linnhe <- read_csv("data/pen_sites_linnhe_2025.csv")
sites_etive <- read_csv("data/SAMS_etive_sites_OS.csv")
linnhe_farms <- read_csv("data/pen_sites_linnhe_2025.csv") |>
  mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1])
etive_order <- c("LY1", "RE8", "RE10", "RE9", "RE7", "RE6", "RE5")


# Influx
mesh_wce <- st_read(glue("{dirs$mesh}/WeStCOMS3_etive28_mesh.gpkg"))
site_vols <- sim.i |>
  mutate(connectivityThresh=30) |>
  select(i, connectivityThresh) |>
  mutate(site_df=list(sites_etive)) |>
  unnest(site_df) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700) %>%
  st_buffer(dist=.$connectivityThresh) |>
  st_intersection(mesh_wce) %>%
  mutate(area=as.numeric(st_area(.))) |>
  st_drop_geometry() |>
  group_by(i, name) |>
  summarise(vol=sum(area*pmin(depth, 20)),
            area=sum(area)) |>
  ungroup() |>
  rename(sim=i)
write_csv(site_vols, "out/field_2025/site_influx_vols.csv")

# hourly IP for each day per pen
library(furrr); library(future)
plan(multisession, workers=20)
c_df <- future_map_dfr(
  dirrf(dirs$out, "connectivity.*csv"),
  ~load_connectivity(.x,
                     source_names=pens_linnhe$pen,
                     dest_names=sites_etive$name,
                     liceScale=1) |>
    calc_influx(destination, value, depthRange) |>
    mutate(stage=if_else(grepl("Immature", .x), "Nauplius", "Copepodid"),
           sim=str_sub(str_split_fixed(.x, "sim_", 2)[,2], 1, 2),
           elapsedHours=str_sub(str_split_fixed(basename(.x), "_", 4)[,4], 1, -5))
) |>
  mutate(time=ymd(start_date) + dhours(as.numeric(elapsedHours)-1)) |>
  rename(name=destination) |>
  inner_join(site_vols, by=join_by(name, sim)) |>
  mutate(influx_m2=influx/area) |>
  select(-area, -vol, -elapsedHours) |>
  complete(name, sim, time, stage, depthRange,
           fill=list(influx=0, influx_m2=0, N_influx=0)) |>
  mutate(station=factor(name, levels=etive_order),
         stage=factor(stage, levels=c("Nauplius", "Copepodid"))) |>
  select(-name)
plan(sequential)
data.table::fwrite(c_df, "out/field_2025/influx_hourly.csv")

plan(multisession, workers=20)
c_df_all <- future_map_dfr(dirrf(dirs$out, "connectivity.*csv"),
    ~load_connectivity(.x,
                       source_names=pens_linnhe$pen,
                       dest_names=sites_etive$name,
                       liceScale=1) |>
      mutate(stage=if_else(grepl("Immature", .x), "Nauplius", "Copepodid"),
             sim=str_sub(str_split_fixed(.x, "sim_", 2)[,2], 1, 2),
             elapsedHours=str_sub(str_split_fixed(basename(.x), "_", 4)[,4], 1, -5))
) |>
  mutate(time=ymd(start_date) + dhours(as.numeric(elapsedHours)-1)) |>
  rename(name=destination, influx=value) |>
  mutate(source=str_split_fixed(source, "_", 2)[,1]) |>
  group_by(name, source, sim, time, stage, depthRange) |>
  summarise(influx=sum(influx)) |>
  ungroup() |>
  inner_join(site_vols, by=join_by(name, sim)) |>
  mutate(influx_m2=influx/area) |>
  complete(name, source, sim, time, stage, depthRange,
           fill=list(influx=0, influx_m2=0)) |>
  mutate(station=factor(name, levels=etive_order),
         stage=factor(stage, levels=c("Nauplius", "Copepodid"))) |>
  select(-name, -vol, -area)

dat_df <- readRDS("data/field_data_2025_TEMP_processed.rds") |>
  mutate(station=factor(Station_simple, levels=levels(Station_simple), labels=etive_order),
         Count=if_else(Stage=="Copepodids", Lepeophtheirus_salmonis, Count),
         stage=if_else(Stage=="Copepodids", "Copepodid", "Nauplius")) |>
  select(Sample_id, Date_collected, station, Depth, Depth_F, hour, Salinity, stage, Count, m3)

site_dates_df <- dat_df |>
  summarise(.by=c(station, Date_collected))

c_df <- data.table::fread("out/field_2025/influx_hourly.csv")
influx_obs_df <- inner_join(
  c_df |>
    # summarise(influx_m2_lnMn=expm1(mean(log1p(influx_m2))),
    #           influx_m2_sqrtMn=mean(sqrt(influx_m2))^2,
    #           influx_m2_rtrtMn=mean(sqrt(sqrt(influx_m2)))^4,
    #           influx_m2=mean(influx_m2),
    #           .by=c(time, stage, depthRange, station)) |>
    mutate(Date_collected=date(time),
           hour=hour(time),
           Depth=case_when(depthRange=="0.5-1.5m" ~ 1,
                           depthRange=="11.5-12.5m" ~ 12,
                           depthRange=="19.5-20.5m" ~ 20)),
  site_dates_df,
  by=join_by(Date_collected, station)
) |>
  left_join(
    dat_df,
    by=join_by(Date_collected, hour, station, Depth, stage)
  )

influx_obs_df |>
  filter(station %in% etive_order[-(2:3)]) |>
  # filter(Count > 0) |>
  ggplot(aes(influx_m2, (Count/m3))) +
  geom_abline(linetype=3, colour="grey") +
  geom_point() +
  stat_smooth(method="lm", se=F) +
  facet_grid(Depth~station)

influx_obs_df |>
  filter(station %in% etive_order[-(2:3)]) |>
  filter(!is.na(Depth_F)) |>
  ggplot(aes(influx_m2, (Count/m3), colour=Depth_F)) +
  geom_point(alpha=0.5) +
  stat_smooth(method="lm", se=F) +
  scale_colour_viridis_d(direction=-1, begin=0.1) +
  facet_wrap(~stage*station, scales="fixed", nrow=2)


influx_obs_df |>
  filter(stage=="Copepodid") |>
  filter(Depth <= 20) |>
  ggplot(aes(hour, colour=Depth)) +
  geom_line(aes(y=influx_m2, group=paste(Depth, sim))) +
  geom_point(aes(y=Count/m3)) +
  scale_colour_viridis_c(direction=-1, begin=0.1) +
  facet_wrap(~station*Date_collected, scales="free", nrow=4)



c_df_all |>
  ggplot(aes(time, source, fill=influx_m2)) +
  geom_raster() +
  scale_fill_viridis_c(option="turbo") +
  facet_grid(depthRange*stage~station)

c_df_all |>
  # filter(stage=="Copepodid") |>
  ggplot(aes(time, source, fill=log(influx_m2))) +
  geom_raster() +
  scale_fill_viridis_c(option="turbo") +
  facet_grid(stage*depthRange~station)

c_df |>
  ggplot(aes(time, depthRange, fill=(influx_m2))) +
  geom_raster() +
  scale_fill_viridis_c(option="turbo") +
  facet_grid(stage~station)

c_df |>
  ggplot(aes(time, influx_m2, colour=depthRange, group=paste(sim, depthRange))) +
  geom_line() +
  scale_colour_viridis_d(end=0.95, direction=-1) +
  facet_grid(stage*depthRange~station)

c_df |>
  filter(stage=="Copepodid") |>
  ggplot(aes(time, influx_m2, colour=depthRange, group=paste(sim, depthRange))) +
  geom_line() +
  scale_colour_viridis_d(end=0.95, direction=-1) +
  facet_grid(sim~station)

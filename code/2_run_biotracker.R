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
theme_set(theme_bw() + theme(panel.grid=element_blank()))

set.seed(11)


# define parameters -------------------------------------------------------

cores_per_sim <- 20
parallel_sims <- 1
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
             jdk="/usr/local/java/jre1.8.0_211/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v1-0-0.jar",
             out=glue("{getwd()}/out/biotracker")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hydro1="E:/hydro/etive28/Archive",
               hydro2="E:/hydro/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-22.0.1/bin/javaw",
               # jdk="C:/Users/sa04ts/.jdks/openjdk-19/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v1-0-0.jar",
               out=glue("D:/EtiveLice/out/biotracker"))
)

post_dir <- glue("{dirs$proj}/data/lit/fit")
egg_post <- list(constant=tibble(b=rnorm(4000, 28.2, 2)),
                 logistic=read_csv(glue("{post_dir}/egg_temp_logistic_post.csv"), show_col_types=F),
                 linear=read_csv(glue("{post_dir}/egg_temp_linear_post.csv"), show_col_types=F)) |>
  map(~.x |> mutate(across(everything(), ~signif(.x, 3))) |>
        unite("all", everything(), sep=",") %>%
        .$all)
mort_post <- list(constant=tibble(b=plogis(rnorm(4000, qlogis(0.1), 0.25))),
                 logistic=read_csv(glue("{post_dir}/mort_sal_logistic_post.csv"), show_col_types=F)) |>
  map(~.x |> mutate(across(everything(), ~signif(.x, 3))) |>
        unite("all", everything(), sep=",") %>%
        .$all)
n_sim <- 6
sim.i <- tibble(D_h=runif(n_sim, 0.05, 0.5),
                D_hVert=runif(n_sim, 0.005, 0.05),
                mortSal_fn=sample(c("constant", "logistic"), n_sim, replace=T),
                eggTemp_fn=sample(c("constant", "logistic"), n_sim, replace=T),
                lightThreshCopepodid=runif(n_sim, 0.00001, 0.0005),
                lightThreshNauplius=runif(n_sim, 0.05, 0.75),
                swimUpSpeedMean=runif(n_sim, 1e-4, 1e-2) * -1,
                swimDownSpeedMean=runif(n_sim, 1e-4, 1e-2),
                passiveSinkRateSal=sample(c(T, F), n_sim, replace=T),
                salinityThreshMin=runif(n_sim, 20, 28),
                salinityThreshMax=pmin(salinityThreshMin + runif(n_sim, 0.1, 6), 32),
                viableDegreeDays=runif(n_sim, 35, 45),
                connectivityThresh=runif(n_sim, 100, 200)) |>
  mutate(across(where(is.numeric), ~signif(.x, 5))) |>
  mutate(i=str_pad(row_number(), 2, "left", "0"),
         outDir=glue("{dirs$out}/sim_{i}/")) |>
  rowwise() |>
  mutate(eggTemp_b=sample(egg_post[[eggTemp_fn]], 1),
         mortSal_b=sample(mort_post[[mortSal_fn]], 1)) |>
  ungroup()
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
       checkOpenBoundaries="true",
       # meshes and environment
       mesh1=glue("{dirs$mesh}/etive28_mesh.nc"),
       mesh1Domain="etive28",
       datadir=glue("{dirs$hydro1}/"),
       mesh2=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       mesh2Domain="westcoms2",
       datadir2=glue("{dirs$hydro2}/"),
       # sites
       sitefile=glue("D:/EtiveLice/data/farm_sites.csv"),
       sitefileEnd=glue("D:/EtiveLice/data/farm_sites.csv"),
       siteDensityPath=glue("D:/EtiveLice/data/lice_daily_2023-01-01_2023-12-31.csv"),
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
       salinityThreshMin=sim.i$salinityThreshMin[.x],
       salinityThreshMax=sim.i$salinityThreshMax[.x],
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
       connectDepth1_max=20,
       recordPsteps= "false",
       recordVertDistr="false",
       recordElemActivity="false"))


# run simulations ---------------------------------------------------------

plan(multisession, workers=parallel_sims)
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
etive_farms <- c("FFMC84", "FFMC32", "APT1", "SAR1", "FFMC27")
# Site conditions
site_env_df <- dirf(glue("{dirs$out}/sim_01"), "siteConditions") |>
  map_dfr(~read_csv(.x, show_col_types=F,
                    col_select=c(site, depth, mesh, u, v, w, uv, salinity, temperature)) |>
            mutate(date=ymd(str_split_fixed(basename(.x), "_", 3)[,2]))) |>
  filter(site %in% etive_farms) |>
  rename(sepaSite=site)
write_csv(site_env_df, "data/sim/inputs/farm_env.csv")

# Influx
mesh_fp <- st_read(glue("{dirs$mesh}/etive28_mesh_footprint.gpkg"))
site_i <- read_csv("D:/EtiveLice/data/farm_sites.csv")
c_long <- map_dfr(dirrf(dirs$out, "connectivity.*csv"),
                  ~load_connectivity(.x, site_i$sepaSite, liceScale=1) |>
                    mutate(sim=str_sub(str_split_fixed(.x, "sim_", 2)[,2], 1, 2)))
site_areas <- sim.i |>
  select(i, connectivityThresh) |>
  mutate(site_df=list(site_i)) |>
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
    rename(sepaSite=destination) |>
    inner_join(site_areas, by=join_by(sepaSite, sim)) |>
    mutate(influx_m2=influx/area) |>
    select(-area),
  c_long |>
    calc_self_infection(source, destination, value, sim, date) |> rename(sepaSite=source) |>
    inner_join(site_areas, by=join_by(sepaSite, sim)) |>
    mutate(self_m2=self/area) |>
    select(-area)
) |>
  reduce(full_join) |>
  complete(sepaSite, sim, date,
           fill=list(influx=0, influx_m2=0, N_influx=0,
                     self=0, self_m2=0))

write_csv(c_daily, "data/sim/inputs/influx.csv")

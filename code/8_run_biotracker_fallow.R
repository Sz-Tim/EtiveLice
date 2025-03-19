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

set.seed(101)


# define parameters -------------------------------------------------------

cores_per_sim <- 15
parallel_sims <- 1
start_date <- "2023-03-01"
end_date <- "2023-03-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))
etive_farms <- c("FFMC84", "FFMC32", "APT1", "SAR1")

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hydro1="/home/sa04ts/hydro/etive28/Archive",
             hydro2="/home/sa04ts/hydro/WeStCOMS2/Archive",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v1-0-0.jar",
             out=glue("{getwd()}/out/fallow")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hydro1="E:/hydro/etive28/Archive",
               hydro2="E:/hydro/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker.jar",
               out=glue("D:/EtiveLice/out/fallow"))
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
n_sim <- 25
swim_mx <- MASS::mvrnorm(n_sim, c(0,0), matrix(c(1, 0.8, 0.8, 1), nrow=2))
light_mx <- MASS::mvrnorm(n_sim, c(0,0), matrix(c(1, 0.8, 0.8, 1), nrow=2))
sim.i <- tibble(D_h=runif(n_sim, 0.05, 0.5),
                D_hVert=runif(n_sim, 0.0005, 0.005),
                mortSal_fn=sample(c("constant", "logistic"), n_sim, replace=T),
                eggTemp_fn=sample(c("constant", "logistic"), n_sim, replace=T),
                lightThreshCopepodid=qunif(pnorm(light_mx[,1]), (2e-6)^0.5, (2e-4)^0.5)^2,
                lightThreshNauplius=qunif(pnorm(light_mx[,2]), (0.05)^0.5, (0.5)^0.5)^2,
                swimUpSpeedMean=-(qunif(pnorm(swim_mx[,1]), (1e-4)^0.5, (2e-2)^0.5))^2,
                swimDownSpeedMean=(qunif(pnorm(swim_mx[,2]), (1e-4)^0.5, (2e-2)^0.5))^2,
                passiveSinkRateSal=sample(c(T, F), n_sim, replace=T),
                salinityThreshMin=runif(n_sim, 20, 28),
                salinityThreshMax=pmin(salinityThreshMin + runif(n_sim, 0.1, 6), 32),
                viableDegreeDays=runif(n_sim, 35, 45),
                connectivityThresh=30) |>
  mutate(across(where(is.numeric), ~signif(.x, 5))) |>
  mutate(i=str_pad(row_number(), 3, "left", "0"),
         outDir=glue("{dirs$out}/sim_{i}/")) |>
  rowwise() |>
  mutate(eggTemp_b=sample(egg_post[[eggTemp_fn]], 1),
         mortSal_b=sample(mort_post[[mortSal_fn]], 1)) |>
  ungroup()
write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))
sim_seq <- 1:nrow(sim.i)
sim_seq <- sim_seq[1]



# fallow Etive ------------------------------------------------------------

boundary_release <- st_read(glue("{dirs$proj}/data/etive28_boundary_centroids.gpkg")) |>
  select(i, geom) |>
  sevcheck::add_lonlat(drop_geom=T) |>
  mutate(i=paste0("elem", i)) |>
  rename(sepaSite=i, easting=lon, northing=lat)
write_csv(boundary_release, glue("{dirs$proj}/data/etive28_boundary_centroids.csv"))
write_csv(boundary_release, glue("D:/EtiveLice/data/etive28_boundary_centroids.csv"))
init_df <- boundary_release |>
  select(sepaSite) |>
  bind_cols(tibble(date=ymd(start_date) - 1 + 1:nDays,
                   init=100000) |>
              mutate(date=format(date, "%Y%m%d")) |>
              pivot_wider(names_from=date, values_from=init))
write_csv(init_df, glue("{dirs$proj}/data/lice_daily_etive28_boundary_centroids.csv"))
write_csv(init_df, glue("D:/EtiveLice/data/lice_daily_etive28_boundary_centroids.csv"))

farms_df <- read_csv(glue("{dirs$proj}/data/farm_sites_2023.csv"))
etive_df <- farms_df |>
  filter(sepaSite %in% etive_farms)
write_csv(etive_df, glue("{dirs$proj}/data/etive_farms_2023.csv"))
write_csv(etive_df, glue("D:/EtiveLice/data/etive_farms_2023.csv"))
linnhe_df <- farms_df |>
  filter(! sepaSite %in% etive_farms)
write_csv(linnhe_df, glue("{dirs$proj}/data/linnhe_farms_2023.csv"))
write_csv(linnhe_df, glue("D:/EtiveLice/data/linnhe_farms_2023.csv"))
init_linnhe <- linnhe_df |>
  select(sepaSite) |>
  bind_cols(tibble(date=ymd(start_date) - 1 + 1:nDays,
                   init=100000) |>
              mutate(date=format(date, "%Y%m%d")) |>
              pivot_wider(names_from=date, values_from=init))
write_csv(init_linnhe, glue("{dirs$proj}/data/lice_daily_linnhe_farms.csv"))
write_csv(init_linnhe, glue("D:/EtiveLice/data/lice_daily_linnhe_farms.csv"))



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
       nparts=100,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh1=glue("{dirs$mesh}/etive28_mesh.nc"),
       mesh1Domain="etive28",
       datadir=glue("{dirs$hydro1}/"),
       mesh2=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       mesh2Domain="westcoms2",
       datadir2=glue("{dirs$hydro2}/"),
       # sites
       sitefile= glue("D:/EtiveLice/data/linnhe_farms_2023.csv"),
       sitefileEnd=glue("D:/EtiveLice/data/etive_farms_2023.csv"),
       siteDensityPath="D:/EtiveLice/data/lice_daily_linnhe_farms.csv",
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
       connectDepth1_max=5,
       connectDepth2_min=5,
       connectDepth2_max=20,
       recordPsteps= "true",
       pstepsInterval=1,
       recordImmature="true",
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

mesh_fp <- st_read(glue("{dirs$mesh}/etive28_mesh_footprint.gpkg"))
mesh <- st_read(glue("{dirs$mesh}/etive28_mesh.gpkg"))
sim.i <- read_csv(glue("{dirs$out}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i)) |> slice_head(n=1)
src <- read_csv("D:/EtiveLice/data/linnhe_farms_2023.csv")
dest <- read_csv("D:/EtiveLice/data/etive_farms_2023.csv")
out_dir <- "D:/EtiveLice/out/fallow/"

# Particle densities
plan(multisession, workers=18)
date_seq <- seq(ymd(start_date), ymd(end_date), by=1) |> str_remove_all("-")
lims_naup <- lims_cop <- tibble(copDens_rtrt=c(0,0))
mesh_psteps <- st_read(glue("{dirs$mesh}/etive28_mesh.gpkg")) |>
  mutate(vol_top20m=area*pmin(depth, 20)) |>
  st_drop_geometry() |>
  select(i, area, vol_top20m)
for(i in 1:length(date_seq)) {
  # Immature: nauplii
  ps_i <- load_psteps_simSets(out_dir,
                              mesh_psteps,
                              sim.i, ncores=1, liceScale=1,
                              stage=paste0("Immature_", date_seq[i]), per_m2=TRUE, trans="4th_rt")
  i_ts <- grep("^t_", names(ps_i), value=T)
  cat("Starting", date_seq[i], "\n")
  for(j in seq_along(i_ts)) {
    ps_j <- ps_i |>
      filter(!is.na(i)) |>
      select("sim", "i", all_of(i_ts[j])) |>
      drop_na() |>
      pivot_wider(names_from=sim, values_from=starts_with("t_"))
    ps_j |>
      saveRDS(glue("{out_dir}/processed/hourly/Immature_{date_seq[i]}_{i_ts[j]}.rds"))
    if(nrow(ps_j) > 0) {
      lims_naup$copDens_rtrt <- range(c(lims_naup$copDens_rtrt, range(ps_j$sim_001)))
    }
  }
  # Mature: copepodids
  ps_i <- load_psteps_simSets(out_dir,
                              mesh_psteps,
                              sim.i, ncores=1, liceScale=1,
                              stage=paste0("Mature_", date_seq[i]), per_m2=TRUE, trans="4th_rt")
  i_ts <- grep("^t_", names(ps_i), value=T)
  for(j in seq_along(i_ts)) {
    ps_j <- ps_i |>
      filter(!is.na(i)) |>
      select("sim", "i", all_of(i_ts[j])) |>
      drop_na() |>
      pivot_wider(names_from=sim, values_from=starts_with("t_"))
    ps_j |>
      saveRDS(glue("{out_dir}/processed/hourly/Mature_{date_seq[i]}_{i_ts[j]}.rds"))
    if(nrow(ps_j) > 0) {
      lims_cop$copDens_rtrt <- range(c(lims_cop$copDens_rtrt, range(ps_j$sim_001)))
    }
  }
  cat("", j)
}
saveRDS(lims_naup, glue("{out_dir}/processed/hourly_pslims_naup.rds"))
saveRDS(lims_cop, glue("{out_dir}/processed/hourly_pslims_cop.rds"))
cat("\n")
plan(sequential)
gc()



# anim --------------------------------------------------------------------
library(tidyverse); library(glue)
library(sf)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(rstan)
library(yardstick)
library(terra)
library(recipes)
library(ggpubr)
library(cowplot)
library(scico)
library(ggdist)
library(ggnewscale)
theme_set(theme_bw() + theme(panel.grid=element_blank()))

# WeStCOMS mesh
mesh_fp <- st_read(glue("{dirs$mesh}/etive28_mesh_footprint.gpkg"))
mesh_sf <- st_read(glue("{dirs$mesh}/etive28_mesh.gpkg")) |> select(i, geom)
site_i <- read_csv("data/etive_farms_2023.csv") |>
  st_as_sf(coords=c("easting", "northing"), crs=27700)

etive_panel <- ggplot() +
  geom_sf(data=mesh_fp, fill="grey", colour="grey30", linewidth=0.2) +
  guides(fill=guide_colourbar(title.position="top", direction="horizontal")) +
  theme_classic() +
  theme(legend.position="bottom",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"))

fig_temp_dir <- "D:/EtiveLice/figs/singleSimTemp/"
timesteps <- str_split_i(dir("D:/EtiveLice/out/fallow/processed/hourly", "_t_"), "_t_", 2) |>
  str_remove(".rds") |>
  unique()
lims_cop <- readRDS("D:/EtiveLice/out/fallow/processed/hourly_pslims_cop.rds") |>
  mutate(mn=pmin(copDens_rtrt, 1)) |>
  mutate(mn_N=pmin(mn^4, 1))
lims_naup <- readRDS("D:/EtiveLice/out/fallow/processed/hourly_pslims_naup.rds") |>
  mutate(mn=pmin(copDens_rtrt, 1)) |>
  mutate(mn_N=pmin(mn^4, 1))
breaks <- list(mn=c(0, 0.01, 0.1, 0.5, 1),
               mn_N=seq(0, 1, by=0.25))

plan(sequential)
library(doFuture)
foreach(i=rev(seq_along(timesteps)),
        .options.future=list(globals=structure(TRUE, add=c("fig_temp_dir")))) %dofuture% {

          timestep <- ymd_hms("2023-03-01 00:00:00") +
            dhours(as.numeric(timesteps[i])-1)

          if(file.exists(glue("{fig_temp_dir}etive_{format(timestep, '%F_%H')}.png"))) {
            next
          }

          ps_naup <- dirf("D:/EtiveLice/out/fallow/processed/hourly",
                          glue("Immature.*_{timesteps[i]}.rds")) |>
            readRDS()
          ps_cop <- dirf("D:/EtiveLice/out/fallow/processed/hourly",
                         glue("Mature.*_{timesteps[i]}.rds")) |>
            readRDS()
          naup_exist <- nrow(ps_naup) > 0
          cop_exist <- nrow(ps_cop) > 0
          if(naup_exist) {
            ps_naup <- ps_naup |> filter(sim_001 > 0) |>
            mutate(mn=pmin(sim_001, lims_naup$mn[2]))
          }
          if(cop_exist) {
            ps_cop <- ps_cop |> filter(sim_001 > 0) |>
              mutate(mn=pmin(sim_001, lims_cop$mn[2]))
          }

          ps_i <- bind_rows(ps_naup, ps_cop) |>
            group_by(i) |>
            summarise(across(where(is.numeric), sum)) |>
            mutate(mn=pmin(sim_001, lims_naup$mn[2] + lims_cop$mn[2]))
          if(nrow(ps_i)==0) {
            next
          }

          if(naup_exist) {
            fig_a <- etive_panel +
              geom_sf(data=mesh_sf |> inner_join(ps_naup), aes(fill=mn), colour=NA) +
              geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
              scale_fill_viridis_c("Nauplii/m2/h", option="turbo", limits=lims_naup$mn,
                                   breaks=breaks$mn^0.25, labels=breaks$mn) +
              ggtitle(format(timestep, "%b-%d %H:%M"))
            ggsave(filename=glue("{fig_temp_dir}etive_naup_{format(timestep, '%F_%H')}.png"),
                   plot=fig_a, width=8, height=5)
            fig_a <- etive_panel +
              geom_sf(data=mesh_sf |> inner_join(ps_naup), aes(fill=mn^4), colour=NA) +
              scale_fill_viridis_c("Nauplii/m2/h", option="turbo", limits=lims_naup$mn_N,
                                   breaks=breaks$mn_N, labels=breaks$mn_N) +
              ggtitle(format(timestep, "%b-%d %H:%M"))
            ggsave(filename=glue("{fig_temp_dir}etive_naup-N_{format(timestep, '%F_%H')}.png"),
                   plot=fig_a, width=8, height=5)
          }
          if(cop_exist) {
            fig_a <- etive_panel +
              geom_sf(data=mesh_sf |> inner_join(ps_cop), aes(fill=mn), colour=NA) +
              geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
              scale_fill_viridis_c("Copepodids/m2/h", option="turbo", limits=lims_cop$mn,
                                   breaks=breaks$mn^0.25, labels=breaks$mn) +
              ggtitle(format(timestep, "%b-%d %H:%M"))
            ggsave(filename=glue("{fig_temp_dir}etive_cop_{format(timestep, '%F_%H')}.png"),
                   plot=fig_a, width=8, height=5)
            fig_a <- etive_panel +
              geom_sf(data=mesh_sf |> inner_join(ps_cop), aes(fill=mn^4), colour=NA) +
              scale_fill_viridis_c("Copepodids/m2/h", option="turbo", limits=lims_cop$mn_N,
                                   breaks=breaks$mn_N, labels=breaks$mn_N) +
              ggtitle(format(timestep, "%b-%d %H:%M"))
            ggsave(filename=glue("{fig_temp_dir}etive_cop-N_{format(timestep, '%F_%H')}.png"),
                   plot=fig_a, width=8, height=5)
          }
          fig_a <- etive_panel +
            geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=mn), colour=NA) +
            geom_sf(data=site_i, colour="violet", shape=1, size=0.5) +
            scale_fill_viridis_c("Lice/m2/h", option="turbo", limits=range(lims_naup$mn + lims_cop$mn),
                                 breaks=breaks$mn^0.25, labels=breaks$mn) +
            ggtitle(format(timestep, "%b-%d %H:%M"))
          ggsave(filename=glue("{fig_temp_dir}etive_all_{format(timestep, '%F_%H')}.png"),
                 plot=fig_a, width=8, height=5)
          fig_a <- etive_panel +
            geom_sf(data=mesh_sf |> inner_join(ps_i), aes(fill=mn^4), colour=NA) +
            scale_fill_viridis_c("Lice/m2/h", option="turbo", limits=range(lims_naup$mn + lims_cop$mn)^4,
                                 breaks=breaks$mn_N, labels=breaks$mn_N) +
            ggtitle(format(timestep, "%b-%d %H:%M"))
          ggsave(filename=glue("{fig_temp_dir}etive_all-N_{format(timestep, '%F_%H')}.png"),
                 plot=fig_a, width=8, height=5)


          gc()
        }

library(av)
sets <- c("etive_all_", "etive_all-N_",
          "etive_cop_", "etive_cop-N_",
          "etive_naup_", "etive_naup-N_")
for(i in sets) {
  dirf(fig_temp_dir, glue("{i}.*png")) |>
    av_encode_video(glue("figs/hourly_anim_{i}_singleSim.mp4"),
                    framerate=12, vfilter="scale=-2:'min(720,ih)'")
}




# Influx
c_long <- map_dfr(dirrf(dirs$out, "connectivity.*csv"),
                  ~load_connectivity(.x,
                                     source_names=src$sepaSite,
                                     dest_names=dest$sepaSite,
                                     liceScale=1) |>
                    mutate(sim=str_sub(str_split_fixed(.x, "sim_", 2)[,2], 1, 3)))

# mean hourly IP for each day
c_daily <- c_long |>
  calc_influx(destination, value, sim, date) |>
  complete(destination, sim, date,
           fill=list(influx=0, N_influx=0)) |>
  mutate(influx_m2=influx/(pi*30^2))


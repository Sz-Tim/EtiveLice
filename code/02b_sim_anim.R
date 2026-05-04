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

set.seed(1001)


# define parameters -------------------------------------------------------

cores_per_sim <- 20
parallel_sims <- 1
start_date <- "2025-03-01"
end_date <- "2025-12-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hf0="/home/sa04ts/hydro/WeStCOMS2/Archive",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v2-2-2.jar",
             dat="/home/sa04ts/EtiveLice/data",
             out=glue("{getwd()}/out/biotracker")),
  windows=list(proj=getwd(),
               mesh="E:/hydro",
               hf0="E:/hydro/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v2-2-2.jar",
               dat="E:/EtiveLice/data",
               out="E:/EtiveLice/out/jpm_anim")
)

sim.i <- read_csv("E:/EtiveLice/out/jpm/sim_i.csv")
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
       nparts=100,
       releaseInterval=1,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh0=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       meshType0="FVCOM",
       hfDir0=glue("{dirs$hf0}/"),
       hfDirPrefix0="netcdf_",
       hfFilePrefix0="westcoms2",
       # sites
       sitefile=glue("{dirs$dat}/pen_sites_widerLinnhe_2025.csv"),
       sitefileEnd=glue("{dirs$dat}/pen_sites_widerLinnhe_2025.csv"),
       siteDensityPath=glue("{dirs$dat}/lice_daily_widerLinnhe_2025-03-01_2025-12-31.csv"),
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
       recordConnectivity="false",
       recordPsteps="true",
       recordImmature="true",
       pstepsInterval=1,
       recordSiteEnv="false",
       recordVertDistr="false",
       recordElemActivity="false"))


# run simulations ---------------------------------------------------------

for(i in sim_seq) {
  setwd(dirs$proj)
  cat("Starting", sim.i$i[i], "\n")
  biotrackR::run_biotracker(
    jdk_path=dirs$jdk,
    jar_path=dirs$jar,
    f_properties=glue::glue("{dirs$out}/sim_{sim.i$i[i]}.properties"),
    sim_dir=glue::glue("{sim.i$outDir[i]}")
  )
}




# process -----------------------------------------------------------------

library(sf)
mesh_sf <- dir(dirs$mesh, "WeStCOMS2_mesh.gpkg", full.names=T) |>
  st_read() |>
  st_crop(c(xmin=100000, xmax=225000, ymin=680000, ymax=810000))
mesh_fp <- dir(dirs$mesh, "WeStCOMS2_meshFootprint.gpkg", full.names=T) |>
  st_read() |>
  st_crop(c(xmin=100000, xmax=225000, ymin=680000, ymax=810000))

# pSteps_df <- dir(dirs$out, "pstepsMature", recursive=T, full.names=T)[-(1:(24*14))][1:(24*7) + 6*30*24] |>
#   map_dfr(~load_psteps(.x, liceScale=1) |>
#             filter(i %in% mesh_sf$i)) |>
#   pivot_longer(starts_with("t_")) |>
#   drop_na() |>
#   mutate(time=ymd(start_date) + hours(str_sub(name, 3, -1)) - hours(1))
site_df <- read_csv(glue("{dirs$dat}/pen_sites_widerLinnhe_2025.csv")) |>
  mutate(sepaSite=str_split_fixed(pen, "_", 2)[,1]) |>
  summarise(easting=mean(easting),
            northing=mean(northing),
            .by=sepaSite)
sim.i <- read_csv(glue("{dirs$out}/sim_i.csv"))


# animation ---------------------------------------------------------------

for(s in 1:nrow(sim.i)) {
  ps_f <- dir(sim.i$outDir[s], "pstepsMature", recursive=T, full.names=T)[-(1:(24*31))]
  nat_breaks <- c(0.001, 0.01, 0.1, 1)
  nat_lims <- c(0.0001, 3)

  for(j in seq_along(ps_f)) {
    ps_j <- load_psteps(ps_f[j], liceScale=1) |>
      filter(i %in% mesh_sf$i) |>
      drop_na() |>
      rename_with(~"value", .cols=2) |>
      mutate(time=ymd(start_date) +
               hours(str_sub(basename(ps_f[j]), 23, -5)) - hours(1),
             time=strftime(time, format="%F %H"))

    ps_sf <- inner_join(mesh_sf, ps_j, by="i")
    p <- ggplot(ps_sf) +
      geom_sf(data=mesh_fp, fill="grey20") +
      geom_sf(aes(fill=pmin(log(value/area)), ), colour=NA) +
      geom_point(data=site_df, aes(easting, northing), colour="white", shape=1) +
      scale_fill_viridis_c(expression("Lice" %.% m^-2),
                           option="turbo",
                           limits=log(nat_lims),
                           breaks=log(nat_breaks),
                           labels=nat_breaks,
                           oob=scales::oob_squish) +
      scale_x_continuous("Longitude", breaks=c(-6, -5), limits=c(120000, 220000),
                         expand=c(0,0), oob=scales::oob_keep) +
      scale_y_continuous("Latitude", breaks=c(56, 56.5, 57), limits=c(685000, 800000),
                         expand=c(0,0), oob=scales::oob_keep) +
      facet_wrap(~time) +
      theme_dark() +
      theme(panel.grid=element_blank(),
            legend.key.height=unit(1, "cm"),
            legend.key.width=unit(0.2, "cm"),
            axis.title=element_blank())
    ggsave(paste0("figs/temp/sim_", s, "_",
                  str_pad(str_sub(basename(ps_f[j]), 23, -5), 4, "left", "0"),
                  ".png"),
           p, width=5, height=4.8)
  }

  anim_f <- dirf("figs/temp", paste0("sim_", s, "_"))
  av::av_encode_video(anim_f, paste0("figs/hourly_anim_sim", s, ".mp4"),
                      framerate=12)
}


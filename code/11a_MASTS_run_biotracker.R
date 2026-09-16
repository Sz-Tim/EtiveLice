# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# MASTS sensitivity analysis talk
# Run simulations


# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR@dev")
library(doFuture)
library(sf)
dirf("code/fn", ".R") |> walk(source)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
set.seed(1003)


# define parameters -------------------------------------------------------

remake_sim_i <- T
os <- get_os()
cores_per_sim <- if_else(os=="linux", 20, 15)
parallel_sims <- if_else(os=="linux", 2, 2)
start_date <- "2025-01-01"
end_date <- "2025-06-30"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

dirs <- switch(
  os,
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/hydro/meshes",
             hf0="/home/sa04ts/hydro/WeStCOMS3/Archive",
             hf2="/home/sa04ts/hydro/WeStCOMS2/SWAN/Archive_daily/",
             jdk="/home/sa04ts/.jdks/jdk-23.0.1/bin/java",
             jar="/home/sa04ts/biotracker/biotracker_v2-2-3.jar",
             out=glue("{getwd()}/out/MASTS_salinity")),
  windows=list(proj="E:/EtiveLice", #getwd(),
               mesh="E:/hydro",
               hf0="E:/hydro/WeStCOMS3/Archive",
               hf2="E:/hydro/WeStCOMS2/SWAN/Archive_daily/",
               jdk="C:/Users/sa04ts/.jdks/openjdk-23.0.2/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/03_packages/biotracker/out/biotracker_v2-2-3.jar",
               out="E:/EtiveLice/out/MASTS_salinity")
)

if(remake_sim_i) {
  post_dir <- glue("{dirs$proj}/data/lit/fit")
  egg_post <- list(constant=read_csv(glue("{post_dir}/egg_temp_constant_post.csv"), show_col_types=F),
                   logistic=read_csv(glue("{post_dir}/egg_temp_logistic_post.csv"), show_col_types=F),
                   linear=read_csv(glue("{post_dir}/egg_temp_linear_post.csv"), show_col_types=F)) |>
    map(~.x |> summarise(across(everything(), mean)) |> mutate(across(everything(), ~signif(.x, 3))) |>
          unite("all", everything(), sep=",") %>%
          .$all)
  mort_post <- list(constant=read_csv(glue("{post_dir}/mort_sal_constant_post.csv"), show_col_types=F),
                    logistic=read_csv(glue("{post_dir}/mort_sal_logistic_post.csv"), show_col_types=F)) |>
    map(~.x |> summarise(across(everything(), mean)) |> mutate(across(everything(), ~signif(.x, 3))) |>
          unite("all", everything(), sep=",") %>%
          .$all)
  sink_post <- read_csv(glue("{post_dir}/sink_sal_linear_post.csv"), show_col_types=F) |>
    summarise(across(everything(), mean)) |>
    mutate(across(everything(), ~.x/1e3)) #mm/s to m/s

  swim_N <- param_seq(min_i=0.05, max_i=5, length.out=11, mean_i=0.5, sd_i=1.25) # mm/s
  swim_C <- param_seq(min_i=0.1, max_i=10, length.out=11, mean_i=1, sd_i=5) # mm/s
  sal_NC <- param_seq(min_i=20, max_i=30, length.out=11, mean_i=25, sd_i=2)

  sim.i <- tibble(
    variableDh=F,
    variableDhV=F,
    D_h=0.1,
    D_hVert=0.001,
    mortSal_fn="logistic",
    eggTemp_fn="logistic",
    swimLightLevel=T,
    lightThreshNauplius=0.4,
    lightThreshCopepodid=2e-5,
    tempPrefNauplius="none",
    swimUpSpeedNaupliusMean=1e-3,
    swimDownSpeedNaupliusMean=1.25e-3,
    swimUpSpeedCopepodidMean=5e-3,
    swimDownSpeedCopepodidMean=1.25e-3,
    passiveSinkRateSal=F,
    salinityThreshNaupliusMax=sal_NC$vals+1e-4,
    salinityThreshNaupliusMin=sal_NC$vals,
    salinityThreshCopepodidMax=sal_NC$vals+1e-4,
    salinityThreshCopepodidMin=sal_NC$vals,
    viableDegreeDays=40,
    maxDepth=5e3,
    connectivityThresh=30,
    eggTemp_b=egg_post$logistic[1],
    mortSal_b=mort_post$logistic[1],
    passiveSinkInt=sink_post$Intercept[1],
    passiveSinkSlope=sink_post$slope[1]
  ) |>
    mutate(across(where(is.numeric), ~signif(.x, 5))) |>
    mutate(i=str_pad(row_number(), nchar(n()), "left", "0"),
           outDir=glue("{dirs$out}/sim_{i}/"))
  write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))
} else {
  sim.i <- read_csv(glue("{dirs$out}/sim_i.csv"))
  if(os=="linux") {
    sim.i$outDir <- str_replace(sim.i$outDir, "E:/EtiveLice/out/MASTS", dirs$out)
  } else {
    sim.i$outDir <- str_replace(sim.i$outDir, "/home/sa04ts/EtiveLice/out/MASTS", dirs$out)
  }
}

sim_seq <- 1:nrow(sim.i)


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
       nparts=25,
       releaseInterval=1,
       checkOpenBoundaries="true",
       # meshes and environment
       mesh0=glue("{dirs$mesh}/WeStCOMS3_mesh.nc"),
       hfDir0=glue("{dirs$hf0}/"),
       hfFilePrefix0="westcoms3",
       hfDir2="",
       hfFilePrefix2="",
       hfDirPrefix2="",
       hfDirSuffix2="",
       # sites
       sitefile=glue("{dirs$proj}/data/farm_sites_widerLinnhe_2024-2025.csv"),
       sitefileEnd=glue("{dirs$proj}/data/farm_sites_widerLinnhe_2024-2025.csv"),
       siteDensityPath=glue("{dirs$proj}/data/lice_daily_widerLinnhe_2025-01-01_2025-12-31.csv"),
       daylightPath=glue("{dirs$proj}/data/day_hours_2025-01-01_2025-12-31.csv"),
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
       swimLightLevel=sim.i$swimLightLevel[.x],
       lightThreshCopepodid=sim.i$lightThreshCopepodid[.x],
       lightThreshNauplius=sim.i$lightThreshNauplius[.x],
       tempPrefNauplius=sim.i$tempPrefNauplius[.x],
       swimUpSpeedCopepodidMean=sim.i$swimUpSpeedCopepodidMean[.x],
       swimUpSpeedCopepodidStd=0,
       swimDownSpeedCopepodidMean=sim.i$swimDownSpeedCopepodidMean[.x],
       swimDownSpeedCopepodidStd=0,
       swimUpSpeedNaupliusMean=sim.i$swimUpSpeedNaupliusMean[.x],
       swimUpSpeedNaupliusStd=0,
       swimDownSpeedNaupliusMean=sim.i$swimDownSpeedNaupliusMean[.x],
       swimDownSpeedNaupliusStd=0,
       passiveSinkingIntercept=sim.i$passiveSinkInt[.x],
       passiveSinkingSlope=sim.i$passiveSinkSlope[.x],
       viableDegreeDays=sim.i$viableDegreeDays[.x],
       connectivityThresh=sim.i$connectivityThresh[.x],
       # recording
       verboseSetUp="true",
       recordConnectivity="true",
       connectImmature="false",
       connectivityInterval=24,
       connectDepth1_min=0,
       connectDepth1_max=5,
       connectDepth2_min=5,
       connectDepth2_max=10,
       connectDepth3_min=10,
       connectDepth3_max=30,
       recordVertDistr="false",
       recordPsteps="true",
       recordImmature="true",
       pstepsInterval=1,
       recordSiteEnv=ifelse(.x==1, "true", "false"),
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
    if(file.exists(glue("{sim.i$outDir[i]}/pstepsMature/pstepsMature_20250630_4344.csv"))) {
      next
    }
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


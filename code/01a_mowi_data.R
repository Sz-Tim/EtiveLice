# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Mowi data cleaning



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR@dev")
library(zoo)
library(janitor)
library(sf)
library(readxl)

mowi_f <- "data/aquaculture/MowiLiceData_LinnheRegion_2022-2025.xlsx"
sepa_key <- c("APT1"="Etive 4",
              "ARDG1"="Ardgour",
              "CAG1"="Camas Glas",
              "CALL1"="Leven",
              "FFMC84"="Etive 3",
              "GORS1"="Gorsten",
              "INV1"="Invasion Bay",
              "KING1"="Kingairloch",
              "MCLN1"="Macleans Nose",
              "SAR1"="Etive 6")


# cleaning ----------------------------------------------------------------

mowi_df <- full_join(
  imap_dfr(
    sepa_key,
    ~read_xlsx(mowi_f, sheet=.x,
               col_types=c(rep("skip", 8), "text", rep("numeric", 3), "skip", rep("numeric", 2)),
               col_names=c("date", "AFg", "AFn", "AM", "Ch", "PA"),
               skip=1) |>
      mutate(sepaSite=.y,
             date=ymd(date),
             AF=AFg + AFn)) |>
    drop_na(date),
  imap_dfr(
    sepa_key,
    ~read_xlsx(mowi_f, sheet=.x,
               col_types=c("skip", "text", rep("numeric", 3), rep("skip", 10)),
               col_names=c("date", "nFish", "nFishSampled", "biomass"),
               skip=1) |>
      mutate(sepaSite=.y,
             date=ymd(date))),
  by=join_by(sepaSite, date)
  ) |>
  complete(sepaSite, date=seq(first(date), last(date), by=1), fill=list(nFish=0, biomass=0))

write_csv(mowi_df, "data/aquaculture/mowi_cleaned.csv")



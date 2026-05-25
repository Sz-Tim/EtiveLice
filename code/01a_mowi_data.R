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



# treatments --------------------------------------------------------------

mowi_f <- "data/aquaculture/Mowi_LiceTreatments_LinnheRegion_2022-2025_COPY.xlsx"

test <- read_xlsx(mowi_f, sheet=sepa_key[1])

mowi_trt <- imap_dfr(
  sepa_key,
  ~read_xlsx(mowi_f, sheet=.x) |>
    mutate(sepaSite=.y,
           date=dmy(Date),
           duration=End-Start)) |>
  select(where(~any(!is.na(.x)))) |>
  rename(PenName=`Pen Name`,
         DoseUnit=`Unit...12`,
         ActiveSubstance=`Active Substance`,
         ActiveSubstanceAmount=`Active substance amount`,
         ActiveSubstanceUnit=`Unit...17`,
         TreatmentType=heading) |>
  select(sepaSite, PenName, date, Method, TreatmentType,
         starts_with("Active"), starts_with("Dose"), duration) |>
  mutate(Type=paste0(Method, "_", str_to_camel(ActiveSubstance))) |>
  arrange(sepaSite, date, PenName, Type) |>
  group_by(sepaSite) |>
  mutate(maxPens=n_distinct(PenName)) |>
  ungroup() |>
  arrange(Method, Type) |>
  mutate(Method=factor(Method),
         Type=factor(Type)) |>
  arrange(sepaSite, date, PenName, Type)
mowi_trt_site <- mowi_trt |>
  summarise(nPens=n_distinct(PenName),
            propPens=nPens/max(maxPens),
            totSubstance=sum(ActiveSubstanceAmount),
            totDose=sum(Dose),
            totDuration=sum(duration),
            .by=c(sepaSite, date, Method, ActiveSubstance, Type, ends_with("Unit"))) |>
  mutate(MethodNum=as.numeric(Method),
         TypeNum=as.numeric(Type))

write_csv(mowi_trt_site, "data/aquaculture/mowi_trt_cleaned.csv")


mowi_trt_site |>
  count(Method, ActiveSubstance, Type)



mowi_trt_site |>
  ggplot(aes(sepaSite, fill=Type)) +
  geom_bar(position="fill", colour="grey30") +
  scale_fill_brewer(type="qual", palette="Paired")





logit_trt_pop_mu_sd <- c(0, 1)
logit_trt_grp_sd <- 1
logit_trt_substance_sd <- 1


# Parameters: logit_psi_pop, sigma_grps, sigma_types
# transformed_parameters: logit_psi_grp, logit_psi_type

# hyper prior: logit(psi[pop]) ~ N(mu, sigma)
trt_pars_pop <- tibble(level="pop",
                       name="pop",
                       logit_psi=rnorm(1, logit_trt_pop_mu_sd[1], logit_trt_pop_mu_sd[2]),
                       psi=boot::inv.logit(logit_psi))
# group hierarchy: logit(psi[grp]) ~ N(logit(psi[pop]), sigma_grps)
trt_pars_grp <- tibble(level="grp",
                       name=unique(mowi_trt_siteSum$Method),
                       logit_psi=rnorm(4, trt_pars_pop$logit_psi[1], logit_trt_grp_sd),
                       psi=boot::inv.logit(logit_psi),
                       Method=name)
# type hierarchy: logit(psi[type]) ~ N(logit(psi[grp]), sigma_types)
trt_pars_substance <- mowi_trt_siteSum |>
  summarise(.by=c(Method, ActiveSubstance, MethodSubstance)) |>
  mutate(logit_psi=rnorm(n(),
                         trt_pars_grp$logit_psi[match(Method, trt_pars_grp$name)],
                         logit_trt_substance_sd),
         psi=boot::inv.logit(logit_psi),
         level="substance") |>
  rename(name=MethodSubstance) |>
  select(level, Method, name, logit_psi, psi)

bind_rows(
  trt_pars_pop,
  trt_pars_grp,
  trt_pars_substance
) |>
  mutate(level=factor(level, levels=c("pop", "grp", "substance"))) |>
  ggplot(aes(psi, Method, size=level)) +
  geom_point(shape=1) +
  scale_size_manual(values=5:2) +
  xlim(0,1) +
  theme_bw()





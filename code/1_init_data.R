# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Initial data



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(zoo)
library(janitor)
library(sf)
library(raster)
library(spaths)
library(readxl)
conflicted::conflicts_prefer(dplyr::select,
                             dplyr::filter,
                             tidyr::extract)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
# source("code/00_fn.R")



# load data ---------------------------------------------------------------

mesh.fp <- st_read("../03_packages/WeStCOMS/data/WeStCOMS2_meshFootprint.gpkg")
farm_bbox <- list(xmin=125000, xmax=225500, ymin=629000, ymax=800000)

# Biomass is managed by SEPA, while sea lice is managed by MSS
# Separate site names are used, and there is not always total overlap
# SEPA sites are directly connected to licences, so they are used as the
# locations, with the nearest sea lice counts applied to each

# Download data
download_lice_counts("2023-01-01", "2023-12-31", "data/ms_sea_lice_2021-present.csv")
download_ms_site_details("data/ms_site_details.csv")
download_sepa_licenses("data/se_licence_conditions.csv")

# Process datasets
sepa.locs <- read_csv("data/se_licence_conditions.csv") |>
  filter(salmon==TRUE & waterType=="Seawater") |>
  filter(between(easting, farm_bbox$xmin, farm_bbox$xmax)) |>
  filter(between(northing, farm_bbox$ymin, farm_bbox$ymax)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  mutate(dist=as.numeric(st_distance(., mesh.fp)[,1])) |>
  filter(dist < 1000) |>
  rename(sepaSite=sepaSiteId) |>
  mutate(operator=if_else(operator=="Wester Ross Salmon Ltd", "Wester Ross Fisheries Ltd", operator)) |>
  st_drop_geometry() |>
  # relocate sites that are outside WeStCOMS2 mesh due to coastal simplification
  mutate(easting=if_else(sepaSite=="RONAR1", 161050, easting),
         northing=if_else(sepaSite=="RONAR1", 855080, northing)) |>
  mutate(easting=if_else(sepaSite=="PLOI1", 206760, easting),
         northing=if_else(sepaSite=="PLOI1", 916230, northing)) |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) |>
  select(sepaSite, operator, easting, northing, geometry)
mss.locs <- read_csv("data/ms_site_details.csv") |>
  filter(aquacultureType=="Fish" & waterType=="Seawater") |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  mutate(dist=as.numeric(st_distance(., mesh.fp)[,1])) |>
  filter(dist < 1000) |>
  rename(mssID=marineScotlandSiteId) |>
  mutate(operator=if_else(operator=="Bakkafrost Scotland", "Bakkafrost Scotland Ltd", operator),
         operator=if_else(operator=="Landcatch Natural Selection Ltd", "Landcatch Ltd", operator)) |>
  select(mssID, operator, easting, northing, geometry)
sites.i <- sepa.locs %>%
  mutate(mssID=mss.locs$mssID[st_nearest_feature(., mss.locs)]) |>
  left_join(mss.locs |>
              st_drop_geometry() |>
              select(mssID, operator) |>
              rename(mss_operator=operator))
sites.i$dist <- apply(st_distance(sites.i, mss.locs), 1, min)
sites_validation <- sites.i |>
  filter(operator == mss_operator) |>
  slice_min(dist, by=mssID) |>
  select(sepaSite, mssID, easting, northing, geometry)

# TODO: update this --- copied from HoloLice....
out.df <- read_csv("data/lice_biomass_2017-01-01_2024-06-30.csv")


# save output -------------------------------------------------------------

linnhe_farms <- c("APT1", "ARDG1", "CALL1", "FFMC19", "FFMC20", "FFMC27",
                  "FFMC32", "FFMC40A", "FFMC40B", "FFMC41", "FFMC56", "FFMC60",
                  "FFMC84", "GORS1", "KING1", "SAR1", "SHUI1")
out.df |>
  filter(sepaSite %in% linnhe_farms) |>
  select(sepaSite, date, weeklyAverageAf, actualBiomassOnSiteTonnes, total_AF) |>
  write_csv(glue("data/lice_biomass_{min(out.df$date)}_{max(out.df$date)}.csv"))
out.df |>
  filter(sepaSite %in% linnhe_farms) |>
  filter(date >= "2023-01-01") |>
  filter(date <= "2023-12-31") |>
  group_by(sepaSite) |>
  filter(any(actualBiomassOnSiteTonnes > 0)) |>
  summarise() |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv("data/farm_sites.csv")
out.df |>
  filter(sepaSite %in% linnhe_farms) |>
  filter(date >= "2023-01-01") |>
  filter(date <= "2023-12-31") |>
  group_by(sepaSite) |>
  mutate(actualBiomassOnSiteTonnes=first(actualBiomassOnSiteTonnes)) |>
  ungroup() |>
  mutate(total_AF=actualBiomassOnSiteTonnes * weeklyAverageAf * 240,
         date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2023-01-01_2023-12-31.csv"))





# literature --------------------------------------------------------------

library(brms)
theme_set(theme_bw())

## EGG PRODUCTION
# Kragesteen 2023: Egg = f(temperature)
egg_df <- read_csv("data/lit/Kragesteen2023_Table1.csv") |>
  filter(temperature > 4)

# linear
egg_temp_linear <- brm(eggs_day_AG ~ temperature, data=egg_df, cores=4,
                       prior=c(prior(normal(0, 10), class=sigma)))
egg_temp_linear |> saveRDS("data/lit/fit/egg_temp_linear_brmsfit.rds")
as_draws_df(egg_temp_linear) |>
  rename(a=b_Intercept,
         b=b_temperature) |>
  select(a, b) |>
  write_csv("data/lit/fit/egg_temp_linear_post.csv")

# logistic
egg_temp_logistic <- brm(bf(eggs_day_AG ~ eggMax/(1+exp(-k*(temperature - temp0))) + b, nl=T) +
                           lf(eggMax ~ 1, k~1, temp0~1, b~1),
                         data=egg_df, cores=4,
                         prior=c(prior(normal(100, 10), nlpar="eggMax"),
                                 prior(normal(0, 1), nlpar="k"),
                                 prior(normal(10, 10), nlpar="temp0"),
                                 prior(normal(10, 10), nlpar="b", lb=0),
                                 prior(normal(0, 10), class=sigma)))
egg_temp_logistic |> saveRDS("data/lit/fit/egg_temp_logistic_brmsfit.rds")
as_draws_df(egg_temp_logistic) |>
  rename(eggMax=b_eggMax_Intercept,
         k=b_k_Intercept,
         temp0=b_temp0_Intercept,
         b=b_b_Intercept) |>
  select(eggMax, k, temp0, b) |>
  write_csv("data/lit/fit/egg_temp_logistic_post.csv")

# quadratic
egg_temp_quadratic <- brm(bf(eggs_day_AG ~ a * (b + temperature)^2, nl=T) +
                        lf(a ~ 1, b ~ 1),
                      prior=c(prior(normal(0.17, 3), nlpar="a"),
                              prior(normal(4.28, 3), nlpar="b"),
                              prior(normal(0, 3), class=sigma)),
                      data=egg_df, cores=4)
egg_temp_quadratic |> saveRDS("data/lit/fit/egg_temp_quadratic_brmsfit.rds")
as_draws_df(egg_temp_quadratic) |>
  rename(a=b_a_Intercept, b=b_b_Intercept) |>
  select(a, b) |>
  write_csv("data/lit/fit/egg_temp_quadratic_post.csv")






egg_post <- read_csv("data/lit/fit/egg_temp_logistic_post.csv")
temp_seq <- seq(3, 20, by=0.1)
plot(NA, NA, xlim=c(3, 20), ylim=c(0, 150),
     xlab="Temperature", ylab="Eggs per gravid female per day")
for(i in 1:300) {
  lines(temp_seq, egg_post$eggMax[i]/(1+exp(-egg_post$k[i]*(temp_seq - egg_post$temp0[i]))) + egg_post$b[i], col=rgb(0,0,0,0.1))
}
lines(temp_seq, 0.17 * (4.28 + temp_seq)^2, col="red")
lines(temp_seq, median(egg_post$eggMax)/(1+exp(-median(egg_post$k)*(temp_seq - median(egg_post$temp0)))) + median(egg_post$b), col="black", lwd=2)
points(egg_df$temperature, egg_df$eggs_day_AG, col="blue")



plot(NA, NA, xlim=c(3, 20), ylim=c(0, 150),
     xlab="Temperature", ylab="Eggs per gravid female per day")
for(i in 1:300) {
  lines(temp_seq, egg_post$a[i] * (egg_post$b[i] + temp_seq)^2, col=rgb(0,0,0,0.2))
}
lines(temp_seq, 0.17 * (4.28 + temp_seq)^2, col="red")
points(egg_df$temperature, egg_df$eggs_day_AG, col="blue")


egg_post <- read_csv("data/lit/fit/egg_temp_linear_post.csv")
plot(NA, NA, xlim=c(3, 20), ylim=c(0, 150),
     xlab="Temperature", ylab="Eggs per gravid female per day")
for(i in 1:300) {
  lines(temp_seq, egg_post$a[i] + egg_post$b[i] * temp_seq, col=rgb(0,0,0,0.1))
}
lines(temp_seq, 0.17 * (4.28 + temp_seq)^2, col="red")
lines(temp_seq, mean(egg_post$a) + mean(egg_post$b) * temp_seq, col="black", lwd=2)
points(egg_df$temperature, egg_df$eggs_day_AG, col="blue")




## MORTALITY
# assume survival times follow an exponential distribution
# hazard rate = hourly mortality rate = ln(2)/LT50
# sal = 5, 9, 12 are listed as LT50: <1h
# "At 12 ppt and below the initial death rate was rapid, with all copepodids in
# 9 and 5 ppt dying within the first 2 h."
sal_df <- read_csv("data/lit/Bricknell2006_Table1.csv") |>
  mutate(mort_h=-log(0.5)/LT50_h) |>
  mutate(mort_h=case_when(salinity==5 ~ 0.95,
                          salinity==9 ~ 0.9,
                          .default=mort_h),
         logit_mort_h=boot::logit(mort_h))

# Logistic
mort_sal_logistic <- brm(bf(mort_h ~ mortMax/(1+exp(-k*(salinity - sal0))) + b, nl=T) +
                      lf(mortMax ~ 1, k~1, sal0~1, b~1),
                    data=sal_df, cores=4,
                    prior=c(prior(normal(1, 0.1), nlpar="mortMax", lb=0, ub=1),
                            prior(normal(0, 3), nlpar="k"),
                            prior(normal(15, 4), nlpar="sal0"),
                            prior(normal(0.01, 0.1), nlpar="b", lb=0, ub=1),
                            prior(normal(0, 0.4), class=sigma)))
mort_sal_logistic |> saveRDS("data/lit/fit/mort_sal_logistic_brmsfit.rds")
as_draws_df(mort_sal_logistic) |>
  rename(mortMax=b_mortMax_Intercept,
         k=b_k_Intercept,
         sal0=b_sal0_Intercept,
         b=b_b_Intercept) |>
  select(mortMax, k, sal0, b) |>
  write_csv("data/lit/fit/mort_sal_logistic_post.csv")





mort_post <- read_csv("data/lit/fit/mort_sal_logistic_post.csv")
sal_seq <- seq(5, 35, by=0.1)

plot(NA, NA, xlim=c(5, 35), ylim=c(0, 1),
     xlab="Salinity", ylab="Hourly mortality")
for(i in 1:300) {
  lines(sal_seq, mort_post$mortMax[i]/(1+exp(-mort_post$k[i]*(sal_seq - mort_post$sal0[i]))) + mort_post$b[i], col=rgb(0,0,0,0.1))
}
lines(sal_seq, median(mort_post$mortMax)/(1+exp(-median(mort_post$k)*(sal_seq - median(mort_post$sal0)))) + median(mort_post$b), col="red", lwd=2)
lines(sal_seq, mean(mort_post$mortMax)/(1+exp(-mean(mort_post$k)*(sal_seq - mean(mort_post$sal0)))) + mean(mort_post$b), col="red3", lwd=2)
points(sal_df$salinity, sal_df$mort_h, col="blue")



# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Fit integrated population model


# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(cmdstanr)
library(ggdist)
library(cowplot)
library(doFuture)
fread <- data.table::fread
dir("code/fn", ".R", full.names=T) |> walk(source)
theme_set(theme_classic())

prior_only <- F
GQ <- F

mod <- c("noHarm", "Harm",
         "Harm_randIPbg", "Harm_ydayIPbg", "noHarm_ydayIPbg")[1]
fishCol <- c("RW_logit", "BSA")[2]
suffix <- paste0(switch(mod,
                        'noHarm'='_noHarm',
                        'Harm'='',
                        'Harm_randIPbg'='_randIPbg',
                        'Harm_ydayIPbg'='_ydayIPbg',
                        'noHarm_ydayIPbg'='_noHarm_ydayIPbg'),
                 ifelse(prior_only, '_PRIORS', ''),
                 "_", fishCol,
                 ifelse(GQ, '', '_noGQ'))

n_chains <- 3
dat_dir <- "data/aquaculture/mowi_stan/"
out_dir <- "out/ipm_fit/"
fig_dir <- "figs/ipm_fit/"
info <- readRDS(glue("{dat_dir}info.rds"))

stages <- c("Ch1", "Ch2", "PA1", "PA2", "Ad")
stageGrps <- c("Ch", "PA", "Ad")
stage_trans <- c("Ch-PA", "PA-Ad")
trt_meth_ii <- read_csv("data/aquaculture/mowi_trt_cleaned.csv") |>
  summarise(.by=c(MethodNum, TypeNum, Method, Type)) |>
  arrange(TypeNum) |>
  mutate(abbr=paste0(str_sub(Method, 1, 1), "_", str_split_i(Type, "_", 2)))
param_key <- tibble(name=c({if(fishCol=="BSA") {
                             paste0("attach_beta[", 1:6, "]")
                           } else {
                             paste0("attach_beta[", 1:5, "]")
                           }},
                           {if(grepl("noHarm", mod)) {
                             paste0("ensWts_p[", 1:info$nSims, "]")
                           } else {
                             c(paste0("ensWts_harm[1,", 1:info$nSims, "]"),
                               paste0("ensWts_harm[2,", 1:info$nSims, "]"),
                               paste0("ensWts_harm[3,", 1:info$nSims, "]"))
                           }},
                           {if(grepl("randIP", mod)) {
                             c(paste0("IP_bg[", 1:info$nFarms, "]"),
                               paste0("IP_bg_m3[", 1:info$nFarms, "]"))
                           } else if(grepl("ydayIP", mod)) {
                             c(paste0("log_IP_bg_m3_coef[1,", 1:info$nFarms, "]"),
                               paste0("log_IP_bg_m3_coef[2,", 1:info$nFarms, "]"),
                               paste0("log_IP_bg_m3_coef[3,", 1:info$nFarms, "]"))
                           } else {
                             c("IP_bg", "IP_bg_m3")
                           }},
                           paste0("surv_beta[1,", 1:3, "]"),
                           paste0("surv_beta[2,", 1:3, "]"),
                           paste0("surv_int_farm_sd[", 1:3, "]"),
                           paste0("mnDaysStage_beta[1,", 1:2, "]"),
                           paste0("mnDaysStage_beta[2,", 1:2, "]"),
                           paste0("detect_p[", 1:2, "]"),
                           "nb_prec",
                           "IP_scale", "IP_halfSat_m3",
                           paste0("trtEff_type[", 1:8, "]")),
                    label=c({if(fishCol=="BSA") {
                              paste0("attach_", c("Int", "Fish", "Sal", "Temp", "UV", "UVsq"))
                            } else {
                              paste0("attach_", c("Fish", "Sal", "Temp", "UV", "UVsq"))
                            }},
                            {if(grepl("noHarm", mod)) {
                              paste0("ensWts_", 1:info$nSims)
                            } else {
                              c(paste0("ensWt_Int_", 1:info$nSims),
                                paste0("ensWt_cos_", 1:info$nSims),
                                paste0("ensWt_sin_", 1:info$nSims))
                            }},
                            {if(grepl("randIP", mod)) {
                              c(paste0("IP_bg[", 1:info$nFarms, "]"),
                                paste0("IP_bg_m3[", 1:info$nFarms, "]"))
                            } else if(grepl("ydayIP", mod)) {
                              c(paste0("logIPbg_Int_", 1:info$nFarms),
                                paste0("logIPbg_cos_", 1:info$nFarms),
                                paste0("logIPbg_sin_", 1:info$nFarms))
                            } else {
                              c("IP_bg", "IP_bg_m3")
                            }},
                            paste0("surv_Int_", stageGrps),
                            paste0("surv_Sal_", stageGrps),
                            paste0("surv_int_farm_sd_", stageGrps),
                            paste0("mnDaysStage_Int_", stageGrps[1:2]),
                            paste0("mnDaysStage_Temp_", stageGrps[1:2]),
                            paste0("p_detect_", stageGrps[-3]),
                            "neg_binom_prec",
                            "IP_scale", "IP_halfSat_m3",
                            paste0("trtEff_", trt_meth_ii$abbr)
                    )) |>
  mutate(label=factor(label, levels=unique(label)))

keep_pars <- c(ifelse(grepl("ydayIP", mod), "log_IP_bg_m3_coef", "IP_bg_m3"),
               ifelse(grepl("noHarm", mod), "ensWts_p", "ensWts_harm"),
               "attach_beta","surv_beta", "surv_int_farm_sd", "mnDaysStage_beta",
               "detect_p", "nb_prec", "trtEff_type",
               "mu", "mu_GQ", "y_pred", "y_bar_GQ",
               "log_lik", "log_lik_GQ"
)
if(!GQ) keep_pars <- grep("_GQ", keep_pars, invert=T, value=T)



# fit ---------------------------------------------------------------------

iter <- 1000
if(grepl("ydayIP", mod)) {
  prior_ls <- list(prior_IP_bg_m3=c(log(0.05), 0.5))
} else {
  prior_ls <- NULL
}
stan_dat <- make_stan_data(dat_dir, priors_only=prior_only, source="real",
                           GQ_start=switch(GQ+1, NULL, "2025-01-01"),
                           fishCol=fishCol,
                           prior_ls=prior_ls)

mod_full <- cmdstan_model(glue("code/stan/tuning_integrated_population_model{str_remove(suffix, '_PRIORS') |> str_remove(paste0('_', fishCol)) |> str_remove('_noGQ')}.stan"))
fit_full <- mod_full$sample(
  data=stan_dat$dat, init=0, seed=101, refresh=max(iter/100, 1),
  adapt_delta=0.95,
  iter_warmup=iter*2, iter_sampling=iter,
  chains=n_chains, parallel_chains=n_chains
)

out_full_df <- fit_full$draws(
  variables=keep_pars,
  format="df") |>
  pivot_longer(-starts_with("."))
saveRDS(out_full_df, glue("{out_dir}posterior{suffix}.rds"))
out_full_sum <- out_full_df |>
  group_by(name) |>
  sevcheck::get_intervals(value, type="qi")
saveRDS(out_full_sum, glue("{out_dir}posterior_summary{suffix}.rds"))



# log likelihood ----------------------------------------------------------

if(GQ) {
  p_loglik <- out_full_sum |>
    filter(grepl("log_lik_GQ", name)) |>
    mutate(stage=str_sub(name, 12, 12)) |>
    ggplot(aes(med)) +
    geom_histogram() +
    facet_grid(.~stage, scales="free")
  ggsave(glue("{fig_dir}/fig_loglik{suffix}.png"),
         p_loglik, width=9, height=3.5)
}

p_loglik <- out_full_sum |>
  filter(grepl("log_lik\\[", name)) |>
  mutate(stage=str_sub(name, 9, 9)) |>
  ggplot(aes(med)) +
  geom_histogram() +
  facet_grid(.~stage, scales="free")
ggsave(glue("{fig_dir}/fig_loglik_fitted{suffix}.png"),
       p_loglik, width=9, height=3.5)



# parameter distributions -------------------------------------------------

if(grepl("noHarm", mod)) {
  p_ensWts <- list(out_full_df, out_full_sum) |>
    map(~.x |> filter(grepl("ensWts", name)) |> inner_join(param_key, by=join_by(name))) |>
    post_summary_plot(scales="fixed", ncol=if_else(info$nSims > 10, 10, 5))
} else {
  p_ensWts <- list(out_full_df, out_full_sum) |>
    map(~.x |> filter(grepl("ensWts", name)) |> inner_join(param_key, by=join_by(name))) |>
    post_summary_ensWt_plot(scales="fixed", ncol=if_else(info$nSims > 10, 10, 5))
}
p_attach <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("attach_beta", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(scales="free") +
  geom_vline(xintercept=0, linetype=3)
p_surv <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("surv_beta", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=3, scales="free") +
  geom_vline(xintercept=0, linetype=3)
p_surv_sd <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("surv_int_farm_sd", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=3, scales="free") +
  geom_vline(xintercept=0, linetype=3)
p_trt <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("trtEff_type", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=4, scales="free_y") +
  xlim(0, 1)
p_pMoltTemp <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("mnDaysStage_beta", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=2, scales="free_y") +
  geom_vline(xintercept=0, linetype=3)
if(grepl("ydayIP", mod)) {
  p_detectp <- list(out_full_df, out_full_sum) |>
    map(~.x |> filter(grepl("detect_p", name)) |> inner_join(param_key, by=join_by(name))) |>
    post_summary_plot(scales="free_y", ncol=1, nrow=2) +
    xlim(0, 1)
  p_else <- plot_grid(
    list(out_full_df, out_full_sum) |>
      map(~.x |> filter(grepl("log_IP_bg_m3_coef", name)) |> inner_join(param_key, by=join_by(name))) |>
      post_summary_IPbg_plot(ncol=5, scales="fixed"),
    list(out_full_df, out_full_sum) |>
      map(~.x |> filter(grepl("nb_prec", name)) |> inner_join(param_key, by=join_by(name)) |>
            mutate(label=str_split_i(label, "\\[", 1))) |>
      post_summary_plot(ncol=1, scales="free"),
    nrow=1, rel_widths=c(1, 0.2))
  p <- plot_grid(p_ensWts, p_attach, p_surv, p_surv_sd, p_trt,
                 plot_grid(p_pMoltTemp, p_detectp, nrow=1, axis="tblr", align="hv", rel_widths=c(2,1)),
                 p_else, nrow=7, align="v", axis="rl", rel_heights=c(2, 1, 2, 1, 2, 2, 2))
} else {
  p_detectp <- list(out_full_df, out_full_sum) |>
    map(~.x |> filter(grepl("detect_p", name)) |> inner_join(param_key, by=join_by(name))) |>
    post_summary_plot(scales="free_y") +
    xlim(0, 1)
  p_else <- list(out_full_df, out_full_sum) |>
    map(~.x |> filter(grepl("IP_bg|nb_prec", name)) |> inner_join(param_key, by=join_by(name)) |>
          mutate(label=str_split_i(label, "\\[", 1))) |>
    post_summary_plot(ncol=5, scales="free")
  p <- plot_grid(p_ensWts, p_attach, p_surv, p_surv_sd, p_trt, p_pMoltTemp,
                 plot_grid(p_detectp, p_else, nrow=1, axis="tblr", align="hv"),
                 nrow=7, align="v", axis="rl", rel_heights=c(2, 1, 2, 1, 2, 2, 1))
}

ggsave(glue("{fig_dir}/fig_pars{suffix}.png"), p, width=10, height=14)



obs_df <- read_csv("data/aquaculture/mowi_cleaned.csv") |>
  filter(!is.na(nFishSampled)) |>
  mutate(farm=paste("Farm", as.numeric(factor(sepaSite))))

if(GQ) {
  mu_draws_df <- take_mu_draws(out_full_df, NULL,
                               stan_dat$dat, ndraws=min(1e2, iter), GQ=TRUE) |>
    drop_na(mu) # some prior draws give NAs because of negbinom constraints
  p <- mu_draws_df |>
    filter(stage=="Ad") |>
    ggplot() +
    geom_line(aes(day, mu, group=as.character(.draw)), alpha=0.1) +
    geom_point(data=obs_df |> filter(between(date, min(mu_draws_df$day), max(mu_draws_df$day))),
               aes(date, AF/nFishSampled), shape=1, colour="steelblue3") +
    labs(x="Date", y="Mean lice per fish (latent)") +
    {if(any((mu_draws_df |> filter(stage=="Ad"))$mu > 15)) scale_y_continuous(limits=c(0, 15), oob=scales::oob_keep)} +
    scale_x_date(date_labels="%b") +
    facet_grid(farm~., scales="free_y")
  ggsave(glue("{fig_dir}/fig_mu_draws_GQ{suffix}.png"), p, width=10, height=15)

  # mu_pred_df <- mu_draws_df |>
  #   summarise(mn=mean(mu),
  #             lo=quantile(mu, probs=0.25),
  #             hi=quantile(mu, probs=0.75),
  #             .by=c(farm, stage, day)) |>
  #   inner_join(obs_df |>
  #                mutate(across(all_of(c("Ch", "PA", "AF")),  ~.x/nFishSampled)) |>
  #                select(farm, date, Ch, PA, AF) |>
  #                rename(day=date, Ad=AF) |>
  #                pivot_longer(3:5, names_to="stage", values_to="obs"),
  #              by=join_by(stage, farm, day))
  # ggplot(mu_pred_df, aes(mn, xmin=lo, xmax=hi, y=obs)) +
  #   geom_point(shape=1) +
  #   geom_linerange() +
  #   geom_abline() +
  #   facet_grid(farm~stage, scales="free")
  # ggplot(mu_pred_df, aes(mn, xmin=lo, xmax=hi, y=obs)) +
  #   geom_point(shape=1) +
  #   geom_linerange() +
  #   geom_abline() +
  #   facet_wrap(~stage, scales="free")
  # ggplot(mu_pred_df |> filter(stage=="Ad"), aes(mn, xmin=lo, xmax=hi, y=obs)) +
  #   geom_point(shape=1) +
  #   geom_linerange() +
  #   geom_abline() +
  #   facet_wrap(~farm, scales="free")
}

mu_draws_df <- take_mu_draws(out_full_df, NULL,
                             stan_dat$dat, ndraws=min(1e2, iter), GQ=F) |>
  drop_na(mu) # some prior draws give NAs because of negbinom constraints
p <- mu_draws_df |>
  filter(stage=="Ad") |>
  ggplot() +
  geom_line(aes(day, mu, group=as.character(.draw)), alpha=0.1) +
  geom_point(data=obs_df, aes(date, AF/nFishSampled), shape=1, colour="steelblue3") +
  labs(x="Date", y="Mean lice per fish (latent)") +
  {if(any((mu_draws_df |> filter(stage=="Ad"))$mu > 15)) scale_y_continuous(limits=c(0, 15), oob=scales::oob_keep)} +
  scale_x_date(date_labels="%b") +
  facet_grid(farm~., scales="free_y")
ggsave(glue("{fig_dir}/fig_mu_draws{suffix}.png"), p, width=10, height=15)

if("y_pred" %in% keep_pars) {
  sampledDays <- readRDS(glue("{dat_dir}/sampledDays.rds")) |>
    as_tibble() |>
    rename(farm=sepaSite) |>
    mutate(farm=as.character(farm),
           sample=row_number()) |>
    inner_join(readRDS(glue("{dat_dir}/nFishSampled_mx.rds")) |>
                 as_tibble() |>
                 mutate(day=row_number()) |>
                 pivot_longer(-day, names_to="farm", values_to="nFishSampled") |>
                 mutate(farm=str_sub(farm, 2, -1)))

  # TODO: These aren't lined up correctly. Farms are missing from the predictions
  y_df <- out_full_df |>
    filter(grepl("y_pred\\[[1-3]", name)) |>
    separate_wider_delim(name, delim=",", names=c("stage", "sample")) |>
    mutate(sample=as.numeric(str_sub(sample, 1, -2)),
           stage=factor(stage, levels=paste0("y_pred[", 1:3), labels=c("Ch", "PA", "Ad"))) |>
    group_by(.chain, .iteration, .draw, stage, sample) |>
    summarise(value=sum(value, na.rm=T)) |>
    group_by(stage, sample) |>
    sevcheck::get_intervals(value, type="qi") |>
    ungroup() |>
    left_join(sampledDays, by=join_by(sample)) |>
    select(-sample) |>
    mutate(type="Fitted") |> rename(y=med) |>
    bind_rows(expand_grid(farm=as.character(1:info$nFarms),
                          day=1:info$nDays,
                          stage=factor(c("Ch", "PA", "Ad"), levels=c("Ch", "PA", "Ad"))) |>
                mutate(y=c(readRDS(glue("{dat_dir}/y_obs.rds")))) |>
                mutate(type="Obs") |>
                inner_join(sampledDays |> select(-sample), by=join_by(farm, day))) |>
    mutate(day=ymd("2025-01-01") + day - 1,
           farm=paste("Farm", farm),
           y_perFish=y/nFishSampled) |>
    filter(nFishSampled > 0)
  y_df |>
    saveRDS(glue("{dat_dir}/y_fitted{suffix}.rds"))

  p <- y_df |>
    ggplot(aes(day, y_perFish)) +
    geom_linerange(aes(ymin=L10/nFishSampled, ymax=L90/nFishSampled, group=type), linewidth=0.25) +
    geom_point(aes(colour=type, shape=type)) +
    scale_colour_manual("", values=c("Obs"="blue", "Fitted"="black")) +
    scale_shape_manual("", values=c("Obs"=4, "Fitted"=1)) +
    labs(x="Date", y="Mean lice per fish (observed) [50% CI]") +
    scale_x_date(date_labels="%b") +
    facet_grid(stage~farm, scales="free_y")
  ggsave(glue("{fig_dir}/fig_y{suffix}.png"), p, width=15, height=7)

  p <- y_df |>
    filter(stage=="Ad") |>
    ggplot(aes(day, y_perFish)) +
    geom_linerange(aes(ymin=L10/nFishSampled, ymax=L90/nFishSampled, group=type), linewidth=0.25) +
    geom_point(aes(colour=type, shape=type)) +
    scale_colour_manual("", values=c("Obs"="blue", "Fitted"="black")) +
    scale_shape_manual("", values=c("Obs"=4, "Fitted"=1)) +
    labs(x="Date", y="Mean adult female lice per fish (observed) [50% CI]") +
    scale_x_date(date_labels="%b") +
    facet_grid(farm~., scales="free_y")
  ggsave(glue("{fig_dir}/fig_y_AF{suffix}.png"), p, width=6, height=15)
}




# FOR A NEW SCRIPT --------------------------------------------------------

# RESULTS
# Performance using 2025 GQ predictions
# - Bayes R2
# - Bayes Factor over priors (?)
# - Binary metrics using 0.5, 1 lpf as thresholds

# Timeseries:
#    - Attachment rate
#    - Mortality rate
#    - Abundances (mu + observed points)
#    - Population stage distributions

# Analysis of posteriors
# . Attachment covariate effects ------------------------------------------
farm_env_avg <- readRDS(glue("{dat_dir}/farm_env_avg.rds"))

draw_sample <- sample.int(1500, 50)

if(fishCol=="RW_logit") {
  attach_mx <- readRDS(glue("{dat_dir}/attach_env_mx_RW.rds"))
  pAttach_post_draws <- glue("{out_dir}/posterior{suffix}.rds") |>
    readRDS() |>
    filter(.draw %in% draw_sample) |>
    filter(grepl("attach_beta", name)) |>
    select(-.iteration, -.chain) |>
    mutate(beta=paste0("b", str_sub(name, -2, -2))) |>
    select(-name) |>
    pivot_wider(names_from="beta", values_from="value")

  pAttach_post <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        RW=seq_quantiles(c(attach_mx[,,1])[c(attach_mx[,,1]) > -6], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,2]), 0.25, 0.75, length.out=3),
        Temp=seq_quantiles(c(attach_mx[,,3]), 0.1, 0.9, length.out=3),
        UV=seq_quantiles(c(attach_mx[,,4]), 0, 0.995,  length.out=20)
      ) |>
        mutate(UV_sq=UV^2)
    )) |>
    unnest(attach_df) |>
    mutate(pAttach=plogis(b1*RW + b2*Sal + b3*Temp + b4*UV + b5*UV_sq),
           salinity=Sal*farm_env_avg$salinity[2] + farm_env_avg$salinity[1],
           UV_raw=UV*farm_env_avg$uv[2] + farm_env_avg$uv[1],
           RW=plogis(RW),
           temperature=Temp*farm_env_avg$temperature[2] + farm_env_avg$temperature[1])
  pAttach_post2 <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        RW=seq_quantiles(c(attach_mx[,,1])[c(attach_mx[,,1]) > -6], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,2]), 0.01, 0.99, length.out=20),
        Temp=seq_quantiles(c(attach_mx[,,3]), 0.1, 0.9, length.out=3),
        UV=seq_quantiles(c(attach_mx[,,4]), 0.25, 0.75,  length.out=3)
      ) |>
        mutate(UV_sq=UV^2)
    )) |>
    unnest(attach_df) |>
    mutate(pAttach=plogis(b1*RW + b2*Sal + b3*Temp + b4*UV + b5*UV_sq),
           salinity=Sal*farm_env_avg$salinity[2] + farm_env_avg$salinity[1],
           UV_raw=UV*farm_env_avg$uv[2] + farm_env_avg$uv[1],
           RW=plogis(RW),
           temperature=Temp*farm_env_avg$temperature[2] + farm_env_avg$temperature[1])


  pA <- pAttach_post |>
    mutate(salinity=paste(round(salinity, 1), "psu"),
           RW=paste("RW:", round(RW, 2))) |>
    ggplot(aes(UV_raw, pAttach,
               group=paste(RW, salinity, temperature, .draw), colour=temperature)) +
    geom_line(alpha=0.2, linewidth=0.2) +
    scale_colour_viridis_c(option="plasma", end=0.9) +
    ylim(0, NA) +
    labs(x="UV (cm/s)", y="Daily Pr(Attachment)") +
    facet_grid(RW~salinity)
  pB <- pAttach_post2 |>
    mutate(UV_raw=paste(round(UV_raw, 1), "cm/s"),
           UV_raw=factor(UV_raw, levels=unique(UV_raw)),
           RW=paste("RW:", round(RW, 2))) |>
    ggplot(aes(salinity, pAttach,
               group=paste(RW, UV_raw, temperature, .draw), colour=temperature)) +
    geom_line(alpha=0.2, linewidth=0.2) +
    scale_colour_viridis_c(option="plasma", end=0.9) +
    ylim(0, NA) +
    labs(x="Salinity (psu)", y="Daily Pr(Attachment)") +
    facet_grid(RW~UV_raw)
} else if(fishCol=="BSA") {
  attach_mx <- readRDS(glue("{dat_dir}/attach_env_mx_BSA.rds"))
  pAttach_post_draws <- glue("{out_dir}/posterior{suffix}.rds") |>
    readRDS() |>
    filter(.draw %in% draw_sample) |>
    filter(grepl("attach_beta", name)) |>
    select(-.iteration, -.chain) |>
    mutate(beta=paste0("b", str_sub(name, -2, -2))) |>
    select(-name) |>
    pivot_wider(names_from="beta", values_from="value")

  pAttach_post <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        Int=1,
        BSA=seq_quantiles(c(attach_mx[,,2])[c(attach_mx[,,2]) > 0], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,3]), 0.25, 0.75, length.out=3),
        Temp=seq_quantiles(c(attach_mx[,,4]), 0.1, 0.9, length.out=3),
        UV=seq_quantiles(c(attach_mx[,,5]), 0, 0.995,  length.out=20)
      ) |>
        mutate(UV_sq=UV^2)
    )) |>
    unnest(attach_df) |>
    mutate(pAttach=plogis(b1*Int + b2*BSA + b3*Sal + b4*Temp + b5*UV + b6*UV_sq),
           salinity=Sal*farm_env_avg$salinity[2] + farm_env_avg$salinity[1],
           UV_raw=UV*farm_env_avg$uv[2] + farm_env_avg$uv[1],
           BSA=BSA*farm_env_avg$BSA[2],
           temperature=Temp*farm_env_avg$temperature[2] + farm_env_avg$temperature[1])
  pAttach_post2 <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        Int=1,
        BSA=seq_quantiles(c(attach_mx[,,2])[c(attach_mx[,,2]) > 0], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,3]), 0.01, 0.99, length.out=20),
        Temp=seq_quantiles(c(attach_mx[,,4]), 0.1, 0.9, length.out=3),
        UV=seq_quantiles(c(attach_mx[,,5]), 0.25, 0.75,  length.out=3)
      ) |>
        mutate(UV_sq=UV^2)
    )) |>
    unnest(attach_df) |>
    mutate(pAttach=plogis(b1*Int + b2*BSA + b3*Sal + b4*Temp + b5*UV + b6*UV_sq),
           salinity=Sal*farm_env_avg$salinity[2] + farm_env_avg$salinity[1],
           UV_raw=UV*farm_env_avg$uv[2] + farm_env_avg$uv[1],
           BSA=BSA*farm_env_avg$BSA[2],
           temperature=Temp*farm_env_avg$temperature[2] + farm_env_avg$temperature[1])


  pA <- pAttach_post |>
    mutate(salinity=paste(round(salinity, 1), "psu"),
           BSA=paste("BSA:", round(BSA, 2))) |>
    ggplot(aes(UV_raw, pAttach,
               group=paste(BSA, salinity, temperature, .draw), colour=temperature)) +
    geom_line(alpha=0.2, linewidth=0.2) +
    scale_colour_viridis_c(option="plasma", end=0.9) +
    ylim(0, NA) +
    labs(x="UV (cm/s)", y="Daily Pr(Attachment)") +
    facet_grid(BSA~salinity)
  pB <- pAttach_post2 |>
    mutate(UV_raw=paste(round(UV_raw, 1), "cm/s"),
           UV_raw=factor(UV_raw, levels=unique(UV_raw)),
           BSA=paste("BSA:", round(BSA, 2))) |>
    ggplot(aes(salinity, pAttach,
               group=paste(BSA, UV_raw, temperature, .draw), colour=temperature)) +
    geom_line(alpha=0.2, linewidth=0.2) +
    scale_colour_viridis_c(option="plasma", end=0.9) +
    ylim(0, NA) +
    labs(x="Salinity (psu)", y="Daily Pr(Attachment)") +
    facet_grid(BSA~UV_raw)
}

plot_grid(pA + theme(legend.position="none"),
          pB + theme(legend.position="none"),
          get_legend(pA),
          nrow=1, ncol=3, rel_widths=c(1,1,0.3),
          align="hv", axis="tblr", labels=c("A", "B", "")) |>
  ggsave(glue("{fig_dir}/postAttachReg{suffix}.png"),
       plot=_, width=12, height=4)


# . Salinity covariate effects --------------------------------------------
S_range <- range(readRDS(glue("{dat_dir}/sal_mx.rds")))
S_df <- tibble(sal=seq_range(S_range, length.out=100))
pSurv_post <- glue("{out_dir}/posterior{suffix}.rds") |>
  readRDS() |>
  filter(.draw %in% draw_sample) |>
  filter(grepl("surv_beta", name)) |>
  select(-.iteration, -.chain) |>
  inner_join(param_key, by=join_by(name)) |>
  select(-name) |>
  pivot_wider(names_from="label", values_from="value") |>
  mutate(S_df=list(S_df)) |>
  unnest(S_df) |>
  mutate(pSurv_Ch=plogis(surv_Int_Ch + surv_Sal_Ch*sal),
         pSurv_PA=plogis(surv_Int_PA + surv_Sal_PA*sal),
         pSurv_Ad=plogis(surv_Int_Ad + surv_Sal_Ad*sal),
         salinity=sal + 30) |>
  select(.draw, salinity, starts_with("pSurv")) |>
  pivot_longer(starts_with("pSurv"), names_to="Stage", values_to="pSurv") |>
  mutate(Stage=factor(str_sub(Stage, -2, -1), levels=c("Ch", "PA", "Ad")))
p <- pSurv_post |>
  ggplot(aes(salinity, pSurv, colour=Stage, group=.draw)) +
  geom_line(alpha=0.2, linewidth=0.2) +
  scale_colour_manual(values=RColorBrewer::brewer.pal(n=4, name="Paired")[c(1,2,4)]) +
  ylim(0, 1) +
  labs(x="Salinity (psu)", y="Daily Pr(Survival)") +
  facet_wrap(~Stage, nrow=1) +
  guides(colour=guide_legend(override.aes=list(alpha=1, linewidth=0.5))) +
  theme(legend.position="bottom")
ggsave(glue("{fig_dir}/postSalReg{suffix}.png"), p,
       width=8, height=4)



# . Development time (TEMP EFFECT IS REVERSED????) -----------------------
T_range <- range(readRDS(glue("{dat_dir}/temp_mx.rds")))
T_z_range <- range(readRDS(glue("{dat_dir}/temp_z_mx.rds")))

T_df <- tibble(temp=seq_range(T_range, length.out=100),
               temp_z=seq_range(T_z_range, length.out=100))

stageDur_df <-  glue("{out_dir}/posterior{suffix}.rds") |>
  readRDS() |>
  filter(.draw %in% draw_sample) |>
  filter(grepl("mnDaysStage_beta", name)) |>
  select(-.iteration, -.chain) |>
  inner_join(param_key, by=join_by(name)) |>
  select(-name) |>
  pivot_wider(names_from="label", values_from="value") |>
  mutate(T_df=list(T_df)) |>
  unnest(T_df) |>
  mutate(mnDays_Ch=(mnDaysStage_Int_Ch + mnDaysStage_Temp_Ch*temp_z)*2,
         mnDays_PA=(mnDaysStage_Int_PA + mnDaysStage_Temp_PA*temp_z)*2,
         mnDays_Ad=mnDays_Ch + mnDays_PA) |>
  select(.draw, temp, starts_with("mnDays_")) |>
  pivot_longer(starts_with("mnDays"), values_to="mnDays",
               names_to="Stage", names_prefix="mnDays_") |>
  mutate(pMolt=1/mnDays,
         GDD=temp*mnDays,
         Stage=factor(Stage, levels=c("Ch", "PA", "Ad")))
p <- stageDur_df |>
  ggplot(aes(temp, mnDays, colour=Stage, group=.draw)) +
  geom_line(alpha=0.2, linewidth=0.2) +
  scale_colour_manual(values=RColorBrewer::brewer.pal(n=4, name="Paired")[c(1,2,4)]) +
  labs(x="Temperature (C)", y="Mean stage duration (d)") +
  facet_wrap(~Stage, nrow=1) +
  guides(colour=guide_legend(override.aes=list(alpha=1, linewidth=0.5))) +
  theme(legend.position="bottom")
ggsave(glue("{fig_dir}/postStageDur{suffix}.png"), p,
       width=8, height=4)
p <- stageDur_df |>
  filter(Stage != "Ad") |>
  ggplot(aes(temp, 1/mnDays, colour=Stage, group=.draw)) +
  geom_line(alpha=0.2, linewidth=0.2) +
  scale_colour_manual(values=RColorBrewer::brewer.pal(n=4, name="Paired")[c(1,2,4)]) +
  labs(x="Temperature (C)", y="P(molt) per day") +
  facet_wrap(~Stage, nrow=1) +
  guides(colour=guide_legend(override.aes=list(alpha=1, linewidth=0.5))) +
  theme(legend.position="bottom")
ggsave(glue("{fig_dir}/postPMolt{suffix}.png"), p,
       width=6, height=4)




## Main manuscript figures:
# - Model diagram
# - Map
# - Mu timeseries
# - Attachment covariate effects
# -






# comparison --------------------------------------------------------------

if(FALSE) {
  library(loo)

  mod_f <- dir("out/ipm_fit", "posterior", full.names=T) |>
    grep("summary|PRIORS", x=_, invert=T, value=T)

  ll_names <- c("log_lik\\[", "log_lik_GQ", "log_lik\\[", "log_lik_GQ", "log_lik\\[", "log_lik_GQ", "log_lik\\[")

  log_lik_ls <- map2(mod_f, ll_names,
                     ~readRDS(.x) |>
                       filter(grepl(.y, name)) |>
                       mutate(name=factor(name)) |>
                       arrange(name, .chain, .iteration))

  log_lik_df <- map_dfr(seq_along(mod_f),
                        ~log_lik_ls[[.x]] |>
                          summarise(nll=-sum(value),
                                    .by=.draw) |>
                          mutate(mod=basename(mod_f[.x]) |>
                                   str_remove(".rds"))) |>
    mutate(clean=factor(mod,
                        levels=paste0("posterior",
                                      c("_randIPbg_BSA", "_noHarm_ydayIPbg_BSA", "_ydayIPbg_BSA", "_BSA",
                                        "", "_noHarm", "_randIPbg_RW_logit")),
                        labels=c("ranIPbg_ydayEns_BSA", "ydayIPbg_fixEns_BSA", "ydayIPbg_ydayEns_BSA", "fixIPbg_ydayEns_BSA",
                                 "fixIPbg_ydayEns_RW", "fixIPbg_fixEns_RW", "ranIPbg_ydayEns_RW")))

  ggplot(log_lik_df, aes(nll, y=clean)) +
    ggdist::stat_halfeye()

  log_lik_df |>
    summarise(mn=mean(nll),
              md=median(nll),
              q05=quantile(nll, probs=0.05),
              q95=quantile(nll, probs=0.95),
              .by=clean) |>
    arrange(md)

  ll_ar <- map(log_lik_ls, ~array(.x$value, dim=c(1000,3,1500)))
  ll_Ch <- map(ll_ar, ~.x[,,seq(1, 1500, by=3)])
  ll_PA <- map(ll_ar, ~.x[,,seq(2, 1500, by=3)])
  ll_Ad <- map(ll_ar, ~.x[,,seq(3, 1500, by=3)])

  loo_ar <- map(ll_ar, ~loo(.x, r_eff=relative_eff(.x)))
  loo_Ch <- map(ll_Ch, ~loo(.x, r_eff=relative_eff(.x)))
  loo_PA <- map(ll_PA, ~loo(.x, r_eff=relative_eff(.x)))
  loo_Ad <- map(ll_Ad, ~loo(.x, r_eff=relative_eff(.x)))

  loo_compare(loo_ar)
  loo_compare(loo_Ch)
  loo_compare(loo_PA)
  loo_compare(loo_Ad)



  plot(loo_Ad[[4]], label_points=T)



}

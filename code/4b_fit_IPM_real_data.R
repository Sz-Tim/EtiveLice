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

n_chains <- 3
dat_dir <- "data/aquaculture/mowi_stan/"
out_dir <- "out/ipm_fit/"
stages <- c("Ch1", "Ch2", "PA1", "PA2", "Ad")
stageGrps <- c("Ch", "PA", "Ad")
stage_trans <- c("Ch1-Ch2", "Ch2-PA1", "PA1-PA2", "PA2-Ad")
trt_meth_ii <- read_csv("data/aquaculture/mowi_trt_cleaned.csv") |>
  summarise(.by=c(MethodNum, TypeNum, Method, Type)) |>
  arrange(TypeNum) |>
  mutate(abbr=paste0(str_sub(Method, 1, 1), "_", str_split_i(Type, "_", 2)))
param_key <- tibble(name=c(paste0("attach_beta[", 1:5, "]"),
                           paste0("ensWts_p[", 1:10, "]"),
                           "IP_bg", "IP_bg_m3",
                           paste0("surv_beta[1,", 1:5, "]"),
                           paste0("surv_beta[2,", 1:5, "]"),
                           paste0("surv_int_farm_sd[", 1:5, "]"),
                           paste0("mnDaysStage_beta[1,", 1:4, "]"),
                           paste0("mnDaysStage_beta[2,", 1:4, "]"),
                           paste0("detect_p[", 1:2, "]"),
                           "nb_prec",
                           "IP_scale", "IP_halfSat_m3",
                           paste0("trtEff_type[", 1:8, "]")),
                    label=c(paste0("attach_", c("RW", "Sal", "Temp", "UV", "UVsq")),
                            paste0("ensWt_", 1:10),
                            "IP_bg", "IP_bg_m3",
                            paste0("surv_Int_", stages),
                            paste0("surv_Sal_", stages),
                            paste0("surv_int_farm_sd_", stages),
                            paste0("mnDaysStage_Int_", stages[1:4]),
                            paste0("mnDaysStage_Temp_", stages[1:4]),
                            paste0("p_detect_", stageGrps[-3]),
                            "neg_binom_prec",
                            "IP_scale", "IP_halfSat_m3",
                            paste0("trtEff_", trt_meth_ii$abbr)
                    )) |>
  mutate(label=factor(label, levels=unique(label)))

keep_pars <- c("IP_bg_m3", "ensWts_p", "attach_beta",
               "surv_beta", "surv_int_farm_sd", "mnDaysStage_beta",
               "detect_p", "nb_prec", "trtEff_type",
               "mu"
)



# fit ---------------------------------------------------------------------

iter <- 1000
stan_dat <- make_stan_data(dat_dir, priors_only=prior_only, GQ_start="2025-01-01", source="real")

mod_full <- cmdstan_model("code/stan/integrated_population_model.stan")
fit_full <- mod_full$sample(
  data=stan_dat$dat, init=0, seed=101, refresh=max(iter/100, 1),
  iter_warmup=iter, iter_sampling=iter,
  chains=n_chains, parallel_chains=n_chains
)

# Sampler warnings:
# pMolt[5][423, 2] is 1.01913
# also for 422, 418, 421; all farm 5 (FFMC84)
# But does this only happen early on? So far only 10 iter

out_full_df <- fit_full$draws(
  variables=keep_pars,
  format="df") |>
  pivot_longer(-starts_with("."))
saveRDS(out_full_df, glue("{out_dir}posterior{ifelse(prior_only, '_PRIORS', '')}.rds"))
out_full_sum <- out_full_df |>
  group_by(name) |>
  sevcheck::get_intervals(value, type="qi")
saveRDS(out_full_sum, glue("{out_dir}posterior_summary{ifelse(prior_only, '_PRIORS', '')}.rds"))



# parameter distributions -------------------------------------------------

info <- readRDS(glue("{dat_dir}info.rds"))
p_ensWts <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("ensWts", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(scales="free_y", ncol=if_else(info$nSims > 10, 10, 5)) +
  xlim(0, 1)
p_attach <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("attach_beta", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(scales="free") +
  geom_vline(xintercept=0, linetype=3)
p_surv <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("surv_beta", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=5, scales="free") +
  geom_vline(xintercept=0, linetype=3)
p_surv_sd <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("surv_int_farm_sd", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=5, scales="free") +
  geom_vline(xintercept=0, linetype=3)
p_trt <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("trtEff_type", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=4, scales="free_y") +
  xlim(0, 1)
p_pMoltTemp <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("mnDaysStage_beta", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=4, scales="free_y") +
  geom_vline(xintercept=0, linetype=3)
p_detectp <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("detect_p", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(scales="free_y") +
  xlim(0, 1)
p_else <- list(out_full_df, out_full_sum) |>
  map(~.x |> filter(grepl("IP_bg|nb_prec", name)) |> inner_join(param_key, by=join_by(name))) |>
  post_summary_plot(ncol=5, scales="free")

p <- plot_grid(p_ensWts, p_attach, p_surv, p_surv_sd, p_trt, p_pMoltTemp,
               plot_grid(p_detectp, p_else, nrow=1, axis="tblr", align="hv"),
               nrow=7, align="v", axis="rl", rel_heights=c(2, 1, 2, 1, 2, 2, 1))
ggsave(glue("{out_dir}/fig_pars{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=14)

mu_draws_df <- take_mu_draws(out_full_df, NULL,
                             stan_dat$dat, ndraws=min(1e2, iter), GQ=TRUE) |>
  drop_na(mu) # some prior draws give NAs because of negbinom constraints
p <- mu_draws_df |>
  filter(stage=="Ad") |>
  ggplot(aes(day, mu, group=as.character(.draw))) +
  geom_line() +
  labs(x="Date", y="Mean lice per fish (latent)") +
  {if(any((mu_draws_df |> filter(stage=="Ad"))$mu > 15)) scale_y_continuous(limits=c(0, 15), oob=scales::oob_keep)} +
  scale_x_date(date_labels="%b") +
  facet_grid(farm~., scales="free_y")
ggsave(glue("{out_dir}/fig_mu_draws_GQ{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=15)

mu_draws_df <- take_mu_draws(out_full_df, NULL,
                             stan_dat$dat, ndraws=min(1e2, iter), GQ=F) |>
  drop_na(mu) # some prior draws give NAs because of negbinom constraints
p <- mu_draws_df |>
  filter(stage=="Ad") |>
  ggplot(aes(day, mu, group=as.character(.draw))) +
  geom_line() +
  labs(x="Date", y="Mean lice per fish (latent)") +
  {if(any((mu_draws_df |> filter(stage=="Ad"))$mu > 15)) scale_y_continuous(limits=c(0, 15), oob=scales::oob_keep)} +
  scale_x_date(date_labels="%b") +
  facet_grid(farm~., scales="free_y")
ggsave(glue("{out_dir}/fig_mu_draws{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=15)

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

  y_df <- out_full_df |>
    filter(grepl("y_pred", name)) |>
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
                mutate(y=c(readRDS(glue("{dat_dir}/y.rds"))[,1,,])) |>
                mutate(type="True") |>
                inner_join(sampledDays |> select(-sample), by=join_by(farm, day))) |>
    mutate(day=ymd("2023-01-01") + day - 1,
           farm=paste("Farm", farm),
           y_perFish=y/nFishSampled)
  y_df |>
    saveRDS(glue("{dat_dir}/y_sim_fitted{ifelse(prior_only, '_PRIORS', '')}.rds"))

  p <- y_df |>
    ggplot(aes(day, y_perFish)) +
    geom_linerange(aes(ymin=L10/nFishSampled, ymax=L90/nFishSampled, group=type), linewidth=0.25) +
    geom_point(aes(colour=type, shape=type)) +
    scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
    scale_shape_manual("", values=c("True"=4, "Fitted"=1)) +
    labs(x="Date", y="Mean lice per fish (observed) [50% CI]") +
    scale_x_date(date_labels="%b") +
    facet_grid(stage~farm, scales="free_y")
  ggsave(glue("{dat_dir}/fig_y{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=15, height=7)

  p <- y_df |>
    filter(stage=="Ad") |>
    ggplot(aes(day, y_perFish)) +
    geom_linerange(aes(ymin=L10/nFishSampled, ymax=L90/nFishSampled, group=type), linewidth=0.25) +
    geom_point(aes(colour=type, shape=type)) +
    scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
    scale_shape_manual("", values=c("True"=4, "Fitted"=1)) +
    labs(x="Date", y="Mean adult female lice per fish (observed) [50% CI]") +
    scale_x_date(date_labels="%b") +
    facet_grid(farm~., scales="free_y")
  ggsave(glue("{dat_dir}/fig_y_AF{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=6, height=15)
}

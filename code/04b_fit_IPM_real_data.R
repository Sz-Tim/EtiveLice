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
refit <- T

mod <- c("noHarm_randIPbg", "noHarm_ydayIPbg",
         "Harm_randIPbg", "Harm_ydayIPbg")[4]
fishCol <- c("RW_logit", "BSA")[2]
suffix <- paste0(switch(mod,
                        'noHarm_randIPbg'='_noHarm_randIPbg',
                        'noHarm_ydayIPbg'='_noHarm_ydayIPbg',
                        'Harm_randIPbg'='_randIPbg',
                        'Harm_ydayIPbg'='_ydayIPbg'),
                 ifelse(prior_only, '_PRIORS', ''),
                 "_", fishCol,
                 ifelse(GQ, '', '_noGQ'))

n_chains <- 6
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
                             paste0("attach_beta[", 1:7, "]")
                           } else {
                             paste0("attach_beta[", 1:6, "]")
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
                           paste0("trtEff_type[", outer(1:8, 1:3, paste, sep=","), "]")
                           # paste0("trtEff_type[", 1:8, "]")
                          ),
                    label=c({if(fishCol=="BSA") {
                              paste0("attach_", c("Int", "Fish", "Sal", "Temp", "UV", "Tempsq", "UVsq"))
                            } else {
                              paste0("attach_", c("Fish", "Sal", "Temp", "UV", "Tempsq", "UVsq"))
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
                            paste0("trtEff_type[", outer(trt_meth_ii$abbr, stageGrps, paste, sep=","), "]")
                            # paste0("trtEff_", trt_meth_ii$abbr)
                    )) |>
  mutate(label=factor(label, levels=unique(label)))

keep_pars <- c(ifelse(grepl("ydayIP", mod), "log_IP_bg_m3_coef", "IP_bg_m3"),
               ifelse(grepl("noHarm", mod), "ensWts_p", "ensWts_harm"),
               "attach_beta","surv_beta", "surv_int_farm_sd", "mnDaysStage_beta",
               "detect_p", "nb_prec", "trtEff_type",
               "mu", "mu_GQ", "y_pred", "y_bar_GQ",
               "N_attach", "N_attach_GQ",
               "log_lik", "log_lik_GQ"
)
if(!GQ) keep_pars <- grep("_GQ", keep_pars, invert=T, value=T)



# fit ---------------------------------------------------------------------

post_size <- 3000
iter <- post_size/n_chains
if(grepl("ydayIP", mod)) {
  prior_ls <- list(prior_IP_bg_m3=c(log(0.05), 0.5))
} else {
  prior_ls <- NULL
}
stan_dat <- make_stan_data(dat_dir, priors_only=prior_only, source="real",
                           GQ_start=switch(GQ+1, NULL, "2025-01-01"),
                           fishCol=fishCol,
                           prior_ls=prior_ls)

if(refit) {
  mod_full <- cmdstan_model(glue("code/stan/tuning_integrated_population_model{str_remove(suffix, '_PRIORS') |> str_remove(paste0('_', fishCol)) |> str_remove('_noGQ')}.stan"))
  fit_full <- mod_full$sample(
    data=stan_dat$dat, init=0, seed=101, refresh=max(iter/100, 1),
    adapt_delta=0.95,
    # iter_warmup=5, iter_sampling=5, chains=2
    iter_warmup=2000, iter_sampling=iter,
    chains=n_chains, parallel_chains=n_chains
  )

  suffix <- paste0(suffix, "_tq_trtStg")

  out_full_df <- fit_full$draws(
    variables=keep_pars,
    format="df") |>
    pivot_longer(-starts_with("."))
  # saveRDS(out_full_df, glue("{out_dir}posterior{suffix}.rds"))
  out_full_sum <- out_full_df |>
    group_by(name) |>
    sevcheck::get_intervals(value, type="qi")
  saveRDS(out_full_sum, glue("{out_dir}posterior_summary{suffix}.rds"))

  walk(keep_pars,
       ~out_full_df |>
         filter(grepl(.x, name)) |>
         saveRDS(glue("{out_dir}{.x}_post{suffix}.rds")))
} else {
  suffix <- paste0(suffix, "_tq")
  out_full_df <- readRDS(glue("{out_dir}posterior{suffix}.rds"))
  out_full_sum <- readRDS(glue("{out_dir}posterior_summary{suffix}.rds"))
}



# log likelihood ----------------------------------------------------------

if(GQ) {
  p_loglik <- out_full_sum |>
    filter(grepl("log_lik_GQ", name)) |>
    mutate(stage=str_sub(name, 12, 12)) |>
    ggplot(aes(med)) +
    geom_histogram() +
    facet_grid(.~stage, scales="free")
  ggsave(glue("{fig_dir}/loglik{suffix}.png"),
         p_loglik, width=9, height=3.5)
}

p_loglik <- out_full_sum |>
  filter(grepl("log_lik\\[", name)) |>
  mutate(stage=str_sub(name, 9, 9)) |>
  ggplot(aes(med)) +
  geom_histogram() +
  facet_grid(.~stage, scales="free")
ggsave(glue("{fig_dir}/loglik_fitted{suffix}.png"),
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
  post_summary_plot(scales="free", nrow=1, ncol=7) +
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
  post_summary_plot(ncol=8, scales="free_y") +
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

ggsave(glue("{fig_dir}/pars{suffix}.png"), p, width=10, height=14)



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
  ggsave(glue("{fig_dir}/mu_draws_GQ{suffix}.png"), p, width=10, height=15)

  mu_pred_df <- mu_draws_df |>
    summarise(mn=mean(mu),
              lo=quantile(mu, probs=0.25),
              hi=quantile(mu, probs=0.75),
              .by=c(farm, stage, day)) |>
    inner_join(obs_df |>
                 mutate(across(all_of(c("Ch", "PA", "AF")),  ~.x/nFishSampled)) |>
                 select(farm, date, Ch, PA, AF) |>
                 rename(day=date, Ad=AF) |>
                 pivot_longer(3:5, names_to="stage", values_to="obs"),
               by=join_by(stage, farm, day))
  p <- ggplot(mu_pred_df, aes(mn, xmin=lo, xmax=hi, y=obs)) +
    geom_point(shape=1) +
    geom_linerange() +
    geom_abline() +
    facet_grid(farm~stage, scales="free")
  ggsave(glue("figs/ipm_fit/mu_scatter_1{suffix}.png"), p)
  p <- ggplot(mu_pred_df, aes(mn, xmin=lo, xmax=hi, y=obs)) +
    geom_point(shape=1) +
    geom_linerange() +
    geom_abline() +
    facet_wrap(~stage, scales="free")
  ggsave(glue("figs/ipm_fit/mu_scatter_2{suffix}.png"), p)
  p <- ggplot(mu_pred_df |> filter(stage=="Ad"), aes(mn, xmin=lo, xmax=hi, y=obs)) +
    geom_point(shape=1) +
    geom_linerange() +
    geom_abline() +
    facet_wrap(~farm, scales="free")
  ggsave(glue("figs/ipm_fit/mu_scatter_3{suffix}.png"), p)
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
ggsave(glue("{fig_dir}/mu_draws{suffix}.png"), p, width=10, height=15)

if("y_pred" %in% keep_pars) {
  sampledDays <- stan_dat$dat$sample_i |>
    as_tibble() |>
    rename(farm=sepaSite) |>
    mutate(farm=as.character(farm),
           sample=row_number()) |>
    inner_join(stan_dat$dat$nFishSampled_mx |>
                 t() |>
                 as_tibble() |>
                 mutate(day=row_number()) |>
                 pivot_longer(starts_with("V"), names_to="farm", values_to="nFishSampled") |>
                 mutate(farm=str_sub(farm, 2, -1)),
               by=join_by(farm, day))

  y_df <- out_full_sum |>
    filter(grepl("y_pred\\[[1-3]", name)) |>
    separate_wider_delim(name, delim=",", names=c("stage", "sample")) |>
    mutate(sample=as.numeric(str_sub(sample, 1, -2)),
           stage=factor(stage, levels=paste0("y_pred[", 1:3), labels=c("Ch", "PA", "Ad"))) |>
    left_join(sampledDays, by=join_by(sample)) |>
    select(-sample) |>
    inner_join(expand_grid(farm=as.character(1:info$nFarms),
                          day=1:info$nDays,
                          stage=factor(c("Ch", "PA", "Ad"), levels=c("Ch", "PA", "Ad"))) |>
                mutate(y=c(readRDS(glue("{dat_dir}/y_obs.rds")))) |>
                inner_join(sampledDays |> select(farm, day), by=join_by(farm, day)),
               by=join_by(farm, day, stage)) |>
    mutate(day=ymd("2025-01-01") + day - 1,
           farm=paste("Farm", farm),
           across(c(any_of(c("mn", "med", "y")), starts_with("L")), ~.x/nFishSampled)) |>
    filter(nFishSampled > 0)
  y_df |>
    saveRDS(glue("{dat_dir}/y_fitted{suffix}.rds"))

  y_post <- out_full_df |>
    filter(.draw %in% sample.int(3000, 200) & grepl("y_pred\\[[1-3]", name)) |>
    separate_wider_delim(name, delim=",", names=c("stage", "sample")) |>
    mutate(sample=as.numeric(str_sub(sample, 1, -2)),
           stage=factor(stage, levels=paste0("y_pred[", 1:3), labels=c("Ch", "PA", "Ad"))) |>
    left_join(sampledDays, by=join_by(sample)) |>
    select(-sample) |>
    mutate(y_rep=value/nFishSampled,
           farm=paste("Farm", farm))
  pA <- y_post |>
    ggplot(aes(y_rep)) +
    geom_density(aes(group=.draw), alpha=0.1, colour="steelblue3") +
    geom_density(data=y_df, aes(y))
  pB <- y_post |>
    ggplot(aes(y_rep)) +
    geom_density(aes(group=.draw), alpha=0.1, colour="steelblue3") +
    geom_density(data=y_df, aes(y)) +
    facet_wrap(~stage, scales="free")
  pC <- y_post |>
    ggplot(aes(y_rep)) +
    geom_density(aes(group=.draw), alpha=0.1, colour="steelblue3") +
    geom_density(data=y_df, aes(y)) +
    facet_wrap(~stage*farm, scales="free", ncol=10)
  ggsave(glue("{fig_dir}/ppcheck_1{suffix}.png"), pA, width=5, height=5)
  ggsave(glue("{fig_dir}/ppcheck_2{suffix}.png"), pB, width=10, height=4)
  ggsave(glue("{fig_dir}/ppcheck_3{suffix}.png"), pC, width=15, height=8)

  pA <- y_post |>
    ggplot(aes(y_rep)) +
    geom_density(aes(group=.draw), alpha=0.1, colour="steelblue3") +
    geom_density(data=y_df, aes(y)) +
    scale_x_continuous(transform="log1p")
  pB <- y_post |>
    ggplot(aes(y_rep)) +
    geom_density(aes(group=.draw), alpha=0.1, colour="steelblue3") +
    geom_density(data=y_df, aes(y)) +
    facet_wrap(~stage, scales="free") +
    scale_x_continuous(transform="log1p")
  pC <- y_post |>
    ggplot(aes(y_rep)) +
    geom_density(aes(group=.draw), alpha=0.1, colour="steelblue3") +
    geom_density(data=y_df, aes(y)) +
    facet_wrap(~stage*farm, scales="free", ncol=10) +
    scale_x_continuous(transform="log1p")
  ggsave(glue("{fig_dir}/ppcheck_1log{suffix}.png"), pA, width=5, height=5)
  ggsave(glue("{fig_dir}/ppcheck_2log{suffix}.png"), pB, width=10, height=4)
  ggsave(glue("{fig_dir}/ppcheck_3log{suffix}.png"), pC, width=15, height=8)


  p <- y_df |>
    ggplot(aes(mn, y)) +
    # geom_abline() +
    geom_point(alpha=0.2, shape=1) +
    facet_wrap(~stage, scales="free")
  ggsave(glue("{fig_dir}/y_scatter{suffix}.png"), p, width=10, height=5)

  p <- y_df |>
    ggplot(aes(day, mn)) +
    geom_linerange(aes(ymin=L10, ymax=L90), linewidth=0.25) +
    geom_point(shape=1) +
    geom_point(aes(y=y), colour="blue", shape=4) +
    labs(x="Date", y="Mean lice per fish (observed) [80% CI]") +
    scale_x_date(date_labels="%b") +
    facet_grid(stage~farm, scales="free_y")
  ggsave(glue("{fig_dir}/y{suffix}.png"), p, width=15, height=7)

  p <- y_df |>
    filter(stage=="Ad") |>
    ggplot(aes(day, mn)) +
    geom_linerange(aes(ymin=L10, ymax=L90), linewidth=0.25) +
    geom_point(shape=1) +
    geom_point(aes(y=y), colour="blue", shape=4) +
    labs(x="Date", y="Mean adult female lice per fish (observed) [80% CI]") +
    scale_x_date(date_labels="%b") +
    facet_grid(farm~., scales="free_y")
  ggsave(glue("{fig_dir}/y_AF{suffix}.png"), p, width=6, height=15)
}





# rhat --------------------------------------------------------------------

library(furrr)
plan(multicore, workers=10)
out_split <- out_full_df |>
  filter(name %in% param_key$name) |>
  group_split(name)
rhats <- future_map_dbl(out_split, ~rstan::Rhat(matrix(.x$value, nrow=iter, ncol=n_chains))) |>
  set_names(map_chr(out_split, ~.x$name[1]))
plan(sequential)
p <- rhats |>
  as_tibble() |>
  mutate(name=names(rhats),
         type=str_split_i(name, "\\[", 1)) |>
  ggplot(aes(value, name)) +
  geom_point() +
  geom_vline(xintercept=1.05) +
  facet_wrap(~type, scales="free_y", space="free_y", ncol=1)
ggsave(glue("{fig_dir}/rhat{suffix}.png"), p, width=5, height=10)



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

draw_sample <- sample.int(1500, 300)


# . IPbg ------------------------------------------------------------------

library(sf)
mesh_dir <- ifelse(sevcheck::get_os()=="windows", "../03_packages/WeStCOMS/data", "~/hydro/meshes")
farm_bbox <- list(xmin=125000, xmax=225500, ymin=690000, ymax=785000)
farm_bbox <- list(xmin=145000, xmax=225500, ymin=720000, ymax=785000)
linnhe_fp <- st_read(glue("{mesh_dir}/WeStCOMS3_mesh_footprint.gpkg")) |>
  st_crop(unlist(farm_bbox))
site_names <- readRDS(glue("{dat_dir}/site_names.rds"))
farm_sites <- read_csv("data/farm_sites_widerLinnhe_2022-2025.csv")

if(grepl("randIPbg", suffix)) {
  farm_IPbg <- out_full_sum |>
    filter(grepl("IP_bg", name)) |>
    mutate(farm=str_sub(name, 10, -2) |> as.numeric(),
           sepaSite=names(site_names)[farm]) |>
    left_join(farm_sites)

  p <- ggplot(linnhe_fp) +
    geom_sf() +
    geom_point(data=farm_IPbg, aes(easting, northing, fill=med), size=5, shape=21) +
    scale_fill_viridis_c(expression("Posterior median lice "%.%" m"^"-3"%.%" h"^"-1"),
                         option="plasma", begin=0.05, end=0.95, limits=c(0, NA),
                         breaks=seq(0, 0.2, by=0.05),
                         labels=c("0", "0.05", "0.10", "0.15", "0.20")) +
    scale_x_continuous(expand=c(0,0), oob=scales::oob_keep, n.breaks=3) +
    scale_y_continuous(expand=c(0,0), oob=scales::oob_keep, n.breaks=4) +
    ggtitle("Background IP") +
    theme(axis.title=element_blank(),
          legend.position="inside",
          # legend.position.inside=c(0.277, 0.915),
          legend.position.inside=c(0.293, 0.915),
          legend.direction="horizontal",
          legend.title.position="top",
          legend.key.height=unit(0.2, "cm"),
          legend.key.width=unit(1, "cm"),
          legend.background=element_rect(colour="grey30", linewidth=0.2, fill="white"))
  ggsave(glue("{fig_dir}/IPbg_md_map{suffix}.png"),
         plot=p, width=5, height=5)
} else {
  yday_df <- as_tibble(make_ydayh_mx()[(1:(366*24))%%(24)==1,]) |>
    set_names(c("yday_Int", "yday_cos", "yday_sin")) |>
    mutate(yday=1:366,
           date=ymd("2024-01-01") + yday-1) |>
    filter(day(date)==1)
  post_df <- readRDS(glue("{out_dir}log_IP_bg_m3_coef_post{suffix}.rds")) |>
    mutate(farm=as.numeric(str_sub(str_split_i(name, ",", 2), 1, -2)),
           param=str_sub(str_split_i(name, ",", 1), -1, -1),
           param=c("Int", "cos", "sin")[as.numeric(param)]) |>
    select(farm, param, .draw, value) |>
    pivot_wider(names_from="param", values_from="value") |>
    mutate(yday_df=list(yday_df)) |>
    unnest(yday_df) |>
    mutate(log_IP_bg_m3=Int*yday_Int + cos*yday_cos + sin*yday_sin) |>
    group_by(yday, .draw) |>
    mutate(IP_bg_m3=exp(log_IP_bg_m3)) |>
    ungroup() |>
    summarise(med=median(IP_bg_m3),
              .by=c(farm, yday)) |>
    mutate(date_std=ymd("2024-01-01") + yday - 1,
           sepaSite=names(site_names)[farm]) |>
    left_join(farm_sites)
  for(m in 1:12) {
    p <- ggplot(linnhe_fp) +
      geom_sf() +
      geom_point(data=post_df |> filter(month(date_std)==m),
                 aes(easting, northing, fill=med), size=5, shape=21) +
      scale_fill_viridis_c(expression("Posterior median lice "%.%" m"^"-3"%.%" h"^"-1"),
                           option="plasma", begin=0.05, end=0.95, limits=c(0, 0.2),
                           breaks=seq(0, 0.2, by=0.05),
                           labels=c("0", "0.05", "0.10", "0.15", "0.20")) +
      scale_x_continuous(expand=c(0,0), oob=scales::oob_keep, n.breaks=3) +
      scale_y_continuous(expand=c(0,0), oob=scales::oob_keep, n.breaks=4) +
      ggtitle(paste0("Background IP: 01-", month.abb[m])) +
      theme(axis.title=element_blank(),
            legend.position="inside",
            # legend.position.inside=c(0.277, 0.915),
            legend.position.inside=c(0.293, 0.915),
            legend.direction="horizontal",
            legend.title.position="top",
            legend.key.height=unit(0.2, "cm"),
            legend.key.width=unit(1, "cm"),
            legend.background=element_rect(colour="grey30", linewidth=0.2, fill="white"))
    ggsave(glue("{fig_dir}/IPbg_md_map-{str_pad(m, 2, 'left', '0')}{suffix}.png"),
           plot=p, width=5, height=5)
  }
}



# . IP Ensemble weights ---------------------------------------------------

if(!grepl("noHarm", suffix)) {
  yday_df <- as_tibble(make_ydayh_mx()[(1:(366*24))%%24==1,]) |>
    set_names(c("yday_Int", "yday_cos", "yday_sin")) |>
    mutate(yday=1:366)
  ensP_post <- readRDS(glue("{out_dir}ensWts_harm_post{suffix}.rds")) |>
    inner_join(param_key, by=join_by(name)) |>
    mutate(sim=paste("Sim", str_split_i(label, "_", 3)),
           param=str_split_i(label, "_", 2)) |>
    select(sim, param, .draw, .chain, .iteration, value) |>
    pivot_wider(names_from="param", values_from="value") |>
    mutate(yday_df=list(yday_df)) |>
    unnest(yday_df) |>
    mutate(logit_p=Int*yday_Int + cos*yday_cos + sin*yday_sin) |>
    group_by(yday, .draw, .chain, .iteration) |>
    mutate(p=make_compositional(logit_p, method="softmax")) |>
    ungroup()
  ensP_sum <- ensP_post |>
    summarise(md=median(p),
              mn=mean(p),
              q05=quantile(p, probs=0.05),
              q10=quantile(p, probs=0.1),
              q25=quantile(p, probs=0.25),
              q75=quantile(p, probs=0.75),
              q90=quantile(p, probs=0.9),
              q95=quantile(p, probs=0.95),
              .by=c(sim, yday)) |>
    mutate(date_std=ymd("2020-01-01") + yday - 1,
           sim=factor(sim, levels=paste("Sim", 1:50)))

  # ggplot(ensP_sum, aes(date_std, md, group=sim)) +
  #   geom_ribbon(aes(ymin=q10, ymax=q90), colour=NA, alpha=0.1) +
  #   geom_ribbon(aes(ymin=q25, ymax=q75), colour=NA, alpha=0.1) +
  #   geom_line() +
  #   scale_y_continuous("Ensemble weight p", limits=c(0, 1),
  #                      breaks=round(seq(0, 1, by=1/info$nSims), 2)) +
  #   scale_x_datetime("Day of year", date_breaks="3 months", date_labels="%b")

  p <- ggplot(ensP_sum, aes(date_std, md, fill=sim)) +
    geom_area(position="fill", colour="grey30", linewidth=0.1) +
    scale_fill_brewer("Simulation", palette="Paired") +
    scale_y_continuous("Posterior median weight", limits=c(0, 1),
                       expand=c(0,0), oob=scales::oob_keep,
                       breaks=round(seq(0, 1, by=1/info$nSims), 2)) +
    scale_x_datetime("Day of year", date_breaks="1 month", date_labels="%b",
                     expand=c(0,0), oob=scales::oob_keep) +
    ggtitle("Dispersal simulation IP") +
    theme(axis.title.x=element_blank(),
          legend.position="none")
  ggsave(glue("{fig_dir}/IPens_md{suffix}.png"),
         plot=p, width=4.5, height=4.5)

  # ggplot(ensP_sum, aes(date_std, md, fill=sim)) +
  #   geom_area(position="fill", colour="grey30", linewidth=0.1) +
  #   scale_fill_brewer("Simulation", palette="Paired") +
  #   scale_y_continuous("Weight in ensemble", limits=c(0, 1),
  #                      breaks=round(seq(0, 1, by=1/info$nSims), 2)) +
  #   scale_x_datetime("Day of year", date_breaks="1 month", date_labels="%b") +
  #   coord_polar()
  #
  # ensP_post |>
  #   filter(.draw %in% draw_sample[1:20]) |>
  #   mutate(date_std=ymd("2020-01-01") + yday - 1,
  #          sim=factor(sim, levels=paste("Sim", 1:50))) |>
  #   ggplot(aes(date_std, p, fill=sim)) +
  #   geom_area(position="fill", colour="grey30", linewidth=0.1) +
  #   scale_fill_brewer("Simulation", palette="Paired") +
  #   scale_y_continuous("Weight in ensemble", limits=c(0, 1),
  #                      breaks=round(seq(0, 1, by=1/info$nSims), 2)) +
  #   scale_x_datetime("Day of year", date_breaks="3 months", date_labels="%b") +
  #   facet_wrap(~.draw, nrow=4)
  #
  #
  # ensP_post |>
  #   filter(.draw %in% draw_sample) |>
  #   mutate(date_std=ymd("2020-01-01") + yday - 1,
  #          sim=factor(sim, levels=paste("Sim", 1:50))) |>
  #   ggplot(aes(date_std, p, group=.draw)) +
  #   geom_line(alpha=0.2) +
  #   scale_y_continuous("Weight in ensemble", limits=c(0, 1),
  #                      breaks=round(seq(0, 1, by=1/info$nSims), 2)) +
  #   scale_x_datetime("Day of year", date_breaks="3 months", date_labels="%b") +
  #   facet_wrap(~sim, nrow=2)
}



# . Attachment covariate effects ------------------------------------------
farm_env_avg <- readRDS(glue("{dat_dir}/farm_env_avg.rds"))

if(fishCol=="RW_logit") {
  attach_mx <- readRDS(glue("{dat_dir}/attach_env_mx_RW.rds"))
  pAttach_post_draws <- glue("{out_dir}attach_beta_post{suffix}.rds") |>
    readRDS() |>
    filter(.draw %in% draw_sample) |>
    select(-.iteration, -.chain) |>
    mutate(beta=paste0("b", str_sub(name, -2, -2))) |>
    select(-name) |>
    pivot_wider(names_from="beta", values_from="value")

  pAttach_post <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        RW=seq_quantiles(c(attach_mx[,,1])[c(attach_mx[,,1]) > -6], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,2]), 0.01, 0.99, length.out=3),
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
        UV=seq_quantiles(c(attach_mx[,,4]), 0.1, 0.9,  length.out=3)
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
    labs(x="UV (cm/s)", y="Pr(Attach to host)") +
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
    labs(x="Salinity (psu)", y="Pr(Attach to host)") +
    facet_grid(RW~UV_raw)
} else if(fishCol=="BSA") {
  attach_mx <- readRDS(glue("{dat_dir}/attach_env_mx_BSA.rds"))
  pAttach_post_draws <- glue("{out_dir}posterior_attach_beta{suffix}.rds") |>
    readRDS() |>
    filter(.draw %in% draw_sample) |>
    select(-.iteration, -.chain) |>
    mutate(beta=paste0("b", str_sub(name, -2, -2))) |>
    select(-name) |>
    pivot_wider(names_from="beta", values_from="value")

  pAttach_post <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        Int=1,
        BSA=seq_quantiles(c(attach_mx[,,2])[c(attach_mx[,,2]) > 0], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,3]), 0.01, 0.99, length.out=3) |> round(),
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
        UV=seq_quantiles(c(attach_mx[,,5]), 0.1, 0.95,  length.out=3) |> round()
      ) |>
        mutate(UV_sq=UV^2)
    )) |>
    unnest(attach_df) |>
    mutate(pAttach=plogis(b1*Int + b2*BSA + b3*Sal + b4*Temp + b5*UV + b6*UV_sq),
           salinity=Sal*farm_env_avg$salinity[2] + farm_env_avg$salinity[1],
           UV_raw=UV*farm_env_avg$uv[2] + farm_env_avg$uv[1],
           BSA=BSA*farm_env_avg$BSA[2],
           temperature=Temp*farm_env_avg$temperature[2] + farm_env_avg$temperature[1])
  pAttach_post3 <- pAttach_post_draws |>
    mutate(attach_df=list(
      expand_grid(
        Int=1,
        BSA=seq_quantiles(c(attach_mx[,,2])[c(attach_mx[,,2]) > 0], 0.25, 0.75, length.out=2),
        Sal=seq_quantiles(c(attach_mx[,,3]), 0.01, 0.99, length.out=3),
        Temp=seq_quantiles(c(attach_mx[,,4]), 0.01, 0.99, length.out=20),
        UV=seq_quantiles(c(attach_mx[,,5]), 0.1, 0.95,  length.out=3) |> round()
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
    scale_colour_viridis_c("Temperature (C)", option="plasma", end=0.9) +
    ylim(0, NA) +
    labs(x="UV (cm/s)", y="Pr(Attach to host)") +
    facet_grid(BSA~salinity)
  pB <- pAttach_post2 |>
    mutate(UV_raw=paste(round(UV_raw, 1), "cm/s"),
           UV_raw=factor(UV_raw, levels=unique(UV_raw)),
           BSA=paste("BSA:", round(BSA, 2))) |>
    ggplot(aes(salinity, pAttach,
               group=paste(BSA, UV_raw, temperature, .draw), colour=temperature)) +
    geom_line(alpha=0.2, linewidth=0.2) +
    scale_colour_viridis_c("Temperature (C)", option="plasma", end=0.9) +
    ylim(0, NA) +
    labs(x="Salinity (psu)", y="Pr(Attach to host)") +
    facet_grid(BSA~UV_raw)
  pC <- pAttach_post3 |>
    mutate(UV_raw=paste(round(UV_raw, 1), "cm/s"),
           UV_raw=factor(UV_raw, levels=unique(UV_raw)),
           BSA=paste("BSA:", round(BSA, 2))) |>
    ggplot(aes(temperature, pAttach,
               group=paste(BSA, UV_raw, salinity, .draw), colour=salinity)) +
    geom_line(alpha=0.2, linewidth=0.2) +
    scale_colour_distiller("Salinity (psu)", palette="Blues") +
    ylim(0, NA) +
    labs(x="Temperature (C)", y="Pr(Attach to host)") +
    facet_grid(BSA~UV_raw)
}

# plot_grid(pA + theme(legend.position="none"),
#           pB + theme(legend.position="none"),
#           get_legend(pA),
#           nrow=1, ncol=3, rel_widths=c(1,1,0.3),
#           align="hv", axis="tblr", labels=c("A", "B", "")) |>
#   ggsave(glue("{fig_dir}/postAttachReg{suffix}.png"),
#        plot=_, width=12, height=4)
ggsave(glue("{fig_dir}/postAttachRegA{suffix}.png"), pA, width=7, height=5)
ggsave(glue("{fig_dir}/postAttachRegB{suffix}.png"), pB, width=7, height=5)
ggsave(glue("{fig_dir}/postAttachRegC{suffix}.png"), pC, width=7, height=5)


# . Salinity covariate effects --------------------------------------------
S_range <- range(readRDS(glue("{dat_dir}/sal_mx.rds")))
S_df <- tibble(sal=seq_range(S_range, length.out=100))
pSurv_post <- glue("{out_dir}surv_beta_post{suffix}.rds") |>
  readRDS() |>
  filter(.draw %in% draw_sample) |>
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

stageDur_df <-  glue("{out_dir}mnDaysStage_beta_post{suffix}.rds") |>
  readRDS() |>
  filter(.draw %in% draw_sample) |>
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
    grep("noGQ", x=_, value=T) |>
    grep("summary|PRIORS", x=_, invert=T, value=T)

  log_lik_ls <- map(mod_f,
                     ~readRDS(.x) |>
                       filter(grepl("log_lik", name)) |>
                       mutate(name=factor(name)) |>
                       arrange(name, .chain, .iteration))
  saveRDS(log_lik_ls, "out/ipm_fit/log_lik_ls.rds")

  log_lik_df <- map_dfr(seq_along(mod_f),
                        ~log_lik_ls[[.x]] |>
                          summarise(nll=-sum(value),
                                    .by=.draw) |>
                          mutate(mod=basename(mod_f[.x]) |>
                                   str_remove(".rds"))) |>
    mutate(clean=factor(mod,
                        levels=paste0("posterior_",
                                      c("noHarm_ydayIPbg_BSA_noGQ",
                                        "noHarm_randIPbg_BSA_noGQ_tq", "noHarm_ydayIPbg_BSA_noGQ_tq",
                                        "randIPbg_BSA_noGQ_tq", "ydayIPbg_BSA_noGQ_tq")),
                        labels=c("noHarm_ydayIPbg_linTemp", "noHarm_fixIPbg", "noHarm_ydayIPbg",
                                 "Harm_fixIPbg", "Harm_ydayIPbg")))
    # mutate(clean=factor(mod,
    #                     levels=paste0("posterior",
    #                                   c("_randIPbg_BSA", "_noHarm_ydayIPbg_BSA", "_ydayIPbg_BSA", "_BSA",
    #                                     "", "_noHarm", "_randIPbg_RW_logit")),
    #                     labels=c("ranIPbg_ydayEns_BSA", "ydayIPbg_fixEns_BSA", "ydayIPbg_ydayEns_BSA", "fixIPbg_ydayEns_BSA",
    #                              "fixIPbg_ydayEns_RW", "fixIPbg_fixEns_RW", "ranIPbg_ydayEns_RW")))
  saveRDS(log_lik_df, "out/ipm_fit/log_lik_df.rds")

  p <- ggplot(log_lik_df, aes(nll, y=clean)) +
    ggdist::stat_halfeye()
  ggsave("figs/ipm_fit/log_lik_distributions.png", p)

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

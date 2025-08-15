# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Fit integrated ensemble


# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(rstan)
library(cowplot)
library(ggdist)
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1) # reduce_sum() or map_rect() in .stan
source("code/fn/helpers.R")


stages <- c("Ch", "PA", "Ad")
stage_trans <- c("Ch-PA", "PA-Ad")
param_key <- tibble(name=c(paste0("attach_beta[", 1:5, "]"),
                           paste0("ensWts_p[", 1:6, "]"),
                           "IP_bg", "IP_bg_m3",
                           paste0("surv_beta[1,", 1:3, "]"),
                           paste0("surv_beta[2,", 1:3, "]"),
                           paste0("mnDaysStage_beta[1,", 1:2, "]"),
                           paste0("mnDaysStage_beta[2,", 1:2, "]"),
                           paste0("thresh_GDD[", 1:2, ",1]"),
                           paste0("thresh_GDD[", 1:2, ",2]"),
                           "lifespan",
                           paste0("detect_p[", 1:2, "]"),
                           "nb_prec",
                           "IP_scale",
                           "treatEfficacy"),
                    label=c(paste0("attach_", c("RW", "Sal", "UV", "UVsq", "Temp")),
                            paste0("ensWt_", 1:6),
                            "IP_bg", "IP_bg_m3",
                            paste0("surv_Int_", stages),
                            paste0("surv_Sal_", stages),
                            paste0("mnDaysStage_Int_", stages[1:2]),
                            paste0("mnDaysStage_Temp_", stages[1:2]),
                            paste0("moltF_GDD_", stage_trans),
                            paste0("moltM_GDD_", stage_trans),
                            "lifespan_GDD",
                            paste0("p_detect_", stages[-3]),
                            "neg_binom_prec",
                            "IP_scale",
                            "treatEfficacy"
                    )) |>
  mutate(label=factor(label, levels=unique(label)))

for(sim in 1:10) {

  dat_dir <- glue("data/sim/sim_{str_pad(sim, 2, 'left', '0')}/")
  stan_dat <- make_stan_data(dat_dir, priors_only=F)

  # IEM: full model ---------------------------------------------------------

  out_full <- stan(file="code/stan/joint_population_model.stan",
                   data=stan_dat$dat, chains=1, cores=1, iter=10, refresh=1,
                   control=list(max_treedepth=20),
                   init=0, seed=101,
                   pars=c("IP_bg", "IP_bg_m3", "IP_scale", "ensWts_p", "attach_beta",
                          "surv_beta", "mnDaysStage_beta",
                          "detect_p", "nb_prec", "treatEfficacy"))#,
                          # "mu", "y_pred"))
  saveRDS(out_full, glue("{dat_dir}/out_full_stanfit.rds"))

  out_full_df <- as.data.frame(out_full) |>
    select(-lp__) |>
    mutate(iter=row_number()) |>
    pivot_longer(-iter)
  out_full_sum <- out_full_df |>
    group_by(name) |>
    sevcheck::get_intervals(value, type="qi")
  dat_full_df <- tibble(
    name=c(paste0("attach_beta[", 1:stan_dat$dat$nAttachCov, "]"),
           "IP_bg", "IP_bg_m3",
           paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
           paste0("surv_beta[1,", 1:stan_dat$dat$nStages, "]"),
           paste0("surv_beta[2,", 1:stan_dat$dat$nStages, "]"),
           paste0("thresh_GDD[", 1:(stan_dat$dat$nStages-1), ",1]"),
           paste0("thresh_GDD[", 1:(stan_dat$dat$nStages-1), ",2]"),
           "lifespan",
           paste0("detect_p[", 1:(stan_dat$dat$nStages), "]"), "nb_prec"),
    value=c(stan_dat$params$attach_beta,
            stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
            stan_dat$params$ensWts_p,
            t(stan_dat$params$surv_beta),
            stan_dat$params$thresh_GDD,
            stan_dat$params$lifespan,
            stan_dat$params$detect_p,
            stan_dat$params$nb_prec)
  )

  attach_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                   ~.x |> filter(grepl("attach_beta", name)) |>
                     inner_join(param_key, by=join_by(name)))
  p_attach <- post_summary_plot(attach_ls, scales="free") +
    geom_vline(xintercept=0, linetype=3)

  ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                   ~.x |> filter(grepl("ensWts", name)) |>
                     inner_join(param_key, by=join_by(name)))
  p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
    xlim(0, 1)

  surv_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                 ~.x |> filter(grepl("surv_beta", name)) |>
                   inner_join(param_key, by=join_by(name)))
  p_surv <- post_summary_plot(surv_ls, ncol=3, scales="free") +
    geom_vline(xintercept=0, linetype=3)

  pMolt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                  ~.x |> filter(grepl("thresh_GDD|lifespan", name)) |>
                    inner_join(param_key, by=join_by(name)) |>
                    filter(grepl("lifespan|moltF", label)))
  p_pMolt <- post_summary_plot(pMolt_ls, ncol=4, scales="free_y")

  detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                    ~.x |> filter(grepl("detect_p", name)) |>
                      inner_join(param_key, by=join_by(name)))
  p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
    xlim(0, 1)

  else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                 ~.x |> filter(grepl("IP_bg|nb_prec", name)) |>
                   inner_join(param_key, by=join_by(name)))
  p_else <- post_summary_plot(else_ls, ncol=3, scales="free")

  plot_grid(p_ensWts, p_attach, p_surv, p_pMolt, p_detectp, p_else,
            nrow=6, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 1, 1, 1))
  ggsave(glue("{dat_dir}/fig_full_pars.png"), width=10, height=3)

  # IEM: full pop model -----------------------------------------------------

  out_full <- stan(file="code/stan/IEM_pop5stage_mnDaysPrior.stan",
                   data=stan_dat$dat, chains=3, cores=3, iter=40, refresh=1,
                   # control=list(max_treedepth=20),
                   init=0, seed=101,
                   pars=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
                          "surv_beta", "pMoltF_beta",
                          "detect_p", "nb_prec"))
  # saveRDS(out_full, glue("{dat_dir}/out_full_stanfit.rds"))


  out_full_df <- as.data.frame(out_full) |>
    select(-lp__) |>
    mutate(iter=row_number()) |>
    pivot_longer(-iter)
  out_full_sum <- out_full_df |>
    group_by(name) |>
    sevcheck::get_intervals(value, type="qi")
  dat_full_df <- tibble(
    name=c(paste0("attach_beta[", 1:5, "]"),
           "IP_bg", "IP_bg_m3",
           paste0("ensWts_p[", 1:6, "]"),
           paste0("surv_beta[1,", 1:3, "]"), paste0("surv_beta[2,", 1:3, "]"),
           paste0("thresh_GDD[", 1:2, ",1]"), paste0("thresh_GDD[", 1:2, ",2]"),
           "lifespan",
           paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
    value=c(stan_dat$params$attach_beta,
            stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
            stan_dat$params$ensWts_p,
            t(stan_dat$params$surv_beta),
            stan_dat$params$thresh_GDD,
            stan_dat$params$lifespan,
            stan_dat$params$detect_p,
            stan_dat$params$nb_prec)
  )

  ggplot(out_full_df) +
    # stat_pointinterval(.width=c(0.5, 0.8, 0.95)) +
    geom_density(aes(x=value)) +
    geom_pointrange(data=out_full_sum,
                    aes(x=mn, xmin=L25, xmax=L75, y=0), linewidth=1.25) +
    geom_linerange(data=out_full_sum,
                   aes(xmin=L10, xmax=L90, y=0), linewidth=0.8) +
    geom_linerange(data=out_full_sum,
                   aes(xmin=L025, xmax=L975, y=0), linewidth=0.5) +
    ylim(0, NA) +
    geom_point(data=dat_full_df,
               aes(x=value, y=0), colour="red", shape=1) +
    facet_wrap(~name, scales="free", ncol=4)
  ggsave(glue("{dat_dir}/fig_full_pars.png"), width=10, height=3)

}






# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Fit integrated ensemble



# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(cmdstanr)
library(ggdist)
library(cowplot)
source("code/fn/helpers.R")
theme_set(theme_classic())

n_chains <- 3
stages <- c("Ch", "PA", "Ad", "Gr")
stage_trans <- c("Ch-PA", "PA-Ad", "Ad-Gr")
param_key <- tibble(name=c(paste0("attach_beta[", 1:5, "]"),
                           paste0("ensWts_p[", 1:6, "]"),
                           "IP_bg", "IP_bg_m3",
                           paste0("surv_beta[1,", 1:4, "]"),
                           paste0("surv_beta[2,", 1:4, "]"),
                           paste0("pMoltF_beta[1,", 1:3, "]"),
                           paste0("pMoltF_beta[2,", 1:3, "]"),
                           paste0("pMoltM_beta[1,", 1:2, "]"),
                           paste0("pMoltM_beta[2,", 1:2, "]"),
                           paste0("thresh_GDD[", 1:3, ",1]"),
                           paste0("thresh_GDD[", 1:2, ",2]"),
                           "lifespan",
                           paste0("detect_p[", 1:4, "]"),
                           "nb_prec"),
                    label=c(paste0("attach_", c("Int", "RW", "Sal", "UV", "UVsq")),
                            paste0("ensWt_", 1:6),
                            "IP_bg", "IP_bg_m3",
                            paste0("surv_Int_", stages),
                            paste0("surv_Sal_", stages),
                            paste0("pMoltF_Int_", stage_trans),
                            paste0("pMoltF_Temp_", stage_trans),
                            paste0("pMoltM_Int_", stage_trans[-3]),
                            paste0("pMoltM_Temp_", stage_trans[-3]),
                            paste0("moltF_GDD_", stage_trans),
                            paste0("moltM_GDD_", stage_trans[-3]),
                            "lifespan_GDD",
                            paste0("p_detect_", stages),
                            "neg_binom_prec"
                    )) |>
  mutate(label=factor(label, levels=unique(label)))


for(sim in 1:30) {

  dat_dir <- glue("data/sim/sim_{str_pad(sim, 2, 'left', '0')}/")
  stan_dat <- make_stan_data(dat_dir)


  # IEM pt 1: influx -> attachment ------------------------------------------

  # mod_pt1 <- cmdstan_model("code/stan/IEM_pt1_ensAttach.stan")
  # fit_pt1 <- mod_pt1$sample(
  #   data=stan_dat$dat, init=0,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_pt1$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_pt1.csv"))
  #
  # out_pt1_df <- fit_pt1$draws(
  #   variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_pt1_df, glue("{dat_dir}posterior_pt1.csv"))
  #
  # out_pt1_sum <- out_pt1_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_pt1_sum, glue("{dat_dir}posterior_summary_pt1.csv"))
  #
  # dat_pt1_df <- tibble(
  #   name=c(paste0("attach_beta[", 1:5, "]"),
  #          "nb_prec", "IP_bg_m3",
  #          paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
  #          "IP_bg"),
  #   value=c(stan_dat$params$attach_beta,
  #           stan_dat$params$nb_prec, stan_dat$params$IP_bg_m3,
  #           stan_dat$params$ensWts_p,
  #           stan_dat$params$IP_bg)
  # )
  # write_csv(dat_pt1_df, glue("{dat_dir}params_pt1.csv"))
  #
  #
  # attach_ls <- map(list(out_pt1_df, out_pt1_sum, dat_pt1_df),
  #                  ~.x |> filter(grepl("attach_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_attach <- post_summary_plot(attach_ls, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # ensWts_ls <- map(list(out_pt1_df, out_pt1_sum, dat_pt1_df),
  #                  ~.x |> filter(grepl("ensWts", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # else_ls <- map(list(out_pt1_df, out_pt1_sum, dat_pt1_df),
  #                  ~.x |> filter(!grepl("ensWts|attach_beta", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=3, scales="free")
  #
  # plot_grid(p_attach, p_ensWts, p_else, nrow=3, align="hv", axis="trbl")
  # ggsave(glue("{dat_dir}/fig_pt1_pars.png"), width=10, height=6)



  # IEM pt 2: attachment -> adults ------------------------------------------

  # mod_pt2 <- cmdstan_model("code/stan/IEM_pt2_onFish_pop.stan")
  # fit_pt2 <- mod_pt2$sample(
  #   data=stan_dat$dat, init=0, seed=101,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_pt2$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, autodiff_calls, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_pt2_pop.csv"))
  # # cohort_N > N > cohort_stage > cohort_surv
  #
  # out_pt2_df <- fit_pt2$draws(
  #   variables=c("surv_beta", "pMoltF_beta", "pMoltM_beta", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_pt2_df, glue("{dat_dir}posterior_pt2_pop.csv"))
  # out_pt2_sum <- out_pt2_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_pt2_sum, glue("{dat_dir}posterior_summary_pt2_pop.csv"))
  #
  # dat_pt2_df <- tibble(
  #   name=c(paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
  #          "nb_prec"),
  #   value=c(t(stan_dat$params$surv_beta),
  #           stan_dat$params$nb_prec)
  # )
  # write_csv(dat_pt2_df, glue("{dat_dir}params_pt2.csv"))
  #
  # surv_ls <- map(list(out_pt2_df, out_pt2_sum, dat_pt2_df),
  #                  ~.x |> filter(grepl("surv_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltInt_ls <- map(list(out_pt2_df, out_pt2_sum, dat_pt2_df),
  #                  ~.x |> filter(grepl("pMolt.*_beta\\[1", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltTemp_ls <- map(list(out_pt2_df, out_pt2_sum, dat_pt2_df),
  #                     ~.x |> filter(grepl("pMolt.*_beta\\[2", name)) |>
  #                       inner_join(param_key, by=join_by(name)))
  # p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # plot_grid(p_surv, p_pMoltInt, p_pMoltTemp, nrow=3, align="hv", axis="trbl",
  #           rel_heights=c(1, 0.6, 0.6))
  #
  # ggsave(glue("{dat_dir}/fig_pt2_pop_pars.png"), width=10, height=10)



  # IEM pt 3: detection -----------------------------------------------------

  # mod_pt3 <- cmdstan_model("code/stan/IEM_pt3_detection.stan")
  # fit_pt3 <- mod_pt3$sample(
  #   data=stan_dat$dat, init=0,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_pt3$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_pt3.csv"))
  #
  # out_pt3_df <- fit_pt3$draws(
  #   variables=c("detect_p", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_pt3_df, glue("{dat_dir}posterior_pt3.csv"))
  #
  # out_pt3_sum <- out_pt3_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_pt3_sum, glue("{dat_dir}posterior_summary_pt3.csv"))
  #
  # dat_pt3_df <- tibble(
  #   name=c(paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
  #   value=c(stan_dat$params$detect_p,
  #           stan_dat$params$nb_prec)
  # )
  # write_csv(dat_pt3_df, glue("{dat_dir}params_pt3.csv"))
  #
  #
  # detectp_ls <- map(list(out_pt3_df, out_pt3_sum, dat_pt3_df),
  #                  ~.x |> filter(grepl("detect_p", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # else_ls <- map(list(out_pt3_df, out_pt3_sum, dat_pt3_df),
  #                ~.x |> filter(!grepl("detect_p", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=1, scales="free")
  #
  # plot_grid(p_detectp, p_else, nrow=1, align="hv", axis="trbl", rel_widths=c(1, 0.25))
  # ggsave(glue("{dat_dir}/fig_pt3_pars.png"), width=10, height=3)






  # IEM: pt12 model pop -----------------------------------------------------

  # mod_pt12 <- cmdstan_model("code/stan/IEM_pt12_pop.stan")
  # fit_pt12 <- mod_pt12$sample(
  #   data=stan_dat$dat, init=0, seed=101,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_pt12$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_pt12_pop.csv"))
  #
  # out_pt12_df <- fit_pt12$draws(
  #   variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
  #               "surv_beta", "pMoltF_beta", "pMoltM_beta",
  #               "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_pt12_df, glue("{dat_dir}posterior_pt12_pop.csv"))
  # out_pt12_sum <- out_pt12_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_pt12_sum, glue("{dat_dir}posterior_summary_pt12_pop.csv"))
  #
  # dat_pt12_df <- tibble(
  #   name=c(paste0("attach_beta[", 1:5, "]"),
  #          "IP_bg", "IP_bg_m3",
  #          paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
  #          paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
  #          paste0("thresh_GDD[", 1:3, ",1]"), paste0("thresh_GDD[", 1:3, ",2]"),
  #          "lifespan",
  #          paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
  #   value=c(stan_dat$params$attach_beta,
  #           stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
  #           stan_dat$params$ensWts_p,
  #           t(stan_dat$params$surv_beta),
  #           stan_dat$params$thresh_GDD,
  #           stan_dat$params$lifespan,
  #           stan_dat$params$detect_p,
  #           stan_dat$params$nb_prec)
  # )
  # write_csv(dat_pt12_df, glue("{dat_dir}params_pt12.csv"))
  #
  # ensWts_ls <- map(list(out_pt12_df, out_pt12_sum, dat_pt12_df),
  #                  ~.x |> filter(grepl("ensWts", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # attach_ls <- map(list(out_pt12_df, out_pt12_sum, dat_pt12_df),
  #                  ~.x |> filter(grepl("attach_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_attach <- post_summary_plot(attach_ls, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # surv_ls <- map(list(out_pt12_df, out_pt12_sum, dat_pt12_df),
  #                ~.x |> filter(grepl("surv_beta", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltInt_ls <- map(list(out_pt12_df, out_pt12_sum, dat_pt12_df),
  #                    ~.x |> filter(grepl("pMolt.*_beta\\[1", name)) |>
  #                      inner_join(param_key, by=join_by(name)))
  # p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltTemp_ls <- map(list(out_pt12_df, out_pt12_sum, dat_pt12_df),
  #                     ~.x |> filter(grepl("pMolt.*_beta\\[2", name)) |>
  #                       inner_join(param_key, by=join_by(name)))
  # p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # else_ls <- map(list(out_pt12_df, out_pt12_sum, dat_pt12_df),
  #                ~.x |> filter(grepl("IP_bg|nb_prec", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=3, scales="free")
  #
  # plot_grid(p_ensWts, p_attach, p_surv, p_pMoltInt, p_pMoltTemp,  p_else,
  #           nrow=6, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 1, 1, 1))
  # ggsave(glue("{dat_dir}/fig_pt12_pop_pars.png"), width=10, height=10)





  # IEM: full model pop PRIORS ----------------------------------------------

  # mod_full <- cmdstan_model("code/stan/IEM_full_pop_PRIORS.stan")
  # fit_full <- mod_full$sample(
  #   data=stan_dat$dat, init=0, seed=101,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  #
  # out_full_df <- fit_full$draws(
  #   variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
  #               "surv_beta", "pMoltF_beta", "pMoltM_beta",
  #               "detect_p", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_full_df, glue("{dat_dir}posterior_full_pop_PRIORS.csv"))
  # out_full_sum <- out_full_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_full_sum, glue("{dat_dir}posterior_summary_full_pop_PRIORS.csv"))
  #
  # dat_full_df <- tibble(
  #   name=c(paste0("attach_beta[", 1:5, "]"),
  #          "IP_bg", "IP_bg_m3",
  #          paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
  #          paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
  #          paste0("thresh_GDD[", 1:3, ",1]"), paste0("thresh_GDD[", 1:3, ",2]"),
  #          "lifespan",
  #          paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
  #   value=c(stan_dat$params$attach_beta,
  #           stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
  #           stan_dat$params$ensWts_p,
  #           t(stan_dat$params$surv_beta),
  #           stan_dat$params$thresh_GDD,
  #           stan_dat$params$lifespan,
  #           stan_dat$params$detect_p,
  #           stan_dat$params$nb_prec)
  # )
  # write_csv(dat_full_df, glue("{dat_dir}params_full.csv"))
  #
  # ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("ensWts", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # attach_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("attach_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_attach <- post_summary_plot(attach_ls, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # surv_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("surv_beta", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltInt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                    ~.x |> filter(grepl("pMolt.*_beta\\[1", name)) |>
  #                      inner_join(param_key, by=join_by(name)))
  # p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltTemp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                     ~.x |> filter(grepl("pMolt.*_beta\\[2", name)) |>
  #                       inner_join(param_key, by=join_by(name)))
  # p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                   ~.x |> filter(grepl("detect_p", name)) |>
  #                     inner_join(param_key, by=join_by(name)))
  # p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("IP_bg|nb_prec", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=3, scales="free")
  #
  # plot_grid(p_ensWts, p_attach, p_surv, p_pMoltInt, p_pMoltTemp, p_detectp, p_else,
  #           nrow=7, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 1, 1, 1, 1))
  # ggsave(glue("{dat_dir}/fig_full_pop_pars_PRIORS.png"), width=10, height=12)



  # IEM: full model pop F ---------------------------------------------------

  # mod_full <- cmdstan_model("code/stan/IEM_full_pop_FEMALE.stan")
  # fit_full <- mod_full$sample(
  #   data=stan_dat$dat, init=0, seed=101,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_full$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_full_pop_FEMALE.csv"))
  #
  # out_full_df <- fit_full$draws(
  #   variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
  #               "surv_beta", "pMoltF_beta",
  #               "detect_p", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_full_df, glue("{dat_dir}posterior_full_pop_FEMALE.csv"))
  # out_full_sum <- out_full_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_full_sum, glue("{dat_dir}posterior_summary_full_pop_FEMALE.csv"))
  #
  # dat_full_df <- tibble(
  #   name=c(paste0("attach_beta[", 1:5, "]"),
  #          "IP_bg", "IP_bg_m3",
  #          paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
  #          paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
  #          paste0("thresh_GDD[", 1:3, ",1]"), paste0("thresh_GDD[", 1:3, ",2]"),
  #          "lifespan",
  #          paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
  #   value=c(stan_dat$params$attach_beta,
  #           stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
  #           stan_dat$params$ensWts_p,
  #           t(stan_dat$params$surv_beta),
  #           stan_dat$params$thresh_GDD,
  #           stan_dat$params$lifespan,
  #           stan_dat$params$detect_p,
  #           stan_dat$params$nb_prec)
  # )
  # write_csv(dat_full_df, glue("{dat_dir}params_full_FEMALE.csv"))
  #
  # ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("ensWts", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # attach_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("attach_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_attach <- post_summary_plot(attach_ls, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # surv_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("surv_beta", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltInt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                    ~.x |> filter(grepl("pMolt.*_beta\\[1", name)) |>
  #                      inner_join(param_key, by=join_by(name)))
  # p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltTemp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                     ~.x |> filter(grepl("pMolt.*_beta\\[2", name)) |>
  #                       inner_join(param_key, by=join_by(name)))
  # p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                   ~.x |> filter(grepl("detect_p", name)) |>
  #                     inner_join(param_key, by=join_by(name)))
  # p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("IP_bg|nb_prec", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=3, scales="free")
  #
  # plot_grid(p_ensWts, p_attach, p_surv, p_pMoltInt, p_pMoltTemp, p_detectp, p_else,
  #           nrow=7, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 1, 1, 1, 1))
  # ggsave(glue("{dat_dir}/fig_full_pop_FEMALE_pars.png"), width=10, height=12)



  # IEM: full model pop -----------------------------------------------------

  # mod_full <- cmdstan_model("code/stan/IEM_full_pop.stan")
  # fit_full <- mod_full$sample(
  #   data=stan_dat$dat, init=0, seed=101,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_full$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_full_pop.csv"))
  #
  # out_full_df <- fit_full$draws(
  #   variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
  #               "surv_beta", "pMoltF_beta", "pMoltM_beta",
  #               "detect_p", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_full_df, glue("{dat_dir}posterior_full_pop.csv"))
  # out_full_sum <- out_full_df |>
  #   group_by(name) |>
  #   sevcheck::get_intervals(value, type="qi")
  # write_csv(out_full_sum, glue("{dat_dir}posterior_summary_full_pop.csv"))
  #
  # dat_full_df <- tibble(
  #   name=c(paste0("attach_beta[", 1:5, "]"),
  #          "IP_bg", "IP_bg_m3",
  #          paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
  #          paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
  #          paste0("thresh_GDD[", 1:3, ",1]"), paste0("thresh_GDD[", 1:3, ",2]"),
  #          "lifespan",
  #          paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
  #   value=c(stan_dat$params$attach_beta,
  #           stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
  #           stan_dat$params$ensWts_p,
  #           t(stan_dat$params$surv_beta),
  #           stan_dat$params$thresh_GDD,
  #           stan_dat$params$lifespan,
  #           stan_dat$params$detect_p,
  #           stan_dat$params$nb_prec)
  # )
  # write_csv(dat_full_df, glue("{dat_dir}params_full.csv"))
  #
  # ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("ensWts", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # attach_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("attach_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_attach <- post_summary_plot(attach_ls, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # surv_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("surv_beta", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltInt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                    ~.x |> filter(grepl("pMolt.*_beta\\[1", name)) |>
  #                      inner_join(param_key, by=join_by(name)))
  # p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltTemp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                     ~.x |> filter(grepl("pMolt.*_beta\\[2", name)) |>
  #                       inner_join(param_key, by=join_by(name)))
  # p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=5, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                   ~.x |> filter(grepl("detect_p", name)) |>
  #                     inner_join(param_key, by=join_by(name)))
  # p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("IP_bg|nb_prec", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=3, scales="free")
  #
  # plot_grid(p_ensWts, p_attach, p_surv, p_pMoltInt, p_pMoltTemp, p_detectp, p_else,
  #           nrow=7, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 1, 1, 1, 1))
  # ggsave(glue("{dat_dir}/fig_full_pop_pars.png"), width=10, height=12)



  # IEM: full model F -------------------------------------------------------

  iter <- 1000
  mod_full <- cmdstan_model("code/stan/IEM_full_FEMALE.stan")
  fit_full <- mod_full$sample(
    data=stan_dat$dat, init=0, seed=101,
    iter_warmup=iter, iter_sampling=iter, refresh=1,
    chains=n_chains, parallel_chains=n_chains
  )
  fit_full$profiles()[[1]] |>
    as_tibble() |>
    select(name, total_time) |>
    mutate(pct=total_time/sum(total_time)*100) |>
    arrange(desc(pct)) |>
    write_csv(glue("{dat_dir}profile_full_FEMALE.csv"))

  out_full_df <- fit_full$draws(
    variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
                "surv_beta", "thresh_GDD", "lifespan",
                "detect_p", "nb_prec"),
    format="df") |>
    pivot_longer(-starts_with("."))
  write_csv(out_full_df, glue("{dat_dir}posterior_full_FEMALE.csv"))
  out_full_sum <- out_full_df |>
    group_by(name) |>
    sevcheck::get_intervals(value, type="qi")
  write_csv(out_full_sum, glue("{dat_dir}posterior_summary_full_FEMALE.csv"))

  dat_full_df <- tibble(
    name=c(paste0("attach_beta[", 1:5, "]"),
           "IP_bg", "IP_bg_m3",
           paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
           paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
           paste0("thresh_GDD[", 1:3, ",1]"), paste0("thresh_GDD[", 1:3, ",2]"),
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
  p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
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
  ggsave(glue("{dat_dir}/fig_full_FEMALE_pars.png"), width=10, height=12)



  # IEM: full model ---------------------------------------------------------

  # iter <- 1000
  # mod_full <- cmdstan_model("code/stan/IEM_full.stan")
  # fit_full <- mod_full$sample(
  #   data=stan_dat$dat, init=0, seed=101,
  #   iter_warmup=iter, iter_sampling=iter, refresh=10,
  #   chains=n_chains, parallel_chains=n_chains
  # )
  # fit_full$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_full.csv"))
  #
  # out_full_df <- fit_full$draws(
  #   variables=c("IP_bg", "IP_bg_m3", "ensWts_p", "attach_beta",
  #               "surv_beta", "thresh_GDD", "lifespan",
  #               "detect_p", "nb_prec"),
  #   format="df") |>
  #   pivot_longer(-starts_with("."))
  # write_csv(out_full_df, glue("{dat_dir}posterior_full.csv"))
  # out_full_sum <- out_full_df |>
  #     group_by(name) |>
  #     sevcheck::get_intervals(value, type="qi")
  # write_csv(out_full_sum, glue("{dat_dir}posterior_summary_full.csv"))
  #
  # dat_full_df <- tibble(
  #   name=c(paste0("attach_beta[", 1:5, "]"),
  #          "IP_bg", "IP_bg_m3",
  #          paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
  #          paste0("surv_beta[1,", 1:4, "]"), paste0("surv_beta[2,", 1:4, "]"),
  #          paste0("thresh_GDD[", 1:3, ",1]"), paste0("thresh_GDD[", 1:3, ",2]"),
  #          "lifespan",
  #          paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec"),
  #   value=c(stan_dat$params$attach_beta,
  #           stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3,
  #           stan_dat$params$ensWts_p,
  #           t(stan_dat$params$surv_beta),
  #           stan_dat$params$thresh_GDD,
  #           stan_dat$params$lifespan,
  #           stan_dat$params$detect_p,
  #           stan_dat$params$nb_prec)
  # )
  #
  # attach_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("attach_beta", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_attach <- post_summary_plot(attach_ls, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                  ~.x |> filter(grepl("ensWts", name)) |>
  #                    inner_join(param_key, by=join_by(name)))
  # p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # surv_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("surv_beta", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_surv <- post_summary_plot(surv_ls, ncol=4, scales="free") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMolt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                    ~.x |> filter(grepl("thresh_GDD|lifespan", name)) |>
  #                      inner_join(param_key, by=join_by(name)))
  # p_pMolt <- post_summary_plot(pMolt_ls, ncol=3, scales="free_y")
  #
  # detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                   ~.x |> filter(grepl("detect_p", name)) |>
  #                     inner_join(param_key, by=join_by(name)))
  # p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
  #   xlim(0, 1)
  #
  # else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                ~.x |> filter(grepl("IP_bg|nb_prec", name)) |>
  #                  inner_join(param_key, by=join_by(name)))
  # p_else <- post_summary_plot(else_ls, ncol=3, scales="free")
  #
  # plot_grid(p_ensWts, p_attach, p_surv, p_pMolt, p_detectp, p_else,
  #           nrow=6, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 2, 1, 1))
  # ggsave(glue("{dat_dir}/fig_full_pars.png"), width=10, height=12)

}






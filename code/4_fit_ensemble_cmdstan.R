# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Fit integrated ensemble

# TODO: Simplify data inputs based on final structure; calculate from real data
# TODO: Confirm all priors -- hard code in make_stan_data()?

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

prior_only <- T
keep_licePreds <- T
refit <- T
n_parallel <- 3

n_chains <- 3
stages <- c("Ch", "PA", "Ad")
stage_trans <- c("Ch-PA", "PA-Ad")
param_key <- tibble(name=c(paste0("attach_beta[", 1:5, "]"),
                           paste0("ensWts_p[", 1:20, "]"),
                           "IP_bg", "IP_bg_m3",
                           paste0("surv_beta[1,", 1:3, "]"),
                           paste0("surv_beta[2,", 1:3, "]"),
                           paste0("surv_int_farm_sd[", 1:3, "]"),
                           paste0("mnDaysStage_beta[1,", 1:2, "]"),
                           paste0("mnDaysStage_beta[2,", 1:2, "]"),
                           paste0("thresh_GDD[", 1:2, ",1]"),
                           paste0("thresh_GDD[", 1:2, ",2]"),
                           "lifespan",
                           paste0("detect_p[", 1:2, "]"),
                           "nb_prec",
                           "IP_scale", "IP_halfSat_m3",
                           "treatEfficacy"),
                    label=c(paste0("attach_", c("RW", "Sal", "Temp", "UV", "UVsq")),
                            paste0("ensWt_", 1:20),
                            "IP_bg", "IP_bg_m3",
                            paste0("surv_Int_", stages),
                            paste0("surv_Sal_", stages),
                            paste0("surv_int_farm_sd_", stages),
                            paste0("mnDaysStage_Int_", stages[1:2]),
                            paste0("mnDaysStage_Temp_", stages[1:2]),
                            paste0("moltF_GDD_", stage_trans),
                            paste0("moltM_GDD_", stage_trans),
                            "lifespan_GDD",
                            paste0("p_detect_", stages[-3]),
                            "neg_binom_prec",
                            "IP_scale", "IP_halfSat_m3",
                            "treatEfficacy"
                    )) |>
  mutate(label=factor(label, levels=unique(label)))

sim_dirs <- paste0(dir("data/sim", "sim_", include.dirs=T, full.names=T), "/")

plan(multicore, workers=n_parallel)

foreach(sim_dir=sim_dirs, .errorhandling="pass") %dofuture% {
  if(!refit & file.exists(glue("{sim_dir}posterior_summary{ifelse(prior_only, '_PRIORS', '')}.rds"))) {
    out_full_df <- readRDS(glue("{sim_dir}posterior{ifelse(prior_only, '_PRIORS', '')}.rds"))
    out_full_sum <- readRDS(glue("{sim_dir}posterior_summary{ifelse(prior_only, '_PRIORS', '')}.rds"))
    dat_full_df <- fread(glue("{sim_dir}stan_params{ifelse(prior_only, '_PRIORS', '')}.csv")) |>
      as_tibble()
  } else {
    stan_dat <- make_stan_data(sim_dir, priors_only=prior_only)

    # IEM: full model pop -----------------------------------------------------

    iter <- 1000
    mod_full <- cmdstan_model("code/stan/joint_population_model.stan")
    fit_full <- mod_full$sample(
      data=stan_dat$dat, init=0, seed=101, refresh=max(iter/100, 1),
      iter_warmup=iter, iter_sampling=iter,
      chains=n_chains, parallel_chains=n_chains
    )

    keep_pars <- c("IP_bg_m3",
                   "ensWts_p", "attach_beta",
                   "surv_beta", "surv_int_farm_sd", "mnDaysStage_beta",
                   "detect_p", "nb_prec", "treatEfficacy")
    if(keep_licePreds) {
      keep <- c(keep_pars, "mu")
    } else {
      keep <- keep_pars
    }
    out_full_df <- fit_full$draws(
      variables=keep,
      format="df") |>
      pivot_longer(-starts_with("."))
    saveRDS(out_full_df, glue("{sim_dir}posterior{ifelse(prior_only, '_PRIORS', '')}.rds"))
    out_full_sum <- out_full_df |>
      group_by(name) |>
      sevcheck::get_intervals(value, type="qi")
    saveRDS(out_full_sum, glue("{sim_dir}posterior_summary{ifelse(prior_only, '_PRIORS', '')}.rds"))

    dat_full_df <- tibble(
      name=c(paste0("attach_beta[", 1:5, "]"),
             "IP_bg_m3",
             paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
             paste0("surv_beta[1,", 1:3, "]"), paste0("surv_beta[2,", 1:3, "]"),
             paste0("surv_int_farm_sd[", 1:3, "]"),
             paste0("thresh_GDD[", 1:2, ",1]"), paste0("thresh_GDD[", 1:2, ",2]"),
             "lifespan",
             paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec",
             "treatEfficacy"),
      value=c(stan_dat$params$attach_beta,
              stan_dat$params$IP_bg_m3,
              stan_dat$params$ensWts_p,
              t(stan_dat$params$surv_beta),
              stan_dat$params$surv_int_farm_sd,
              stan_dat$params$thresh_GDD,
              stan_dat$params$lifespan,
              stan_dat$params$detect_p,
              stan_dat$params$nb_prec,
              stan_dat$params$treat_efficacy)
    )
    write_csv(dat_full_df, glue("{sim_dir}stan_params{ifelse(prior_only, '_PRIORS', '')}.csv"))

  }

  info <- readRDS(glue("{sim_dir}info.rds"))

  ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                   ~.x |> filter(grepl("ensWts", name)) |>
                     inner_join(param_key, by=join_by(name)))
  p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y", ncol=10) +
    xlim(0, 1)

  attach_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                   ~.x |> filter(grepl("attach_beta", name)) |>
                     inner_join(param_key, by=join_by(name)))
  p_attach <- post_summary_plot(attach_ls, scales="free") +
    geom_vline(xintercept=0, linetype=3)

  surv_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                 ~.x |> filter(grepl("surv_beta", name)) |>
                   inner_join(param_key, by=join_by(name)))
  p_surv <- post_summary_plot(surv_ls, ncol=3, scales="free") +
    geom_vline(xintercept=0, linetype=3)

  surv_sd_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                 ~.x |> filter(grepl("surv_int_farm_sd", name)) |>
                   inner_join(param_key, by=join_by(name)))
  p_surv_sd <- post_summary_plot(surv_sd_ls, ncol=3, scales="free") +
    geom_vline(xintercept=0, linetype=3)

  # pMoltInt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                    ~.x |> filter(grepl("mnDaysStage_beta\\[1", name)) |>
  #                      inner_join(param_key, by=join_by(name)))
  # p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=2, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)
  #
  # pMoltTemp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
  #                     ~.x |> filter(grepl("mnDaysStage_beta\\[2", name)) |>
  #                       inner_join(param_key, by=join_by(name)))
  # p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=2, scales="free_y") +
  #   geom_vline(xintercept=0, linetype=3)

  pMoltTemp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                      ~.x |> filter(grepl("mnDaysStage_beta", name)) |>
                        inner_join(param_key, by=join_by(name)))
  p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=4, scales="free_y") +
    geom_vline(xintercept=0, linetype=3)

  detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                    ~.x |> filter(grepl("detect_p", name)) |>
                      inner_join(param_key, by=join_by(name)))
  p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
    xlim(0, 1)

  else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                 ~.x |> filter(grepl("IP_bg|nb_prec|IP_halfSat|treat", name)) |>
                   inner_join(param_key, by=join_by(name)))
  p_else <- post_summary_plot(else_ls, ncol=5, scales="free")

  p <- plot_grid(p_ensWts, p_attach, p_surv, p_surv_sd, p_pMoltTemp, p_detectp, p_else,
            nrow=7, align="hv", axis="trbl", rel_heights=c(2, 1, 2, 1, 1, 1, 1))
  ggsave(glue("{sim_dir}/fig_pars{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=12)

  if(keep_licePreds) {
    draws <- sample.int(max(out_full_df$.draw), 100)
    mu_draws_df <- out_full_df |>
      filter(.draw %in% draws & grepl("mu", name)) |>
      separate_wider_delim(name, delim=",", names=c("farm", "stage5", "day")) |>
      mutate(farm=str_sub(farm, 4, -1),
             day=as.numeric(str_sub(day, 1, -2)),
             stage=case_when(stage5 %in% 1:2 ~ "Ch",
                             stage5 %in% 3:4 ~ "PA",
                             stage5 == 5 ~ "Ad"),
             stage=factor(stage, levels=c("Ch", "PA", "Ad"))) |>
      group_by(.chain, .iteration, .draw, farm, stage, day) |>
      summarise(value=sum(value)) |>
      ungroup() |>
      mutate(type="Fitted",
             value=pmax(value, 0)) |>
      rename(mu=value) |>
      bind_rows(expand_grid(farm=as.character(1:info$nFarms),
                            day=1:info$nDays,
                            stage=factor(c("Ch", "PA", "Ad"), levels=c("Ch", "PA", "Ad"))) |>
                  mutate(mu=c(readRDS(glue("{sim_dir}/mu.rds"))[,1,,])) |>
                  mutate(type="True")) |>
      mutate(day=ymd("2023-01-01") + day - 1,
             farm=paste("Farm", farm))
    p <- mu_draws_df |>
      filter(stage=="Ad") |>
      ggplot(aes(day, mu, group=as.character(.draw), colour=type, alpha=type)) +
      geom_line() +
      scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
      scale_alpha_manual("", values=c("True"=1, "Fitted"=0.1)) +
      labs(x="Date", y="Mean lice per fish (latent) [50% CI]") +
      {if(any((mu_draws_df |> filter(stage=="Ad"))$mu > 15)) scale_y_continuous(limits=c(0, 15), oob=scales::oob_keep)} +
      scale_x_date(date_labels="%b") +
      facet_grid(farm~., scales="free_y")
    ggsave(glue("{sim_dir}/fig_mu_draws{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=15)

    lice_thresh <- 0.5
    mu_true <- mu_draws_df |>
      filter(is.na(.draw)) |>
      select(farm, stage, day, mu) |>
      rename(mu_true=mu)
    mu_draw_metrics <- mu_draws_df |>
      filter(!is.na(.draw)) |>
      filter(stage=="Ad") |>
      left_join(mu_true, by=join_by(farm, stage, day)) |>
      group_by(.draw) |>
      summarise(TP=sum(mu_true > lice_thresh & mu > lice_thresh),
                TN=sum(mu_true < lice_thresh & mu < lice_thresh),
                FP=sum(mu_true < lice_thresh & mu > lice_thresh),
                FN=sum(mu_true > lice_thresh & mu < lice_thresh),
                N=n()) |>
      mutate(accuracy=(TP+TN)/(TP+TN+FP+FN),
             PPV=TP/(TP+FP),
             NPV=TN/(TN+FN),
             TPR=TP/(TP+FN),
             TNR=TN/(TN+FP),
             S=(TP+FN)/N,
             P=(TP+FP)/N,
             F1=(2*PPV*TPR)/(PPV+TPR),
             MCC=(TP/N - S*P)/sqrt(P*S*(1-S)*(1-P)),
             TSS=(TP*TN - FP*FN)/((TP+FN)*(TN+FP)),
             accBal=(TPR+TNR)/2) |>
      ungroup() |>
      select(-S, -P, -N)

    mu_draw_metrics_wk <- mu_draws_df |>
      filter(!is.na(.draw)) |>
      filter(stage=="Ad") |>
      left_join(mu_true, by=join_by(farm, stage, day)) |>
      mutate(wk=floor_date(day, "week")) |>
      group_by(.draw, farm, stage, wk) |>
      summarise(mu=max(mu), mu_true=max(mu_true)) |>
      group_by(.draw) |>
      summarise(TP=sum(mu_true > lice_thresh & mu > lice_thresh),
                TN=sum(mu_true < lice_thresh & mu < lice_thresh),
                FP=sum(mu_true < lice_thresh & mu > lice_thresh),
                FN=sum(mu_true > lice_thresh & mu < lice_thresh),
                N=n()) |>
      mutate(accuracy=(TP+TN)/(TP+TN+FP+FN),
             PPV=TP/(TP+FP),
             NPV=TN/(TN+FN),
             TPR=TP/(TP+FN),
             TNR=TN/(TN+FP),
             S=(TP+FN)/N,
             P=(TP+FP)/N,
             F1=(2*PPV*TPR)/(PPV+TPR),
             MCC=(TP/N - S*P)/sqrt(P*S*(1-S)*(1-P)),
             TSS=(TP*TN - FP*FN)/((TP+FN)*(TN+FP)),
             accBal=(TPR+TNR)/2) |>
      ungroup() |>
      select(-S, -P, -N)

    write_csv(mu_draw_metrics, glue("{sim_dir}/mu_fitted_metrics_daily{ifelse(prior_only, '_PRIORS', '')}.csv"))
    write_csv(mu_draw_metrics_wk, glue("{sim_dir}/mu_fitted_metrics_weekly{ifelse(prior_only, '_PRIORS', '')}.csv"))


    if("y_pred" %in% keep_pars) {
      sampledDays <- readRDS(glue("{sim_dir}/sampledDays.rds")) |>
        as_tibble() |>
        rename(farm=sepaSite) |>
        mutate(farm=as.character(farm),
               sample=row_number()) |>
        inner_join(readRDS(glue("{sim_dir}/nFishSampled_mx.rds")) |>
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
                    mutate(y=c(readRDS(glue("{sim_dir}/y.rds"))[,1,,])) |>
                    mutate(type="True") |>
                    inner_join(sampledDays |> select(-sample), by=join_by(farm, day))) |>
        mutate(day=ymd("2023-01-01") + day - 1,
               farm=paste("Farm", farm),
               y_perFish=y/nFishSampled)
      y_df |>
        saveRDS(glue("{sim_dir}/y_sim_fitted{ifelse(prior_only, '_PRIORS', '')}.rds"))

      p <- y_df |>
        ggplot(aes(day, y_perFish)) +
        geom_linerange(aes(ymin=L10/nFishSampled, ymax=L90/nFishSampled, group=type), linewidth=0.25) +
        geom_point(aes(colour=type, shape=type)) +
        scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
        scale_shape_manual("", values=c("True"=4, "Fitted"=1)) +
        labs(x="Date", y="Mean lice per fish (observed) [50% CI]") +
        scale_x_date(date_labels="%b") +
        facet_grid(stage~farm, scales="free_y")
      ggsave(glue("{sim_dir}/fig_y{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=15, height=7)

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
      ggsave(glue("{sim_dir}/fig_y_AF{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=6, height=15)
    }
  }
}

plan(sequential)

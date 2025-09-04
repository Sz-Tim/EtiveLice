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
dir("code/fn", ".R", full.names=T) |> walk(source)
theme_set(theme_classic())

prior_only <- FALSE
keep_licePreds <- FALSE

n_chains <- 3
stages <- c("Ch", "PA", "Ad")
stage_trans <- c("Ch-PA", "PA-Ad")
param_key <- tibble(name=c(paste0("attach_beta[", 1:5, "]"),
                           paste0("ensWts_p[", 1:20, "]"),
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
                            paste0("ensWt_", 1:20),
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

sim_dirs <- paste0(dir("data/sim", "sim_", include.dirs=T), "/")
for(sim_dir in sim_dirs) {

  dat_dir <- paste0("data/sims/", sim_dir)
  stan_dat <- make_stan_data(dat_dir, priors_only=prior_only)

  # IEM: full model pop -----------------------------------------------------

  iter <- 1000
  mod_full <- cmdstan_model("code/stan/joint_population_model.stan")
  fit_full <- mod_full$sample(
    data=stan_dat$dat, init=0, seed=101, refresh=max(iter/100, 1),
    iter_warmup=iter, iter_sampling=iter,
    chains=n_chains, parallel_chains=n_chains
  )
  # fit_full$profiles()[[1]] |>
  #   as_tibble() |>
  #   select(name, total_time) |>
  #   mutate(pct=total_time/sum(total_time)*100) |>
  #   arrange(desc(pct)) |>
  #   write_csv(glue("{dat_dir}profile_pop5stage.csv"))

  keep <- c("IP_bg", "IP_bg_m3", "IP_scale", "ensWts_p", "attach_beta",
            "surv_beta", "mnDaysStage_beta",
            "detect_p", "nb_prec", "treatEfficacy")
  if(keep_licePreds) keep <- c(keep, "mu", "y_pred")
  out_full_df <- fit_full$draws(
    variables=keep,
    format="df") |>
    pivot_longer(-starts_with("."))
  write_csv(out_full_df, glue("{dat_dir}posterior_pop5stage{ifelse(prior_only, '_PRIORS', '')}.csv"))
  out_full_sum <- out_full_df |>
    group_by(name) |>
    sevcheck::get_intervals(value, type="qi")
  write_csv(out_full_sum, glue("{dat_dir}posterior_summary_pop_5stage{ifelse(prior_only, '_PRIORS', '')}.csv"))

  dat_full_df <- tibble(
    name=c(paste0("attach_beta[", 1:5, "]"),
           "IP_bg", "IP_bg_m3", "IP_scale",
           paste0("ensWts_p[", 1:stan_dat$dat$nSims, "]"),
           paste0("surv_beta[1,", 1:3, "]"), paste0("surv_beta[2,", 1:3, "]"),
           paste0("thresh_GDD[", 1:2, ",1]"), paste0("thresh_GDD[", 1:2, ",2]"),
           "lifespan",
           paste0("detect_p[", 1:stan_dat$dat$nStages, "]"), "nb_prec",
           "treatEfficacy"),
    value=c(stan_dat$params$attach_beta,
            stan_dat$params$IP_bg, stan_dat$params$IP_bg_m3, stan_dat$params$IP_scale,
            stan_dat$params$ensWts_p,
            t(stan_dat$params$surv_beta),
            stan_dat$params$thresh_GDD,
            stan_dat$params$lifespan,
            stan_dat$params$detect_p,
            stan_dat$params$nb_prec,
            stan_dat$params$treat_efficacy)
  )
  write_csv(dat_full_df, glue("{dat_dir}params_pop_5stage{ifelse(prior_only, '_PRIORS', '')}.csv"))

  ensWts_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                   ~.x |> filter(grepl("ensWts", name)) |>
                     inner_join(param_key, by=join_by(name)))
  p_ensWts <- post_summary_plot(ensWts_ls, scales="free_y", ncol=9) +
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

  pMoltInt_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                     ~.x |> filter(grepl("mnDaysStage_beta\\[1", name)) |>
                       inner_join(param_key, by=join_by(name)))
  p_pMoltInt <- post_summary_plot(pMoltInt_ls, ncol=2, scales="free_y") +
    geom_vline(xintercept=0, linetype=3)

  pMoltTemp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                      ~.x |> filter(grepl("mnDaysStage_beta\\[2", name)) |>
                        inner_join(param_key, by=join_by(name)))
  p_pMoltTemp <- post_summary_plot(pMoltTemp_ls, ncol=2, scales="free_y") +
    geom_vline(xintercept=0, linetype=3)

  detectp_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                    ~.x |> filter(grepl("detect_p", name)) |>
                      inner_join(param_key, by=join_by(name)))
  p_detectp <- post_summary_plot(detectp_ls, scales="free_y") +
    xlim(0, 1)

  else_ls <- map(list(out_full_df, out_full_sum, dat_full_df),
                 ~.x |> filter(grepl("IP_bg|nb_prec|IP_scale|treat", name)) |>
                   inner_join(param_key, by=join_by(name)))
  p_else <- post_summary_plot(else_ls, ncol=5, scales="free")

  p <- plot_grid(p_ensWts, p_attach, p_surv, p_pMoltInt, p_pMoltTemp, p_detectp, p_else,
            nrow=7, align="hv", axis="trbl", rel_heights=c(1, 1, 2, 1, 1, 1, 1))
  ggsave(glue("{dat_dir}/fig_pars_pop_5stage{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=12)

  if(keep_licePreds) {
    mu_df <- out_full_df |>
      filter(grepl("mu", name)) |>
      separate_wider_delim(name, delim=",", names=c("farm", "stage5", "day")) |>
      mutate(farm=str_sub(farm, 4, -1),
             day=as.numeric(str_sub(day, 1, -2)),
             stage=case_when(stage5 %in% 1:2 ~ "Ch",
                             stage5 %in% 3:4 ~ "PA",
                             stage5 == 5 ~ "Ad"),
             stage=factor(stage, levels=c("Ch", "PA", "Ad"))) |>
      group_by(.chain, .iteration, .draw, farm, stage, day) |>
      summarise(value=sum(value)) |>
      group_by(farm, stage, day) |>
      sevcheck::get_intervals(value, type="qi") |>
      ungroup() |>
      mutate(type="Fitted") |> rename(mu=med) |>
      bind_rows(expand_grid(farm=as.character(1:5),
                            day=1:stan_dat$dat$nDays,
                            stage=factor(c("Ch", "PA", "Ad"), levels=c("Ch", "PA", "Ad"))) |>
                  mutate(mu=c(readRDS(glue("{dat_dir}/mu.rds"))[,1,,])) |>
                  mutate(type="True")) |>
      mutate(day=ymd("2023-01-01") + day - 1,
             farm=paste("Farm", farm))
    mu_df |>
      saveRDS(glue("{dat_dir}/mu_sim_fitted_pop_5stage{ifelse(prior_only, '_PRIORS', '')}.rds"))

    p <- mu_df |>
      ggplot(aes(day, mu)) +
      geom_ribbon(aes(ymin=L025, ymax=L975, group=type), colour=NA, fill="grey80") +
      geom_line(aes(colour=type)) +
      scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
      labs(x="Date", y="Mean lice per fish (latent)") +
      scale_x_date(date_labels="%b") +
      facet_grid(stage~farm, scales="free_y")
    ggsave(glue("{dat_dir}/fig_pop_5stage_mu{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=7)

    p <- mu_df |>
      filter(stage=="Ad") |>
      ggplot(aes(day, mu)) +
      geom_ribbon(aes(ymin=L025, ymax=L975, group=type), colour=NA, fill="grey80") +
      geom_line(aes(colour=type)) +
      scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
      labs(x="Date", y="Mean adult female lice per fish (latent)") +
      scale_x_date(date_labels="%b") +
      facet_grid(farm~., scales="free_y")
    ggsave(glue("{dat_dir}/fig_pop_5stage_mu_AF{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=6, height=8)


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
      bind_rows(expand_grid(farm=as.character(1:5),
                            day=1:stan_dat$dat$nDays,
                            stage=factor(c("Ch", "PA", "Ad"), levels=c("Ch", "PA", "Ad"))) |>
                  mutate(y=c(readRDS(glue("{dat_dir}/y.rds"))[,1,,])) |>
                  mutate(type="True") |>
                  inner_join(sampledDays |> select(-sample), by=join_by(farm, day))) |>
      mutate(day=ymd("2023-01-01") + day - 1,
             farm=paste("Farm", farm),
             y_perFish=y/nFishSampled)
    y_df |>
      saveRDS(glue("{dat_dir}/y_sim_fitted_pop_5stage{ifelse(prior_only, '_PRIORS', '')}.rds"))

    p <- y_df |>
      ggplot(aes(day, y_perFish)) +
      geom_linerange(aes(ymin=L025/nFishSampled, ymax=L975/nFishSampled, group=type), linewidth=0.25) +
      geom_point(aes(colour=type, shape=type)) +
      scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
      scale_shape_manual("", values=c("True"=4, "Fitted"=1)) +
      labs(x="Date", y="Mean lice per fish (observed)") +
      scale_x_date(date_labels="%b") +
      facet_grid(stage~farm, scales="free_y")
    ggsave(glue("{dat_dir}/fig_pop_5stage_y{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=10, height=7)

    p <- y_df |>
      filter(stage=="Ad") |>
      ggplot(aes(day, y_perFish)) +
      geom_linerange(aes(ymin=L025/nFishSampled, ymax=L975/nFishSampled, group=type), linewidth=0.25) +
      geom_point(aes(colour=type, shape=type)) +
      scale_colour_manual("", values=c("True"="blue", "Fitted"="black")) +
      scale_shape_manual("", values=c("True"=4, "Fitted"=1)) +
      labs(x="Date", y="Mean adult female lice per fish (observed)") +
      scale_x_date(date_labels="%b") +
      facet_grid(farm~., scales="free_y")
    ggsave(glue("{dat_dir}/fig_pop_5stage_y_AF{ifelse(prior_only, '_PRIORS', '')}.png"), p, width=6, height=8)
  }
}


# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Helper functions


seq_range <- function(x, ...) {
  x_range <- range(x)
  seq(x_range[1], x_range[2], ...)
}

seq_quantiles <- function(x, lo=0.025, hi=0.975, ...) {
  x_q <- quantile(x, probs=c(lo, hi))
  seq(x_q[1], x_q[2], ...)
}



make_compositional <- function(x, method="softmax") {
  switch(method,
         "proportional" = x/sum(x),
         "softmax" = exp(x)/sum(exp(x)),
         paste0("Error: Method must be 'proportional' or 'softmax', but was '", method, "'")
  )
}



pow <- function(x, y) {
  x^y
}

multiply <- function(x, y) {
  x * y
}

runif_minmax <- function(minmax, n) {
  runif(n, minmax[1], minmax[2])
}

qunif_minmax <- function(minmax, p) {
  qunif(p, minmax[1], minmax[2])
}







render_qmd <- function(input_file, output_path, file_ext, ...) {
  # Extract just the input file name (without the file-extension)
  file_name <- xfun::sans_ext(input_file)

  # render the input document and output file will be in the
  # current working directory.
  quarto::quarto_render(input = input_file, output_format = file_ext, ...)

  # name of the rendered output file
  output_name <- paste0(file_name, ".", file_ext)

  # move the file to the output path
  fs::file_move(paste0(output_name), paste0(output_path, "sim_", str_sub(output_path, -3, -2), ".", file_ext))

  msg <- paste0(paste0(output_name, collapse = " and "), " moved to ", output_path)
  message(msg)
}







post_summary_plot <- function(df_ls, ncol=6, nrow=1, scales="fixed") {
  ggplot(data=df_ls[[1]]) +
    geom_pointrange(data=df_ls[[2]],
                    aes(x=mn, xmin=L25, xmax=L75, y=0), linewidth=1.25) +
    geom_linerange(data=df_ls[[2]],
                   aes(xmin=L10, xmax=L90, y=0), linewidth=0.8) +
    geom_linerange(data=df_ls[[2]],
                   aes(xmin=L025, xmax=L975, y=0), linewidth=0.5) +
    geom_density(aes(value)) +
    ylim(0, NA) +
    geom_point(data=df_ls[[3]],
               aes(x=value, y=0), colour="red", shape=1, size=2) +
    facet_wrap(~label, ncol=ncol, scales=scales) +
    theme(axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}




verification_plot <- function(df, ncol=6, nrow=1, scales="fixed") {
  ggplot(df, aes(y=sim)) +
    geom_pointrange(aes(x=mn, xmin=L25, xmax=L75, colour=fit), linewidth=1.25, shape=1) +
    geom_linerange(aes(xmin=L10, xmax=L90, colour=fit), linewidth=0.8) +
    geom_errorbarh(aes(xmin=L025, xmax=L975, colour=fit), linewidth=0.5, height=0.2) +
    geom_point(aes(x=value), colour="red", shape="|", size=4, stroke=2) +
    scale_colour_manual("Fit", values=c("#377eb8", "#4daf4a")) +
    facet_wrap(~label, ncol=ncol, scales=scales) +
    theme(axis.title=element_blank())
}



pr_in_CI_plot <- function(df) {
  # df |>
  #   group_by(fit, sim) |>
  #   summarise(inCI=mean(inCI)) |>
  #   ggplot(aes(x=inCI)) +
  #   geom_histogram(breaks=seq(-0.05, 1.05, by=0.1)) +
  #   facet_grid(fit~.)
  # ggplot(df, aes(y=sim, fill=inCI)) +
  #   geom_bar(position="fill") +
  #   scale_fill_manual(values=c("grey", "green3"), guide="none") +
  #   scale_x_continuous("Proportion of parameters in 95% CIs") +
  #   facet_grid(.~fit) +
  #   theme(axis.title.y=element_blank())
  ggplot(df, aes(y=fit, fill=inCI)) +
    geom_bar(position="fill", colour="grey30") +
    scale_fill_manual(values=c("FALSE"="grey", "TRUE"="green3"), guide="none") +
    scale_x_continuous("Proportion of parameters in 95% CIs") +
    # facet_grid(.~fit) +
    theme(axis.title.y=element_blank())
}


scatter_post_mean_plot <- function(df) {
  lims <- range(c(df$mn, df$value, df$L25, df$L75))
  ggplot(df, aes(mn, value, colour=fit)) +
    geom_abline(linetype=3) +
    geom_point(shape=1, alpha=0.75) +
    geom_linerange(aes(xmin=L25, xmax=L75)) +
    scale_colour_manual("Fit", values=c("#377eb8", "#4daf4a")) +
    facet_grid(.~fit) +
    coord_fixed(xlim=lims, ylim=lims) +
    labs(x="Posterior mean & 50% CI", y="True value")
}




take_mu_draws <- function(out_full_df, f_mu_true=NULL, dat, ndraws=100, GQ=T) {
  draws <- sample.int(max(out_full_df$.draw), ndraws)
  if(GQ) {
    mu_ii <- (1:dat$nDays_GQ)+dat$nDays
  } else {
    mu_ii <- 1:dat$nDays
  }
  mu_true_df <- expand_grid(farm=as.character(1:dat$nFarms),
                            day=1:ifelse(GQ, dat$nDays_GQ, dat$nDays),
                            stage=factor(c("Ch", "PA", "Ad"), levels=c("Ch", "PA", "Ad"))) |>
    mutate(mu=c(readRDS(f_mu_true)[,1,mu_ii,])) |>
    mutate(type="True")
  mu_draws_df <- out_full_df |>
    filter(.draw %in% draws & grepl(ifelse(GQ, "mu_GQ", "mu(?!_GQ)"), name, perl=T)) |>
    separate_wider_delim(name, delim=",", names=c("farm", "stage5", "day")) |>
    mutate(farm=str_sub(farm, ifelse(GQ, 7, 4), -1),
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
    bind_rows(mu_true_df) |>
    mutate(day=ymd("2023-01-01") + ifelse(GQ, dat$nDays, 0) + day - 1,
           farm=paste("Farm", farm))
}



calc_mu_metrics <- function(mu_draws_df, lice_thresh) {
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
  return(list(daily=mu_draw_metrics,
              weekly=mu_draw_metrics_wk))
}




pivot_spp <- function(df, spp_names=c("Lernaeocera_branchialis",
                                      "Lepeophtheirus_salmonis",
                                      "Caligus_elongatus"),
                      spp_levels=c("Nauplii",
                                   "Lepeophtheirus salmonis",
                                   "Caligus elongatus",
                                   "Lernaeocera branchialis")) {
  df |>
    filter(Stage=="Copepodids") |>
    select(-Count) |>
    pivot_longer(any_of(spp_names), names_to="Species", values_to="Count") |>
    drop_na(Count) |>
    mutate(Species_clean=str_replace(Species, "_", " "),
           Species_clean=factor(Species_clean,
                                levels=spp_levels))
}



add_gam_preds <- function(df, mod_fit, suffix, se=T, exclude='t2(Depth,Station_simple,Month)') {
  preds <- predict(mod_fit, newdata=df, se.fit=se, exclude=exclude)
  if(se) {
    pred_df <- tibble(predMn=c(preds$fit),
                      predSE=c(preds$se.fit)) |>
      mutate(predLo1=predMn - 1*predSE,
             predHi1=predMn + 1*predSE,
             predLo2=predMn - 2*predSE,
             predHi2=predMn + 2*predSE) |>
      rename_with(~paste0(.x, suffix, recycle0=T), everything())
  } else {
    pred_df <- tibble(predMn=c(preds)) |>
      rename_with(~paste0(.x, suffix, recycle0=T), everything())
  }
  bind_cols(df, pred_df)
}


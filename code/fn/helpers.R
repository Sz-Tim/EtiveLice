# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Helper functions


seq_range <- function(x, ...) {
  x_min <- min(x)
  x_max <- max(x)
  seq(x_min, x_max, ...)
}



make_compositional <- function(x, method="softmax") {
  switch(method,
         "proportional" = x/sum(x),
         "softmax" = exp(x)/sum(exp(x)),
         paste0("Error: Method must be 'proportional' or 'softmax', but was '", method, "'")
  )
}





make_stan_data <- function(dat_dir, dateRange=NULL, source="sim", priors_only=FALSE,
                           force_detect_priors=FALSE) {
  library(tidyverse)
  library(glue)

  if(source=="sim") {
    info <- readRDS(glue("{dat_dir}info.rds"))
    params <- readRDS(glue("{dat_dir}params.rds"))
    if(!is.null(dateRange)) {
      dates <- seq_range(info$dateRange, by=1) |>
        between(dateRange[1], dateRange[2]) |>
        which()
      info$nDays <- length(dates)
    } else {
      dates <- 1:info$nDays
    }

    nDays <- info$nDays
    nFarms <- info$nFarms
    nSims <- info$nSims
    nStages <- info$nStages
    nAttachCov <- length(params$attach_beta)
    nSurvCov <- nrow(params$surv_beta)

    stan_dat <- list(
      nDays=nDays,
      nFarms=nFarms,
      nSims=nSims,
      nStages=nStages,
      nPens=info$nPens,
      IP_volume=pi*30^2*20,
      y_bar_minimum=1e-5,
      nAttachCov=nAttachCov,
      nSurvCov=nSurvCov,
      y=readRDS(glue("{dat_dir}y.rds"))[,,dates,],
      IP_mx=readRDS(glue("{dat_dir}IP_mx.rds"))[,dates,],
      attach_env_mx=readRDS(glue("{dat_dir}attach_env_mx.rds"))[,dates,],
      sal_mx=readRDS(glue("{dat_dir}sal_mx.rds"))[,dates,],
      temp_mx=readRDS(glue("{dat_dir}temp_mx.rds"))[dates,],
      temp_z_mx=readRDS(glue("{dat_dir}temp_z_mx.rds"))[dates,],
      nFish_mx=readRDS(glue("{dat_dir}nFish_mx.rds"))[dates,],
      sample_i=readRDS(glue("{dat_dir}sampledDays.rds")),
      nFishSampled_mx=readRDS(glue("{dat_dir}nFishSampled_mx.rds"))[dates,],
      y_attach=readRDS(glue("{dat_dir}y_attach.rds")),
      N_attach=readRDS(glue("{dat_dir}N_attach.rds")),
      ensIP=readRDS(glue("{dat_dir}ensIP.rds")),
      mu_true=readRDS(glue("{dat_dir}mu.rds")),
      detect_p=params$detect_p,
      prior_attach_beta=cbind(c(-2.5, rep(0.25, nAttachCov-2), -0.01),
                              c(1, rep(0.25, nAttachCov-2), 0.02)),
      prior_surv_beta=array(c(rep(c(4, rep(0.2, nSurvCov-1)), nStages),
                              rep(c(0.5, rep(0.1, nSurvCov-1)), nStages)),
                            dim=c(nSurvCov, nStages, 2),
                            dimnames=list(c("int", "temp"),
                                          c("Ch", "Pr", "Ad"),
                                          c("mu", "sd"))),
      prior_thresh_GDD_F=cbind(c(130, 325),
                               c(10, 30)),
      prior_thresh_GDD_M=cbind(c(130, 325),
                               c(10, 30)),
      prior_lifespan=c(1500, 20),
      prior_pMolt_F=array(c(qlogis(1/15), 0.5, qlogis(1/20), 0.5,
                            rep(c(0.25, 0.25), 2)),
                          dim=c(2, nStages-1, 2),
                          dimnames=list(c("int", "temp"),
                                        c("Pr", "Ad"),
                                        c("mu", "sd"))),
      prior_pMolt_M=array(c(qlogis(1/15), 0.5, qlogis(1/20), 0.5,
                            rep(c(0.25, 0.25), 2)),
                          dim=c(2, nStages-1, 2),
                          dimnames=list(c("int", "temp"),
                                        c("Pr", "Ad"),
                                        c("mu", "sd"))),
      prior_mnDaysStage_F=array(c(params$mnDaysStageCh[,1], params$mnDaysStagePA[,1],
                                  params$mnDaysStageCh[,2], params$mnDaysStagePA[,2]),
                                dim=c(2, nStages-1, 2),
                                dimnames=list(c("int", "temp"),
                                              c("Pr", "Ad"),
                                              c("mu", "sd"))),
      prior_logit_detect_p=cbind(c(-1, 1),
                                 c(0.5, 0.5))
    )

    if(!is.null(dateRange)) {
      inRng <- stan_dat$sampledDays > min(dates) & stan_dat$sampledDays < max(dates)
      stan_dat$sampledDays <- stan_dat$sample_i[,colSums(inRng)==info$nFarms]
    }
    stan_dat$nSamples <- nrow(stan_dat$sample_i)
    stan_dat$sample_ii <- stan_dat$sample_i |>
      as_tibble() |>
      mutate(index=row_number()) |>
      group_by(sepaSite) |>
      summarise(start=min(index),
                end=max(index)) |>
      select(-sepaSite) |>
      as.matrix()
  }
  stan_dat$sample_prior_only <- as.numeric(priors_only)

  if(force_detect_priors) {
    stan_dat$prior_logit_detect_p <- cbind(qlogis(params$detect_p[1:2]),
                                           c(1e-1, 1e-1))
  }

  return(list(dat=stan_dat, params=params))
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
    geom_bar(position="fill") +
    scale_fill_manual(values=c("grey", "green3"), guide="none") +
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
    labs(x="Posterior mean", y="True value")
}

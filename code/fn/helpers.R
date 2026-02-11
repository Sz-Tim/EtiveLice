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

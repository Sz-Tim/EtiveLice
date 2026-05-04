# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# 2025 field campaign



# setup -------------------------------------------------------------------

library(tidyverse)
library(ggdist)
library(cowplot)
library(brms)
library(plotly)
library(deSolve)
theme_set(ggthemes::theme_few(11) +
            theme(panel.grid.major=element_line(colour="grey95", linewidth=0.1)))

proj_dir <- rprojroot::find_rstudio_root_file()
source(paste0(proj_dir, "/code/fn/helpers.R"))


save_figs <- T


# data munging ------------------------------------------------------------

pump_df <- tibble(depth=c(1,1,1,20,12),
                  lpm=c(386.5335235, 383.8504037, 375.8113731, 320.1323251, 339.4637935)) |>
  mutate(m3pm=lpm/1000)
pump_lm <- lm(m3pm ~ depth, data=pump_df)
calc_m3_pumped <- function(depth_m, time_min, pump_lm) {
  predict(pump_lm, newdata=data.frame(depth=depth_m)) * time_min
}

wc_df <- read_csv(paste0(proj_dir, "/data/westcoms_salinity_profiles.csv")) |>
  mutate(Depth=-depth)
wc_interp <- wc_df |>
  group_by(Station_code, Date_collected, hour) |>
  summarise(sal_approx=list(approxfun(Depth, sal, rule=2)))

spp_df <- read_csv(paste0(proj_dir, "/data/field_data_2025_sppIDs_TEMP.csv")) |>
  select(-Copepodids) |>
  mutate(Stage="Copepodids")

dat_full_df <- read_csv(paste0(proj_dir, "/data/field_data_2025_TEMP.csv")) |>
  mutate(Date_collected=dmy(Date_collected),
         Month=factor(month(Date_collected, label=T), levels=month.abb),
         Start=if_else(is.na(Start), End - duration(30, "minutes"), Start),
         hour=hour(Start),
         Station_simple=fct_collapse(Station_code, `LY1(S)`=c("LY1", "LY1_S")),
         Bout=paste0(Station_code, "_", Date_collected, "_", Start),
         Bout_simple=paste0(Station_simple, "_", Date_collected, "_", Start),
         Salinity=if_else(is.na(Salinity_CTD), Salinity_idronaut, Salinity_CTD)) |>
  pivot_longer(c("Copepodids", "Nauplii"), names_to="Stage", values_to="Count") |>
  mutate(Detected=Count > 0) |>
  drop_na(Count) |>
  group_by(Bout, Stage) |>
  mutate(Any=any(Detected)) |>
  ungroup() |>
  left_join(spp_df, by=join_by(Sample_id, Stage)) |>
  mutate(across(any_of(names(spp_df)[2:4]),
                ~case_when(is.na(.x) & Stage=="Copepodids" & Count==0 ~ 0,
                           .default=.x))) |>
  mutate(m3=calc_m3_pumped(Depth, Duration, pump_lm),
         Station_code=factor(Station_code,
                             levels=c("LY1", "LY1_S", "RE8", "RE10", "RE9", "RE7", "RE6", "RE5")),
         Station_simple=factor(Station_simple,
                               levels=c("LY1(S)", "RE8", "RE10", "RE9", "RE7", "RE6", "RE5")),
         Stage=factor(Stage, levels=c("Nauplii", "Copepodids")),
         Depth_F=factor(Depth, levels=1:40)) |>
  left_join(wc_interp) |>
  mutate(Salinity_WC=NA_real_)

for(i in 1:nrow(dat_full_df)) {
  if(!is.null(dat_full_df$sal_approx[[i]])) {
    dat_full_df$Salinity_WC[i] <- dat_full_df$sal_approx[[i]](dat_full_df$Depth[i])
  }
}

saveRDS(dat_full_df, paste0(proj_dir, "/data/field_data_2025_TEMP_processed.rds"))

dat_main_df <- dat_full_df |>
  filter(Date_collected < "2025-10-22")
dat_12h_df <- dat_full_df |>
  filter(Date_collected == "2025-10-22")

scale_col_depth <- scale_colour_manual("Depth (m)",
                                       values=viridis::viridis(40, direction=-1, end=0.95) |>
                                         set_names(levels(dat_full_df$Depth_F)))
scale_fill_depth <- scale_fill_manual("Depth (m)",
                                      values=viridis::viridis(40, direction=-1, end=0.95) |>
                                        set_names(levels(dat_full_df$Depth_F)))
spp_levels <- c("Nauplii",
                "Lepeophtheirus salmonis",
                "Caligus elongatus",
                "Lernaeocera branchialis")



# Fig: salinity x lice_m3 -------------------------------------------------

dat_full_df |>
  filter(!is.na(Salinity)) |>
  ggplot(aes(Salinity, Count/m3, colour=Depth_F)) +
  geom_point(size=2, shape=1, stroke=1) +
  scale_col_depth +
  labs(x="Salinity (psu)", y=expression("Lice "%.%" m"^-3)) +
  facet_wrap(~Stage) +
  theme(legend.position="inside",
        legend.position.inside=c(0.07, 0.67))
if(save_figs) {
  ggsave(paste0(proj_dir, "/figs/field_ms/lice_salinity_Naup_Cop.png"), width=7, height=3.5)
}

dat_full_df |>
  filter(!is.na(Salinity)) |>
  pivot_spp() |>
  bind_rows(dat_full_df |> filter(Stage=="Nauplii") |> mutate(Species_clean="Nauplii")) |>
  mutate(Species_clean=factor(Species_clean, levels=spp_levels)) |>
  ggplot(aes(Salinity, Count/m3, colour=Depth_F)) +
  geom_point(size=2, shape=1, stroke=1) +
  scale_col_depth +
  labs(x="Salinity (psu)", y=expression("Lice "%.%" m"^-3)) +
  facet_wrap(~Species_clean, nrow=1) +
  theme(legend.position="inside",
        legend.position.inside=c(0.05, 0.8),
        legend.key.height=unit(0.12, "cm"),
        legend.background=element_blank())
if(save_figs) {
  ggsave(paste0(proj_dir, "/figs/field_ms/lice_salinity_Spp.png"), width=10.5, height=3.5)
}




# Fig: 12h samples --------------------------------------------------------

dat_12h_df |>
  pivot_spp() |>
  bind_rows(dat_12h_df |> filter(Stage=="Nauplii") |> mutate(Species_clean="Nauplii")) |>
  mutate(Species_clean=factor(Species_clean, levels=spp_levels)) |>
  mutate(hour_dec=hour + minute(Start)/60) |>
  ggplot(aes(hour_dec, Depth, fill=Count/m3)) +
  # geom_tile(height=2, width=0.6, colour=NA, position=position_nudge(x=0.3)) +
  geom_raster(position=position_nudge(x=0.3)) +
  scale_fill_viridis_c(expression("Lice "%.%" m"^-3), option="plasma") +
  scale_y_reverse("Depth (m)", limits=c(21, 0), breaks=c(1, 12, 20), oob=scales::oob_keep) +
  scale_x_continuous("Hour", minor_breaks=7:18, breaks=7:18) +
  facet_wrap(~Species_clean, nrow=1) +
  theme(legend.position="bottom",
        legend.title.position="top",
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1, "cm"),
        panel.grid.minor.x=element_line(colour="grey95", linewidth=0.1))
if(save_figs) {
  ggsave(paste0(proj_dir, "/figs/field_ms/lice_12h.png"), width=10.5, height=4)
}



# Models: Pres ~ Depth ----------------------------------------------------

n_df <- dat_main_df |>
  filter(Stage=="Nauplii") |>
  mutate(Detected=as.numeric(Detected))
Lb_df <- dat_main_df |>
  pivot_spp() |>
  filter(Species=="Lernaeocera_branchialis") |>
  mutate(Detected=as.numeric(Count > 0))
Ls_df <- dat_main_df |>
  pivot_spp() |>
  filter(Species=="Lepeophtheirus_salmonis") |>
  mutate(Detected=as.numeric(Count > 0))
Ce_df <- dat_main_df |>
  pivot_spp() |>
  filter(Species=="Caligus_elongatus") |>
  mutate(Detected=as.numeric(Count > 0))

# Depth
library(mgcv)
det_formula <- formula(Detected ~ s(Depth, bs="tp", k=5, m=2) +
                         t2(Depth, Station_simple, Month, k=5, bs=c("tp", "re", "re"), m=2, full=TRUE) +
                         offset(log(m3)))
pDet_out_n <- gam(det_formula, data=n_df, family="binomial", method="REML")
pDet_out_Lb <- gam(det_formula, data=Lb_df, family="binomial", method="REML")
pDet_out_Ls <- gam(det_formula, data=Ls_df, family="binomial", method="REML")
pDet_out_Ce <- gam(det_formula, data=Ce_df, family="binomial", method="REML")



# . 0-20m -----------------------------------------------------------------
pred_df_FE <- expand_grid(Depth=seq(0, 20, by=0.1),
                       Month=unique(Ls_df$Month)[1],
                       Station_simple=levels(Ls_df$Station_simple)[1],
                       m3=1) |>
  add_gam_preds(pDet_out_n, "_Naup") |>
  add_gam_preds(pDet_out_Ls, "_Ls") |>
  add_gam_preds(pDet_out_Ce, "_Ce") |>
  add_gam_preds(pDet_out_Lb, "_Lb") |>
  mutate(across(starts_with("pred"), boot::inv.logit)) |>
  select(-starts_with("predSE")) |>
  pivot_longer(starts_with("pred"),
               names_to=c("metric", "model"),
               names_pattern="pred(.*)_(.*)",
               values_to="pred") |>
  pivot_wider(names_from=metric, values_from=pred) |>
  mutate(model=factor(model, levels=c("Naup", "Ls", "Ce", "Lb"),
                      labels=spp_levels))
pred_df_RE <- expand_grid(Depth=seq(0, 20, by=0.1),
                          Month=unique(Ls_df$Month),
                          Station_simple=levels(Ls_df$Station_simple),
                          m3=1) |>
  add_gam_preds(pDet_out_n, "_Naup", se=F, exclude=NULL) |>
  add_gam_preds(pDet_out_Ls, "_Ls", se=F, exclude=NULL) |>
  add_gam_preds(pDet_out_Ce, "_Ce", se=F, exclude=NULL) |>
  add_gam_preds(pDet_out_Lb, "_Lb", se=F, exclude=NULL) |>
  mutate(across(starts_with("pred"), boot::inv.logit)) |>
  pivot_longer(starts_with("pred"),
               names_to=c("metric", "model"),
               names_pattern="pred(.*)_(.*)",
               values_to="pred") |>
  pivot_wider(names_from=metric, values_from=pred) |>
  mutate(model=factor(model, levels=c("Ls", "Ce", "Lb", "Naup"),
                      labels=levels(pred_df_FE$model)))
x_adj <- rep(0:9, times=4)*0.00375
y_adj <- rep(1:4, each=10)*0.275
obs_df <- bind_rows(n_df |> mutate(model="Naup"),
                    Ls_df |> mutate(model="Ls"),
                    Ce_df |> mutate(model="Ce"),
                    Lb_df |> mutate(model="Lb")) |>
  select(Station_simple, Month, Depth, model, Detected) |>
  mutate(model=factor(model, levels=c("Naup", "Ls", "Ce", "Lb"),
                      labels=levels(pred_df_FE$model))) |>
  group_by(model, Depth, Detected) |>
  mutate(row_id=row_number(),
         x_adj_i=x_adj[row_id],
         y_adj_i=y_adj[row_id],
         y_pos=Depth + (y_adj_i - mean(unique(y_adj_i))),
         x_pos=Detected*0.225 + x_adj_i*if_else(Detected==1, 1, -1) - 0.0025) |>
  ungroup()
obs_counts <- obs_df |>
  count(model, Depth, Detected)
axis_lines <- tibble(y1=unique(obs_df$Depth),
                     y2=unique(obs_df$Depth),
                     x1=0,
                     x2=0.22) |>
  mutate(grp=row_number())

pred_df_FE |>
  ggplot() +
  geom_segment(data=axis_lines, aes(x1, y1, xend=x2, yend=y2, group=grp),
               colour="grey85", linewidth=0.1) +
  geom_vline(xintercept=c(0, 0.22), colour="grey30", linewidth=0.25) +
  geom_hline(yintercept=0, colour="grey30", linewidth=0.9) +
  geom_point(data=obs_df |> filter(Depth <= 20),
             aes(x_pos, y_pos, colour=factor(Detected)), shape=1, size=0.5) +
  scale_colour_manual(values=c("grey30", "red3"), guide="none") +
  geom_ribbon(aes(Mn, xmin=Lo2, xmax=Hi2, Depth), alpha=0.2, colour=NA) +
  geom_ribbon(aes(Mn, xmin=Lo1, xmax=Hi1, Depth), alpha=0.2, colour=NA) +
  geom_path(aes(Mn, Depth), linewidth=0.5) +
  scale_y_reverse("Depth (m)", breaks=c(1, 3, 12, 20),
                  limits=c(0, 21),
                  expand=c(0, 0),
                  oob=scales::oob_keep) +
  scale_x_continuous(expression("Pr(presence) per m"^3),
                     breaks=seq(0, 0.2, by=0.05),
                     labels=scales::label_percent()
                     ) +
  theme(legend.key.spacing.y=unit(0.2, "cm"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="grey85", linewidth=0.1),
        panel.spacing.x=unit(0.7, "cm"),
        panel.border=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.25),
        axis.line.x.bottom=element_line(),
        axis.line.y=element_line(),
        axis.line.x.top=element_line()) +
  facet_wrap(~model, nrow=1)

if(save_figs) {
  ggsave(paste0(proj_dir, "/figs/field_ms/prDet_m3.png"), width=10, height=3.5)
}


x_adj <- rep(0:9, times=4)
x_scale <- c(`Nauplii`=0.001,
             `Lepeoptheirus salmonis`=0.001,
             `Caligus elongatus`=0.0005,
             `Lernaeocera branchialis`=0.00375)
x_lim <- c(`Nauplii`=0.046,
           `Lepeophtheirus salmonis`=0.07,
           `Caligus elongatus`=0.022,
           `Lernaeocera branchialis`=0.225)
y_adj <- rep(1:4, each=10)*0.275
obs_df <- bind_rows(n_df |> mutate(model="Naup"),
                    Ls_df |> mutate(model="Ls"),
                    Ce_df |> mutate(model="Ce"),
                    Lb_df |> mutate(model="Lb")) |>
  select(Station_simple, Month, Depth, model, Detected) |>
  mutate(model=factor(model, levels=c("Naup", "Ls", "Ce", "Lb"),
                      labels=levels(pred_df_FE$model))) |>
  group_by(model, Depth, Detected) |>
  mutate(row_id=row_number(),
         x_adj_i=x_adj[row_id]*x_scale[model],
         y_adj_i=y_adj[row_id],
         y_pos=Depth + (y_adj_i - mean(unique(y_adj_i))),
         x_pos=Detected*x_lim[model] + x_adj_i*if_else(Detected==1, 1, -1) +
           x_scale[model]*if_else(Detected==1, 1, -1)*0.5) |>
  ungroup() |>
  filter(Depth <= 20)
axis_lines <- tibble(model=names(x_lim),
                     x2=x_lim) |>
  mutate(dat=list(tibble(y1=unique(obs_df$Depth),
                         y2=unique(obs_df$Depth),
                         x1=0))) |>
  unnest(dat) |>
  mutate(grp=row_number(),
         model=factor(model, levels=levels(pred_df_FE$model)))
dummy_lims <- expand_grid(Depth=range(pred_df_FE$Depth),
                          model=levels(pred_df_FE$model),
                          row_id=1:40,
                          Detected=0:1) |>
  mutate(x_adj_i=x_adj[row_id]*x_scale[model],
         y_adj_i=y_adj[row_id],
         y_pos=Depth + (y_adj_i - mean(unique(y_adj_i))),
         x_pos=Detected*x_lim[model] + x_adj_i*if_else(Detected==1, 1, -1) +
           x_scale[model]*if_else(Detected==1, 1, -1)*0.5) |>
  mutate(model=factor(model, levels=levels(pred_df_FE$model)))

pred_df_FE |>
  ggplot() +
  geom_point(data=dummy_lims, aes(x_pos, y_pos), colour="white") +
  geom_segment(data=axis_lines, aes(x1, y1, xend=x2, yend=y2, group=grp),
               colour="grey85", linewidth=0.1) +
  geom_vline(xintercept=c(0), colour="grey30", linewidth=0.25) +
  geom_vline(data=axis_lines, aes(xintercept=x2), colour="grey30", linewidth=0.25) +
  geom_hline(yintercept=0, colour="grey30", linewidth=0.9) +
  geom_point(data=obs_df,
             aes(x_pos, y_pos, colour=factor(Detected)), shape=1, size=0.5) +
  scale_colour_manual(values=c("grey30", "red3"), guide="none") +
  geom_ribbon(aes(Mn, xmin=Lo2, xmax=Hi2, Depth), alpha=0.2, colour=NA) +
  geom_ribbon(aes(Mn, xmin=Lo1, xmax=Hi1, Depth), alpha=0.2, colour=NA) +
  geom_path(aes(Mn, Depth), linewidth=0.5) +
  scale_y_reverse("Depth (m)", breaks=c(1, 3, 12, 20),
                  limits=c(0, 21),
                  expand=c(0, 0),
                  oob=scales::oob_keep) +
  scale_x_continuous(expression("Pr(presence) per m"^3),
                     labels=scales::label_percent(accuracy=1)
  ) +
  theme(legend.key.spacing.y=unit(0.2, "cm"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.spacing.x=unit(0.7, "cm"),
        panel.border=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.25),
        axis.line.x.bottom=element_line(),
        axis.line.y=element_line(),
        axis.line.x.top=element_line()) +
  facet_wrap(~model, nrow=1, scales="free_x")

if(save_figs) {
  ggsave(paste0(proj_dir, "/figs/field_ms/prDet_m3_free.png"), width=10, height=3.5)
}



# . 0-40m -----------------------------------------------------------------

pred_df_FE <- expand_grid(Depth=seq(0, 40, by=0.1),
                          Month=unique(Ls_df$Month)[1],
                          Station_simple=levels(Ls_df$Station_simple)[1],
                          m3=1) |>
  add_gam_preds(pDet_out_n, "_Naup") |>
  add_gam_preds(pDet_out_Ls, "_Ls") |>
  add_gam_preds(pDet_out_Ce, "_Ce") |>
  add_gam_preds(pDet_out_Lb, "_Lb") |>
  mutate(across(starts_with("pred"), boot::inv.logit)) |>
  select(-starts_with("predSE")) |>
  pivot_longer(starts_with("pred"),
               names_to=c("metric", "model"),
               names_pattern="pred(.*)_(.*)",
               values_to="pred") |>
  pivot_wider(names_from=metric, values_from=pred) |>
  mutate(model=factor(model, levels=c("Naup", "Ls", "Ce", "Lb"),
                      labels=spp_levels))
pred_df_RE <- expand_grid(Depth=seq(0, 40, by=0.1),
                          Month=unique(Ls_df$Month),
                          Station_simple=levels(Ls_df$Station_simple),
                          m3=1) |>
  add_gam_preds(pDet_out_n, "_Naup", se=F, exclude=NULL) |>
  add_gam_preds(pDet_out_Ls, "_Ls", se=F, exclude=NULL) |>
  add_gam_preds(pDet_out_Ce, "_Ce", se=F, exclude=NULL) |>
  add_gam_preds(pDet_out_Lb, "_Lb", se=F, exclude=NULL) |>
  mutate(across(starts_with("pred"), boot::inv.logit)) |>
  pivot_longer(starts_with("pred"),
               names_to=c("metric", "model"),
               names_pattern="pred(.*)_(.*)",
               values_to="pred") |>
  pivot_wider(names_from=metric, values_from=pred) |>
  mutate(model=factor(model, levels=c("Naup", "Ls", "Ce", "Lb"),
                      labels=levels(pred_df_FE$model)))
x_adj <- rep(0:9, times=4)*0.01
y_adj <- rep(1:4, each=10)*0.4
obs_df <- bind_rows(n_df |> mutate(model="Naup"),
                    Ls_df |> mutate(model="Ls"),
                    Ce_df |> mutate(model="Ce"),
                    Lb_df |> mutate(model="Lb")) |>
  select(Station_simple, Month, Depth, model, Detected) |>
  mutate(model=factor(model, levels=c("Naup", "Ls", "Ce", "Lb"),
                      labels=levels(pred_df_FE$model))) |>
  group_by(model, Depth, Detected) |>
  mutate(row_id=row_number(),
         x_adj_i=x_adj[row_id],
         y_adj_i=y_adj[row_id],
         y_pos=Depth + (y_adj_i - mean(unique(y_adj_i))),
         x_pos=Detected*0.41 + x_adj_i*if_else(Detected==1, 1, -1) - 0.005) |>
  ungroup()
obs_counts <- obs_df |>
  count(model, Depth, Detected)
axis_lines <- tibble(y1=unique(obs_df$Depth),
                     y2=unique(obs_df$Depth),
                     x1=0,
                     x2=0.4) |>
  mutate(grp=row_number())

pred_df_FE |>
  mutate(across(all_of(c("Mn", "Hi1", "Hi2")), ~pmin(.x, 0.4))) |>
  ggplot() +
  geom_segment(data=axis_lines, aes(x1, y1, xend=x2, yend=y2, group=grp),
               colour="grey85", linewidth=0.1) +
  geom_vline(xintercept=c(0, 0.4), colour="grey30", linewidth=0.25) +
  geom_hline(yintercept=0, colour="grey30", linewidth=0.9) +
  geom_point(data=obs_df, aes(x_pos, y_pos, colour=factor(Detected)), shape=1, size=0.5) +
  scale_colour_manual(values=c("grey", "red3"), guide="none") +
  geom_ribbon(aes(Mn, xmin=Lo2, xmax=Hi2, Depth), alpha=0.2, colour=NA) +
  geom_ribbon(aes(Mn, xmin=Lo1, xmax=Hi1, Depth), alpha=0.2, colour=NA) +
  geom_path(aes(Mn, Depth), linewidth=0.5) +
  scale_y_reverse("Depth (m)", breaks=c(1, 3, 12, 20, 30, 40),
                  limits=c(0, 41),
                  expand=c(0, 0),
                  oob=scales::oob_keep) +
  scale_x_continuous(expression("Pr(presence) per m"^3),
                     breaks=seq(0, 0.4, by=0.1),
                     labels=scales::label_percent()) +
  theme(legend.key.spacing.y=unit(0.2, "cm"),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="grey85", linewidth=0.1),
        panel.spacing.x=unit(0.7, "cm"),
        panel.border=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.25),
        axis.line.x.bottom=element_line(),
        axis.line.y=element_line(),
        axis.line.x.top=element_line()) +
  facet_wrap(~model, nrow=1)

if(save_figs) {
  ggsave(paste0(proj_dir, "/figs/field_ms/prDet_m3_40m.png"), width=10, height=3.5)
}

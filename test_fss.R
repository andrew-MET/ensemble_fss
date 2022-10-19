library(harp)
library(dplyr)
library(meteogrid)
library(here)
library(purrr)
library(forcats) # only used for plotting

source(here("binary_prob.R"))
source(here("nbhd_upscale.R"))
source(here("fss.R"))
source(here("ens_read_and_fbs.R"))
source(here("dfss.R"))

# FSS settings
accum_hours   <- 1
thresholds    <- c(0.9, 0.95, 0.99) #c(0.1, 1)
quantile_thresh <- TRUE
nbhd_radius   <- seq(0, 100, 2)
groupings     <- "lead_time"

# fcst file settings
fcst_date_times <- 2022081500#seq_dates(2019020100, 2019020200, "1d")
fcst_lead_times <- seq(3, 48, 3)
fcst_param      <- "Pcp"
fcst_model      <- "meps"
fcst_dir        <- ""
fcst_template   <- "meps_lagged_6h_subset"
fcst_fmt_opts   <- netcdf_opts("met_norway_eps")
fcst_members    <- seq(1, 6)

# obs file settings
obs_param     <- "precipitation_amount"
obs_dir       <- "/lustre/storeB/immutable/archive/projects/metproduction/yr_short"
obs_template  <- "{YYYY}/{MM}/{DD}/met_analysis_1_0km_nordic_{YYYY}{MM}{DD}T{HH}Z.nc"
obs_fmt_opts  <- netcdf_opts(proj4_var = "projection_lcc")

# Get the verification domain
member_col <- paste0(fcst_model[1], "_mbr000")
dom <- read_forecast(
  date_times       = fcst_date_times[1],
  fcst_model       = fcst_model[1],
  parameter        = fcst_param,
  members          = 0,
  lead_time        = 0,
  file_path        = fcst_dir,
  file_template    = fcst_template,
  file_format_opts = fcst_fmt_opts,
  return_data      = TRUE
)  %>%
  pluck(1) %>%
  pull({{member_col}}) %>%
  pluck(1) %>%
  as.geodomain() %>%
  subgrid(72, 750, 80, 900)

# Read data and compute unsummarized Fractions Brier Score and ref
arg_df  <- expand.grid(
  fcst_lead_times = fcst_lead_times,
  fcst_date_times = fcst_date_times
)

fbs_all <- purrr::map2_dfr(
  arg_df$fcst_date_times,
  arg_df$fcst_lead_times,
  ens_read_and_fbs,
  fcst_model,
  fcst_param,
  fcst_members,
  fcst_dir,
  fcst_template,
  fcst_fmt_opts,
  accum_hours,
  obs_param,
  obs_dir,
  obs_template,
  obs_fmt_opts,
  dom,
  thresholds,
  nbhd_radius,
  num_cores = length(nbhd_radius),
  quantile_thresholds = quantile_thresh
)

# Compute the FSS - always group by fcst_model, threshold and nbhd_length
# and user supplied group
fixed_groups <- c("fcst_model", "threshold", "nbhd_length")
thresh_col <- "threshold"
thresh_label <- "Threshold [mm]"
forcats_fun <- fct_rev
if (quantile_thresh) {
  fixed_groups <- c("fcst_model", "quantile", "nbhd_length")
  thresh_col <- "quantile"
  thresh_label <- "Quantile"
  forcats_fun <- fct_inorder
}

groupings <- union(
  groupings,
  fixed_groups
)

fss_df <- fbs_all %>%
  dplyr::group_by(!!!rlang::syms(groupings)) %>%
  dplyr::summarize(fss = 1 - (sum(fbs) / sum(fbs_ref)))

# Plot FSS as a function of lead time - line colour is nbhd_length in km,
# and one panel for each threshold
ggplot(fss_df, aes(x = lead_time, y = fss, colour = fct_inorder(factor(nbhd_length * 2.5)))) +
  geom_line() +
  facet_grid(rows = vars(.data[[thresh_col]]), cols = vars(fcst_model)) +
  labs(
    x      = "Lead Time [h]",
    y      = "Fractions Skill Score",
    colour = "Neighbourhood\nlength [km]"
  )

# Plot FSS as a heat map - one panel for each lead time
ggplot(fss_df, aes(x = fct_inorder(factor(nbhd_length * 2.5)), y = forcats_fun(factor(.data[[thresh_col]])), fill = fss)) +
  geom_raster() +
  scale_fill_gradientn(colours = c("steelblue4", "white", "indianred4"), limits = c(0, 1)) +
  facet_wrap(vars(lead_time)) +
  #facet_grid(rows = vars(lead_time), cols = vars(fcst_model)) +
  labs(
    x    = "Neighbourhood Length [km]",
    y    = thresh_label,
    fill = "Fractions\nSkill Score"
  ) +
  coord_fixed(expand = FALSE) #+
  #geom_text(aes(label = format(round(fss, digits = 2), nsmall = 2)))

# Plot dfss and efss
breaks <- seq(0, 1, 0.05)
break_labels <- paste(
  sprintf(breaks[1:(length(breaks) - 1)], fmt = "%#.2f"),
  sprintf(breaks[2:length(breaks)], fmt = "%#.2f"),
  sep = "-"
)
cols <- colorRampPalette(c("steelblue4", "white", "indianred4"))(length(break_labels))
middle_contour <- which.min(abs(breaks - 0.5)) - 1
contour_cols <- rep("grey70", length(breaks))
contour_cols[middle_contour] <- "grey20"

ggplot(
  tidyr::pivot_longer(fbs_all, c(efss_mean, dfss_mean)),
  aes(y = nbhd_length * 2.5, x = lead_time, z = value)
) +
  geom_contour_filled(breaks = breaks) +
  scale_fill_manual(values = cols, labels = break_labels) +
  geom_contour(
    aes(colour = factor(after_stat(level))),
    breaks = breaks
  ) +
  scale_colour_manual(
    values = contour_cols,
    guide = "none"
  ) +
  labs(
    y = "Neighbourhood Length [km]",
    x = "Lead Time [h]",
    fill = "FSS"
  ) +
  coord_cartesian(expand = FALSE) +
  facet_wrap(vars(name)) +
  scale_x_continuous(breaks = seq(0, 180, 6)) +
  theme_bw()



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
thresholds    <- c(0.9, 0.95, 0.99, 0.999) #c(0.1, 1)
quantile_thresh <- TRUE
nbhd_radius   <- seq(0, 100, 2)
groupings     <- "lead_time"

# fcst file settings
fcst_date_times <- 2022110600#seq_dates(2019020100, 2019020200, "1d")
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

# Parallel settings
num_cores <- 20
if (num_cores == "auto") {
  num_cores <- length(nbhd_radius)
}

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

start_time <- Sys.time()

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
  num_cores = num_cores,
  quantile_thresholds = quantile_thresh
)

print(Sys.time() - start_time)

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
  fbs_all$quantile <- fct_inorder(fbs_all$quantile)
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
  ) +
  scale_x_continuous(breaks = seq(0, 180, 6)) +
  theme_bw()

# Plot FSS as a heat map - one panel for each lead time
ggplot(fss_df, aes(x = nbhd_length * 2.5, y = .data[[thresh_col]], fill = fss)) +
  geom_raster() +
  scale_fill_gradientn(colours = c("steelblue4", "white", "darkolivegreen4"), limits = c(0, 1)) +
  facet_wrap(vars(fct_inorder(paste0("Lead Time = ",lead_time, "h")))) +
  #facet_grid(rows = vars(lead_time), cols = vars(fcst_model)) +
  labs(
    x    = "Neighbourhood Length [km]",
    y    = thresh_label,
    fill = "Fractions\nSkill Score"
  ) +
  coord_cartesian(expand = FALSE) #+
  #geom_text(aes(label = format(round(fss, digits = 2), nsmall = 2)))

# Plot dfss and efss
breaks <- seq(0, 1, 0.05)
break_labels <- paste(
  sprintf(breaks[1:(length(breaks) - 1)], fmt = "%#.2f"),
  sprintf(breaks[2:length(breaks)], fmt = "%#.2f"),
  sep = "-"
)
cols <- colorRampPalette(c("steelblue4", "white", "darkolivegreen4"))(length(break_labels))

ggplot(
  tidyr::pivot_longer(fbs_all, c(efss_mean, dfss_mean)),
  aes(y = nbhd_length * 2.5, x = lead_time, z = value)
) +
  geom_contour_filled(breaks = breaks) +
  scale_fill_manual(values = cols, labels = break_labels, drop = FALSE) +
  geom_contour(breaks = breaks, colour = "grey70") +
  geom_contour(breaks = 0.5, colour =" grey20") +
  labs(
    y = "Neighbourhood Length [km]",
    x = "Lead Time [h]",
    fill = "Fractions\nSkill Score"
  ) +
  coord_cartesian(expand = FALSE) +
  facet_grid(cols = vars(name), rows = vars(fct_inorder(.data[[thresh_col]]))) +
  scale_x_continuous(breaks = seq(0, 180, 6)) +
  theme_bw()

# Extract the FSS = 0.5 contour and plot to compare dfss and efss

skillful_line <- function(lt, nb, z) {

  x   <- unique(lt)
  y   <- unique(nb)
  z   <- t(matrix(z, ncol = length(x), nrow = length(y)))

  res <- contourLines(x, y, z, levels = 0.5)
  res <- lapply(
    seq_along(res),
    function(i) dplyr::mutate(data.frame(res[[i]]), group = i)
  )

  res <- dplyr::bind_rows(res)

  if (nrow(res) == 0) return(NULL)

  dplyr::transmute(
    res,
    lead_time = x, nbhd_length = y, group = group
  )

}

skillful_scale <- dplyr::group_by(
  tidyr::pivot_longer(fbs_all, c(efss_mean, dfss_mean)),
  .data[[thresh_col]], .data[["name"]]
) %>%
  dplyr::summarise(
    skill_line = list(
      skillful_line(
        .data[["lead_time"]], .data[["nbhd_length"]], .data[["value"]]
      )
    )
  ) %>%
  tidyr::unnest(skill_line)

ggplot(
  skillful_scale,
  aes(lead_time, nbhd_length * 2.5, colour = name, group = paste(name, group))
) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(.data[[thresh_col]]), ncol = 1) +
  scale_colour_manual(NULL, values = c(dfss_mean = "steelblue", efss_mean = "darkorange3")) +
  scale_x_continuous(breaks = seq(0, 180, 6)) +
  labs(
    x = "Lead Time [h]",
    y = "Neighbourhood Length [km]"
  ) +
  theme_bw()

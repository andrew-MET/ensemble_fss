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

# FSS settings
accum_hours   <- 3
thresholds    <- c(0.1, 1, 2, 4)
nbhd_radius   <- seq(0, 10)
groupings     <- "lead_time"

# fcst file settings
fcst_date_times <- seq_dates(2019020100, 2019020200, "1d")
fcst_lead_times <- seq(3, 6, 3)
fcst_param      <- "Pcp"
fcst_model      <- "meps"
fcst_dir        <- ""
fcst_template   <- "meps_subset"
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
  subgrid(52, 750, 30, 900)

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
  nbhd_radius
)

# Compute the FSS - always group by fcst_model, threshold and nbhd_length
# and user supplied group
groupings <- union(
  groupings,
  c("fcst_model", "threshold", "nbhd_length")
)

fss_df <- fbs_all %>%
  dplyr::group_by(!!!rlang::syms(groupings)) %>%
  dplyr::summarize(fss = 1 - (sum(fbs) / sum(fbs_ref)))

# Plot FSS as a function of lead time - line colour is nbhd_length in km,
# and one panel for each threshold
ggplot(fss_df, aes(x = lead_time, y = fss, colour = fct_inorder(factor(nbhd_length * 2.5)))) +
  geom_line() +
  facet_grid(rows = vars(threshold), cols = vars(fcst_model)) +
  labs(
    x      = "Lead Time [h]",
    y      = "Fractions Skill Score",
    colour = "Neighbourhood\nlength [km]"
  )

# Plot FSS as a heat map - one panel for each lead time
ggplot(fss_df, aes(x = fct_inorder(factor(nbhd_length * 2.5)), y = fct_rev(factor(threshold)), fill = fss)) +
  geom_raster() +
  geom_text(aes(label = format(round(fss, digits = 2), nsmall = 2))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  facet_wrap(vars(lead_time)) +
  #facet_grid(rows = vars(lead_time), cols = vars(fcst_model)) +
  labs(
    x    = "Neighbourhood Length [km]",
    y    = "Threshold [mm]",
    fill = "Fractions\nSkill Score"
  ) +
  coord_fixed(expand = FALSE)






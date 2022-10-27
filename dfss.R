# Some tests for spatial spread using dFSS
# The scale at which FSS begins to exceed 0.5 for each unique pair of members

#library(harpMET)
library(harpIO)
library(harpVis)
library(dplyr)
library(here)

source(here("binary_prob.R"))
source(here("nbhd_upscale.R"))
source(here("fss.R"))

# Read some random members precipitation forecasts plus the control
members <- c(0, sample(seq(1, 29), 6))
# fcst <- read_meps(
#   2022081518, "Pcp", seq(6, 36), members = members
# ) %>%
#   harpIO::accumulate(1)
#
# # Plot for a single lead time
# map_path <- get_map(dom = get_domain(fcst, meps_mbr000), poly = FALSE)
# ggplot(filter(harpPoint::gather_members(fcst)$meps, lead_time == 18)) +
#   geom_georaster(aes(geofield = forecast)) +
#   geom_path(aes(x, y), map_path, colour = "seagreen3") +
#   facet_wrap(vars(member), nrow = 2) +
#   scale_fill_viridis_c(
#     NULL,
#     option = "C",
#     limits = c(1, NA),
#     na.value = "transparent",
#     trans = "log",
#     breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64),
#     direction = 1
#   ) +
#   coord_equal(expand = FALSE) +
#   theme_harp_map() +
#   theme(panel.background = element_rect(fill = "grey30"))
#
#

# For each member combination at each lead time find the
# scale at which the fss goes above 0.5

# Function to find the radius for a pair of members
fss_pair <- function(mbr1, mbr2, threshold = 0.1, test_radii = seq(0, 10)) {

  for (radius in test_radii) {
    pair_fss <- fss(mbr1, mbr2, thresh = threshold, radius = radius)
    if (pair_fss >= 0.5) break
  }

  cover1 <- paste0(sum(binary_prob(mbr1, threshold)),  "/",  length(mbr1))
  cover2 <- paste0(sum(binary_prob(mbr2, threshold)),  "/",  length(mbr2))

  list(rad = radius, cvr1 = cover1, cvr2 = cover2, fss = pair_fss)
}

fss_pair_prob <- function(mbr1, mbr2, radius) {
  mbr1 <- nbhd_upscale(mbr1, radius)
  mbr2 <- nbhd_upscale(mbr2, radius)
  data.frame(fbs = fbs(mbr1, mbr2), fbs_ref = fbs_ref(mbr1, mbr2))
}

# Function to create unique pairs of members
unique_pairs <- function(x) {
  x <- unique(x)
  g <- function(i) {
    z <- setdiff(x, x[seq_len(i)])
    if (length(z)) {
      cbind(x[i], z, deparse.level = 0)
    }
  }
  do.call(rbind, lapply(seq_along(x), g))
}

# Function to get dFSS by comparing all members against each other
dfss_row <- function(fcst_row, threshold = 0.1, test_radii = seq(0, 10), truncate = FALSE, num_cores = 1) {

  # Takes 1 data frame with one row - looping over rows is done by the calling function

  # Generate all of the member combinations
  members <- grep("_mbr[[:digit:]]", colnames(fcst_row), value = TRUE)
  member_pairs <- unique_pairs(members)
  pair_fss <- function(i) {
    res <- fss_pair(
      fcst_row[[member_pairs[i, 1]]][[1]],
      fcst_row[[member_pairs[i, 2]]][[1]],
      threshold  = threshold,
      test_radii = test_radii
    )
    tibble::tibble(
      fcdate     = fcst_row[["fcdate"]][1],
      lead_time  = fcst_row[["lead_time"]][1],
      mbrA       = member_pairs[i, 1],
      mbrB       = member_pairs[i, 2],
      mbrA_cover = res[["cvr1"]],
      mbrB_cover = res[["cvr2"]],
      fss        = res[["fss"]],
      scale      = res[["rad"]] * 2 + 1,
      truncated  = res[["fss"]] < 0.5
    )
  }

  pair_fss_prob <- function(i) {
    res <- fss_pair_prob(
      fcst_row[[member_pairs[i, 1]]][[1]],
      fcst_row[[member_pairs[i, 2]]][[1]],
      test_radii
    )
  }

  if (num_cores > 1) {
    available_cores <- parallel::detectCores()
    if (num_cores > available_cores) {
      warning(num_cores, "requested, but only", available_cores, "found.")
      num_cores <- available_cores
    }
  }

  if (num_cores > 1) {
    res <- do.call(
      rbind,
      parallel::mclapply(1:nrow(member_pairs), pair_fss_prob, mc.cores = num_cores)
    )
  } else {
    res <- do.call(rbind, lapply(1:nrow(member_pairs), pair_fss_prob))
  }

  data.frame(
    dfss_mean = 1 - (sum(res[["fbs"]]) / sum(res[["fbs_ref"]])),
    dfss_sd   = sd(1 - (res[["fbs"]] / res[["fbs_ref"]]))
  )
}

efss_row <- function(
  fcst_row, obs_col = "obs", threshold = 0.1, test_radii = seq(0, 10), num_cores = 1
) {

  members <- grep("_mbr[[:digit:]]", colnames(fcst_row), value = TRUE)
  pair_fss <- function(i) {
    res <- fss_pair(
      fcst_row[[members[i]]][[1]],
      fcst_row[[obs_col]][[1]],
      threshold  = threshold,
      test_radii = test_radii
    )
    tibble::tibble(
      fcdate     = fcst_row[["fcdate"]][1],
      lead_time  = fcst_row[["lead_time"]][1],
      mbr        = members[i],
      obs        = obs_col,
      mbr_cover  = res[["cvr1"]],
      obs_cover  = res[["cvr2"]],
      fss        = res[["fss"]],
      scale      = res[["rad"]] * 2 + 1,
      truncated  = res[["fss"]] < 0.5
    )
  }

  pair_fss_prob <- function(i) {
    res <- fss_pair_prob(
      fcst_row[[members[i]]][[1]],
      fcst_row[[obs_col]][[1]],
      test_radii
    )
  }

  if (num_cores > 1) {
    do.call(
      rbind,
      parallel::mclapply(seq_along(members), pair_fss_prob, mc.cores = num_cores)
    )
  }

  do.call(rbind, lapply(seq_along(members), pair_fss_prob))

}

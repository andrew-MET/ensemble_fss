# Some tests for spatial spread using dFSS
# The scale at which FSS begins to exceed 0.5 for each unique pair of members

library(harpMET)
library(harpIO)
library(harpVis)
library(dplyr)
library(here)

source(here("binary_prob.R"))
source(here("nbhd_upscale.R"))
source(here("fss.R"))

# Read some random members precipitation forecasts plus the control
members <- c(0, sample(seq(1, 29), 6))
fcst <- read_meps(
  2022081518, "Pcp", seq(6, 36), members = members
) %>%
  accumulate(1)

# Plot for a single lead time
map_path <- get_map(dom = get_domain(fcst, meps_mbr000), poly = FALSE)
ggplot(filter(harpPoint::gather_members(fcst)$meps, lead_time == 18)) +
  geom_georaster(aes(geofield = forecast)) +
  geom_path(aes(x, y), map_path, colour = "seagreen3") +
  facet_wrap(vars(member), nrow = 2) +
  scale_fill_viridis_c(
    NULL,
    option = "C",
    limits = c(1, NA),
    na.value = "transparent",
    trans = "log",
    breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64),
    direction = 1
  ) +
  coord_equal(expand = FALSE) +
  theme_harp_map() +
  theme(panel.background = element_rect(fill = "grey30"))



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

  if (pair_fss < 0.5) {
    radius <- radius + 1000 + pair_fss
  }

  list(radius, cover1, cover2)
}

# Function to get the FSS scale for each pair of members
fss_all_pairs <- function(.fcst, threshold = 0.1, test_radii = seq(0, 10)) {

  out <- list()

  # Generate all of the member combinations
  members <- grep("_mbr[[:digit:]]", colnames(.fcst), value = TRUE)
  for (i in 1:nrow(.fcst)) {
    pair_counter <- 0
    inrow <- list()
    message("Forecast : ", .fcst$fcdate[[i]], " | Lead Time : ", .fcst$lead_time[[i]])
    for (j in 1:(length(members) - 1)) {
      for (k in (j + 1):length(members)) {
        pair_counter <- pair_counter + 1
        res <- fss_pair(
          .fcst[[members[j]]][[i]],
          .fcst[[members[k]]][[i]],
          threshold,
          test_radii
        )
        inrow[[pair_counter]] <- res[[1]]
        message(
          members[j], " vs ", members[k], " :",  inrow[[pair_counter]],
          " ", res[[2]], " ", res[[3]]
        )
      }
    }
    out[[i]] <- unlist(inrow)
  }
  out
}

truncate_scales <- function(scls) {
  scls[scls > 1000] <- floor(scls[scls > 1000] - 1000)
  scls
}


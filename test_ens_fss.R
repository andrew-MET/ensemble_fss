library(harp)
library(harpMET)
library(here)
library(dplyr)

source(here("binary_prob.R"))
source(here("nbhd_upscale.R"))
source(here("fss.R"))
source(here("dfss_methods.R"))
source(here("accumulate.R"))
source(here("get_prob.R"))
source(here("ens_fss.R"))

fcst <- read_meps(
  2022081600, "Pcp", lead_time = seq(6, 12, 3), members = seq(0, 5)
) %>%
  accumulate(3) %>%
  .[[1]]

obs <- read_met_analysis(seq_dates(2022081606, 2022081612), "Pcp") %>%
  accumulate(met_analysis, 3) %>%
  .[[1]]

fcst <- join_to_fcst(fcst, obs)

rm(obs)
gc()



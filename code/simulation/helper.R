library(tidyverse)
# library(IMR)
library(kableExtra)
library(magrittr)
library(scales)
library(RSSthemes)

#------------------------------------------------------------------------------


#  source("./code/simulation/config_default.R")
if (file.exists("./code/simulation/config.R")) {
  source("./code/simulation/config.R")
}
source("./code/simulation/simulation_function.R")

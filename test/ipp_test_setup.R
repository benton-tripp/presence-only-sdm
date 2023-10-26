
library(sf)
library(terra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(data.table)
library(knitr)
library(purrr)
library(caret)
library(spatstat)
library(tidyr)

# Set seed for splitting and modeling
set.seed(19)

source("test/ipp_metrics.R")
source("test/get_object.R")


# Load the dataset saved in part 2 of the study 
df <- readRDS("artifacts/final_data/final_data.rds") %>% setDT()

# Define some global variables that will be referenced throughout the modeling 
states <- sort(unique(df$state))
species <- sort(unique(df$common.name))
spec.state <- expand.grid(species=species, 
                          state=states, 
                          stringsAsFactors=F)
source("test/ipp_get_var_imp.R")

# Load the dataset saved in part 2 of the study 
df <- readRDS("artifacts/final_data/final_data.rds") %>% setDT()

# Define some global variables that will be referenced throughout the modeling 
states <- sort(unique(df$state))
species <- sort(unique(df$common.name))
spec.state <- expand.grid(species=species, 
                          state=states, 
                          stringsAsFactors=F)

# Rasters
rasters <- states %>% 
  set_names() %>%
  purrr::map(~{
    r <- rast(file.path("artifacts/preprocessed_rasters_updated",
                        paste0(.x, ".tif")))
    r <- names(r) %>%
      set_names() %>%
      purrr::map(~terra::project(r[[.x]], 
                                 crs(paste0("EPSG:", 4326)))) %>%
      rast()
  })


### Getting the Data - Bird Species Observations and Environmental Rasters


# SETUP -------------------------------------------------------------------

# Load necessary libraries
library(sf)
library(terra)
library(stars)
library(auk)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(purrr)
library(stringr)
library(fs)

setwd("C:/users/bento/gis630_data_preprocessing")

# Hard coded path for data saved in external hard drive
ext.data.path <- "D:/AvianAnalyticsData"

# GET US STATE BOUNDARY DATA ----------------------------------------------

# There are four regions explored in this analysis from 2016 through 2019 (corresponding to
# the 4 boundary datasets):
#   
# -   Colorado, 2016-2019
# -   North Carolina, 2016-2019
# -   Oregon, 2016-2019
# -   Vermont, 2016-2019
# 
# For simplicity, only data for the four states that are being used as observation area
# (https://hub.arcgis.com/datasets/1612d351695b467eba75fdf82c10884f/explore?filters=eyJTVEFURV9BQkJSIjpbIkNPIiwiVlQiLCJOQyIsIk9SIl19&location=48.814319%2C163.610769%2C2.35)
#  needs to be downloaded (as a .shp file). By default, the zipped files should be saved in
#  *US_State_Boundaries.zip*. Extract the files into a folder of the same name within the
#  data directory, and use the following code to split the states into separate
#  files:


# Define the state abbreviations
states <- c("CO", "NC", "OR", "VT")

# Set the output directory
state.output.dir <- "data/US_State_Boundaries"

output.paths <- file.path(state.output.dir, 
                         paste0(states, "_State_Boundaries.shp"))
if (!all(file.exists(output.paths))) {
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(state.output.dir)) {
    dir.create(state.output.dir)
  }
  
  # Set the path for the shapefile
  shapefile.path <- file.path(ext.data.path,
                              "US_State_Boundaries/US_State_Boundaries.shp")
  
  # Load the shapefile into an sf object
  gdf <- st_read(shapefile.path)
  
  # Transform the coordinate reference system to EPSG:5070
  gdf <- st_transform(gdf, crs = 5070)
  
  
  # Loop through state abbreviations and save individual shapefiles
  for (sa in states) {
    tryCatch({
      state.gdf <- gdf[gdf$STATE_ABBR == sa, ]
      output.path <- file.path(state.output.dir, paste0(sa, "_State_Boundaries.shp"))
      if (!file.exists(output.path)) st_write(state.gdf, output.path, quiet = T)
    }, error = function(e) {
      message(paste0('Failed to save shapefile for state ', sa, ': ', e$message))
    })
  }
}



# GET SPECIES DATA --------------------------------------------------------

# Species included are:
#   
# -   Belted Kingfisher
# -   Cedar Waxwing
# -   Downy Woodpecker
# -   Ruddy Duck
# -   Sanderling
# -   Sandhill Crane
# -   Sharp-shinned Hawk
# -   Wild Turkey

# Species Observation Data Pre-Processing Steps:
# 1.  Download each of the regions from the eBird Download Page 
# (https://ebird.org/data/download), filtered by date range (you will 
# need to request access annually). They will initially be downloaded 
# as compressed folders, so the contents will need to be extracted. 
# There is also an eBird API, but the use of the API is beyond the 
# scope of this project.
# 2.  Using the R `auk` package (https://cornelllabofornithology.github.io/auk/), 
# the contents can be filtered and saved within the project data directory:

auk_set_ebd_path("data", overwrite = T)

# Specify where your eBird datasets were downloaded to
ebird.download.dir <- "data/ebird_downloads"

if (!dir.exists(ebird.download.dir)) {
  dir.create(ebird.download.dir)
}

for (sa in states) {
  dname <- file.path(ebird.download.dir, sa)
  if (!dir.exists(dname)) dir.create(dname)
}

# Define species
species <- c("Sandhill Crane", "Sharp-shinned Hawk",
             "Wild Turkey", "Downy Woodpecker",
             "Sanderling", "Cedar Waxwing", 
             "Belted Kingfisher", "Ruddy Duck")

# Parse eBird downloads
for (dir in list.dirs(ebird.download.dir)[-1]) {
  # List .txt files that start with "ebd_"
  in.file <- list.files(path = dir, pattern = "^ebd_.*\\.txt$", full.names = T)[[1]]
  # Use regular expression to extract state abbreviation
  state.abbreviation <- sub(".*_US-([A-Z]{2})_.*", "\\1", in.file)
  out.file <- paste0("data/", state.abbreviation, ".txt")
  # Read in the filtered data using the `auk` library, saving to `out.file`
  auk_ebd(in.file) %>%
    auk_species(species = species) %>% 
    auk_filter(file = out.file, overwrite=T, execute=T)
}

# Show examples (NC data)
# 
# Import kable and leaflet
# library(kable)
# library(leaflet)
# 
# readr::read_delim("data/NC.txt", delim = "\t") %>%
#   group_by(`COMMON NAME`) %>%
#   summarize(count = n()) %>%
#   kable()
#   
# 
# # Read in sample data
# nc.shc <- readr::read_delim("data/NC.txt", delim = "\t") %>%
#   dplyr::filter(`COMMON NAME` == "Sandhill Crane")
# 
# leaflet(nc.shc) %>%
#   addTiles() %>%  # This adds the OpenStreetMap tiles
#   addMarkers(~LONGITUDE, ~LATITUDE) %>%
#   setView(lng = mean(nc.shc$LONGITUDE), lat = mean(nc.shc$LATITUDE), zoom = 7) %>%
#   addProviderTiles(providers$Stamen.Toner)


# Some basic pre-processing
preprocess.obs <- function(data.path) {
  data <- readr::read_delim(data.path, delim = "\t")
  names(data) <- gsub(" ", "_", tolower(names(data)))
  data <- data %>%
    filter(approved == 1) %>%
    dplyr::select(common_name, observation_count, latitude, longitude) %>%
    group_by(common_name, latitude, longitude) %>%
    summarize(observation_count = sum(as.numeric(observation_count), na.rm=T)) %>%
    ungroup() %>%
    as.data.table()
  return(data)
}

obs <- purrr::map(states, function(.x) {
  out <- preprocess.obs(paste0("data/", .x, ".txt"))
  fwrite(out, paste0("data/observations_", .x, ".csv"))
  out
})

names(obs) <- states



# RASTER PREPROCESSING ----------------------------------------------------
# General Raster Pre-Processing Steps
# 
# 1.  Load the environmental raster datasets & State Boundary data
# 2.  Reproject to CRS EPSG:5070
# 3.  Resample each to 5000 x 5000 meters
# 4.  Mask each raster using State Boundary (ensure oceans and great 
#     lakes are set to `nodata`)
# 5.  Output updated raster data


# `terra` equivalent of base R `any()` function
terra.any <- function(r) {
  freqs <- freq(r)
  any(freqs[,1] == 1)
}

# Function to fix NA values in a raster (e.g., NULL values are equal to 999, 
# but should be NA)
fix.raster.na <- function(data.path, 
                          na.val, 
                          raster.name,
                          states=c("CO", "NC", "OR", "VT")) {
  
  for (state in states) {
    
    # Set path of the input raster
    input.raster.path <- file.path(data.path, raster.name)
    
    # Create a temporary output path
    temp.output.path <- file.path(data.path, paste0("temp_", raster.name))
    
    # Read the input raster
    cat(paste0("Loading raster from ", input.raster.path, "...\n"))
    input.raster <- rast(input.raster.path)
    
    # Identify cells with value `na.val`
    msk <- input.raster == na.val
    cat("Checking for improperly formatted NA values...\n")
    if (terra.any(msk)) {
      cat("Updating NA values...\n")
      # Set those cells to NA
      input.raster[msk] <- NA
      
      # Save the output raster to the temporary file path
      cat("Saving updated raster...\n")
      writeRaster(input.raster, temp.output.path)
      
      # Remove temporary prefix and overwrite the original raster
      file.rename(temp.output.path, input.raster.path)
    }
  }
}

# Recursive function to list directories containing .adf files
# Usage:
# raster.directories <- list.raster.dirs.recursive("your_start_directory_path")
list.rasters.recursive <- function(directory.path) {
  
  # List all directories recursively
  all.dirs <- list.dirs(directory.path, recursive = T, full.names = T)
  
  # Filter directories that contain .adf files
  raster.dirs <- all.dirs[sapply(all.dirs, function(dir.path) {
    any(grepl("\\.adf$", list.files(dir.path, full.names = F)))
  })]
  
  return(raster.dirs)
}


general.raster.preprocessing <- function(
    raster.name,
    out.raster.name,
    data.path = "data",
    state.boundary.path = "data/US_State_Boundaries",
    states = c("NC", "CO", "OR", "VT"),
    state.file.suffix = "_State_Boundaries.shp",
    out.path = "../gis630/data",
    crs = 5070,
    wildcard = "*.tif",
    resolution = 5000,
    agg="mean",
    recursive.adf.path=F,
    crop.by.state=T
) {
  # Example usage:
  # general.raster.preprocessing(
  #   raster.name = "canopy/nlcd_2016_treecanopy_2019_08_31",
  #   out.raster.name = "canopy",
  #   out.path = "data/canopy"
  # )
  
  # Make sure output directory exists
  if (!dir.exists(out.path)) dir.create(out.path)
  
  if (!recursive.adf.path) {
    # List all the raster files in the workspace directory
    rasters <- list.files(path = file.path(data.path, raster.name), 
                          pattern = wildcard, 
                          full.names = T)
  } else {
    # List directories containing .adf files
    rasters <- list.rasters.recursive(data.path)
  }
  
  cat("Rasters:", paste(rasters, collapse=", "), "\n")
  
  # Loop through each raster file in the list
  for (i in 1:length(rasters)) {
    
    raster.path <- rasters[[i]]
    
    if (crop.by.state) {
      
      # Read raster
      cat(paste0("[", i, "/", length(rasters), "]:"), "Reading", raster.path, "\n")
      raster <- rast(raster.path)
      
      # Loop through each state
      for (state in states) {
        
        state.out.file <- file.path(out.path, 
                                    paste0(out.raster.name, "_", 
                                           state, ".tif"))
        if (!file.exists(state.out.file)) {
          # Construct the file path to the state boundary shapefile
          state.path <- file.path(state.boundary.path, 
                                  paste0(state, state.file.suffix))
          
          # Load state shape and reproject it to match raster's CRS
          state.shape <- vect(state.path)
          reprojected.shape <- project(state.shape, crs(raster))
          
          # Crop raster by reprojected state shape
          cat("\tCropping raster for", state, "\n")
          # masked.raster <- terra::mask(raster, reprojected.shape)
          masked.raster <- terra::crop(raster, reprojected.shape, mask=T)
          
          # Reproject
          cat("\tReprojecting raster for", state, "\n")
          reprojected.raster <- project(masked.raster, 
                                        crs(paste0("EPSG:", crs)))
          
          # Clean up
          rm(masked.raster)
          gc()
          
          
          ### TODO: FILL NUMERIC VARIATIONS OF NA VALUES WITH NA
          # Check each:
          # DEM - 
          # Land Cover - 
          # Canopy - 
          # Urban Imperviousness -
          # Weather - 
          # Soil - 
          # Hydrography - 
          # NDVI - 
  
          # Resample
          current.res <- terra::res(reprojected.raster)
          cat("\tCurrent Resolution:", current.res, "\n")
          cat("\tTarget Factor:", resolution/terra::res(reprojected.raster)[1], "\n")
          if (any(current.res != resolution)) {
            cat("\tResampling raster for", state, "\n")
            resampled.raster <- terra::aggregate(
              reprojected.raster, 
              fact=c(resolution/current.res[1], resolution/current.res[2]), 
              fun = agg, 
              expand = T) %>% suppressWarnings()
          } else {
            resampled.raster <- reprojected.raster
          }
          
          # Clean up
          rm(reprojected.raster)
          gc()
          
          # Save raster
          cat("\tSaving to", state.out.file, "\n")
          writeRaster(resampled.raster, state.out.file)
          
          # Clean up
          rm(resampled.raster)
          gc()
        }
      }
    } else {
      path.parts <- stringr::str_split(raster.path, "\\.")[[1]]
      if (length(path.parts) > 1) path.parts <- path.parts[[1]]
      path.parts <- stringr::str_split(path.parts, "/")[[1]]
      r.name <- path.parts[[length(path.parts)]]
      out.file <- file.path(out.path, 
                            paste0(out.raster.name, "_", 
                                   r.name, ".tif"))
      if (!file.exists(out.file)) {
        
        # Read raster
        cat(paste0("[", i, "/", length(rasters), "]:"), "Reading", raster.path, "\n")
        raster <- rast(raster.path)
        
        # Reproject
        cat("\tReprojecting raster...\n")
        reprojected.raster <- project(raster, 
                                      crs(paste0("EPSG:", crs)))
        
        # Clean up
        rm(raster)
        gc()
        
        # Resample
        current.res <- terra::res(reprojected.raster)
        cat("\tCurrent Resolution:", current.res, "\n")
        cat("\tTarget Factor:", resolution/terra::res(reprojected.raster)[1], "\n")
        if (any(current.res != resolution)) {
          cat("\tResampling raster...\n")
          resampled.raster <- terra::aggregate(
            reprojected.raster, 
            fact=c(resolution/current.res[1], resolution/current.res[2]), 
            fun = agg, 
            expand = T) %>% suppressWarnings()
        } else {
          resampled.raster <- reprojected.raster
        }
        
        # Clean up
        rm(reprojected.raster)
        gc()
        
        # Save raster
        cat("\tSaving to", out.file, "\n")
        writeRaster(resampled.raster, out.file)
        
        # Clean up
        rm(resampled.raster)
        gc()
      }
    }
    
    # Clean up
    rm(raster)
    gc()
  }
  cat("Finished general pre-processing.\n")
}


# DIGITAL ELEVATION MODEL (DEM) -------------------------------------------


# Combine DEM rasters (Source data has the US split into parts)
# 
# By default, `terra` does most of its work on the disk (as 
# opposed to RAM). However, the temp data for this particular
# process was using ~70GB of storage. This function is set up 
# to reduce this overhead as much as possible.
combine.dem.parts <- function(dem.dir = "data/dem/raw_dem", 
                              output.dir = "data/dem", 
                              output.raster.name = "mosaic_dem") {
  
  if (!dir.exists(output.dir)) dir.create(output.dir)
  
  out.name <- file.path(output.dir, paste0(output.raster.name, ".tif"))
  
  # Get list of rasters recursively
  rasters <- list.rasters.recursive(dem.dir)
  
  # Read the first raster
  cat("1: Reading first raster from", rasters[1], "out of", 
      length(rasters), "rasters. \nWriting to", out.name, ".\n\n")
  
  # Save the combined raster
  writeRaster(rast(rasters[1]), 
              out.name,
              gdal=c("COMPRESS=DEFLATE", "TFW=YES"),
              verbose=T,
              overwrite=T)
  cat("----------------------------\n")
  
  # If there are more rasters, mosaic them iteratively
  if (length(rasters) > 1) {
    for (i in 2:length(rasters)) {
      cat(paste0(i, ":"), "Reading raster from", rasters[i], "out of", 
          length(rasters), "rasters.\nAdding [", i, "/", length(rasters), 
          "] rasters to mosaic.\n\n")
      
      # Save the combined raster
      writeRaster(mosaic(rast(out.name), rast(rasters[i]), fun=mean), 
                  out.name,
                  gdal=c("COMPRESS=DEFLATE", "TFW=YES"),
                  verbose=T,
                  overwrite=T)
      cat("----------------------------\n")
      
      # Clean up
      gc()
    }
  }
  
  cat("Mosaic completed successfully.\n")
}

dem.path <- file.path(ext.data.path, "dem")

if (!file.exists(file.path(dem.path, "mosaic_dem.tif"))) {
  # Combine rasters 
  combine.dem.parts(dem.dir=file.path(dem.path, "raw_dem"),
                    output.dir=dem.path)
}

if (!all(file.exists(paste0("data/dem/dem_", states, ".tif")))) {
  # General pre-processing on each of the DEM parts
  general.raster.preprocessing(
    data.path = ext.data.path,
    raster.name = "dem", 
    out.raster.name = "dem",
    out.path = "data/dem",
    resolution = 5000,
    wildcard="\\.tif$"
  )
}



# URBAN IMPERVIOUSNESS ----------------------------------------------------

# Use the combined raster for the general raster pre-processing
general.raster.preprocessing(
  data.path = ext.data.path,
  raster.name="urban_imperviousness", 
  out.raster.name="urban_imperviousness",
  out.path="data/urban_imperviousness",
  resolution=5000,
  wildcard="\\.tif$"
)


# LAND COVER --------------------------------------------------------------

# Using the combined raster, apply "general raster preprocessing" 
# (resample, reproject, mask)
general.raster.preprocessing(
  data.path = ext.data.path,
  raster.name="land_cover/nlcd_2019_land_cover_l48_20210604", 
  out.raster.name="land_cover",
  out.path="data/land_cover",
  resolution=5000,
  wildcard="\\.tif$",
  agg="modal"
)


# CANOPY ------------------------------------------------------------------

# Using the combined raster, apply "general raster preprocessing" 
# (resample, reproject, mask)
general.raster.preprocessing(
  raster.name="canopy/nlcd_2016_treecanopy_2019_08_31", 
  out.raster.name="canopy",
  out.path="data/canopy",
  resolution=5000
)


# WEATHER -----------------------------------------------------------------

get.weather.data <- function(data.path) {
  # Downloads and processes weather raster data for specified variables and 
  # years at a 4km resolution 
  # and 30-year monthly normals at an 800m resolution. The function downloads,
  #  aggregates, and resamples 
  # the rasters before trimming them to the North Carolina boundary.
  
  if (!dir.exists(data.path)) {
    stop(paste("Data path '", data.path, "' not found."))
  }
  
  vars <- c("ppt", "tmax", "tmin")
  
  # Setup for yearly 4km resolution 
  yrs <- c(2017, 2018, 2019)
  pairs <- expand.grid(vars, yrs)
  
  # Setup for 30 year normal monthly 800m resolution
  mnths <- sprintf("%02d", 1:12)
  norm.pairs <- expand.grid(vars, mnths)
  
  ### Get Raster Data ######
  cat("Getting explanatory Weather Rasters...\n")
  
  # Data documentation https://www.prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf
  out.path <- file.path(data.path, "weather/")
  dir.create(out.path, recursive = TRUE, showWarnings = FALSE)
  
  # 4km yearly data (for 2017-2019)
  for (i in 1:nrow(pairs)) {
    v <- pairs[i, 1]
    y <- pairs[i, 2]
    
    dwnld.out <- file.path(out.path, paste0(v, "_", y, ".zip"))
    dwnld.path <- file.path(out.path, paste0(v, "_", y))
    
    if (!dir.exists(dwnld.path)) {
      url <- paste0("https://services.nacse.org/prism/data/public/4km/", v, "/", y)
      cat(paste("Downloading weather data from", url, "...\n"))
      download.file(url, dwnld.out)
      cat(paste("Saved", v, "/", y, "to", dwnld.out, "\n"))
      unzip(dwnld.out, exdir = dwnld.path)
      cat(paste("Extracted", v, "/", y, "from", dwnld.out, "to", dwnld.path, "\n"))
      file.remove(dwnld.out)
    }
  }
  
  # 800m monthly data (30 year normals)
  for (i in 1:nrow(norm.pairs)) {
    v <- norm.pairs[i, 1]
    m <- norm.pairs[i, 2]
    
    dwnld.out <- file.path(out.path, paste0(v, "_", m, ".zip"))
    dwnld.path <- file.path(out.path, paste0(v, "_", m))
    
    if (!dir.exists(dwnld.path)) {
      url <- paste0("https://services.nacse.org/prism/data/public/normals/800m/", 
                    v, "/", m)
      cat(paste("Downloading weather data from", url, "...\n"))
      download.file(url, dwnld.out)
      cat(paste("Saved", v, "/", m, "to", dwnld.out, "\n"))
      unzip(dwnld.out, exdir = dwnld.path)
      cat(paste("Extracted", v, "/", m, "from", dwnld.out, "to", dwnld.path, "\n"))
      file.remove(dwnld.out)
    }
  }
  
  cat("Finished getting weather data.\n")
}

get.weather.data("data")



get.rasters.from.dir <- function(var, periods, out.path) {
  rasters <- list()
  for(p in periods) {
    dir.path <- file.path(out.path, paste0(var, "_", p))
    raster.file <- list.files(path = dir.path, pattern = paste0(var, "\\.bil$"),
                              full.names = T)[1]
    rasters[[p]] <- rast(raster.file)
  }
  return(rasters)
}

data.path <- "data/weather"
wspace <- "data/weather/aggregated"
if(!dir.exists(wspace)) {
  dir.create(wspace)
}
max.temp.data <- "max_temp.tif"
min.temp.data <- "min_temp.tif"
avg.prcp.data <- "avg_prcp.tif"
states <- c("CO", "NC", "OR", "VT")
vars <- c("ppt", "tmax", "tmin")
yrs <- 2017:2019
pairs <- expand.grid(vars, yrs)

mnths <- sprintf("%02d", 1:12)
norm.pairs <- expand.grid(vars, mnths)

for(v in vars) {
  cat(paste0("Checking ", v, "...\n"))
  out.path <- file.path(wspace, v)
  if(!dir.exists(out.path)) {
    dir.create(out.path)
  }
  raster.name <- ifelse(v == "tmax", max.temp.data, 
                        ifelse(v == "tmin", min.temp.data, avg.prcp.data))
  rasters <- get.rasters.from.dir(v, yrs, data.path)
  agg.func <- switch(v,
                     tmax = "max",
                     tmin = "min",
                     ppt = "mean")
  
  agg.rasters <- cellStats(do.call(stack, rasters), fun = agg.func)
  writeRaster(agg.rasters, file = file.path(out.path, raster.name), 
              format = "GTiff", overwrite = T)
  cat(paste0("Finished aggregating ", v, " for all years.\n"))
  
  rasters <- get.rasters.from.dir(v, mnths, data.path)
  norm.raster.name <- paste0(v, "_30yr_800m.tif")
  agg.rasters <- cellStats(do.call(stack, rasters), fun = agg.func)
  writeRaster(agg.rasters, file = file.path(out.path, norm.raster.name), 
              format = "GTiff", overwrite = T)
  
  initial.weight.4km <- 3.0
  initial.weight.800m <- 3.0 / 30.0
  total.weight <- initial.weight.4km + initial.weight.800m
  normalized.weight.4km <- initial.weight.4km / total.weight
  normalized.weight.800m <- initial.weight.800m / total.weight
  
  combined.raster <- rast(file.path(out.path, raster.name)) * 
    normalized.weight.4km + 
    rast(file.path(out.path, norm.raster.name)) * 
    normalized.weight.800m
  writeRaster(combined.raster, 
              file = file.path(out.path, paste0("aggregated_", raster.name)), 
              format = "GTiff", overwrite = T)
  
  general.raster.preprocessing(
    raster.name = paste0("weather/aggregated/", v),
    out.raster.name = v,
    out.path = "data/weather",
    wildcard = "agg*"
  )
}

cat("Finished weather data pre-processing.\n")


# HYDROGRAPHY (WATER BODIES AND COASTS) ------------------------------------

convert.hydro.gdb <- function(data.path = "data/hydrography/",
                              gdb.name = "hydrusm010g.gdb_nt00897/hydrusm010g.gdb",
                              data.sets = c("Coastline", "Waterbody"),
                              output.path = "data/hydrography/") {
  
  # Iterate over each dataset and convert
  for (ds in data.sets) {
    gdb.path <- file.path(data.path, gdb.name, ds)
    
    # Check if the path exists before trying to read it
    if (dir.exists(gdb.path)) {
      # Read the geodatabase layer using st_read
      data <- st_read(dsn = file.path(data.path, gdb.name), layer = ds, quiet = T)
      
      # Write the shapefile
      st_write(obj = data, dsn = file.path(output.path, paste0(ds, ".shp")), 
               quiet = T)
      
      cat(sprintf("Converted %s from geodatabase to shapefile.\n", ds))
    } else {
      cat(sprintf("Cannot find the dataset %s at %s\n", ds, gdb.path))
    }
  }
}

convert.hydro.gdb()


data.path <- "data"
hyd.data.path <- file.path(data.path, "hydrography")
states <- c("CO", "NC", "OR", "VT")
resolution <- 5000  # in the same units as your CRS

update.array <- function(arr) {
  arr.updated <- matrix(as.numeric(NA), nrow = nrow(arr), ncol = ncol(arr))
  
  # Perform the dilation operation multiple times
  dilated1 <- terra::focal(arr, w = matrix(1,3,3), fun = max)
  dilated2 <- terra::focal(dilated1, w = matrix(1,3,3), fun = max)
  dilated3 <- terra::focal(dilated2, w = matrix(1,3,3), fun = max)
  dilated4 <- terra::focal(dilated3, w = matrix(1,3,3), fun = max)
  
  arr.updated[dilated1 == 1 & arr != 1] <- 0.8
  arr.updated[dilated2 == 1 & dilated1 != 1] <- 0.6
  arr.updated[dilated3 == 1 & dilated2 != 1] <- 0.4
  arr.updated[dilated4 == 1 & dilated3 != 1] <- 0.2
  arr.updated[arr == -1] <- -1
  arr.updated[arr == 1] <- 1
  
  return(arr.updated)
}

cat("Reading waterbody data...\n")
gdf <- st_read(file.path(hyd.data.path, 'Waterbody.shp'))

# Filter out records where Feature == "Lake Dry"
gdf <- gdf %>% filter(Feature != "Lake Dry")

# Overwrite the original shapefile
st_write(gdf, file.path(hyd.data.path, 'Waterbody.shp'))

gdf$raster.value <- 0

for (state in states) {
  cat(paste("Reading", state, "boundary data...\n"))
  state.boundary <- st_read(file.path(data.path, paste0(state, "_State_Boundaries.shp")))
  
  state.boundary$raster.value <- 1
  
  cat("Reprojecting...\n")
  gdf <- st_transform(gdf, crs = 5070)
  state.boundary <- st_transform(state.boundary, crs = 5070)
  
  cat("Combining US states...\n")
  state.boundary <- st_union(state.boundary)
  
  bounds <- st_bbox(state.boundary)
  
  xs <- seq(from = bounds["xmin"], to = bounds["xmax"], by = resolution)
  ys <- seq(from = bounds["ymin"], to = bounds["ymax"], by = resolution)
  grid.points <- expand.grid(xs, ys) %>%
    st_as_sf(coords = c("Var1", "Var2"), crs = 5070) %>%
    rename(X = Var1, Y = Var2)
  
  grid.points$is.within.bounds <- as.integer(st_within(grid.points, state.boundary))
  grid.points$is.water <- as.integer(st_within(grid.points, gdf))
  grid.points[is.na(grid.points$is.water), "is.water"] <- -1
  grid.points <- grid.points %>%
    select(-X, -Y) %>%
    st_as_sf()
  
  raster <- st_rasterize(grid.points)
  raster.updated <- update.array(raster)
  
  out.fn <- file.path(hyd.data.path, paste0(state, "_Waterbody.tif"))
  terra::writeRaster(raster.updated, out.fn, overwrite = TRUE)
  
  cat(paste("Writing final raster for", state, "...\n"))
}


# Set file paths
input.path <- "data/hydrography/Coastline.shp"
output.path <- "data/hydrography/Coastline.Buffer.shp"

# Read the shapefile
coastline <- st_read(input.path)

# Create buffer
coastline.buffer <- st_buffer(coastline, dist = 7500)

# Merge all buffered features into a single feature
coastline.buffer.merged <- st_union(coastline.buffer)

# Write the buffered feature back to a shapefile
st_write(coastline.buffer.merged, output.path)

cat("Buffer operation completed.\n")


# Define your desired resolution
resolution <- 5000

# Load the coastline buffer shapefile
coastline.buffer <- st_read('data/hydrography/Coastline_Buffer.shp')

# Reproject the coastline buffer to the desired CRS
coastline.buffer <- st_transform(coastline.buffer, crs = 5070)

# Determine the extent of the coastline buffer
ext <- st_bbox(coastline.buffer)

# Create a grid of points within this extent
grid <- expand.grid(x = seq(ext$xmin, ext$xmax, by = resolution), 
                    y = seq(ext$ymin, ext$ymax, by = resolution))

# Convert to sf object
points <- st_as_sf(grid, coords = c("x", "y"), crs = 5070)

# Assign each point a value based on whether it falls within the coastline buffer polygon or not
points$is.buffer <- as.integer(st_within(points, coastline.buffer))

# Convert to stars object
raster.grid <- st_rasterize(points)

# Define the output filename
out.fn <- 'data/hydrography/Coastline_Buffer_Raster.tif'

# Save to raster file
write_stars(raster.grid, out.fn, driver = "GTiff", overwrite = T)

# Message
cat("Raster file created successfully!\n")

general.raster.preprocessing(
  raster.name="hydrography", 
  out.raster.name="coastline",
  out.path="data/hydrography",
  wildcard="Coastline_Buffer*",
  raster.scale=5000
)



# SOIL --------------------------------------------------------------------

# Open the data source
vector.ds <- "data/soil/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp"
raster.fn <- "data/soil/raster/gsmsoilmu_a_us.tif"

shp.to.raster <- function(vector.ds, raster.fn, akey) {
  vec <- vect(vector.ds)
  pixel.size <- 0.01  # Define this according to your needs
  
  ext <- ext(vec)
  size.x <- ceiling((ext[2] - ext[1]) / pixel.size)
  size.y <- ceiling((ext[4] - ext[3]) / pixel.size)
  
  # Create the raster dataset
  rast.ds <- rast(nrows=size.y, ncols=size.x, xmin=ext[1], xmax=ext[2], 
                  ymin=ext[3], ymax=ext[4])
  proj <- crs(vec)
  crs(rast.ds) <- proj
  
  # Rasterize
  rast.ds <- rasterize(vec, rast.ds, field=akey, fun=max) # assuming we want to retain the max value where polygons overlap
  writeRaster(rast.ds, raster.fn, overwrite=T)
}

# Load the shapefile into an sf object
gdf <- st_read("data/soil/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp")

# Create a tibble that maps each unique MUSYM to a unique integer
musym.codes <- tibble(MUSYM = unique(gdf$MUSYM), 
                      MUSYM_CODE = 1:length(unique(gdf$MUSYM)))

# Merge the MUSYM_CODE field into the sf object
gdf <- left_join(gdf, musym.codes, by = "MUSYM")

# Save the sf object back to a shapefile
st_write(gdf, "data/soil/codes/gsmsoilmu_a_us.shp", delete_dsn = T)

vector.ds <- "data/soil/codes/gsmsoilmu_a_us.shp"
raster.fn <- "data/soil/raster/gsmsoilmu_a_us.tif"
shp.to.raster(vector.ds, raster.fn, akey="MUSYM_CODE")

general.raster.preprocessing(
  raster.name="soil", 
  out.raster.name="soil",
  out.path="data/soil",
  resolution=5000
)


# VEGETATION INDEX --------------------------------------------------------

# Iterate through each season
for (season in c("Spring", "Summer", "Fall", "Winter")) {
  cat(paste0("Applying raster preprocessing for ", season, "...\n"))
  tryCatch({
    # Using the combined raster, apply "general raster preprocessing" 
    # (resample, reproject, mask)
    general.raster.preprocessing(
      data.path=ext.data.path,
      raster.name=paste0("NDVI/US_eVSH_NDVI-", season, "-2021"), 
      out.raster.name=paste0(season, "_NDVI"),
      out.path="data/NDVI",
      wildcard="*1KM\\.VI_NDVI.*\\.tif$",
      resolution=5000)
    cat("-----------------\n")
  }, error = function(e) {
    cat(paste0("An error occurred while processing ", season, ": ", e$message, "\n"))
  })
}
cat("Finished NDVI pre-processing.\n")


# FINAL PRE-PROCESSING STEPS ----------------------------------------------

# Final Pre-Processing Steps (Ensuring Shape/Size Conformity)


get.file.paths.list <- function(
    base.path = "data",
    states = c("CO", "NC", "OR", "VT"),
    categories = c(
      "canopy/canopy_",
      "dem/dem_",
      "hydrography/coastline_",
      "hydrography/Waterbody",
      "land_cover/land_cover_",
      "NDVI/Fall_NDVI_",
      "NDVI/Spring_NDVI_",
      "NDVI/Summer_NDVI_",
      "NDVI/Winter_NDVI_",
      "soil/soil_",
      "urban_imperviousness/urban_imperviousness_",
      "weather/ppt_",
      "weather/tmax_",
      "weather/tmin_"
    ), verbose = F) {
  
  file.paths.list <- list()
  
  for(category in categories) {
    category.name <- str_split(category, "/")[[1]][1]
    filename.prefix <- str_split(category, "/")[[1]][2]
    
    for(state in states) {
      if(filename.prefix == "Waterbody") {
        file.path <- paste0(base.path, "/", category.name, "/", state,
                            "_", filename.prefix, ".tif")
      } else {
        file.path <- paste0(base.path, "/", category, state, ".tif")
      }
      inner.key <- gsub("_$", "", filename.prefix)
      file.paths.list[[category.name]][[inner.key]][[state]] <- file.path
    }
  }
  
  if(verbose) {
    walk(file.paths.list, function(category) {
      cat(names(category), ":\n")
      walk(category, function(subcategory) {
        cat("  ", names(subcategory), ":\n")
        walk(subcategory, function(path, state) {
          cat("    ", state, ": ", path, "\n")
        })
      })
    })
  }
  
  return(file.paths.list)
}

file.paths.list <- get.file.paths.list(verbose=T)

for(category in names(file.paths.list)) {
  cat(category, ":\n")
  
  for(subcategory in names(file.paths.list[[category]])) {
    cat("  ", subcategory, ":\n")
    
    for(state in names(file.paths.list[[category]][[subcategory]])) {
      file.path <- file.paths.list[[category]][[subcategory]][[state]]
      tryCatch({
        raster <- rast(file.path)
        cat("    ", state, " dimensions: ", dim(raster)[2:1], "\n")
      }, error = function(e) {
        cat("    ", state, ": Could not read file. ", conditionMessage(e), "\n")
      })
    }
  }
}


dst.crs <- "EPSG:5070"
data.path <- "data"

resample.raster <- function(src.path, dst.path, dst.crs, target, 
                            resampling = "nearest") {
  
  # Reading the target raster
  tar <- rast(target)
  
  # Create the masking polygon from the target raster
  ext <- ext(tar)
  mask.geom <- st_as_sfc(st_bbox(c(ext[1], ext[3], ext[2], 
                                   ext[4]), crs = crs(tar)))
  
  # Read the source raster
  src <- rast(src.path)
  
  # Reprojecting the raster
  reprojected.raster <- project(src, dst.crs, method = resampling)
  
  # Aligning the raster to match target dimensions
  aligned.raster <- align(reprojected.raster, tar)
  
  # Applying the mask to the raster
  masked.raster <- mask(aligned.raster, mask.geom)
  
  # Write to the destination file
  writeRaster(masked.raster, dst.path)
  
}

target.raster.list <- file.paths.list[["NDVI"]][["Fall_NDVI"]]
outpath <- file.path(data.path, "rasters")
dir.create(outpath, showWarnings = FALSE)

walk2(file.paths.list, names(file.paths.list), function(categories, category.name) {
  walk2(categories, names(categories), function(subcategories, subcategory.name) {
    walk2(subcategories, names(subcategories), function(file.path, state) {
      raster.name <- tools::file_path_sans_ext(basename(file.path))
      cat("Resampling and transforming", raster.name, "...\n")
      destination.raster <- file.path(outpath, raster.name, '.tif')
      if(file.path == target.raster.list[[state]]) {
        file.copy(file.path, destination.raster)
      } else {
        if(grepl("land_cover", raster.name) || grepl("soil", raster.name)) {
          resampling.method <- "nearest"
        } else {
          resampling.method <- "bilinear"
        }
        resample.raster(file.path, destination.raster, dst.crs,
                        target.raster.list[[state]], resampling.method)
      }
    })
  })
})

cat("Fixing land cover NA values...\n")
fix.lc.na(data.path = "data/rasters", na.val = 0)
cat("Finished.\n")



folder.path <- 'data/rasters'

# Get the list of files in the folder
file.list <- fs::dir_ls(folder.path, regexp = "\\.tif$")

# Loop through all the .tif files in the folder
for (file.path in file.list) {
  
  # Open the raster file
  raster.obj <- terra::rast(file.path)
  
  # Print the shape of the raster file
  file.name <- fs::path_file(file.path)
  cat(sprintf("Shape of %s: %s\n", file.name, 
              paste0(dim(raster.obj), collapse = " x ")))
}



data.path <- 'data/'
states <- c("CO", "NC", "OR", "VT")

for (state in states) {
  
  input.name <- paste0("observations_", state, ".csv")
  output.file <- paste0("observations_", state, ".shp")
  
  cat(sprintf("Reading %s data...\n", state))
  
  bird.data.path <- fs::path_join(c(data.path, input.name))
  df <- readr::read_csv(bird.data.path)
  
  # Convert the bird data to a sf object
  cat(sprintf("Converting %s data to sf...\n", state))
  geo.df <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
  
  # Define target CRS
  target.crs <- "EPSG:5070"
  cat(sprintf("Updating CRS for %s...\n", state))
  
  # If the bird data is not in the correct CRS, convert it
  if (st_crs(geo.df) != target.crs) {
    geo.df <- st_transform(geo.df, target.crs)
  }
  
  cat(sprintf("Getting centroid of %s data...\n", state))
  # Get centroid of the entire sightings
  centroid <- st_centroid(st_union(geo.df))
  
  cat(sprintf("Calculating distances from centroid for %s...\n", state))
  # Assign each point a distance attribute, being the distance from the centroid
  geo.df$distance <- st_distance(geo.df, centroid)
  
  cat(sprintf("Sorting by distance from centroid in %s data...\n", state))
  # Sort the sf object by the distance attribute
  geo.df <- geo.df %>% arrange(distance)
  
  cat(sprintf("Saving %s data to %s...\n", state, output.file))
  # Drop the distance attribute (we don't need it anymore)
  geo.df$distance <- NULL
  
  st_write(geo.df, fs::path_join(c(data.path, output.file)))
  
  cat("--------------\n")
}


data.path <- 'data/'

states <- c("CO", "NC", "OR", "VT")

for (state in states) {
  
  input.name <- paste0("observations_", state, ".shp")
  out.name <- paste0("all_data_", state, ".shp")
  output.no.null <- paste0("all_data_no_null_", state, ".shp")
  csvfile <- fs::path_join(c(data.path, paste0("final_data_", state, ".csv")))
  
  raster.path <- fs::path_join(c(data.path, "rasters"))
  
  # Create a named list to store raster paths
  raster.paths <- list()
  
  # Loop through all the files in the folder
  for (file.name in fs::dir_ls(raster.path, regexp = paste0('.*', state, '.*\\.tif$'))) {
    raster.name <- fs::path_file(file.name) %>%
      fs::path_ext_remove() %>%
      stringr::str_remove(paste0("(_", state, "|", state, "_)"))
    raster.paths[[raster.name]] <- file.name
  }
  
  cat(sprintf("Extracting points to values for %s...\n", state))
  
  # Load observations shapefile
  obs <- st_read(fs::path_join(c(data.path, input.name)))
  
  # Extract raster values to the points
  for (r_name in names(raster.paths)) {
    r <- terra::rast(raster.paths[[r_name]])
    obs[[r_name]] <- terra::extract(r, obs)
  }
  
  cat(sprintf("Selecting non-null points for %s...\n", state))
  
  # Select points that have non-null raster values for all raster layers
  non_null_condition <- !rowSums(is.na(obs[names(raster.paths)])) 
  obs_no_null <- obs[non_null_condition,]
  
  # Write to a shapefile
  st_write(obs_no_null, fs::path_join(c(data.path, output.no.null)))
  
  # Save as CSV
  cat(sprintf("Saving %s final data to csv...\n", state))
  write_csv(as.data.frame(obs_no_null), csvfile)
  
  # Fix names
  cat("Fixing column names...\n")
  df <- read_csv(csvfile)
  df <- df %>%
    rename(
      common_name = common_nam,
      observations = observatio,
      ndvi_spring = Spring_NDV,
      ndvi_summer = Summer_NDV,
      ndvi_winter = Winter_NDV,
      ndvi_fall = Fall_NDVI,
      urban_imperviousness = urban_impe
    )
  
  write_csv(df, csvfile)
  cat(sprintf('%s column names: %s\n', state, paste(names(df), collapse = ', ')))
  cat("Finished.\n")
  cat("--------------\n")
}


for (state in states) {
  
  csvfile <- fs::path_join(c(data.path, paste0("final_data_", state, ".csv")))
  
  # Load the CSV file
  df <- read_csv(csvfile)
  
  # Perform one-hot encoding on the "common_name" field
  df.encoded <- df %>%
    mutate(common.name = paste0("species_", tolower(common.name))) %>%
    mutate(common.name = gsub(" ", "_", common.name)) %>%
    mutate(common.name = gsub("-", "_", common.name)) %>%
    pivot_wider(names_from = common.name, values_from = common.name, 
                values_fill = list(common.name = 0), values_fn = length) %>%
    bind_cols(df, .)
  
  # Save the encoded DataFrame to a new CSV file
  encoded.csvfile <- fs::path_join(c(data.path, paste0("encoded_data_", state, ".csv")))
  write_csv(df.encoded, encoded.csvfile)
  
  cat(sprintf("Processed and saved one-hot encoded data for %s.\n", state))
}





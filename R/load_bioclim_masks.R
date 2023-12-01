cwd <- getwd()
setwd("C:/Users/bento/gis630")

load.bioclim.masks <- function() {
  # Load libraries
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
  library(tidyr)
  
  # Set seed for splitting
  set.seed(19)
  
  # Load dataset containing observation points, 
  # but we can remove the "absence" points
  df <- setDT(readRDS("artifacts/final_data/final_data.rds"))[
    observations > 0, .(common.name, state, lon, lat, geometry)]
  
  # Define some global variables that will be referenced throughout the modeling 
  states <- sort(unique(df$state))
  species <- sort(unique(df$common.name))
  spec.state <- expand.grid(species=species, 
                            state=states, 
                            stringsAsFactors=F)
  
  # Load Rasters
  r.list <- states %>%
    set_names() %>%
    purrr::map(~rast(file.path("artifacts/final_rasters", 
                               paste0(.x, ".tif"))))
  
  # Pre-processed rasters
  # Define min/max scaling function for rasters
  min.max.scale.raster <- function(r, x=NULL, na.rm=T) {
    # If training data is provided
    if (!is.null(x)) {
      min.x <- min(x, na.rm=na.rm)
      max.x <- max(x, na.rm=na.rm)
      # Or just scale based on full raster values
    } else {
      min.x <- min(values(r), na.rm=na.rm)
      max.x <- max(values(r), na.rm=na.rm)
    }
    (r - min.x) / (max.x - min.x)
  }
  
  # Get model or other object from cache if it has been saved before
  get.object <- function(obj, file.name, obj.path, 
                         read.func=readRDS, save.func=saveRDS, ...) {
    f.path <- file.path(obj.path, file.name)
    if (!dir.exists(obj.path)) {
      dir.create(obj.path)
    }
    # Object cache check
    if (file.exists(f.path)) {
      obj <- read.func(f.path)
    } else {
      save.func(obj, f.path, ...)
    }
    obj
  }
  
  # Get training/testing spatstat data, and 
  # pre-process (scale) rasters 
  r.list <- states %>%
    set_names() %>%
    purrr::map(function(st) {
      get.object(obj = {
        # Get raster
        r <- r.list[[st]]
        # Compute filtered & pre-processed rasters
        names(r) %>%
          set_names() %>%
          purrr::map(function(.x) {
            cat("Pre-processing", .x, "data in", st, "\n")
            r.layer <- r[[.x]]
            if (length(unique(values(r.layer))) > 2) {
              # Scale the data
              return(min.max.scale.raster(r.layer))
            } else {
              return(r.layer)
            }
          }) %>%
          terra::rast()
      },
      file.name=paste0(st, ".tif"),
      obj.path="artifacts/preprocessed_rasters_updated",
      read.func=terra::rast,
      save.func=terra::writeRaster)
    })
  
  
  stratified.split.idx <- function(df, p=0.7, lat.lon.bins=25) {
    # Cut along lat/lon values to create grids (lat.bin & lon.bin)
    # lat.lon.bins is the number of divisions you want
    df$lat.bin <- cut(df$lat, breaks=lat.lon.bins, labels = F)
    df$lon.bin <- cut(df$lon, breaks=lat.lon.bins, labels = F)
    
    # Create a new variable combining the stratification variables
    df %>%
      # stratify on lat/lon bins, species, state
      mutate(strata = paste(lat.bin, lon.bin, common.name, state)) %>%
      pull(strata) %>%
      # Create the data partitions
      createDataPartition(., p = p, list = F) %>%
      suppressWarnings()
  }
  
  prepare.data <- function(df, p=.7, lat.lon.bins=25) {
    train.index <- stratified.split.idx(df, p=p, lat.lon.bins = lat.lon.bins)
    df.train <- df[train.index, ]
    df.test <- df[-train.index, ]
    
    list(train = df.train, 
         test = df.test,
         index = train.index)
  }
  
  # Get train/test indices
  train.test <- prepare.data(df, .7)
  
  # Split datatset
  df.train <- df[train.test$index,] 
  df.test <- df[-train.test$index,]
  
  bioclim.suitability <- function(sf.subset, r) {
    
    # Remove spatial layers (lat, lon, etc.)
    r <- r[[names(r)[!names(r) %in% 
                       c("lat", "lat.sqrd", "lon", 
                         "lon.lat.interaction", "lon.sqrd")]]]
    
    # Extract values from rasters at presence locations
    values.at.presences <- terra::extract(r, st_transform(sf.subset, crs=crs(r)))
    
    # Determine suitability masks based on unique values
    suitable.masks <- names(r) %>%
      set_names() %>%
      purrr::map(function(name) {
        unique.vals <- unique(values.at.presences[[name]])
        if (length(unique.vals) == 2) {
          # Binary method
          mode.val <- as.numeric(names(sort(table(values.at.presences[[name]]), 
                                            decreasing = T)[1]))
          return(r[[name]] == mode.val)
        } else {
          # 10th and 90th percentile method
          lower.bound <- quantile(values.at.presences[[name]], 0.1, na.rm = T)
          upper.bound <- quantile(values.at.presences[[name]], 0.9, na.rm = T)
          return((r[[name]] >= lower.bound) & (r[[name]] <= upper.bound))
        }
      }) %>% rast()
    
    # Combine individual masks to get areas that have the highest suitablity
    sum(suitable.masks)
  }
  
  bcr.dir <- "artifacts/bioclim_suitable_rast"
  if (!dir.exists(bcr.dir )) dir.create(bcr.dir )
  bioclim.list <- purrr::map(1:nrow(spec.state), function(i) {
    st=spec.state[i,]$state
    spec=spec.state[i,]$species
    sf.subset <- sf::st_as_sf(df[state == st & common.name==spec], 
                              crs=4326, sf_column_name = "geometry")
    f <- file.path(bcr.dir, paste0(st, "_", spec, ".tif"))
    if (!file.exists(f)) {
      suitable.areas <- bioclim.suitability(sf.subset, r.list[[st]])
      terra::writeRaster(suitable.areas, f)
    } else {
      suitable.areas <- rast(f)
    }
    list(
      bioclim=suitable.areas,
      state=st,
      species=spec,
      df=sf.subset
    )
  })
  
  # Load 5k resolution mask rasters
  mask.dir <- "artifacts/masks_5k"
  masks <- states %>%
    set_names() %>%
    purrr::map(function(st) {
      species %>%
        set_names() %>%
        purrr::map(function(spec) {
          fname <- file.path(mask.dir, paste0(st, "_", spec, ".tif"))
          rast(fname)
        }) %>% rast()
    })
  
  # Load the original sample regions for pseudo-absence points
  pa.regions.dir <- "artifacts/pseudo_absence_regions"
  pa.regions <- states %>%
    set_names() %>%
    purrr::map(function(st) {
      species %>%
        set_names(species) %>%
        purrr::map(function(spec) {
          fname <- file.path(
            pa.regions.dir,
            gsub(" |\\-", "_", paste0(paste(st, tolower(spec), sep="_"), ".tif"))
          )
          rast(fname)
        }) %>%
        rast()
    })
  
  purrr::walk(1:nrow(spec.state), function(i) {
    st <- spec.state[i,]$state
    spec <- spec.state[i,]$species
    msk <- masks[[st]][[spec]]
    pa.region <- pa.regions[[st]][[spec]]
    idx <- purrr::detect_index(bioclim.list, ~.x$state == st & 
                                 .x$species == spec)
    bioclim.list[[idx]]$mask <<- msk
    bioclim.list[[idx]]$pseudo.absence.region <<- pa.region
    qnt.5 <- quantile(values(bioclim.list[[idx]]$bioclim), .5, na.rm=T)
    bioclim.list[[idx]]$qnt.5 <<- qnt.5
    bioclim.list[[idx]]$final.pseudo.absence.region <<-
      bioclim.list[[idx]]$bioclim <= qnt.5 & terra::lapp(pa.region, function(.x) {
        ifelse(is.na(.x), 0, 1)}) == 1
  })
  
  return(bioclim.list)
}

bioclim.list <- load.bioclim.masks()


setwd(cwd)

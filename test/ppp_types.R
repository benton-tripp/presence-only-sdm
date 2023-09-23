
# Minimal example (test script)

# Load libraries
library(spatstat)
library(dplyr)
library(terra)
library(sf)
library(data.table)
library(purrr)
library(ggplot2)

# Load pre-processed (training) data
df.train <- readRDS("test/df.train.rds")
# Load pre-processed rasters (covariates)
r <- rast("test/NC_covariates.tif")

# Load feature selection results from EDA
fs.df <- readRDS("test/fs.df.rds") %>%
  filter(!(var1 %in% c("lon.sqrd", "lat.sqrd")) & 
           !(var2 %in% c("lon.sqrd", "lat.sqrd"))) %>%
  head(25)

# Specify the state and bird species
st <- "NC"
spec <- "Sharp-shinned Hawk"

# filter by state & species
ipp.df <- df.train %>% 
  filter(state == st & common.name == spec) %>%
  # select location point
  select(lon, lat, presence) %>% 
  # Get the unique points, since I am not accounting for the temporal nature of the data
  unique() 

# Get raster ranges
get.range <- function(r, .crs=4326) {
  e <- r %>% 
    terra::ext() %>%
    terra::project(x=., 
                   from=crs(r, proj=T), 
                   to=st_crs(.crs, parameters=T)$Wkt)
  list(
    lon=c(e$xmin[[1]], e$xmax[[1]]),
    lat=c(e$ymin[[1]], e$ymax[[1]])
  )
}

# Ensure counter-clockwise direction
is.ccw <- function(p) {
  tryCatch({owin(poly=p); T}, error=function(e) F)
}

# Check that all polygons were traversed in the right direction
ensure.ccw <- function(polygon.list) {
  lapply(polygon.list, function(p) {
    # Check the first polygon (outer boundary)
    if (!is.ccw(p)) {
      p$x <- rev(p$x)
      p$y <- rev(p$y)
    }
    p
  })
}

# Function to convert polygon to required list format
convert.to.list.format <- function(sf.object) {
  polygons.list <- lapply(1:nrow(sf.object), function(i) {
    sfc <- st_geometry(sf.object[i,])
    if (class(sfc)[[1]] == "sfc_MULTIPOLYGON") {
      sfc <- st_cast(sfc, "POLYGON")
    }
    lapply(sfc, function(poly) {
      list(x = st_coordinates(poly)[,1], 
           y = st_coordinates(poly)[,2])
    })
  })
  
  # If the object has only one row, we unlist one level to adhere 
  # to the described format for `ppp` objects
  if (length(polygons.list) == 1) {
    polygons.list <- polygons.list[[1]]
  }
  
  if (length(polygons.list) > 1) {
    # Ensure counter-clockwise
    polygons.list <- ensure.ccw(polygons.list)
  }
  
  return(polygons.list)
}

# Get the first layer, set it to either NA or TRUE
r.poly <- terra::project(x=sum(r, na.rm=F), 
                         y=st_crs(4326, parameters=T)$Wkt) %>% 
  lapp(function(z) ifelse(is.na(z), NA, T)) %>%
  terra::as.polygons() %>%
  # Convert to polygon
  st_as_sf()

# Convert polygon to list
r.poly.list <- convert.to.list.format(r.poly)

# Convert your point dataframe to an sf object
ipp.sf <- st_as_sf(ipp.df, coords = c("lon", "lat"), crs = 4326)

# Get indices of points that are within the polygon
valid.pts <- sapply(st_intersects(ipp.sf, r.poly), function(x) length(x) > 0)

# Filter out invalid points
ipp.df <- filter(ipp.df, valid.pts)
# Subset df by presence
ipp.presence <- filter(ipp.df, presence == 1)
ipp.dummies <- filter(ipp.df, presence == 0)

# Convert the data to a ppp objects
locations <- ppp(ipp.presence$lon, 
                 ipp.presence$lat, 
                 poly=r.poly.list) 
dummy.locs <- ppp(ipp.dummies$lon, 
                  ipp.dummies$lat, 
                  poly=r.poly.list) 

# Create Quadscheme
Q <- quadscheme(locations, dummy.locs)
# Quadrature scheme (Berman-Turner)
# 789 data points, 858 dummy points
#   (provided manually)
#   Total weight 12.5619752974165


# Get list of distinct variable names
r.fs <- c(fs.df$var1, fs.df$var2) %>% 
  unique() %>% 
  sort()

# Define min/max scaling function for rasters
min.max.scale <- function(r, na.rm=T) {
  min.r <- min(values(r), na.rm=na.rm)
  max.r <- max(values(r), na.rm=na.rm)
  (r - min.r) / (max.r- min.r)
}


# Load/compute filtered & pre-processed rasters
covariates <- r.fs %>%
  set_names() %>%
  purrr::map(function(.x) {
    out <- r[[.x]] %>% 
      # Scale the data
      min.max.scale() %>%
      # Reproject to same CRS as point data
      terra::project(x=., 
                     y=st_crs(4326, parameters=T)$Wkt) %>%
      as.data.frame(xy=T) %>%
      filter(!is.na(!!.x)) %>%
      as.im()
  })


# Create formula
.f <- paste(fs.df$variable, collapse=" + ") %>% 
  paste("~", .) %>% 
  as.formula()

# GLM Model
fit.glm <- ppm(Q=Q, trend=.f, covariates=covariates, 
               rbord=.05, method="mpl")

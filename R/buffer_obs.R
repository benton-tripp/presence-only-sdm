# SETUP -------------------------------------------------------------------

# Load libraries
library(sf)
library(terra)
library(dplyr)
library(data.table)


# LOAD DATA ---------------------------------------------------------------

# Get raster data
states <- c("CO", "NC", "OR", "VT")
r.files <- paste0("data/final_rasters/", states, ".tif")
r.list <- purrr::map(r.files, rast)
names(r.list) <- states

# Get observation data
data.files <- list.files("data/final", full.names = T)
df <- purrr::map_df(data.files, readRDS) 


# Function to compute bounding box from centroid
compute.bbox <- function(x, y, half.res.x, half.res.y) {
  c(x - half.res.x, x + half.res.x, y - half.res.y, y + half.res.y)
}

# Function to generate a single POLYGON from the bounding box coordinates
make.polygon <- function(xmin, ymin, xmax, ymax) {
  m <- matrix(c(xmin, ymin, 
                xmax, ymin, 
                xmax, ymax, 
                xmin, ymax, 
                xmin, ymin), 
              ncol = 2, byrow = T)
  st_polygon(list(m))
}

get.grid.geoms <- function(r, .crs=NULL) {
  
  # Calculate the centroids of each cell
  centroids <- terra::xyFromCell(r, seq_len(ncell(r)))
  # Get resolution / 2 for x & y
  half.res.x <- res(r)[1] / 2
  half.res.y <- res(r)[2] / 2
  # Compute bounding box for each centroid
  bboxes <- t(apply(centroids, 1, function(pt) {
    compute.bbox(pt[1], pt[2], half.res.x, half.res.y)
  }))
  # Create dataframe
  dt <- as.data.table(bboxes)
  colnames(dt) <- c("xmin", "xmax", "ymin", "ymax")
  
  dt[, `:=` (
    # Add centroid lat/lon values to dataframe
    lon=centroids[, 1],
    lat=centroids[, 2],
    # Add i (row) and j (column) indices
    i=rowFromCell(r, 1:ncell(r)),
    j=colFromCell(r, 1:ncell(r))
  )]
  
  # Convert i and j to cell number
  dt[, cell := cellFromRowCol(r, i, j)]
  
  # Extract cell #'s of raster r where any values are NA
  r.na <- which(apply(as.data.frame(values(is.na(r))), MARGIN = 1, FUN = any))
  
  # Remove those cells from the data
  dt <- dt[-r.na]
  
  # Convert centroids to spatial points
  dt[, centroid := purrr::map2(lon, lat, ~st_point(cbind(.x, .y)))]
  
  # Create bounding box polygons for all rows and assign to geometry column
  dt[, bbox := purrr::pmap(.l=list(xmin, ymin, xmax, ymax), 
                           .f=make.polygon)]
  # Make spatial frame
  df <- st_sf(dt, 
              bbox = st_sfc(dt$bbox, crs=st_crs(r)),
              centroid = st_sfc(dt$centroid, crs=st_crs(r)))
  
  # Update CRS
  if (!is.null(.crs)) {
    df <- st_transform(df, st_crs(.crs)) %>%
      st_set_geometry("centroid") %>%
      st_transform(st_crs(.crs)) %>%
      st_set_geometry("bbox")
  }
  
  df %>% select(i, j, cell, bbox, centroid)
}

if (!file.exists("test/grid.df.rds")) {
  # Get grid data
  .grids <- purrr::map(r.list, ~get.grid.geoms(.x, .crs=4326))
  names(.grids) <- states
  
  # Add to spatial dataframe
  grid.df <- purrr::map_df(states, ~mutate(.grids[[.x]], state=.x))
  saveRDS(grid.df, "test/grid.df.rds")
} else {
  grid.df <- readRDS("test/grid.df.rds")
}

# Get observation data
obs.df <- df %>%
  select(state, common.name, observation.point=geometry)

# Confirm matching CRS
obs.df <- st_transform(obs.df, st_crs(grid.df))

if (!file.exists("test/joined.df.rds")) {
  # Join the observation points to grid cells based on spatial intersection
  joined.df <- st_join(obs.df, select(grid.df, -state), left = T)
  saveRDS(joined.df, "test/joined.df.rds")
} else {
  joined.df <- readRDS("test/joined.df.rds")
}

# Extract all i and j values from the joined.df that intersected with obs.df
intersected.grid <- unique(joined.df[, c("i", "j", "common.name", "state")]) %>% 
  as.data.frame()


# BUFFERING RASTER DATA ---------------------------------------------------

mask.update <- function(i, mask.raster, obs.df, obs.field="observation.point",
                        dist=10000, u="m") {
  obs.pt <- st_transform(obs.df[i, "observation.point"], st_crs(mask.raster))
  # Create a buffer around the point, ensuring correct units
  buf <- st_buffer(obs.pt, dist=units::as_units(paste(dist, u)))
  terra::rasterize(buf, mask.raster, update=T, field=1)
}

# For each observation point, you can now create a distance 
# raster and then mask cells within the buffer distance
get.buffered.zones <- function(r, obs.df, obs.field="observation.point",
                               dist=10000, u="m") {
  # Create an empty raster with the same extent and resolution as r
  mask.raster <- terra::rast(r)
  # Recursively update mask.raster with additional buffered regions
  purrr::walk(1:nrow(obs.df), function(i) {
    mask.raster <<- mask.update(i, mask.raster, obs.df, 
                                obs.field=obs.field, dist=dist, u=u)
  }, .progress=T)
  mask.raster
}

# ~1-2 minutes per 1k records
# obs.df.test <- obs.df.vt %>% slice(sample(1:nrow(.), 1e2)) 
# mask.test <- get.buffered.zones(r.list$VT, obs.df.test)
# plot(mask.test)

gc()

# Get masks by state, species
masks <- purrr::map(states, function(state) {
    specs <- sort(unique(obs.df$common.name))
    spec.masks <- purrr::map(specs, function(spec, st=state) {
      fname <- file.path("artifacts", "masks", paste0(st, "_", spec, ".tif"))
      if (file.exists(fname)) {
        cat("Reading", spec, st, "mask from", fname, "\n")
        r.mask <- rast(fname)
      } else {
        cat("Computing", spec, st, "mask, and saving to", fname, "\n")
        r.mask <- get.buffered.zones(r.list[[st]], 
                           filter(obs.df, state == st & common.name == spec))
        terra::writeRaster(r.mask, fname, overwrite=T)
      }
      r.mask
    })
    names(spec.masks) <- specs
    specs
  })
names(masks) <- states


# BUFFERING VECTOR DATA ---------------------------------------------------------------
# 
# cell.buffer <- function(pt, bboxes, dist=10000, u="m") {
#   
#   # Create a buffer around the point, ensuring correct units
#   buf <- st_buffer(pt, dist=units::as_units(paste(dist, u)))
#   # sqrt(st_area(buf) / pi) should be ~ 10000
#   
#   # Fix CRS to Web Mercator
#   buf <- st_transform(buf, 3857) # EPSG:3857 is Web Mercator
#   bboxes <- st_transform(bboxes, 3857)
#   
#   # Identify which bboxes intersect with the buffer
#   intersects.indices <- which(st_intersects(bboxes, buf, sparse=F))
#   # Select cells that intersect with the buffer
#   buf.cells <- bboxes[intersects.indices]
#   
#   # Union the selected cells into one polygon; Update CRS
#   out.poly <- st_union(buf.cells) %>%
#     st_transform(st_crs(pt))
#   
#   out.poly
# }
# 
# 
# obs.buffer.states <- purrr::map(unique(joined.df$state), function(st) {
#   obs.buffer.species <- purrr::map(unique(joined.df$common.name), function(spec) {
#     f <- paste0("tests/obs.buffer.", spec, ".", st, ".rds")
#     cat("Getting buffers for", spec, "at", st, "\n")
#     if (!file.exists(f)) {
#       spec.state.df <- joined.df %>% filter(state == st & common.name == spec)
#       state.grid.df <- grid.df %>% filter(state == st)
#       obs.buffer <- purrr::map(1:nrow(spec.state.df),
#                                function(i) {
#                                  pt <- spec.state.df[i, "observation.point"]
#                                  cell.buffer(pt, state.grid.df$bbox)
#                                }, .progress=T)
#       saveRDS(obs.buffer, f)
#     } else {
#       obs.buffer <- readRDS(f)
#     }
#     obs.buffer
#   })
#   names(obs.buffer.species) <- unique(joined.df$common.name)
# })
# names(obs.buffer.states) <- unique(joined.df$state)
# 
# 
# 

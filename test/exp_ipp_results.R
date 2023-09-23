ipp.models <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  tibble(
    common.name=spec,
    state=st,
    glm.path=file.path("artifacts/models/ipp_glm_mpl", 
                       paste0(spec, "_", st, "_ipp_glm_mpl.rds")),
    gam.path=file.path("artifacts/models/ipp_gam_mpl", 
                       paste0(spec, "_", st, "_ipp_gam_mpl.rds"))
  )
})

f <- readRDS(ipp.models[1, ]$glm.path)$fit

# Diagnostics
d <- diagnose.ppm(f, plot.it=F)

p.y <- ggplot(data=d$Ycts)



suppressPackageStartupMessages(
  suppressMessages(suppressWarnings(source("test/ipp_test_setup.R"))))

spec <- "Sanderling"
st <- "CO"
covariates.keep <- 50
Q <- readRDS(file.path("artifacts", "train_spatstat_Q_2",
                       paste0(st, "_", spec, "_Q.rds")))
Q.test <- readRDS(file.path("artifacts", "test_spatstat_Q_2",
                            paste0(st, "_", spec, "_Q.rds")))

fs.df <- var.imp %>% 
  filter(state == st & common.name == spec) %>%
  mutate(var1 = purrr::map_chr(variable, 
                               ~ stringr::str_split(.x, "\\:")[[1]][1]),
         var2 = purrr::map_chr(variable, ~ {
           split_result <- stringr::str_split(.x, "\\:")[[1]]
           if(length(split_result) > 1) split_result[2] else NA_character_
         })) %>%
  mutate(variable = ifelse(is.na(var2), 
                           var1, 
                           paste(var1, var2, sep = ":"))) %>%
  # Keep only pre-determined # of variables/interactions
  head(covariates.keep)

if (nrow(fs.df) > 0) {
  covariates <- c(fs.df$var1, fs.df$var2) %>% 
    unique() %>% 
    sort()
  
  # Load/compute filtered & pre-processed rasters
  covariates <- covariates %>%
    set_names() %>%
    purrr::map(function(.x) {
      file <- file.path("artifacts/spatstat_imgs", 
                        paste0(st, "_", .x, ".rds"))
      readRDS(file)
    })
  
  if (length(covariates) < covariates.keep) {
    covariates.keep <- length(covariates) # in case there are fewer available
  } else if (length(covariates) > covariates.keep) {
    covariates <- covariates[1:covariates.keep]
  }
  
  # No interactions
  .f <- paste(names(covariates), collapse = " + ") %>% 
    paste("~", .) %>% 
    as.formula()
}

fit.glm <- ppm(Q=Q, trend=.f, covariates=covariates, 
               rbord=.05, method="mpl") 
s <- summary(fit.glm)
imp.vars <-rownames(s$coefs.SE.CI[s$coefs.SE.CI$Ztest == "***",])[-1]
covariates <- covariates[imp.vars]
covariates.keep <- length(imp.vars)
.f <- paste(names(covariates), collapse = " + ") %>% 
  paste("~", .) %>% 
  as.formula()
fit.glm <- ppm(Q=Q, trend=.f, covariates=covariates, 
               rbord=.05, method="mpl") 

fit.glm$internal$glmfit$converged

inten <- predict(fit.glm)
pixarea <- with(inten, xstep * ystep)
prob <- inten * pixarea

locations.train <- data.table::rbindlist(
  list(
    data.table(x=Q$data$x, y=Q$data$y, obs=T), 
    data.table(x=Q$dummy$x, y=Q$dummy$y, obs=F)
  )
) 

covariates.rasters <- names(covariates) %>%
  set_names() %>%
  purrr::map(function(x) {
    r <- terra::rast(covariates[[x]])
    set.crs(r, st_crs(4326, parameters=T)$Wkt)
    names(r) <- x
    r
  })
purrr::walk2(names(covariates), covariates.rasters, function(n, r) {
  locations.train[, (n) := 
                    terra::extract(r, cbind(locations.train$x, 
                                            locations.train$y))]
})
glm.pred <- spatstat.model::predict.ppm(
  fit.glm, 
  locations=locations.train[, .(x,y)],
  type="trend", se=T)
glm.train <- cbind(
  locations.train,
  data.table(
    est=glm.pred$estimate,
    se=glm.pred$se)
)
glm.train[, `:=` (
  p.obs = est * pixarea
)]
optimal.threshold <- optimize.f1(glm.train)
cm <- get.acc(glm.train, optimal.threshold)
acc <- tibble(
  common.name=spec,
  state=st,
  covariate.count=covariates.keep,
  optimal.threshold=optimal.threshold 
) %>%
  cbind(as.list(c(cm$overall, cm$byClass)) %>% 
          as_tibble()) %>%
  select(common.name:Accuracy, Sensitivity, Specificity, F1)
cat("Train Results:\n Species:", spec, "\n",
    "State:", st, "\n",
    "Covariates:", covariates.keep, "\n",
    "Optimal Threshold:", optimal.threshold, "\n",
    "Accuracy:", acc$Accuracy, "\n",
    "F1:", acc$F1, "\n",
    "Sensitivity (TP Rate):", acc$Sensitivity, "\n",
    "Specificity (TN Rate):", acc$Specificity, "\n")

  
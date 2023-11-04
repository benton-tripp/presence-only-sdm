# Load variable importance from fitted LASSO models
lasso.model.path="artifacts/models/lasso_3_fs"

# Define min/max scaling function for rasters
min.max.scale <- function(x, na.rm=T) {
  min.x <- min(x, na.rm=na.rm)
  max.x <- max(x, na.rm=na.rm)
  (x - min.x) / (max.x - min.x)
}

get.var.imp <- function(st, spec, dir="artifacts/models/lasso_3_fs") {
  files <- list.files(dir) %>%
    .[grepl(paste(tolower(gsub(" ", "_", spec)), st, sep="_"), .)] %>%
    file.path(dir, .)
  if (length(files) == 0) return(
    tibble(
      common.name=character(0),
      state=character(0),
      variable=character(0),
      importance=numeric(0),
      wt=numeric(0),
      weighted.imp=numeric(0)
    )
  )
  var.imp <- purrr::map_df(files, ~{
    fit <- readRDS(.x)
    coef.df <- coef(fit$finalModel, s = fit$bestTune$lambda) %>%
      as.matrix() %>%
      as.data.frame()
    # Remove the intercept
    coef.df <- coef.df[-1, , drop = F]
    # Create a data frame of variable importance
    var.importance <- tibble(
      common.name = spec,
      state = st,
      variable = rownames(coef.df),
      importance = abs(coef.df[,1])
    ) %>%
      # Rank variables by importance
      arrange(state, common.name, -importance, variable) %>%
      # Only look at variables where imp. > 0
      filter(importance > 0)
  }) %>%
    mutate(n=n()) %>%
    group_by(common.name, state, variable, n, .drop=F) %>%
    summarize(importance=median(importance), 
              wt=n()) %>%
    ungroup() %>%
    mutate(wt=wt/n) %>%
    dplyr::select(-n) %>%
    mutate(weighted.imp=min.max.scale(wt*importance)) %>%
    arrange(-weighted.imp, variable) 
}

var.imp <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  get.var.imp(st, spec)
})
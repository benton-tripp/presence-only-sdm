# Setup -------------------------------------------------------------------

if (!file.exists("artifacts/summary_img.RData")) {
  
  library(data.table)
  library(dplyr)
  library(purrr)
  library(terra)
  library(caret)
  library(sf)
  library(flextable)
  library(officer)
  library(spatstat)
  library(htmltools)
  library(htmlwidgets)
  
  # Load pre-processed data
  source("R/load_preprocessed_data.R")
  # Load other utility function
  source("R/get_objects.R")
  
  model.scores <- purrr::map_df(c("ipp_glm_mpl_2", "maxent", 
                                  "logreg", "knn", "tree", 
                                  "randomforest", "xgboost"), function(m) {
                                    m.type <- case_when(
                                      m == "ipp_glm_mpl_2" ~ "ipp",
                                      m == "maxent" ~ "maxent",
                                      m == "logreg" ~ "logistic regression",
                                      m == "knn" ~ "knn",
                                      m == "tree" ~ "classification tree",
                                      m == "randomforest" ~ "random forest",
                                      m == "xgboost" ~ "xgboost",
                                      T ~ "")
                                    purrr::map_df(1:nrow(spec.state), function(i) {
                                      spec <- spec.state[i,]$species
                                      st <- spec.state[i,]$state
                                      results.path <- file.path(
                                        paste0("artifacts/test_results/", m),
                                        paste0(spec, "_", st, "_", m, ".rds"))
                                      readRDS(results.path)$test.accuracy
                                    }) %>% mutate(model.type=m.type)
                                  })
  
  
  reshape.data <- function(data, metric) {
    # Ensure that the metric is one of the expected columns
    if (!metric %in% c("Accuracy", "Sensitivity", 
                       "Specificity", "F1")) {
      stop("Invalid metric specified.")
    }
    # Reshape the data
    data %>%
      select(common.name, state, model.type, !!sym(metric)) %>%
      mutate(model.type = tolower(gsub(" ", ".", model.type))) %>%
      tidyr::pivot_wider(names_from = model.type, 
                         values_from = !!sym(metric),
                         names_sep=".",
                         names_prefix=paste0(tolower(metric), "."))
  }
  
  metrics <- c("Accuracy", "Sensitivity", "Specificity", "F1") %>%
    set_names() %>%
    purrr::map(~reshape.data(model.scores, .x))
  
  
  model.scores %>% setDT()
  model.scores[, `:=` (
    model.type = factor(model.type,
                        levels=unique(model.scores$model.type)),
    common.name = factor(common.name,
                         levels=unique(model.scores$common.name)),
    state = factor(state, 
                   levels=unique(model.scores$state))
  )]
  
  plt.data <- model.scores %>%
    melt(.,
         id.vars = c("common.name", "state", "model.type"), 
         measure.vars = c("Accuracy", "Sensitivity", 
                          "Specificity", "F1"),
         variable.name = "metric", 
         value.name = "metric.value")
  
  plt.1 <- ggplot(plt.data, 
                  aes(x = metric, y = metric.value, fill = model.type)) +
    geom_boxplot(outlier.shape = "o") +
    stat_summary(aes(group = model.type), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray"),
          plot.background = element_rect(fill = "white"),
          strip.placement = "inside") + 
    coord_flip() +
    facet_wrap(~metric, ncol=2, scales="free",
               strip.position="top")
  
  plt.2 <- states %>%
    set_names() %>%
    purrr::map(~ggplot(plt.data[state == .x], 
                  aes(x = metric, y = metric.value, fill = model.type)) +
    geom_boxplot(outlier.shape = "o") +
    stat_summary(aes(group = model.type), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray"),
          plot.background = element_rect(fill = "white"),
          strip.placement = "inside") + 
    coord_flip() +
    facet_wrap(~metric, ncol=2, scales="free",
               strip.position="top") +
      labs(title=paste("Model Metrics for", .x))
    )
  
  
  plt.3 <- species %>%
    set_names() %>%
    purrr::map(~ggplot(plt.data[common.name==.x], 
                  aes(x = metric, y = metric.value, fill = model.type)) +
    geom_boxplot(outlier.shape = "o") +
    stat_summary(aes(group = model.type), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray"),
          plot.background = element_rect(fill = "white"),
          strip.placement = "inside") + 
    coord_flip() +
    facet_wrap(~metric, ncol=2, scales="free",
               strip.position="top") +
      labs(title=paste("Model Metrics for", .x))
    )
  
  # The Kruskal-Wallis Rank Sum Test is a non-parametric alternative 
  # to one-way ANOVA. It's used to compare more than two groups and 
  # does not assume a normal distribution. The null hypothesis for this 
  # test is that all groups have the same median, with a threshold of 0.05.
  kwt <- kruskal.test(Accuracy ~ model.type, data = model.scores)
  kwt.spec <- kruskal.test(Specificity ~ model.type, data = model.scores)
  kwt.sens <- kruskal.test(Sensitivity ~ model.type, data = model.scores)
  kwt.f1 <- kruskal.test(F1 ~ model.type, data = model.scores)
  # Non-Parametric Pairwise Comparisons: To identify which specific 
  # model types differ from each other, conduct post-hoc pairwise 
  # comparisons using the non-parametric Dunn's test.
  dunn.res <- dunn.test::dunn.test(model.scores$Accuracy, 
                                   model.scores$model.type, 
                                   method="bonferroni")
  dunn.res.spec <- dunn.test::dunn.test(model.scores$Specificity, 
                                        model.scores$model.type, 
                                        method="bonferroni")
  
  dunn.res.sens <- dunn.test::dunn.test(model.scores$Sensitivity, 
                                        model.scores$model.type, 
                                        method="bonferroni")
  
  dunn.res.f1 <- dunn.test::dunn.test(model.scores$F1, 
                                      model.scores$model.type, 
                                      method="bonferroni")
  # Effect Size: To understand the magnitude of differences, not just 
  # the statistical significance, calculate Cliff's Delta for each pair 
  # of model types. This metric quantifies the difference in the likelihood 
  # that a randomly selected case from one group (model-type 1) scores higher 
  # than a randomly chosen case from another group (model-type 2), and vice 
  # versa. Essentially, Cliff's Delta is a measure of the 'success rate 
  # difference' between two groups, providing insight into the practical
  #  significance of their performance disparity.
  get.trtmnt <- function(comp) {
    stringr::str_split(comp, " - ")[[1]][[1]]
  }
  
  get.cntrl <- function(comp) {
    stringr::str_split(comp, " - ")[[1]][[2]]
  }
  
  
  get.sig.eff.size <- function(model.scores, dunn.res, meas="Accuracy") {
    pairs <- expand.grid(
      treatment=unique(model.scores$model.type),
      control=unique(model.scores$model.type)
    ) %>% 
      filter(treatment != control)
    
    eff.size <- purrr::map_df(1:nrow(pairs), ~{
      treatment <- pairs[.x,]$treatment
      control <- pairs[.x,]$control
      t.acc <- model.scores[[meas]][model.scores$model.type == treatment]
      c.acc <- model.scores[[meas]][model.scores$model.type == control]
      cd <- effsize::cliff.delta(t.acc, c.acc)
      data.table(
        treatment=treatment,
        control=control,
        delta=cd$estimate,
        lower=cd$conf.int[[1]],
        upper=cd$conf.int[[2]],
        variance=cd$var,
        magnitude=cd$magnitude
      )
    })
    
    eff.m <- reshape2::acast(eff.size, treatment ~ control, 
                             value.var='delta', margins=F)
    
    sig.efsz <- data.table(
      comp=dunn.res$comparisons,
      p.adj=dunn.res$P.adjusted
    ) %>%
      .[, .(treatment=sapply(comp, get.trtmnt), 
            control=sapply(comp, get.cntrl),
            p.adj)] %>%
      eff.size[., on=.(treatment, control)] %>%
      .[p.adj <= 0.05]
    
    .SDcols <- names(sig.efsz)[c(3:6,8)]
    sig.efsz <- cbind(sig.efsz[, 1:2], 
                      sig.efsz[, lapply(.SD, signif, 3), 
                               .SDcols=.SDcols]) %>%
      setnames(old=names(.), new=c("model.1", "model.2", 
                                   names(.)[3:ncol(.)]))
    
    list(
      mat=eff.m,
      significant.eff.size=sig.efsz,
      eff.size=eff.size,
      dunn.res=dunn.res
    )
    
  }
  
  
  res.acc <- get.sig.eff.size(model.scores, dunn.res, "Accuracy")
  res.sens <- get.sig.eff.size(model.scores, dunn.res.sens, "Sensitivity")
  res.spec <- get.sig.eff.size(model.scores, dunn.res.spec, "Specificity")
  res.f1 <- get.sig.eff.size(model.scores, dunn.res.f1, "F1")
  
  save.image(file="artifacts/summary_img.RData")
} else {
  read.image(file="artifacts/summary_img.RData")
}


# Tables -------------------------------------------------------------------

get.mean <- function(x, s=4, na.rm=T) signif(mean(x, na.rm=na.rm), digits=s)

acc.metrics.tbl <- model.scores %>%
  .[, .(acc=get.mean(Accuracy),
        sens=get.mean(Sensitivity),
        spec=get.mean(Specificity),
        f1=get.mean(F1)),
    by=model.type] %>%
  flextable(cwidth=c(1.5, rep(.75, 4)))

acc.tbl <- res.acc$significant.eff.size %>% 
  .[, p.adj := signif(p.adj, digits=2)] %>% 
  flextable(cwidth=c(1.25, 1.25, rep(.75, 4), .9))
sens.tbl <- res.sens$significant.eff.size %>% 
  .[, p.adj := signif(p.adj, digits=2)] %>% 
  flextable(cwidth=c(1.25, 1.25, rep(.75, 4), 1.05))
spec.tbl <- res.spec$significant.eff.size %>% 
  .[, p.adj := signif(p.adj, digits=2)] %>% 
  flextable(cwidth=c(1.25, 1.25, rep(.75, 4), .9))
f1.tbl <- res.f1$significant.eff.size %>% 
  .[, p.adj := signif(p.adj, digits=2)] %>% 
  flextable(cwidth=c(1.25, 1.25, rep(.75, 4), .9))

# Summary Plots -------------------------------------------------------------------


plt.1
plt.2
plt.3

record.corr.plt <- function(m) {
  corrplot::corrplot(m, diag = F, is.corr=F, tl.col="black")
  recordPlot()
}

svg(filename="docs/summary_plots/acc_mat.svg")
acc.mat.plt <- record.corr.plt(res.acc$mat)
dev.off()

svg(filename="docs/summary_plots/sens_mat.svg")
sens.mat.plt <- record.corr.plt(res.sens$mat)
dev.off()

svg(filename="docs/summary_plots/spec_mat.svg")
spec.mat.plt <- record.corr.plt(res.spec$mat)
dev.off()

svg(filename="docs/summary_plots/f1_mat.svg")
f1.mat.plt <- record.corr.plt(res.f1$mat)
dev.off()


# Prediction Plots --------------------------------------------------------

get.ipp.pred.map <- function(fit, locations) {
  # Covariates
  covariates <- fit$covariates
  covariates.rasters <- names(covariates) %>%
    set_names() %>%
    purrr::map(function(x) {
      r <- terra::rast(covariates[[x]])
      set.crs(r, st_crs(4326, parameters=T)$Wkt)
      names(r) <- x
      r
    })
  
  purrr::walk2(names(covariates), covariates.rasters, function(n, r) {
    locations[, (n) := 
                terra::extract(r, cbind(locations$x, 
                                        locations$y))]
  })
  
  EW <- spatstat.model::predict.ppm(
    fit, 
    locations=locations,
    type="trend")
  
  pred <- locations[, .(x, y, trend=EW)] %>%
    terra::rast()
}

get.maxent.pred.map <- function(fit, locations, rasts) {
  covariates <- attr(fit, "presence") %>% names()  
  covariates.rasters <- covariates %>%
    set_names() %>%
    purrr::map(~rasts[[.x]])
  purrr::walk2(covariates, covariates.rasters, function(n, r) {
    locations[, (n) := 
                terra::extract(r, cbind(locations$x, 
                                        locations$y))]
  })
  pred <- predict(fit, locations)
  
  pred <- locations[, .(x, y, prob=pred)] %>%
    terra::rast()
}

get.pred.map <- function(fit, locations, rasts) {
  covariates <- attr(fit$terms, "dataClasses") %>% 
    names() %>%
    .[. %in% names(rasts)]
  covariates.rasters <- covariates %>%
    set_names() %>%
    purrr::map(~rasts[[.x]])
  purrr::walk2(covariates, covariates.rasters, function(n, r) {
    locations[, (n) := 
                terra::extract(r, cbind(locations$x, 
                                        locations$y))]
  })
  pred <- predict(fit, locations, type="prob")$Presence
  
  pred <- locations[, .(x, y, prob=pred)] %>%
    terra::rast()
}

pred.rast.path <- "artifacts/prediction_rasters"
pred.rasters <- purrr::map(c("ipp_glm_mpl_2", "maxent", 
                       "logreg", "knn", "tree", 
                       "randomforest", "xgboost"), function(m) {
                         m.type <- case_when(
                           m == "ipp_glm_mpl_2" ~ "ipp",
                           m == "maxent" ~ "maxent",
                           m == "logreg" ~ "logistic regression",
                           m == "knn" ~ "knn",
                           m == "tree" ~ "classification tree",
                           m == "randomforest" ~ "random forest",
                           m == "xgboost" ~ "xgboost",
                           T ~ "")
                         species %>%
                           set_names %>%
                           purrr::map(function(spec) {
                             states %>%
                               set_names %>%
                               purrr::map(function(st) {
                                 cat(m.type, st, spec, "\n")
                                 r.path <- file.path(pred.rast.path, 
                                                     paste0(spec, "_", st, "_", 
                                                            m.type, ".tif"))
                                 if (!file.exists(r.path)) {
                                   # Regular expression pattern to match the filename
                                   if (m.type == "ipp") {
                                     pattern <- paste0("^", spec, "_", st, 
                                                       "_([0-9]{1,2})_ipp_glm_mpl_2\\.rds$")
                                   } else {
                                     pattern <- paste0("^", spec, "_", st, "_", m, "\\.rds$")
                                   }
                                   mpath <- list.files(paste0("artifacts/models/", m), 
                                                       full.names = T,
                                                       pattern=pattern)
                                   
                                   fit <- tryCatch({readRDS(mpath)}, 
                                                   error=function(e) {
                                                     cat("Failed to load model for", 
                                                         spec, st, "\n")
                                                     NULL})
                                   rasts <- rasters[[st]]
                                   r <- rasts[[1]]
                                   locations <- terra::as.data.frame(r, xy = T) %>%
                                     as.data.table() %>%
                                     .[, .(x, y)]
                                   
                                   if (m.type == "ipp" & !is.null(fit)) {
                                     # Intensity
                                     pred.raster <- get.ipp.pred.map(fit, locations)
                                   } else if (m.type == "maxent" & !is.null(fit)) {
                                     pred.raster <- get.maxent.pred.map(fit, locations, rasts)
                                   } else if (!is.null(fit)) {
                                     pred.raster <- get.pred.map(fit, locations, rasts)
                                   } else {
                                     pred.raster <- NULL
                                   }
                                   if (!is.null(pred.raster)) {
                                     terra::writeRaster(pred.raster, r.path)
                                   }
                                 } else {
                                   pred.raster <- terra::rast(r.path)
                                 }
                                 pred.raster
                               })
                           })
                       })
names(pred.rasters) <- c("ipp", "maxent", "logreg", "knn", "tree", "rf", "xgb")

pred.rasters$rf$Sanderling$NC %>% plot()

# Presence/Absence Plots ---------------------------------------------------------


# Covariate Plots ----------------------------------------------------------------


# Save Outputs -------------------------------------------------------------------

doc <- read_docx() %>%
  body_add_par("Species Distribution Modeling Benchmark Study Summary and Results") %>%
  body_add_par("\n\n") %>%
  body_add_par("Overview", style="heading 1") %>%
  body_add_par("\n\n") %>%
  body_add_par("Data Overview", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("Pseudo-Absence Selection", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("Results", style="heading 1") %>%
  body_add_par("\n\n") %>%
  body_add_par("Appendix", style="heading 1") %>%
  body_add_par("\n\n") %>%
  # Overall Metrics
  body_add_par("Overall Metric Summaries", style="heading 2") %>%
  body_add_break() %>%
  body_add_par("Average Metrics by Model Type", 
                        style="Table Caption") %>%
  body_add_flextable(acc.metrics.tbl) %>% 
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type", 
                        style="Image Caption") %>%
  body_add_gg(plt.1) %>%
  body_add_break() %>%
  # State Summary Plots
  body_add_par("Metric Summary Plots by State", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("Distribution of Metrics by Model Type - CO", 
                        style="Image Caption") %>%
  body_add_gg(plt.2$CO) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - NC", 
                        style="Image Caption") %>%
  body_add_gg(plt.2$NC) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - OR", 
                        style="Image Caption") %>%
  body_add_gg(plt.2$OR) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - VT", 
                        style="Image Caption") %>%
  body_add_gg(plt.2$VT) %>%
  body_add_break() %>%
  # Species Summary Plots
  body_add_par("Metric Summary Plots by Species", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("Distribution of Metrics by Model Type - Belted Kingfisher",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Belted Kingfisher`) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Cedar Waxwing",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Cedar Waxwing`) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Downy Woodpecker",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Downy Woodpecker`) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Ruddy Duck",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Ruddy Duck`) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Sanderling",
                        style="Image Caption") %>%
  body_add_gg(plt.3$Sanderling) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Sandhill Crane",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Sandhill Crane`) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Sharp-shinned Hawk",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Sharp-shinned Hawk`) %>%
  body_add_break() %>%
  body_add_par("Distribution of Metrics by Model Type - Wild Turkey",
                        style="Image Caption") %>%
  body_add_gg(plt.3$`Wild Turkey`) %>%
  body_add_break() %>%
  # Accuracy
  body_add_par("Model Accuracy Evaluation", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("Accuracy by Model Type - Significant Pairwise Comparisons and Effect Size", style="Table Caption") %>%
  body_add_flextable(acc.tbl) %>%
  body_add_par("\n\n") %>%
  body_add_par("Accuracy by Model Type - Pairwise Comparisons and Effect Size", 
               style="Image Caption") %>%
  body_add_img("docs/summary_plots/acc_mat.svg", 
                        width=5.5, height=5.5) %>%
  body_add_break() %>%
  # Sensitivity
  body_add_par("Model Sensitivity Evaluation", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("Sensitivity by Model Type - Significant Pairwise Comparisons and Effect Size",
               style="Table Caption") %>%
  body_add_flextable(sens.tbl) %>%
  body_add_par("\n\n") %>%
  body_add_par("Sensitivity by Model Type - Pairwise Comparisons and Effect Size", 
               style="Image Caption") %>%
  body_add_img("docs/summary_plots/sens_mat.svg", 
                        width=5.5, height=5.5) %>%
  body_add_break() %>%
  # Specificity
  body_add_par("Model Specificity Evaluation", style="heading 2") %>%
  body_add_par("\n\n") %>%
  # (N/A - there are 0 significant observations)
  # body_add_par("Specificity by Model Type - Significant Pairwise Comparisons and Effect Size", style="Table Caption") %>%
  # body_add_flextable(spec.tbl) %>%
  body_add_par("Specificity by Model Type - Pairwise Comparisons and Effect Size", 
               style="Image Caption") %>%
  body_add_img("docs/summary_plots/spec_mat.svg", 
                        width=5.5, height=5.5) %>%
  body_add_break() %>%
  # F1 Score
  body_add_par("Model F1 Score Evaluation", style="heading 2") %>%
  body_add_par("\n\n") %>%
  body_add_par("F1 Score by Model Type - Significant Pairwise Comparisons and Effect Size",
               style="Table Caption") %>%
  body_add_flextable(f1.tbl) %>%
  body_add_par("\n\n") %>%
  body_add_par("F1 Score by Model Type - Pairwise Comparisons and Effect Size", 
               style="Image Caption") %>%
  body_add_img("docs/summary_plots/f1_mat.svg", 
                        width=5.5, height=5.5) %>%
  body_add_par("\n\n") 
  
  
print(doc, target="docs/FinalSummary.docx")







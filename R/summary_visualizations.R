# Setup -------------------------------------------------------------------

if (!file.exists("artifacts/summary_img.RData")) {
  
  library(data.table)
  library(dplyr)
  library(purrr)
  library(terra)
  library(caret)
  library(sf)
  
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
  
  plt.2 <- ggplot(plt.data, 
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
    facet_wrap(~metric + state, ncol=2, scales="free",
               strip.position="top")
  
  
  plt.3 <- ggplot(plt.data, 
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
    facet_wrap(~metric + common.name, ncol=2, scales="free",
               strip.position="top")
  
  # The Kruskal-Wallis Rank Sum Test is a non-parametric alternative to one-way ANOVA. 
  # It's used to compare more than two groups and does not assume a normal distribution.
  # The null hypothesis for this test is that all groups have the same median, with a
  # threshold of 0.05.
  kwt <- kruskal.test(Accuracy ~ model.type, data = model.scores)
  # Non-Parametric Pairwise Comparisons: To identify which specific model types differ 
  # from each other, conduct post-hoc pairwise comparisons using the 
  # non-parametric Dunn's test.
  dunn.res <- dunn.test::dunn.test(model.scores$Accuracy, 
                                   model.scores$model.type, 
                                   method="bonferroni")
  # Effect Size: To understand the magnitude of differences, not just the statistical 
  # significance, calculate Cliff's Delta for each pair of model types. This metric 
  # quantifies the difference in the likelihood that a randomly selected case from 
  # one group (model-type 1) scores higher than a randomly chosen case from another 
  #  group (model-type 2), and vice versa. Essentially, Cliff's Delta is a measure 
  #  of the 'success rate difference' between two groups, providing insight into
  #   the practical significance of their performance disparity.
  pairs <- expand.grid(
    treatment=unique(model.scores$model.type),
    control=unique(model.scores$model.type)
  ) %>% 
    filter(treatment != control)
  
  eff.size <- purrr::map_df(1:nrow(pairs), ~{
    treatment <- pairs[.x,]$treatment
    control <- pairs[.x,]$control
    t.acc <- model.scores$Accuracy[model.scores$model.type == treatment]
    c.acc <- model.scores$Accuracy[model.scores$model.type == control]
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
  
  eff.m <- reshape2::acast(eff.size, treatment ~ control, value.var='delta', margins=F)
  plt.4 <- corrplot::corrplot(eff.m, diag = F)
  
  get.trtmnt <- function(comp) {
    stringr::str_split(comp, " - ")[[1]][[1]]
  }
  
  get.cntrl <- function(comp) {
    stringr::str_split(comp, " - ")[[1]][[2]]
  }
  
  sig.efsz <- data.table(
    comp=dunn.res$comparisons,
    p.adj=dunn.res$P.adjusted
  ) %>%
    .[, .(treatment=sapply(comp, get.trtmnt), 
          control=sapply(comp, get.cntrl),
          p.adj)] %>%
    eff.size[., on=.(treatment, control)] %>%
    .[p.adj <= 0.05]
  
  save.image(file="artifacts/summary_img.RData")
} else {
  read.image(file="artifacts/summary_img.RData")
}


# Table -------------------------------------------------------------------

get.mean <- function(x, s=4, na.rm=T) signif(mean(x, na.rm=na.rm), digits=s)

model.scores %>%
  .[, .(acc=get.mean(Accuracy),
        sens=get.mean(Sensitivity),
        spec=get.mean(Specificity),
        f1=get.mean(F1)),
    by=model.type] %>%
  flextable::flextable()
  
corrplot::corrplot(eff.m, diag = F, is.corr=F, tl.col="black")

.SDcols <- names(sig.efsz)[c(3:6,8)]
cbind(sig.efsz[, 1:2], sig.efsz[, lapply(.SD, signif, 3), .SDcols=.SDcols]) %>%
  setnames(old=names(.), new=c("model.1", "model.2", names(.)[3:ncol(.)])) %>%
  flextable::flextable()


kwt.spec <- kruskal.test(Specificity ~ model.type, data = model.scores)
dunn.res.spec <- dunn.test::dunn.test(model.scores$Specificity, 
                                 model.scores$model.type, 
                                 method="bonferroni")
kwt.sens <- kruskal.test(Sensitivity ~ model.type, data = model.scores)

dunn.res.sens <- dunn.test::dunn.test(model.scores$Sensitivity, 
                                 model.scores$model.type, 
                                 method="bonferroni")

dunn.res.f1 <- dunn.test::dunn.test(model.scores$F1, 
                                      model.scores$model.type, 
                                      method="bonferroni")

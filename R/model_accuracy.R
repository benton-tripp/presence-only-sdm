
get.f1 <- function(locations.df, thresh=0.5, method="F1") {
  metrics <- locations.df[, .(pred=factor(ifelse(p.obs > thresh, T, F), 
                                          levels=c("FALSE", "TRUE")), 
                              obs=factor(obs, levels=c("FALSE", "TRUE")))] %>% 
    confusionMatrix(reference=.$pred, data=.$obs, 
                    positive = "TRUE", 
                    mode="everything") %>%
    as.list() %>%
    .[["byClass"]] %>%
    as.list()
  if (method == "F1") {
    return(metrics$F1)
  } else if (method == "SS") {
    # Penalize heavily if either is near zero or NA
    if (is.na(metrics$Sensitivity) || is.na(metrics$Specificity)) {
      return(-1e6) 
    } else if (metrics$Sensitivity < .05 || metrics$Specificity < .05) {
      return(-1e6) 
    }
    
    # Combine Specificity and Sensitivity (e.g., geometric mean)
    combined.metric = sqrt(metrics$Sensitivity * metrics$Specificity)
    return(combined.metric)
  }
}

optimize.f1 <- function(locations.df) {
  objective.fn <- function(thresh) {
    # Negative because we want to maximize
    -get.f1(locations.df=locations.df, thresh=thresh, method="SS")  
  }
  opt.result <- optimize(objective.fn, lower=0.05, upper=0.95)
  return(opt.result$minimum)  # Return the optimal threshold
}

get.acc <- function(test, thresh) {
  df <- test[, .(pred=factor(ifelse(p.obs > thresh, T, F),
                             levels=c(F,T)), 
                 obs=factor(obs, levels=c(F,T)))]
  cm <- confusionMatrix(df$pred, df$obs, positive = "TRUE", mode="everything")
}

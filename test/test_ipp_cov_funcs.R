

# This script loops through each function in r.funcs, applies it to the sample data (sample_x and sample_y), and prints the results. This way, you can see if there are any errors or unexpected values being produced by any of the functions.

# Generate some sample data for testing purposes
sample_x <- df$lon[1:1000]
sample_y <- df$lat[1:1000]
sample_data <- data.frame(x = sample_x, y = sample_y)

# Testing each function in r.funcs
for(func_name in names(r.funcs)) {
  cat("Testing function for:", func_name, "\n")
  
  # Try to evaluate the function and print its results
  # We use `purrr::map2` to loop through each pair of x and y values
  results <- purrr::map2_dbl(sample_x, sample_y, 
                             ~r.funcs[[func_name]](.x, .y, r.ipp, func_name))
  
  cat("Results for", func_name, ":\n")
  cat(sum(is.na(results)), "NA values\n")
}

# git clone git@github.com:spatstat/spatstat.model.git
devtools::install("spatstat.model")
library(spatstat.model)


evalfxy <- function(f, x, y, extra) {
  if(length(extra) == 0)
    return(f(x,y))
  ## extra arguments must be matched explicitly by name
  ok <- names(extra) %in% names(formals(f))
  z <- do.call(f, append(list(x,y), extra[ok]))
  return(z)
}


vf <- lapply(covariates[isfun], evalfxy, x=x, y=y,
             extra=append(covfunargs, markinfo))


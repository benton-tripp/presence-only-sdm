



# Adapted from `spatstat.explore::lohboot`
# https://github.com/spatstat/spatstat.explore/blob/2757eb3d5801c2a4e6cd0beb091c5b10f311d0dc/R/lohboot.R


# X: 
#   A point pattern (object of class "ppp").
# 
# confidence: 
#   Confidence level, as a fraction between 0 and 1.
# nx,ny: 
#   Integers. If block=TRUE, divide the window into nx*ny rectangles.
# nsim:
#   Number of bootstrap simulations.
# type: 
#   Integer. Type of quantiles. Argument passed to `quantile.default` controlling the way the
#   quantiles are calculated. (R `stats` generic function `quantile` has 9 types, also
#   defaulting to 7).

blok.boot.loh <- function(X,
                          r,
                          .formula,
                          r.funcs,
                          covfunargs,
                          method="logi",
                          confidence=0.95,
                          nx = 4, ny = nx,
                          nsim=200,
                          type=7,
                          ...) {
  ## Correct block bootstrap as described by Loh.
  stopifnot(is.ppp(X))
  
  # validate confidence level
  stopifnot(confidence > 0.5 && confidence < 1)
  alpha <- 1 - confidence
  probs <- 1-alpha
  rank <- nsim * probs
    
  ## TODO: Set up for ppm model
  ## Fit ppm model
  f <- ppm(Q=.formula, covariates=r.funcs, covfunargs=covfunargs, method=method)
  coefs <- f$coef
  
  ## parse edge correction info
  correction <- f$correction
  
  switch(correction,
         none = {
           ckey <- clab <- "un"
           cadj <- "uncorrected"
         },
         border = {
           ckey <- "border"
           clab <- "bord"
           cadj <- "border-corrected"
         },
         translate = {
           ckey <- clab <- "trans"
           cadj <- "translation-corrected"
         },
         isotropic = {
           ckey <- clab <- "iso"
           cadj <- "Ripley isotropic corrected"
         })

  # first n columns are the local stats for the n points of X
  n <- npoints(X)
  y <- pred$v 
  nr <- nrow(y)

  W <- spatstat.geom::Window(X)
  ## Divides window into rectangular quadrats and returns as a tesselation
  grid.tess <- spatstat.geom::quadrats(boundingbox(W), nx=nx, ny=ny)
  ## Classify points of X into grid tiles
  block.idx <- spatstat.geom::tileindex(X$x, X$y, grid.tess)
  ## Use only 'full' blocks
  if (!is.rectangle(W)) {
    blocks <- tiles(grid.tess)
    fullblocks <- sapply(blocks, is.subset.owin, B = W)
    if (sum(fullblocks) < 2) {
      stop("Not enough blocks are fully contained in the window", call.=F)
    }
    warning(paste("For non-rectangular windows,", 
                  "only blocks fully contained in the window are used:",
                  paste(sum(fullblocks), "were used and",
                        sum(!fullblocks),
                        "were ignored.")
    ), call.=F)
    ## adjust classification of points of X
    indexmap <- cumsum(fullblocks)
    indexmap[!fullblocks] <- NA 
    block.idx <- indexmap[block.idx]
    ## adjust total number of points 
    n <- sum(!is.na(block.idx))
    block.factor <- factor(block.idx, levels=unique(indexmap[!is.na(indexmap)]))
  } else {
    block.factor <- factor(block.idx)
  }
  nmarks <- length(levels(block.factor))
  
  ## Average the local function values in each block
  ymarks <- by(t(y), block.factor, colSums, na.rm=T, simplify=F)
  ## Ensure empty data yield zero
  if(any(isempty <- sapply(ymarks, is.null))) {
    ymarks[isempty] <- rep(list(numeric(nr)), sum(isempty))
  }
  ymarks <- as.matrix(do.call(cbind, ymarks)) * nmarks/n
  ## average all the marks
  ymean <- .rowMeans(ymarks, na.rm=T, nr, nmarks)
  ## Average the marks in each block
  ystar <- matrix(NA, nrow=nr, ncol=nsim)
  for(i in 1:nsim) {
    ## resample nblocks blocks with replacement
    ind <- sample( nmarks , replace=T)
    ## average their local function values
    ystar[,i] <-  .rowMeans(ymarks[,ind], nr, nmarks, na.rm=T)
  }
  
  ## compute quantiles of deviation
  ydif <- sweep(ystar, 1, ymean)
  ydev <- apply(abs(ydif), 2, max, na.rm=T)
  crit <- quantile(ydev, probs=probs, na.rm=T, type=type)
  hilo <- rbind(ymean - crit, ymean + crit)
  
  ## create fv object
  df <- data.frame(r=f$r,
                   theo=theo,
                   ymean,
                   lo=hilo[1L,],
                   hi=hilo[2L,])
  
  colnames(df)[3L] <- ckey
  CIlevel <- paste(100 * confidence, "%% confidence", sep="")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            paste(cadj, "estimate of %s"),
            paste("lower", CIlevel, "limit for %s"),
            paste("upper", CIlevel, "limit for %s"))
} 

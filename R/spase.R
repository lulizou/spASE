#' Main function for detecting spatial allele-specific expression (ASE) patterns
#'
#' Fits a binomial Gaussian process GAM model for each gene,
#' extracts approximate goodness-of-fit (chi-squared) test p-values, and stores
#' predicted estimates and standard errors for visualization.
#'
#' @description This function fits a spatial binomial model, estimates a
#' smoothed probability surface, and stores information for plotting results.
#' Input must include either a list of two gene x pixel count data matrices or
#' a data frame containing a sparse representation of that, and a data frame or
#' matrix of pixel IDs and spatial coordinates. Optionally, single cell types
#' may be included for each pixel.
#'
#' @param data a data object containing gene expression counts for both alleles
#' (types). Accepted format is a
#'  \code{SpatialExperiment} object with exactly two assays, one for each
#'  allele. The assays can be sparse matrices and should be named by allele/type.
#'  By default, the probability of the first assay will be modeled as the
#'  binomial proportion p. Rows should be genes and columns should be pixel IDs.
#'  spatialCoords should optionally contain the cell types corresponding to each
#'  pixel.
#'
#' @param predict (optional) a \code{data.frame} or \code{matrix} containing the
#' coordinates at which to evaluate the model. If \code{NULL} then a default
#' set of 400 coordinates is automatically generated from the range of the data.
#'
#' @param model a string specifying which model to do, either 'all', 'cell.type',
#' or 'cell.type.interaction'. 'all' is the default and just uses all data to fit
#' the spatial model. 'cell.type' includes cell type as a covariate; each cell
#' type gets its own mean (intercept), 'cell.type.interaction' includes
#' interaction terms between cell type and the spatial spline, which allows each
#' cell type to have its own intercept and spatial fit. Typically challenging
#' to fit the latter two models with very sparse data.
#'
#' @param min.pixels numeric specifying the minimum number of pixels a gene
#' should be present on to fit
#'
#' @param min.k numeric, sets the minimum $k$ value for fitting the GAM model.
#' Default is \code{100}.
#'
#' @param max.k numeric, sets the maximum $k$ value for fitting the GAM model.
#' Default is \code{300}.
#'
#' @param cores number of cores to use for parallelization. Default is 1.
#'
#' @param verbose whether or not to print a lot of status messages. Default is
#' FALSE.
#'
#' @return A list containing following output:
#' \itemize{ \item{\code{summaryStats}}{ a data frame containing a summary of
#' the results including approximate goodness of fit (chi-squared) p-values and
#' effect sizes of cell types if applicable}
#' \item{\code{fittedVals}}{ a matrix containing the fitted values for each
#' gene} \item{\code{seVals}}{ a matrix containing the standard errors of the
#' fitted values for each gene}}
#'
#' @import SpatialExperiment
#' @import doSNOW
#' @import mgcv
#'
#' @export

spase <- function(se,
                  predict = NULL, model = 'all',
                  min.pixels = 100,
                  min.k = 100, max.k = 300, cores = 1, verbose = F) {

  if(!is(se, 'SpatialExperiment')) {
    stop('input data must be SpatialExperiment with 2 assays')
  }
  if (length(assays(se)) != 2) {
    stop('SpatialExperiment must have exactly 2 assays, one for each allele/type')
  }
  se <- se[rowSums(assays(se)[[1]]+assays(se)[[2]]!=0)>min.pixels,]
  if (verbose) {
    message(paste('filtered to', length(se), 'genes based on min.pixels >=',
                  min.pixels))
  }
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  genes <- rownames(se)
  alleles <- names(assays(se))
  colnames(spatialCoords(se)) <- c('pixel','x','y')
  pb <- txtProgressBar(max = length(genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  result <- foreach(
    i = 1:length(genes),
    .combine = function(...) {
      mapply('rbind', ..., SIMPLIFY = FALSE)
    },
    .multicombine = TRUE,
    .init = list(data.frame(), matrix(), matrix()),
    .options.snow = opts) %dopar% {
      # make gene data frame
      geneDf <- data.frame(pixel = colnames(se),
                           a1 = unname(assays(se)[[alleles[1]]][i,]),
                           a2 = unname(assays(se)[[alleles[2]]][i,]))
      geneDf <- geneDf[which(geneDf$a1+geneDf$a2 > 0),]
      geneDf <- merge(geneDf, spatialCoords(se), by='pixel', all.x=T)
      REML <- r <- 1:10*10
      for (j in 1:length(r)) {
        g <- gam(cbind(a1, a2) ~ s(x, y, bs = 'gp', m = c(2,r[j]),
                                   k = min.k, data = geneDf,
                                   family = 'binomial', method = 'REML'))
        REML[j] <- mt$gcv.ubre
      }
      REML.min <- r[which(REML == min(REML))]
      g <- gam(cbind(a1, a2) ~ s(x, y, bs = 'gp', m = c(2,REML.min),
                                 k = min.k, data = geneDf,
                                 family = 'binomial', method = 'REML'))
      k.ind <- k.check(g)[[3]] # get the k-index
      k.current <- min.k
      while (k.ind < 1) { # if k-index less than 1, try higher k value
        k.current <- k.current+100
        if (k.current > (nrow(geneDf)-2)) {
          k.current <- k.current-100
          break
        }
        if (k.current < max.k) {
          break
        }
        try_g <- gam(cbind(a1, a2) ~ s(x, y, bs = 'gp', m = c(2,REML.min),
                                       k = k.current, data = geneDf,
                                       family = 'binomial', method = 'REML'))
        k.ind.new <- k.check(try_g)[[3]]
        if ((k.ind.new > k.ind) & (summary(try_g)$edf > summary(g)$edf+1)) {
          k.ind <- k.ind.new
          g <- try_g
        } else {
          break
        }
      }

    }
}

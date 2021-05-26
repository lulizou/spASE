#' Main function for detecting ASE in single-cell data
#'
#' Fits a beta-binomial model for each gene and returns a results data frame.
#'
#' @description This function fits a beta-binomial model and  returns estimates
#' of p and confidence intervals
#'
#' @param matrix1 a matrix of counts where  the rows are genes and the columns
#' are cells for allele 1. This one is the one that gets its probability
#' modeled. Must  have  row  names and column names to identify genes and cells.
#'
#' @param matrix2 a matrix of counts where the rows are genes and the columns
#' are cells for allele 2.
#'
#' @param min.cells numeric specifying the minimum number of cells a gene
#' should be present on to fit. Default is 10 cells.
#'
#' @param cores number of cores to use for parallelization. Default is 1.
#'
#' @param verbose whether or not to print a lot of status messages. Default is
#' FALSE.
#'
#' @return A data frame of results containing the gene names, estimated p and
#' standard errors on the logit scale.
#'
#' @import aod
#'
#' @export

scase <- function(matrix1, matrix2,
                  min.cells = 10, cores = 1, verbose = F) {

  if(!is(matrix1, 'matrix') | !is(matrix2, 'matrix')) {
    stop('input matrix1 and matrix2 must be a matrix')
  }
  if (ncol(matrix1)!=ncol(matrix2)) {
    stop('matrices dont have same number of columns')
  }
  if (nrow(matrix1) != nrow(matrix2)) {
    stop('matrices dont have same number of rows')
  }
  if (is.null(rownames(matrix1))) {
    warning('input matrix1 doesnt have rownames')
  }
  if (is.null(colnames(matrix1))) {
    warning('input matrix1 doesnt have colnames')
  }
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  genes <- rownames(matrix1)
  pb <- txtProgressBar(max = length(genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  result <- foreach(
    i = 1:length(genes),
    .combine = function(...) {
      mapply('rbind', ..., SIMPLIFY = FALSE)
    },
    .options.snow = opts) %dopar% {
      y <- matrix1[i,]
      total <- y + matrix2[i,]
      y <- y[total > 0]; total <- total[total > 0]
      n <- length(y)
      if (all(y==total)) {
        data.frame(gene = genes[i], logit.p = NA, logit.p.sd = NA,
                   phi = NA, phi.sd = NA, flag = 'monoallelic1')
      } else if (sum(y)==0) {
        data.frame(gene = genes[i], logit.p = NA, logit.p.sd = NA,
                   phi = NA, phi.sd = NA, flag = 'monoallelic2')
      }
      conv_issue <-  FALSE
      tryCatch({
        fit <- betabin(cbind(Y,total-Y)~1, ~1, data = data.frame(total=total, Y=Y))
      }, warning = function(w) { conv_issue <<- TRUE},
      error = function(e) { conv_issue <<- TRUE}
      )
      if (conv_issue) {
        data.frame(gene = genes[i], logit.p = NA, logit.p.sd = NA,
                   phi = NA, phi.sd = NA, flag = 'convergence')
      }
      fitsum <- summary(fit)
      data.frame(gene = genes[i], logit.p = fitsum@Coef$Estimate,
                 logit.p.sd = fitsum@Coef$`Std. Error`,
                 phi = fitsum@Phi$Estimate, phi.sd = fitsum@Phi$`Std. Error`,
                 flag = '')
    }
  return(result)
}

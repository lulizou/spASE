#' Main function for estimating ASE in single-cell data
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
#' @param genes which genes to fit the model on, should be  a  subset of  the
#' rownames of matrix1, default is all of the genes
#'
#' @param add.var amount of additional variance  to add, default is 0. If
#' non-zero, addds additional columns to result data  frame suffixed with  .adj
#' to indicate calculated with the adjusted standard error
#'
#' @param verbose whether or not to print a lot of status messages. Default is
#' FALSE.
#'
#' @return A data frame of results containing the gene names, estimated p and
#' standard errors on the logit scale.
#'
#' @importFrom aod betabin
#' @import foreach
#' @import parallel
#' @import doSNOW
#'
#' @export
#'
#'

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

scase <- function(matrix1, matrix2, covariates=NULL,
                  min.cells = 10, cores = 1, genes = NULL, add.var=0,
                  verbose = F) {

  if(!(is(matrix1, 'matrix') | is(matrix1, 'sparseMatrix')) &
     !(is(matrix2, 'matrix') | is(matrix2, 'sparseMatrix'))) {
    stop('input matrix1 and matrix2 must be a matrix or sparseMatrix')
  }
  if (ncol(matrix1)!=ncol(matrix2)) {
    stop('matrices dont have same number of columns')
  }
  if (nrow(matrix1) != nrow(matrix2)) {
    stop('matrices dont have same number of rows')
  }
  if (is.null(rownames(matrix1))) {
    stop('input matrix1 doesnt have rownames')
  }
  if (is.null(colnames(matrix1))) {
    stop('input matrix1 doesnt have colnames')
  }
  if (!is.null(covariates)) {
    colnames(covariates)[1] <- 'cell'
    message(paste('assuming covariates first column is cell, using',
                  colnames(covariates[,-1]), 'as baseline covariates'))
    if (any(!is.factor(covariates[,-1]))) {
      stop('cannot handle non-factor covariates rn; either convert all covariates
           to factor or remove non-factor columns')
    }
  }
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  # Smart-seq uses NAs in the matrices which are different from 0's;
  # make sure don't use cells for which gene counts are NA in one allele
  matrix1[is.na(matrix2)] <- NA
  matrix2[is.na(matrix1)] <- NA
  numcells <- rowSums((matrix1>0)|(matrix2>0), na.rm=T)
  remove.idx <- which(numcells < min.cells)
  matrix1 <- matrix1[-remove.idx,]; matrix2 <- matrix2[-remove.idx,]
  if (is.null(genes)) {
    genes <- rownames(matrix1)
    message(paste(nrow(matrix1), 'genes pass min threshold of', min.cells, 'cells'))
  } else {
    message('using user-supplied genes')
  }
  pb <- txtProgressBar(max = length(genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  result <- foreach(
    i = 1:length(genes),
    .combine = rbind,
    .options.snow = opts) %dopar% {
      y <- matrix1[genes[i],]
      idx <- !is.na(y)
      y <- y[idx]
      total <- y + matrix2[genes[i],idx]
      y <- y[total > 0]; total <- total[total > 0]
      n <- sum(total)
      if (is.null(covariates)) {
        if (all(y==total)) {
          return(data.frame(gene = genes[i], totalUMI = n, totalCells = length(y),
                            logit.p = NA, logit.p.sd = NA,
                            phi = NA, phi.sd = NA, flag = 'monoallelic1'))
        } else if (sum(y)==0) {
          return(data.frame(gene = genes[i], totalUMI = n, totalCells = length(y),
                            logit.p = NA, logit.p.sd = NA,
                            phi = NA, phi.sd = NA, flag = 'monoallelic2'))
        }
        fit <- betabinase(y, total)
        res  <- fit$result
        if (!is.null(res)){
          return(data.frame(gene = genes[i], totalUMI = n, totalCells = length(y),
                            logit.p = unname(res@param['(Intercept)']),
                            logit.p.sd = sqrt(res@varparam[1,1]),
                            phi = unname(res@param['phi.(Intercept)']),
                            phi.sd = sqrt(res@varparam[2,2]),
                            flag = ''))
        } else {
          return(data.frame(gene = genes[i], totalUMI = n, totalCells = length(y),
                            logit.p = NA, logit.p.sd = NA,
                            phi = NA, phi.sd = NA, flag = 'conv'))
        }
      } else {
        if (all(y==total)) {
          return(data.frame(gene = genes[i], factor = NA, lvl = NA, totalUMI = n,
                            totalCells = length(y),logit.p = NA, logit.p.sd = NA,
                            phi = NA, phi.sd = NA, flag = 'monoallelic1'))
        } else if (sum(y)==0) {
          return(data.frame(gene = genes[i], factor = NA, lvl = NA, totalUMI = n,
                            totalCells = length(y),logit.p = NA, logit.p.sd = NA,
                            phi = NA, phi.sd = NA, flag = 'monoallelic2'))
        }
        covari <- left_join(data.frame(cell=names(y)), covariates, by='cell')
        print(covari)
        baseline.covari <- covari[,-1,drop=F]
        bcov.names <- colnames(baseline.covari)
        # fit one model to each cell type and return one entry per gene per level
        dfres <- NULL
        for (b in bcov.names) {
          lvls <- table(baseline.covari[,b])
          lvls <- names(lvls[lvls>min.cells])
          for (l in lvls) {
            ii <- which(baseline.covari[,b]==l)
            n <- sum(y[ii])
            ncell <- length(ii)
            fit <- betabinase(y[ii], total[ii])
            res  <- fit$result
            if (!is.null(res)) {
              if (is.null(dfres)) {
                dfres <- data.frame(gene = genes[i], factor = b, lvl = l, totalUMI = n,
                                    totalCells = ncell,
                                    logit.p = unname(res@param['(Intercept)']),
                                    logit.p.sd = sqrt(res@varparam[1,1]),
                                    phi = unname(res@param['phi.(Intercept)']),
                                    phi.sd = sqrt(res@varparam[2,2]), flag = '')
              } else {
                dfres <- rbind(dfres,
                               data.frame(gene = genes[i], factor = b, lvl = l, totalUMI = n,
                                          totalCells = ncell,
                                          logit.p = unname(res@param['(Intercept)']),
                                          logit.p.sd = sqrt(res@varparam[1,1]),
                                          phi = unname(res@param['phi.(Intercept)']),
                                          phi.sd = sqrt(res@varparam[2,2]), flag = ''))
              }
            }
          }
        }
        return(dfres)
      }
    }
  stopCluster(cl)
  result$p <- expit(result$logit.p)
  result$ci.low <- expit(result$logit.p - 2*result$logit.p.sd)
  result$ci.high <- expit(result$logit.p + 2*result$logit.p.sd)
  result$z <- result$logit.p/result$logit.p.sd
  result$pval <- pnorm(abs(result$z), lower.tail=F)
  result$qval <- p.adjust(result$pval, method='BH')
  if (add.var > 0) {
    result$logit.p.sd.adj <- sqrt(result$logit.p.sd^2 + add.var)
    result$ci.low.adj <- expit(result$logit.p-2*result$logit.p.sd.adj)
    result$ci.high.adj <- expit(result$logit.p+2*result$logit.p.sd.adj)
    result$z.adj <- result$logit.p/result$logit.p.sd.adj
    result$pval.adj <- pnorm(abs(result$z.adj), lower.tail=F)
    result$qval.adj <- p.adjust(result$pval.adj, method='BH')
  }
  return(result)
}

#' Main function for detecting spatial allele-specific expression (ASE) patterns
#'
#' Fits a beta-binomial thin plate regression spline for each gene.
#'
#' @description This function fits a spatial beta-binomial model, estimates a
#' 2D smoothed probability function, and stores information for plotting results.
#' Uses a likelihood ratio test to test for a significant spatial fit over
#' optional baseline covariates such as cell type.
#'
#' @param matrix1 matrix or sparseMatrix of allele counts for allele to be
#' modeled. Row names should be genes and column names should be pixel IDs.
#'
#' @param matrix2 matrix or sparseMatrix of allele counts for second allele.
#'
#' @param covariates data frame that matches pixel IDs to spatial coordinates
#' and other covariates. First column assumed to be pixel IDs, next two columns
#' assumed to be x,y coordinates. Any additional columns assumed to be
#' covariates that will be included in the baseline model.
#'
#' @param method string specifying which method to use. Default is
#' `betabinomial`. Other options are `quasibinomial` and `apeglm`.
#'
#' @param df integer, sets the number of degrees of freedom to use for the
#' smoothing spline. Default is 5. Usually want to increase this or test multiple.
#'
#' @param genes optional vector of genes to run on. Good for testing things out.
#'
#' @param min.pixels integer specifying the minimum number of pixels a gene
#' should be present on to fit. Default is 100.
#'
#' @param min.pixels.per.factor integer specifying the minimum number of pixels
#' a gene should be present on if there are factor covariates e.g. cell type.
#' For example, the default of 15 requires that at least 15 pixels are that cell
#' type, otherwise that cell type will not be included in the model for that gene.
#'
#' @param min.umi integer specifying the minimum number of UMI a gene should
#' have to fit. Default is 500.
#'
#' @param cores number of cores to use for parallelization. Default is 1.
#'
#' @param verbose whether or not to print a lot of status messages. Default is
#' FALSE.
#'
#' @return If method is betabinomial or quasibinomial, a list containing the
#' following output:
#' \itemize{ \item{\code{results}}{ a data frame containing a summary of
#' the results including p-values and q-values for goodness of spatial fit over
#' baseline covariates provided }
#' \item{\code{fits}}{ a list of beta-binomial model fit objects for each gene}}.
#'
#'
#'
#' @import aod
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @import Matrix
#' @importFrom dplyr left_join
#'
#' @export

spase <- function(matrix1, matrix2, covariates, method = 'betabinomial', df = 5,
                  genes = NULL, min.pixels = 100, min.pixels.per.factor = 15,
                  min.umi = 500, cores = 1, verbose = F) {

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
    warning('Warning: input matrix1 doesnt have rownames, will label genes
            using row index')
  }
  if (is.null(colnames(matrix1))) {
    warning('Warning: input matrix1 doesnt have colnames')
  }
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  numpixels <- rowSums((matrix1>0)|(matrix2>0), na.rm=T)
  numumi <- rowSums(matrix1)+rowSums(matrix2)
  remove.idx <- which((numpixels < min.pixels)|(numumi < min.umi))
  if (length(remove.idx)>0) {
    matrix1 <- matrix1[-remove.idx,]; matrix2 <- matrix2[-remove.idx,]
  }
  message(paste(nrow(matrix1), 'genes pass min threshold of', min.pixels,
                'pixels,', min.umi, 'UMI'))
  if (is.null(genes)) {
    if (is.null(rownames(matrix1))) {
      genes <- seq(1,nrow(matrix1))
    } else {
      genes <- rownames(matrix1)
    }
  } else {
    if (!all(genes %in% rownames(matrix1))) {
      stop('not all specified genes are in the matrix after thresholding')
    }
  }
  baseline.cov <- ''
  num.cov <- ncol(covariates)
  mydim <- 2
  if (num.cov<2) {
    stop('covariates should contain 1) pixel names, 2) and possibly 3), spatial coordinates')
  } else if (num.cov==2) {
    message(paste('found 2 columns in covariates; going to assume that first column is pixel names and second column is 1D coordinates'))
    colnames(covariates) <- c('pixel','x1')
    coord.mean <- mean(covariates$x1)
    coord.sd <- sd(covariates$x1)
    covariates$x1 <- (covariates$x1-coord.mean)/coord.sd
    mydim <- 1
  } else if (num.cov==3) {
    message(paste('found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates'))
    coord.mean <- apply(covariates[,2:3], 2, mean)
    coord.sd <- sd(as.matrix(covariates[,2:3]))
    covariates[,2:3] <- as.data.frame(sweep(covariates[,2:3],2,coord.mean,'-')/coord.sd)
    colnames(covariates)[1:3] <- c('pixel','x1','x2')
  } else if (num.cov>3) {
    coord.mean <- apply(covariates[,2:3], 2, mean)
    coord.sd <- sd(as.matrix(covariates[,2:3]))
    covariates[,2:3] <- as.data.frame(sweep(covariates[,2:3],2,coord.mean,'-')/coord.sd)
    colnames(covariates)[1:3] <- c('pixel','x1','x2')
    baseline.cov <- colnames(covariates)[-c(1:3)]
    message(paste('using', paste(baseline.cov, collapse=','),
                  'as baseline covariates'))
    baseline.cov.classes <- sapply(covariates, 'class')[-c(1:3)]
  }
  # get the spline matrix for all the coordinates
  if (mydim==1) {
    sm <- getSplineMatrix(covariates[,c('x1'),drop=F], df=df)
  } else if (mydim==2) {
    sm <- getSplineMatrix(covariates[,c('x1','x2')], df=df)
  }
  spline.mean <- sm[['mean']]
  spline.sd <- sm[['sd']]
  covariates <- cbind(covariates, sm[['X']])
  spline.idx <- (num.cov+1):ncol(covariates)
  if (length(genes)==0) {
    stop('no genes passed minimum thresholds. try lowering them')
  }
  pb <- txtProgressBar(max = length(genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  result <- foreach(
    i = 1:length(genes),
    .combine = function(l1, l2) {
      result <- rbind(l1$result, l2$result)
      fits <- c(l1$fits, l2$fits)
      return(list(df=df, result=result, fits=fits))
    },
    .init = list(df=df, data.frame(), list()),
    .packages = c('aod'),
    .options.snow = opts) %dopar% {
      y <- matrix1[genes[i],]
      total <- y + matrix2[genes[i],]
      present.idx <- which(!is.na(total) & total>0)
      y <- y[present.idx]; total <- total[present.idx]
      n <- sum(total)
      nspots <- length(total)
      covari <- left_join(data.frame(pixel=names(y)), covariates, by='pixel')
      # accounting for pixels w/o coordinates:
      if (mydim > 1) {
        remove.idx <- which(is.na(covari$x1)|is.na(covari$x2))
      } else {
        remove.idx <- which(is.na(covari$x1))
      }
      if (length(remove.idx)>0) {
        warning(paste(length(remove.idx), 'pixels did not have matching coordinates,
                      removing them, but you should probably double-check'))
        y <- y[-remove.idx]; total<-total[-remove.idx]
        covari <- covari[-remove.idx,]
      }
      el <- list()
      if (all(y==total)) {
        el[[genes[i]]] <- NA
        return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                      phi.baseline = NA, phi.full = NA,
                                      chisq.p = NA, flag = 'monoallelic1'), fits=el))
      } else if (sum(y)==0) {
        el[[genes[i]]] <- NA
        return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                      phi.baseline = NA, phi.full = NA,
                                      chisq.p = NA, flag = 'monoallelic2'), fits=el))
      }
      if (baseline.cov=='') {
        if (method == 'betabinomial') {
          baseline.model <- betabinase(y, total)
        } else if (method == 'quasibinomial') {
          baseline.model <- quasibinase(y, total)
        }
        mm <- NULL
        if (is.null(baseline.model$result)) {
          el[[genes[i]]] <- NA
          return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                        phi.baseline = NA, phi.full = NA,
                                        chisq.p = NA, flag = 'conv'), fits=el))
        }
      } else {
        # if covariate is a factor, make sure each level has enough data to fit
        # otherwise, remove the level. if no levels then baseline is the
        # intercept model.
        for (j in 1:length(baseline.cov)) {
          if (baseline.cov.classes[j]=='factor') {
            count.lvl <- table(covari[,baseline.cov[j]])
            not.enough <- which(count.lvl < min.pixels.per.factor)
            present.lvls <- names(count.lvl)[-not.enough]
            if (length(present.lvls)>1) {
              keep.idx <- which(covari[,baseline.cov[j]] %in% present.lvls)
              y <- y[keep.idx]; total <- total[keep.idx]
              covari <- covari[keep.idx,]
            }
          }
        }
        if (length(y) > min.pixels) {
          baseline.covari <- covari[,-c(1:3,spline.idx),drop=F]
          if (ncol(baseline.covari)>1) {
            mm <- model.matrix(~.-1, data=baseline.covari,
                               contrasts.arg=lapply(baseline.covari,contrasts,
                                                    contrasts=FALSE))
          } else {
            cov.name <- colnames(baseline.covari)[1]
            mm <- model.matrix(~.-1, data=baseline.covari,
                               contrasts.arg=list(cov.name=
                                                    diag(nlevels(baseline.covari[,1]))))
          }
          empties <- which(colSums(mm)==0)
          if (length(empties)>0) {
            mm <- mm[,-empties] # in case of empty levels
          }
          if (method == 'betabinomial') {
            baseline.model <- betabinase(y, total, mm)
          } else if (method == 'quasibinomial') {
            baseline.model <- quasibinase(y, total, mm)
          }
          if (is.null(baseline.model$result)) {
            el[[genes[i]]] <- NA
            return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                          phi.baseline = NA, phi.full = NA,
                                          chisq.p = NA, flag = 'conv'), fits=el))
          }
        } else {
          el[[genes[i]]] <- NA
          return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                        phi.baseline = NA, phi.full = NA,
                                        chisq.p = NA, flag = 'min pixels'), fits=el))
        }
      }
      not.intercept <- which(colSums(covari[,spline.idx]==1) != nrow(covari))
      if (is.null(mm)) {
        mm <- as.matrix(covari[,spline.idx])
      } else {
        mm <- as.matrix(cbind(data.frame(covari[,spline.idx][,not.intercept]),mm))
      }
      if (method == 'betabinomial') {
        smooth.model <- betabinase(y, total, mm)
      } else if (method == 'quasibinomial') {
        smooth.model <- quasibinase(y, total, mm)
      }
      if (is.null(smooth.model$result)) {
        el[[genes[i]]] <- NA
        return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                      phi.baseline = NA, phi.full = NA,
                                      chisq.p = NA, flag = 'conv'), fits=el))
      }
      if (method=='betabinomial') {
        lrt.stat <- 2*(smooth.model$result@logL-baseline.model$result@logL)
        df.diff <- baseline.model$result@df.residual-smooth.model$result@df.residual
      } else if (method == 'quasibinomial') {
        lrt.stat <- 2*(logLik(smooth.model$result)-logLik(baseline.model$result))
        df.diff <- baseline.model$result$df.residual-smooth.model$result$df.residual
      }
      pval <- 1 - pchisq(abs(lrt.stat), df.diff)
      el[[genes[i]]] <- smooth.model$result
      return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                    phi.baseline = baseline.model$result@random.param[[1]],
                                    phi.full = smooth.model$result@random.param[[1]],
                                    chisq.p = pval,flag = ''), fits=el))
    }
  result$result$qval <- p.adjust(result$result$chisq.p, method='BH')
  result$coord.mean <- coord.mean
  result$coord.sd <- coord.sd
  result$spline.mean <- spline.mean
  result$spline.sd <- spline.sd
  stopCluster(cl)
  return(result)

}

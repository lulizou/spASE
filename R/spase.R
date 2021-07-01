#' Main function for detecting spatial allele-specific expression (ASE) patterns
#'
#' Fits a binomial GAM for each gene.
#'
#' @description This function fits a spatial binomial model, estimates a
#' smoothed probability surface, and stores information for plotting results.
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
#' @return A list containing following output:
#' \itemize{ \item{\code{results}}{ a data frame containing a summary of
#' the results including p-values and q-values for goodness of spatial fit over
#' baseline covariates provided }
#' \item{\code{fits}}{ a list of beta-binomial model fit objects for each gene}}
#'
#' @import aod
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @importFrom mgcv smoothCon s
#' @import Matrix
#' @importFrom dplyr left_join
#'
#' @export

spase <- function(matrix1, matrix2, covariates, df = 5,
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
    warning('input matrix1 doesnt have rownames')
  }
  if (is.null(colnames(matrix1))) {
    warning('input matrix1 doesnt have colnames')
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
    genes <- rownames(matrix1)
  } else {
    if (!all(genes %in% rownames(matrix1))) {
      stop('not all specified genes are in the matrix after thresholding')
    }
  }
  baseline.cov <- ''
  num.cov <- ncol(covariates)
  if (num.cov<3) {
    stop('covariates should contain 1) pixel names, 2) and 3), spatial coordinates')
  } else if (num.cov>3) {
    baseline.cov <- colnames(covariates)[-c(1:3)]
    message(paste('using', paste(baseline.cov, collapse=','),
                  'as baseline covariates'))
    baseline.cov.classes <- sapply(covariates, 'class')[-c(1:3)]
  }
  colnames(covariates)[1:3] <- c('pixel','x1','x2')
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
      y <- y[total > 0]; total <- total[total > 0]
      n <- sum(total)
      nspots <- length(total)
      covari <- left_join(data.frame(pixel=names(y)), covariates, by='pixel')
      # accounting for pixels w/o coordinates:
      remove.idx <- which(is.na(covari$x1)|is.na(covari$x2))
      if (length(remove.idx)>0) {
        warning(paste(length(remove.idx), 'pixels did not have matching coordinates'))
        y <- y[-remove.idx]; total<-total[-remove.idx]
        covari <- covari[-remove.idx,]
      }
      el <- list()
      if (all(y==total)) {
        el[[genes[i]]] <- NA
        return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                      chisq.p = NA, flag = 'monoallelic1'), fits=el))
      } else if (sum(y)==0) {
        el[[genes[i]]] <- NA
        return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                      chisq.p = NA, flag = 'monoallelic2'), fits=el))
      }
      if (baseline.cov=='') {
        baseline.model <- betabinase(y, total)
        mm <- NULL
        if (is.null(baseline.model$result)) {
          el[[genes[i]]] <- NA
          return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
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
          baseline.covari <- covari[,-c(1:3),drop=F]
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
          mm <- mm[,-which(colSums(mm)==0)] # in case of empty levels
          baseline.model <- betabinase(y, total, mm)
          if (is.null(baseline.model$result)) {
            el[[genes[i]]] <- NA
            return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                          chisq.p = NA, flag = 'conv'), fits=el))
          }
        } else {
          el[[genes[i]]] <- NA
          return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                        chisq.p = NA, flag = 'min pixels'), fits=el))
        }
      }
      sm <- smoothCon(s(x1,x2,k=df,fx=T,bs='tp'),data=covari)[[1]]
      not.intercept <- which(colSums(sm$X==1) != nrow(sm$X))
      if (is.null(mm)) {
        mm <- as.matrix(data.frame(sm$X))
      } else {
        mm <- as.matrix(cbind(data.frame(sm$X[,not.intercept]),mm))
      }
      smooth.model <- betabinase(y, total, mm)
      if (is.null(smooth.model$result)) {
        el[[genes[i]]] <- NA
        return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                      chisq.p = NA, flag = 'conv'), fits=el))
      }
      lrt.stat <- 2*(smooth.model$result@logL-baseline.model$result@logL)
      df.diff <- df.residual(baseline.model$result)-
        df.residual(smooth.model$result)
      pval <- 1 - pchisq(abs(lrt.stat), df.diff)
      el[[genes[i]]] <- smooth.model$result
      return(list(result=data.frame(gene = genes[i], totalUMI = n, totalSpots = nspots,
                                    chisq.p = pval,flag = ''), fits=el))
    }

  result$result$qval <- p.adjust(result$result$chisq.p, method='BH')
  stopCluster(cl)
  return(result)
}

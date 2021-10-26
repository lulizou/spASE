#' Helper function to fit apeglm model while handling common warnings/errors
#'
#' Just wraps apeglm to return either the fit or warning description.
#'
#' @description This function fits a beta-binomial model with shrinkage
#' from the package apeglm
#' and wraps it so that warnings and/or errors are output.
#'
#' @param y vector of counts to be modeled
#'
#' @param total vector of total counts
#'
#' @param mod.mat model matrix that contains the design matrix of
#' interest.
#'
#'
#' @return A list the result, warnings, and errors
#'
#' @importFrom apeglm apeglm bbEstDisp
#'
#' @export

apeglmase <- function(y, total, mod.mat) {
  res <- list(result = NULL, warning = NULL, error = NULL)
  if (length(y) != length(total)) {
    stop('y and total dont have the same length')
  }
  withCallingHandlers({
    tryCatch({
      ase.cts <- as.matrix(y)
      cts <- as.matrix(total)
      if (is.null(mod.mat)) {
        X <- matrix(rep(1,ncol(cts)))
        coef <- 1
      } else {
        X <- mod.mat
        coef <- ncols(mod.mat)
      }
      theta.hat.0 <- 100 # rough estimate of dispersion
      param <- cbind(theta.hat.0, cts)
      fit.mle <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                        no.shrink=TRUE, log.link=FALSE, method="betabinC")
      theta.hat <- bbEstDisp(success=ase.cts, size=cts, x=X,
                             beta=fit.mle$map, minDisp=1, maxDisp=500)
      param <- cbind(theta.hat, cts)
      fit.mle <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                        no.shrink=TRUE, log.link=FALSE, method="betabinCR")
      theta.hat <- bbEstDisp(success=ase.cts, size=cts, x=X,
                             beta=fit.mle$map, minDisp=1, maxDisp=500)
      mle <- cbind(fit.mle$map, fit.mle$sd)
      param <- cbind(theta.hat, cts)
      fit.map <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                            coef=coef, mle=mle, log.link=FALSE, method = "betabinCR")
    })
  },
  warning = function(w) {
    res$warning <<- as.character(w$message)
    invokeRestart('muffleWarning')
  },
  error =  function(e) {
    res$error <<- as.character(e$message)
  })
  if (is.null(res$warning) & is.null(res$error)) {
    res$result <- fit.map
  }
  return(res)
}

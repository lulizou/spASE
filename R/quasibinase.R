#' Helper function to fit quasibinomial model while handling common warnings/errors
#'
#' Just wraps glm to return either the fit or warning description. Warnings
#' and errors occur when the model doesn't converge. These will be flagged for
#' removing downstream.
#'
#' @description This function fits a quasibinomial model from the package aod
#' and wraps it so that warnings and/or errors are output.
#'
#' @param y vector of counts to be modeled
#'
#' @param total vector of total counts
#'
#' @param mod.mat optional model matrix that contains the design matrix of
#' interest. Used for adding the thin plate regression spline basis functions
#' in the `spase` function and for adding cell type as a covariate. Note that
#' we assume there is an intercept column in mod.mat.
#'
#'
#' @return A list the result, warnings, and errors
#'
#' @importFrom aod betabin
#'
#' @export

quasibinase <- function(y, total, mod.mat=NULL) {
  res <- list(result = NULL, warning = NULL, error = NULL)
  if (length(y) != length(total)) {
    stop('y and total dont have the same length')
  }
  withCallingHandlers({
    tryCatch({
      if (is.null(mod.mat)) {
        fit <- glm(cbind(y, total-y)~1,
                       data = data.frame(y=y, total=total),
                   family = quasibinomial)
      } else {
        fit <- glm(cbind(y,total-y)~mod.mat-1,
                       data = data.frame(y=y, total=total),
                       family = quasibinomial)
      }
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
    res$result <- fit
  }
  return(res)
}

#' expit function
#'
#'
#' @description literally just the expit function
#'
#' @param x the number to expit
#' @return the expited version of x
#'
#' @export

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

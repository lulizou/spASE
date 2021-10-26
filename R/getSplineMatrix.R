#' Helper function to get normalized spline basis function matrix
#'
#' @description Takes in coordinates and degrees of freedom generates
#' the normalized spline basis function matrix using the `mgcv` `smoothCon`
#' function.
#'
#' @param coords data frame or matrix of coordinates (either 1D or 2D)
#'
#' @param df degrees of freedom to use for constructing the thin plate
#' regression spline. Default is 5.
#'
#'
#' @return normalized spline basis function matrix
#'
#' @importFrom mgcv smoothCon s
#'
#' @export

getSplineMatrix <- function(coords, df = 5) {
  if (ncol(coords) == 1) {
    colnames(coords) <- 'X'
    sm <- smoothCon(s(X,k=df,fx=T,bs='tp'),data=coords)[[1]]
  }
  if (ncol(coords) == 2) {
    colnames(coords) <- c('X','Y')
    sm <- smoothCon(s(X,Y,k=df,fx=T,bs='tp'),data=coords)[[1]]
  }
  if (ncol(coords) > 2) {
    stop('doesnt work for more than 2 coordinate columns')
  }
  mm <- as.matrix(data.frame(sm$X))
  X2 <- cbind(mm[,(df - 2):df], mm[,1:(df-3)]) #move intercept to front
  X2_mean <- apply(X2[,2:df],2, mean)
  X2_sd <- apply(X2[,2:df],2, sd)
  X2[,2:df] <- sweep(X2[,2:df], 2, X2_mean,  '-')
  X2[,2:df] <- sweep(X2[,2:df], 2, X2_sd, '/') #standardize
  colnames(X2) <- paste0('spline',seq(1,df))
  return(list(X=X2, mean=X2_mean, sd=X2_sd))
}

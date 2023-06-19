#' Main function for detecting spatial allele-specific expression (ASE) patterns
#' considering cell type mixtures.
#'
#' Fits the intercept-only model (only cell type effects) and the spatial
#' model (non-parametric thin plate smoothing splines)
#'
#' @description This function fits the spatial beta-binomial model with
#' cell type mixtures.
#'
#' @param myRCTD RCTD object to run on
#'
#' @param coords the coordinates of the pixels
#'
#' @param cell_types the cell types to test; if left as NULL, tries to test all
#'
#' @param df integer, sets the number of degrees of freedom to use for the
#' smoothing spline. Default is 5. Usually want to increase this or test multiple.
#'
#' @param genes optional vector of genes to run on. Good for testing things out.
#'
#' @param cores number of cores to use for parallelization. Default is 1.
#'
#'
#' @return If method is betabinomial or quasibinomial, a list containing the
#' following output:
#' \itemize{ \item{\code{results}}{ a data frame containing a summary of
#' the results including p-values and q-values for goodness of spatial fit over
#' baseline covariates provided }
#' \item{\code{fits}}{ a list of beta-binomial model fit objects for each gene}}.
#'
#'
#' @export

cside_spase <- function(myRCTD, covariates, df = 5,
                  min.pixels = 100, min.pixels.per.factor = 15,
                  min.umi = 500, cores = 1, verbose = F) {


}

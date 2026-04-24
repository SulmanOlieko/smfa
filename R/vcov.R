################################################################################
#                                                                              #
# R functions for the smfa package                                     #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Variance - Covariance Matrix of estimates                                    #
# Models:                                                                      #
#           -Stochastic metafrontier analysis                                  #
#           -Latent class stochastic metafrontier analysis                     #
#           -Sample selection correction for stochastic metafrontier model     #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Compute variance-covariance matrix of stochastic metafrontier models
#'
#' \code{\link{vcov}} computes the variance-covariance matrix of the maximum
#' likelihood (ML) coefficients from stochastic metafrontier models estimated with
#' \code{\link{smfa}}.
#'
#' @details The variance-covariance matrix is obtained by the inversion of the
#' negative Hessian matrix. Depending on the distribution and the
#' \code{'hessianType'} option, the analytical/numeric Hessian or the bhhh
#' Hessian is evaluated.
#'
#' @param object A stochastic metafrontier model returned
#' by \code{\link{smfa}}.
#' @param ... Currently ignored
#'
#' @name vcov
#'
#' @return The variance-covariance matrix of the maximum likelihood
#' coefficients is returned.
#'
#' @seealso \code{\link{smfa}}, for the stochastic metafrontier analysis model
#' fitting function using cross-sectional or pooled data.
#'
#' @keywords methods vcov
#'
# variance covariance matrix for smfa ----------
#' @rdname vcov
#' @aliases vcov.smfa
#' @importFrom stats vcov
#' @export
vcov.smfa <- function(object, ...) {
  if (!is.null(object$metaFrontierVcov)) {
    return(object$metaFrontierVcov)
  }
  # LP method: no metafrontier vcov; return a named list of group vcov's
  lapply(object$groupModels, function(m) m$invHessian)
}

################################################################################
#                                                                              #
# R functions for the smfa package                                     #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Coefficients extraction                                                      #
# Models:                                                                      #
#           -Stochastic metafrontier analysis                                  #
#           -Latent class stochastic metafrontier analysis                     #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract coefficients of stochastic metafrontier models
#'
#' @description
#' From an object of class \code{'summary.smfa'},
#' \code{\link{coef}} extracts the coefficients,
#' their standard errors, z-values, and (asymptotic) P-values.
#'
#' From on object of class \code{'smfa'}, it extracts
#' only the estimated coefficients.
#'
#' @name coef
#'
#' @param object A stochastic metafrontier model returned by \code{\link{smfa}},
#' or an object of class \code{'summary.smfa'}.
#' @param ... Currently ignored.
#'
#' @return For objects of class \code{'summary.smfa'},
#' \code{\link{coef}} returns a matrix with four columns. Namely, the
#' estimated coefficients, their standard errors, z-values,
#' and (asymptotic) P-values.
#'
#' For objects of class \code{'smfa'}, \code{\link{coef}}
#' returns a numeric vector of the estimated coefficients.
#'
#' @seealso \code{\link{smfa}}, for the stochastic metafrontier analysis model
#' fitting function using cross-sectional or pooled data.
#'
#' @keywords methods coefficients
#'

# coefficients from smfa ----------
#' @rdname coef
#' @aliases coef.smfa
#' @importFrom stats coef
#' @export
coef.smfa <- function(object, ...) {
  if (!is.null(object$metaFrontierParam)) {
    return(object$metaFrontierParam)
  }
  # LP method: no metafrontier params; return a list of group params
  lapply(object$groupModels, function(m) m$mlParam[seq_len(m$nXvar)])
}

# coefficients from summary.smfa ----------
#' @rdname coef
#' @aliases coef.summary.smfa
#' @export
coef.summary.smfa <- function(object, ...) {
  object$metaRes
}

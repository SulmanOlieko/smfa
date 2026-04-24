################################################################################
#                                                                              #
# R functions for the smfa package                                     #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Fitted values of models                                                      #
# Models:                                                                      #
#           -Stochastic metafrontier analysis                                  #
#           -Latent class stochastic metafrontier analysis                     #
#           -Sample selection correction for stochastic metafrontier model     #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract fitted values of stochastic metafrontier models
#'
#' \code{\link{fitted}} returns the fitted frontier values from stochastic
#' metafrontier models estimated with \code{\link{smfa}}.
#'
#' @param object A stochastic metafrontier model returned
#' by \code{\link{smfa}}.
#' @param ... Currently ignored.
#'
#' @name fitted
#'
#' @return A vector of fitted values is returned.
#'
#' @note The fitted values are ordered in the same way as the corresponding
#' observations in the dataset used for the estimation.
#'
#' @seealso \code{\link{smfa}}, for the stochastic metafrontier analysis model
#' fitting function using cross-sectional or pooled data.
#'
#' @keywords methods fitted
#'
# fitted values for smfa (returns metafrontier fitted values) ----------
#' @rdname fitted
#' @aliases fitted.smfa
#' @importFrom stats fitted
#' @export
fitted.smfa <- function(object, ...) {
  object$dataTable$.mf_yhat_meta
}

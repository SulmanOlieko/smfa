################################################################################
#                                                                              #
# R functions for the smfa package                                     #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Information Criteria extraction                                              #
# Models:                                                                      #
#           -Stochastic metafrontier analysis                                  #
#           -Latent class stochastic metafrontier analysis                     #
#           -Sample selection correction for stochastic metafrontier model     #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract information criteria of stochastic metafrontier models
#'
#' \code{\link{ic}} returns information criterion from stochastic
#' metafrontier models estimated with \code{\link{smfa}}.
#'
#' The different information criteria are computed as follows: \itemize{ \item
#' AIC: \eqn{-2 \log{LL} + 2 * K} \item BIC: \eqn{-2 \log{LL} + \log{N} * K}
#' \item HQIC: \eqn{-2 \log{LL} + 2 \log{\left[\log{N}\right]} * K} } where
#' \eqn{LL} is the maximum likelihood value, \eqn{K} the number of parameters
#' estimated and \eqn{N} the number of observations.
#'
#' @name ic
#'
#' @param object A stochastic metafrontier model returned
#' by \code{\link{smfa}}.
#' @param IC Character string. Information criterion measure. Three criteria
#' are available: \itemize{ \item \code{'AIC'} for Akaike information criterion
#' (default) \item \code{'BIC'} for Bayesian information criterion \item
#' \code{'HQIC'} for Hannan-Quinn information criterion }.
#' @param ... Currently ignored.
#'
#' @return \code{\link{ic}} returns a data frame with the values of the 
#' information criteria (AIC, BIC and HQIC) of the maximum likelihood coefficients. 
#' If the \code{IC} argument is provided, it returns only the requested criterion 
#' as a numeric value.
#'
#' @seealso \code{\link{smfa}}, for the stochastic metafrontier analysis model
#' fitting function using cross-sectional or pooled data.
#'
#' @keywords methods AIC BIC HQIC
#'
#' @rdname ic
#' @aliases ic.smfa
#' @importFrom sfaR ic
#' @export
ic.smfa <- function(object, IC = NULL, ...) {
  aic <- -2 * object$mlLoglik + 2 * object$nParm
  bic <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
  hqic <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) * object$nParm
  
  if (is.null(IC)) {
    res <- data.frame(AIC = aic, BIC = bic, HQIC = hqic)
    return(res)
  } else {
    if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
      stop("Unknown information criteria: ", paste(IC), call. = FALSE)
    }
    obj <- switch(IC,
      AIC = aic,
      BIC = bic,
      HQIC = hqic
    )
    return(obj)
  }
}

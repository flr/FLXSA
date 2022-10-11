# xsa.sa.R - DESC
# /xsa.sa.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# xsa.sa {{{

#' An FLXSA estimator module for the mse package
#'
#' The `mse` package can use different stock assessment methods as modules in
#' the estimation step of a management procedure. This functions provides such
#' a module for FLXSA.
#' As a flag for convergence, the returned tracking FLQuant contains the number
#' of maximum number of iterations (maxit).
#' @param stk An FLStock.
#' @param idx An FLIndices.
#'
#' @return A list containing the updated FLStock, and the tracking FLQuant.
#' @examples
#' data(ple4)
#' data(ple4.index)
#' xsa.sa(stk=ple4, idx=ple4.index, args=list(ay=2018), tracking=FLQuant())

xsa.sa <- function(stk, idx, args, tracking, ...) {

	args0 <- list(...)
	
  args0$stock <- stk
	args0$indices <- idx
	
  if(is.null(args0$control)) args0$control <- FLXSA.control()
	
	fit <- do.call('FLXSA', args0)

	stk <- stk + fit

  track(tracking, "conv.est", ac(args$ay)) <- fit@control@maxit
	
  list(stk = stk, tracking = tracking)
} # }}}

# xsa.sa.R - DESC
# /xsa.sa.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# xsa.sa {{{

xsa.sa <- function(stk, idx, ...){
	args <- list(...)
	args$stock <- stk
	args$indices <- idx
	if(is.null(args$control)) args$control <- FLXSA.control()
	tracking <- args$tracking
	args$tracking <- NULL
	fit <- do.call('FLXSA', args)
	stk <- stk + fit
	tracking["convergence",ac(range(stk)["maxyear"]+1)] <- fit@control@maxit
	list(stk = stk, tracking = tracking)
} # }}}

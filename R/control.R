# FLXSA.control.R - DESC
# /FLXSA.control.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# FLXSA.control class {{{

#' A class for FLXSA control options 
#'
#' Runs of the FLXSA method require a number of parameter to be set. Most of them
#' change the behaviour of the solution-searching algorithm
#' Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque eleifend
#' odio ac rutrum luctus. Aenean placerat porttitor commodo. Pellentesque eget porta
#' libero. Pellentesque molestie mi sed orci feugiat, non mollis enim tristique. 
#' Suspendisse eu sapien vitae arcu lobortis ultrices vitae ac velit. Curabitur id 
#' 
#' @name FLXSA.control
#' @rdname FLXSA.control
#' @docType class
#' @aliases FLXSA.control-class
#'
#' @slot x An object of class FLXSA. If provided, the 'FLXSA.control' is
#' initialized with the corresponding values of an XSA analysis stored in the
#' object. This is useful for getting the same initial sloteters for
#' successive analyses. Specifying one or more of the other arguments
#' supersedes default values, or values obtained from this FLXSA object
#' @slot tol The covergence tolerance, i.e. difference between two successive
#' iterations must be lower, to declare convergence of the model.
#' @slot maxit The maximum number of iterations allowed
#' @slot min.nse The minimum value of SE permitted in estimate of N hat
#' @slot fse User set SE of F when shrinking to mean F
#' @slot rage The oldest age for which the two parameter model is used for
#' catchability at age. Note that this value should be one less than the value
#' used in the executable version of XSA
#' @slot qage The age after which catchability is no longer estimated. q at
#' older ages set to the value at this age
#' @slot shk.n If \code{TRUE}, shrinkage to mean N
#' @slot shk.f If \code{TRUE}, shrinkage to mean F
#' @slot shk.yrs The number of years to be used for shrinkage to F for
#' terminal year
#' @slot shk.ages The number of ages to be used for shrinkage to F for
#' terminal age
#' @slot window The time window to consider in the model
#' @slot tsrange The number of years to be used in the time series weighting
#' @slot tspower The power to be used in the time series taper weighting
#' @slot vpa If \code{FALSE} use cohort analysis, otherwise, use VPA
#'
#' @section Validity:
#'
#'   \describe{
#'     \item{VALIDITY}{Neque porro quisquam est qui dolorem ipsum.}
#' }
#'
#' You can inspect the class validity function by using
#'    \code{getValidity(getClassDef('FLCatch'))}
#'
#' @section Accessors:
#' All slots in the class have accessor and replacement methods defined that
#' allow retrieving and substituting individual slots.
#'
#' The values passed for replacement need to be of the class of that slot.
#' A numeric vector can also be used when replacing FLQuant slots, and the
#' vector will be used to substitute the values in the slot, but not its other
#' attributes.
#'
#' @section Constructor:
#' @author Laurence Kell & Philippe Grosjean
#' @seealso \link{FLComp}
#' @keywords classes
#' @examples
#'
#' 	# To create a new FLXSA.control object with default parameters:
#' 	my.xsa.control <- FLXSA.control()
#' 	my.xsa.control
#' 	# Same, but changing values of some parameters
#' 	my.xsa.control <- FLXSA.control(maxit=50, shk.f=FALSE)
#' 	my.xsa.control
#'

setClass("FLXSA.control",
	slots=c(
		tol          ="numeric",
		maxit        ="integer",
		min.nse      ="numeric",
		fse          ="numeric",
		rage         ="integer",
		qage         ="integer",
		shk.n        ="logical",
		shk.f        ="logical",
		shk.yrs      ="integer",
		shk.ages     ="integer",
		window       ="integer",
  	tsrange      ="integer",
  	tspower      ="integer",
  	vpa          ="logical"),
  prototype=prototype(
    tol          =as.double(10e-10),
  	maxit        =as.integer(30),
  	min.nse      =as.double(0.3),
  	fse          =as.double(0.5),
  	rage         =as.integer(-1),
  	qage         =as.integer(10),
  	shk.n        =TRUE,
  	shk.f        =TRUE,
  	shk.yrs      =as.integer(5),
  	shk.ages     =as.integer(5),
  	window       =as.integer(100),
  	tsrange      =as.integer(20),
  	tspower      =as.integer(3),
  	vpa          =FALSE)
)

setValidity("FLXSA.control", function(object) {

	if (object@tol <= 0)
		return("value of tol must be > 0")

	if (object@maxit <= 0)
		return("value of maxit must be > 0")

  if (object@min.nse < 0)
		return("value of min.nse must be > 0")

  if (object@fse < 0)
		return("value of fse must be > 0")

  if (object@rage < -1)
		return("value of rage must be >= -1")

  if (object@qage < 0)
		return("value of qage must be >= 0")

  if (object@shk.yrs < 0)
		return("value of shk.yrs must be >= 0")

  if (object@shk.ages < 0)
		return("value of shk.ages must be >= 0")

  if (object@window <= 5)
		return("value of window must be >= 5")                               

	# Everything is fine
	return(TRUE)
})

# }}}

# FLXSA.control() {{{

#' @rdname FLXSA.control
#' @aliases FLXSA.control

FLXSA.control <- function(x=NULL, tol=10e-10, maxit=30, min.nse=0.3,  
  fse=0.5, rage=0, qage=10, shk.n=TRUE, shk.f=TRUE, shk.yrs=5, shk.ages=5,
  window=100, tsrange=20, tspower=3, vpa=FALSE) {

	  if (is.null(x)){
		  res <- new("FLXSA.control", tol=tol, maxit=as.integer(maxit),
        min.nse=min.nse, fse=fse, rage=as.integer(rage), qage=as.integer(qage),
        shk.n=as.logical(shk.n)[1], shk.f=as.logical(shk.f)[1],
        shk.yrs=as.integer(shk.yrs), shk.ages=as.integer(shk.ages),
        window=as.integer(window), tsrange=as.integer(tsrange),
        tspower=as.integer(tspower), vpa=as.logical(vpa)[1])
	  } else {	# We reuse an FLXSA.control object embedded in an FLXSA object
        if (!is.null(x) & !(is(x, "FLXSA") | is(x, "FLXSA.control")))
    		  	stop("FLXSA must be an 'FLXSA' or an 'FLXSA.control' object!")

        if (is(x, "FLXSA"))
    		   res <- x@control
        if (is(x, "FLXSA.control"))
    		   res <- x
        if (is.null(x))
    		   res <- new("FLXSA.control")
        		    
 		   # ... and possibly redefine some of its parameters
       if (!missing(tol))
       		res@tol <- tol
       if (!missing(maxit))
    			res@maxit <- as.integer(maxit)
    	 if (!missing(min.nse))
    			res@min.nse <- min.nse
   		 if (!missing(fse))
    			res@fse <- fse
     	 if (!missing(rage))
    			res@rage <- as.integer(rage)
       if (!missing(qage))
       		res@qage <- as.integer(qage)
       if (!missing(shk.n))
       		res@shk.n <- as.logical(shk.n)[1]
       if (!missing(shk.f))
       		res@shk.f <- as.logical(shk.f)[1]
       if (!missing(shk.yrs))
        	res@shk.yrs <- as.integer(shk.yrs)
       if (!missing(shk.ages))
       		res@shk.ages <- as.integer(shk.ages)
       if (!missing(window))
       		res@window <- as.integer(window)
       if (!missing(tsrange))
       		res@tsrange <- as.integer(tsrange)
       if (!missing(tspower))
      		res@tspower <- as.integer(tspower)
       if (!missing(vpa))
      		res@vpa <- as.logical(vpa)[1]

  		# Verify that this object is valid
 	  	test <- validObject(res)

		  if (!test)
        stop("Invalid object:", test)
	    }

	return(res)
  }
# }}}

# is.FLXSA.control {{{
# Test if an object is of FLXSA.control class
is.FLXSA.control <- function(x)
	return(inherits(x, "FLXSA.control"))

## show (a replacement for print of S3 classes)
setMethod("show", signature(object="FLXSA.control"),
	function(object){
      n.<-slotNames(object)
	   for (i in 1:length(n.))
         cat(n.[i],"\t\t",slot(object,n.[i]),"\n")
	}
) # }}}

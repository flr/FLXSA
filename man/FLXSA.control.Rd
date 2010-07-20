\name{FLXSA.control}
\alias{FLXSA.control}
\title{ Create a new FLXSA.control object }
\description{
  This function creates a new control object required by XSA model (extended survivors analysis).
}
\usage{
FLXSA.control(x = NULL, tol = 1e-09, maxit = 30, min.nse = 0.3, fse = 0.5,
    rage = 0, qage = 10, shk.n = TRUE, shk.f = TRUE, shk.yrs = 5,
    shk.ages = 5, window = 100, tsrange = 20, tspower = 3, vpa = FALSE)
}
\arguments{
  \item{x}{ An object of class FLXSA. If provided, the 'FLXSA.control' is initialized
    with the corresponding values of an XSA analysis stored in the object. This is useful
    for getting the same initial parameters for successive analyses. Specifying one or
    more of the other arguments supersedes default values, or values obtained from this
    FLXSA object}
  \item{tol}{ The covergence tolerance, i.e. difference between two successive iterations
    must be lower, to declare convergence of the model.}
  \item{maxit}{ The maximum number of iterations allowed }
  \item{min.nse}{ The minimum value of SE permitted in estimate of N hat }
  \item{fse}{ User set SE of F when shrinking to mean F }
  \item{rage}{ The oldest age for which the two parameter model is used for catchability at age. Note that this value should be one less than the value used in the executable version of XSA }
  \item{qage}{ The age after which catchability is no longer estimated. q at older ages set
    to the value at this age }
  \item{shk.n}{ If \code{TRUE}, shrinkage to mean N }
  \item{shk.f}{ If \code{TRUE}, shrinkage to mean F }
  \item{shk.yrs}{ The number of years to be used for shrinkage to F for terminal year }
  \item{shk.ages}{ The number of ages to be used for shrinkage to F for terminal age }
  \item{window}{ The time window to consider in the model }
  \item{tsrange}{ The number of years to be used in the time series weighting }
  \item{tspower}{ The power to be used in the time series taper weighting }
  \item{vpa}{ If \code{FALSE} use cohort analysis, otherwise, use VPA }
}
\details{
  See \code{\link{FLXSA}} for more details about the method and some bibliographical references.
}
\value{
  An \code{FLXSA.control} object with same slots as the arguments of the function (except \code{FLXSA})
}

\author{ Laurence Kell & Philippe Grosjean }

\seealso{ \code{\link{FLXSA}} }
\examples{
	# To create a new FLXSA.control object with default parameters:
	my.xsa.control <- FLXSA.control()
	my.xsa.control
	# Same, but changing values of some parameters
	my.xsa.control <- FLXSA.control(maxit=50, shk.f=FALSE)
	my.xsa.control
}
\keyword{ classes }

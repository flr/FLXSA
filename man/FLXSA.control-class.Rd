\name{FLXSA.control-class}
\docType{class}
\alias{FLXSA.control-class}
\alias{assess,FLXSA.control-method}
\alias{show,FLXSA.control-method}

\title{Class "FLXSA.control"}
\description{control object for the XSA (eXtended Survivors Analysis) model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FLXSA.control", ...)} or through the
creator function \code{\link{FLXSA.control}}.
}
\section{Slots}{
	 \describe{
    \item{\code{tol}:}{ The covergence tolerance, i.e. difference between two successive
        iterations must be lower to declare convergence of the model. Class \code{"numeric"}}
    \item{\code{maxit}:}{The maximum number of iterations allowed. Class \code{"integer"}}
    \item{\code{min.nse}:}{The minimum value of SE permitted in estimate of N hat. Class 
        \code{"numeric"}}
    \item{\code{fse}:}{User set SE of F when shrinking to mean F. Class \code{"numeric"}}
    \item{\code{rage}:}{The oldest age for which the two parameter model is used for
        catchability at age. Class \code{"integer"}}
    \item{\code{qage}:}{The age after which catchability is no longer estimated. Class 
        \code{"integer"}}
    \item{\code{shk.n}:}{If \code{TRUE}, shrinkage to mean N. Class \code{"logical"}}
    \item{\code{shk.f}:}{If \code{TRUE}, shrinkage to mean F. Class \code{"logical"}}
    \item{\code{shk.yrs}:}{The number of years to be used for shrinkage to F for terminal year.
        Class \code{"integer"}}
    \item{\code{shk.ages}:}{The number of ages to be used for shrinkage to F for terminal age. 
        Class \code{"integer"}}
    \item{\code{window}:}{The time window to consider in the model fitting. Class
        \code{"integer"}}
    \item{\code{tsrange}:}{The number of years to be used in the time series weighting. 
        Class \code{"integer"} ~~ }
    \item{\code{tspower}:}{The power to be used in the time series taper weighting. Class
        \code{"integer"}}
    \item{\code{vpa}:}{If \code{FALSE} use cohort analysis, otherwise, use VPA. Class
        \code{"logical"}}
  }
}
\section{Methods}{
  \describe{
    \item{assess}{\code{signature(object = "FLXSA.control")}: Runs FLXSA according to the
        options defined in the object.}
	 }
}
\author{ Laurence Kell & Philippe Grosjean }
\seealso{ \code{\link{FLXSA}}, \code{\link{FLXSA-class}}
}
\keyword{classes}

\name{FLXSA-class}
\docType{class}
\alias{FLXSA-class}
\alias{show,FLXSA-method}

\title{Class FLXSA}
\description{
    A class for the results of an XSA analysis.
    }
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("FLXSA", ...)} and are output by
    calls to \code{\link{FLXSA}}.
}
\section{Slots}{
	 \describe{

    \item{\code{survivors}:}{Object of class \code{"FLQuant"}. Numbers-at-age in the final
        age of a cohort}
    \item{\code{se.int}:}{Object of class \code{"FLQuant"}. Internal standard errors of
        numbers-at-age in the final year }
    \item{\code{se.ext}:}{Object of class \code{"FLQuant"}. External standard errors of
        numbers-at-age in the final year }
    \item{\code{n.fshk}:}{Object of class \code{"FLQuant"}. Numbers-at-age in the final
        age of a cohort derived from shrinkage to the mean F }
    \item{\code{n.nshk}:}{Object of class \code{"FLQuant"}. Numbers-at-age in the final
        age of a cohort derived from shrinkage to the mean N }
    \item{\code{var.fshk}:}{Object of class \code{"FLQuant"}. Numbers-at-age in the final age of a cohort derived from shrinkage to the mean F }
    \item{\code{var.nshk}:}{Object of class \code{"FLQuant"}. Variance assumed for numbers-at-age in the final age of a cohort derived from shrinkage to the mean N }
    \item{\code{q.hat}:}{Object of class \code{"list"}. Estimates of q-at-age by cpue index }
    \item{\code{q2.hat}:}{Object of class \code{"list"}. Estimates of exponent in cpue power
        model }
    \item{\code{control}:}{Object of class \code{"FLXSA.control"}. FLXSA.control object used. }
    \item{\code{catch.n}:}{Object of class \code{"FLQuant"}. Estimated catch numbers-at-age,
        not estimated }
    \item{\code{stock.n}:}{Object of class \code{"FLQuant"}. Estimated stock numbers-at-age }
    \item{\code{harvest}:}{Object of class \code{"FLQuant"}. Estimated catch F-at-age }
    \item{\code{index.name}:}{Object of class \code{"character"}. Names of cpue indices }
    \item{\code{index}:}{Object of class \code{"list"}. cpue indices, adjusted to start of year }
    \item{\code{index.res}:}{Object of class \code{"list"}. cpue index residuals }
    \item{\code{index.hat}:}{Object of class \code{"list"}. predicted cpue indices }
    \item{\code{index.var}:}{Object of class \code{"list"}. Estimated variance of cpue index
        residuals }
    \item{\code{call}:}{A \code{"character"} representation of the call that created the object}
    \item{\code{desc}:}{A description on the object, of class \code{"character"}}
    \item{\code{diagnostics}:}{A data frame with various diagnostics from a run of XSA.}
  }
}
\section{Extends}{
Class \code{"FLAssess"}, directly.
}
\section{Methods}{
    All methods are inherited from FLAssess.
}
\references{ Darby, C. D., and Flatman, S. 1994. Virtual Population Analysis: version 3.1 (Windows/Dos) user guide. Info. Tech. Ser., MAFF Direct. Fish. Res., Lowestoft, (1): 85pp.

Shepherd, J.G. 1992. Extended survivors analysis: an improved method for the analysis of catch-at-age data and catch-per-unit-effort data. Working paper No. 11 ICES Multi-species Assessment Working Group, June 1992, Copenhagen, Denmark. 22pp. (mimeo).

Shepherd, J.G. 1994. Prediction of yearclass strength by calibration regression analysis of multiple recruit index series. ICES J. Mar. Sci. In Prep. 
}
\author{Laurence Kell}
\seealso{
	\code{\link{FLXSA}}
}
\examples{
}
\keyword{classes}

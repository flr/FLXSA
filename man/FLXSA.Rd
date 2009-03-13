\name{FLXSA}
\alias{FLXSA}
\title{ Create a new FLXSA object -run an XSA analysis- }
\description{
  This function runs an XSA (extended survivor analysis) and creates an FLXSA object used to analyse its results.
}
\usage{
 FLXSA(stock, indices, control = FLXSA.control(), desc, diag.flag=TRUE)
}
\arguments{
  \item{stock}{ An FLStock object to be used for the analysis }
  \item{indices}{ An FLIndices object holding the indices of abundance to consider in the model }
  \item{control}{ An \code{FLXSA.control} object giving parameters of the model (see \code{\link{FLXSA.control}}) }
  \item{desc}{ A short description of this analysis }
  \item{diag.flag}{ If TRUE returns all diagnostics, if FALSEonly returns stock.n, harvest and control}
}

\details{

  Virtual population analysis and cohort analysis are essentially accountancy methods whereby a stock's historical population structure may be reconstructed from total catch data given a particular level of natural mortality. Firstly, however, numbers at age in the last year and age have to be found since both methods iterate backwards down a cohort. The main problem in many sequential age based assessment methods is therefore to estimate these terminal population numbers. In XSA these are found from the relationship between catch per unit effort (CPUE), abundance and year class strength.

  Estimates of the catchability for the oldest age in an assessment, tuned by the ad hoc or XSA procedures, are directly dependent on the terminal population or F values used to initialise the underlying VPA. Catchability at the oldest age is therefore under-determined and cannot be utilised without additional information. Within the ad hoc tuning procedures the additional information is obtained by making the assumption that the exploitation pattern on the oldest ages was constant during the assessment time series. F on the oldest age in the final year is estimated as a proportion of an average of the F for preceding ages in the same year. XSA uses an alternative approach by making the assumption that fleet catchability is constant (independent of age) above a certain age. The age (constant for all fleets) is user-defined. For each fleet, the catchability value estimated at the specified age, is used to derive population abundance estimates for all subsequent ages in the fleet data set.

}

\value{

  An \code{FLXSA} object is returned, whith slots:
  \item{n      }{An FLQuant with the number of individuals at age}
  \item{f      }{An FLQuant with the fishing mortality}
  \item{swt    }{An FLQuant with the stock weight}
  \item{mat    }{An FLQuant with the maturity indices}
  \item{qres   }{A list with residuals for q}
  \item{cpue   }{A list with the various cpues}
  \item{wts    }{A list with the various weights}
  \item{control}{The \code{FLXSA.control} object that was used for this analysis}
  \item{call   }{A copy of the call to run this analysis}
  \item{desc   }{A description of the analysis}

}

\references{ Darby, C. D., and Flatman, S. 1994. Virtual Population Analysis: version 3.1 (Windows/Dos) user guide. Info. Tech. Ser., MAFF Direct. Fish. Res., Lowestoft, (1): 85pp.

Shepherd, J.G. 1992. Extended survivors analysis: an improved method for the analysis of catch-at-age data and catch-per-unit-effort data. Working paper No. 11 ICES Multi-species Assessment Working Group, June 1992, Copenhagen, Denmark. 22pp. (mimeo).

Shepherd, J.G. 1994. Prediction of yearclass strength by calibration regression analysis of multiple recruit index series. ICES J. Mar. Sci. In Prep. 
}

\author{ Laurence Kell and Philippe Grosjean }

\note{ See \code{\link{update}} to learn how to update stock data according to an XSA analysis }

\seealso{
    \code{\link{FLXSA.control}}, \code{\link[FLCore]{FLStock-class}}}

\examples{

#TO DO...

}

\keyword{ classes }

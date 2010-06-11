\name{data}
\alias{cod4}
\alias{cod4.indices}
\alias{cod4.stock}
\alias{cod4.f.ref}
\alias{cod4.n.ref}
\alias{her4}
\alias{her4.indices}
\alias{her4.stock}
\alias{her4.f.ref}
\alias{her4.n.ref}
\alias{ple7a}
\alias{ple7a.indices}
\alias{ple7a.stock}
\alias{ple7a.f.ref}
\alias{ple7a.n.ref}

\title{FLXSA data sets}
\description{
    Data sets for testing and validation of the FLXSA package
}
\details{
    The three datasets available in this package (cod4, her4 and ple7a) each consist
    of four objects.
    \describe{
	\item{*.indices}{An FLIndices object holding various indices of abundance for each stock.}
	\item{*.stock}{An stock object with the necessary catch estimates, biological parameters
        and previous stock estimates.}
	\item{*.xsa.control}{The FLXSA.control object for the reference run.}
	\item{*.xsa.ref}{The FLXSA object obtained from the reference run.}
}}
\seealso{
    \code{\link{FLXSA-class}},\code{\link{FLXSA.control-class}}
}
\author{Laurie Kell}
\keyword{data}

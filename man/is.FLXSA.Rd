\name{is.FLXSA}
\alias{is.FLXSA}
\alias{is.FLXSA.control}
\title{is.FLXSA}
\description{
    These two functions return \code{code} if objects are of class FLXSA and FLXSA.control, respectively.
}
\usage{
    is.FLXSA(x)
    is.FLXSA.control(x)
}
\arguments{
  \item{x}{An object to be tested}
}
\value{
    is.FLXSA returns \code{TRUE} if its argument is of class \code{FLXSA} (that is, has
        "FLXSA" amongst its classes) and FALSE otherwise.
    is.FLXSA.control returns \code{TRUE} if its argument is of class \code{FLXSA.control}
    (that is, has "FLXSA.control" amongst its classes) and FALSE otherwise.

}
\examples{
xsa <- FLXSA.control()
is.FLXSA.control(xsa)
}
\keyword{attribute}

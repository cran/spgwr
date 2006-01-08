\name{gwr.gauss}
\alias{gwr.gauss}
\title{GWR Gaussian weights function}
\description{
The function returns a matrix of weights using
the Gaussian scheme:

\deqn{w(g) = e^{{-(d/h)}^2}}

where \eqn{d} are the distances between the observations
and \eqn{h} is the bandwidth.

}
\usage{
gwr.gauss(dist2, bandwidth)
}
\arguments{
  \item{dist2}{matrix of squared distances between observations}
  \item{bandwidth}{bandwidth}
}
\value{
matrix of weights.
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2000,
Quantitative Geography, London: Sage; C. Brunsdon, A.Stewart Fotheringham
and M.E. Charlton, 1996, "Geographically Weighted Regression: A Method for
Exploring Spatial Nonstationarity", Geographical Analysis, 28(4), 281-298;
\url{http://ncg.nuim.ie/ncg/GWR/index.htm}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{gwr.sel}}, \code{\link{gwr}}}

\examples{
plot(seq(-10,10,0.1), gwr.gauss(seq(-10,10,0.1)^2, 3.5), type="l")
}


\keyword{spatial}

\name{gwr.bisquare}
\alias{gwr.bisquare}
\title{GWR bisquare weights function}
\description{
The function returns a vector of weights using
the bisquare scheme:

\deqn{w_{ij}(g) = (1 - (d_{ij}^2/d^2))^2 }
if \eqn{d_{ij} <= d} else \eqn{w_{ij}(g) = 0}, where \eqn{d_{ij}} 
are the distances between the observations
and \eqn{d} is the distance at which weights are set to zero.

}
\usage{
gwr.bisquare(dist2, d)
}
\arguments{
  \item{dist2}{vector of squared distances between observations}
  \item{d}{distance at which weights are set to zero}
}
\value{
matrix of weights.
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2000,
Quantitative Geography, London: Sage; C. Brunsdon, A.Stewart Fotheringham
and M.E. Charlton, 1996, "Geographically Weighted Regression: A Method for
Exploring Spatial Nonstationarity", Geographical Analysis, 28(4), 281-298;
\url{http://gwr.nuim.ie/}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{gwr.sel}}, \code{\link{gwr}}}

\examples{
plot(seq(-10,10,0.1), gwr.bisquare(seq(-10,10,0.1)^2, 6.0), type="l")
}
\keyword{spatial}

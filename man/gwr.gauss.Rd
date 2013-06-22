\name{gwr.gauss}
\alias{gwr.gauss}
\alias{gwr.Gauss}
\title{GWR Gaussian weights function}
\description{
The gwr.gauss function returns a vector of weights using
the Gaussian scheme:

\deqn{w(g) = e^{{-(d/h)}^2}}

where \eqn{d} are the distances between the observations
and \eqn{h} is the bandwidth.

The default (from release 0.5) gwr.Gauss function returns a vector of 
weights using the Gaussian scheme:

\deqn{w(g) = e^{-(1/2) {{(d/h)}^2}}}


}
\usage{
gwr.gauss(dist2, bandwidth)
gwr.Gauss(dist2, bandwidth)
}
\arguments{
  \item{dist2}{vector of squared distances between observations and fit point}
  \item{bandwidth}{bandwidth}
}
\value{
vector of weights.
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2000,
Quantitative Geography, London: Sage; C. Brunsdon, A.Stewart Fotheringham
and M.E. Charlton, 1996, "Geographically Weighted Regression: A Method for
Exploring Spatial Nonstationarity", Geographical Analysis, 28(4), 281-298;
\url{http://gwr.nuim.ie/}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{gwr.sel}}, \code{\link{gwr}}}

\examples{
plot(seq(-10,10,0.1), gwr.Gauss(seq(-10,10,0.1)^2, 3.5), type="l")
}


\keyword{spatial}

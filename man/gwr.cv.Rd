\name{gwr.sel}
\alias{gwr.sel}
\alias{gwr.cv.f}
\alias{gwr.cv.adapt.f}
\alias{gwr.aic.f}
\alias{gwr.aic.adapt.f}
\title{Crossvalidation of bandwidth for geographically weighted regression}
\description{
The function finds a bandwidth for a given geographically weighted regression
by optimzing a selected function. For cross-validation, this scores the root mean square
prediction error for the geographically weighted regressions, choosing the
bandwidth minimizing this quantity. 
}
\usage{
gwr.sel(formula, data=list(), coords, adapt=FALSE, gweight=gwr.gauss, method = "cv", verbose = TRUE, longlat=FALSE)
}
\arguments{
  \item{formula}{regression model formula as in \code{lm}}
  \item{data}{model data frame as in \code{lm}, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame object as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations}
  \item{adapt}{either TRUE: find the proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours), or FALSE --- find global bandwidth}
  \item{gweight}{geographical weighting function, at present only
  \code{gwr.gauss()} or \code{gwr.bisquare()}}
  \item{method}{default "cv" for drop-1 cross-validation, or "aic" for AIC optimisation (depends on assumptions about AIC degrees of freedom)}
  \item{verbose}{if TRUE (default), reports the progress of search for bandwidth}
  \item{longlat}{if TRUE, use distances on an ellipse with WGS84 parameters}
}
\value{
returns the cross-validation bandwidth.
}
\note{Use of method="aic" results in the creation of an n by n matrix, and should not be chosen when n is large.}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley;
\url{http://www.nuim.ie/ncg/GWR/index.htm}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\seealso{\code{\link{gwr.bisquare}}, \code{\link{gwr.gauss}}}

\examples{
data(columbus)
gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y))
gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y), gweight=gwr.bisquare)
}
\keyword{spatial}

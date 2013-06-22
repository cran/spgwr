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
gwr.sel(formula, data=list(), coords, adapt=FALSE, gweight=gwr.Gauss,
 method = "cv", verbose = TRUE, longlat=NULL, RMSE=FALSE, weights,
 tol=.Machine$double.eps^0.25, show.error.messages = FALSE)
}
\arguments{
  \item{formula}{regression model formula as in \code{lm}}
  \item{data}{model data frame as in \code{lm}, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame object as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations}
  \item{adapt}{either TRUE: find the proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours), or FALSE --- find global bandwidth}
  \item{gweight}{geographical weighting function, at present 
  \code{gwr.Gauss()} default, or \code{gwr.gauss()}, the previous default or \code{gwr.bisquare()}}
  \item{method}{default "cv" for drop-1 cross-validation, or "aic" for AIC optimisation (depends on assumptions about AIC degrees of freedom)}
  \item{verbose}{if TRUE (default), reports the progress of search for bandwidth}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself}
  \item{RMSE}{default FALSE to correspond with CV scores in newer references (sum of squared CV errors), if TRUE the previous behaviour of scoring by LOO CV RMSE}
  \item{weights}{case weights used as in weighted least squares, beware of scaling issues --- only used with the cross-validation method, probably unsafe}
  \item{tol}{the desired accuracy to be passed to \code{optimize}}
  \item{show.error.messages}{default FALSE; may be set to TRUE to see error messages if \code{gwr.sel} returns without a value}
}
\details{If the regression contains little pattern, the bandwidth will converge to the upper bound of the line search, which is the diagonal of the bounding box of the data point coordinates for \dQuote{adapt=FALSE}, and 1 for \dQuote{adapt=TRUE}; see the simulation block in the examples below.}
\value{
returns the cross-validation bandwidth.
}
\note{Use of method="aic" results in the creation of an n by n matrix, and should not be chosen when n is large.}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley; Paez A, Farber S, Wheeler D, 2011, "A simulation-based study of geographically weighted regression as a method for investigating spatially varying relationships", Environment and
Planning A 43(12) 2992-3010;
\url{http://gwr.nuim.ie/}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\seealso{\code{\link{gwr.bisquare}}, \code{\link{gwr.gauss}}}

\examples{
data(columbus)
gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y))
gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y), gweight=gwr.bisquare)
\dontrun{
data(georgia)
set.seed(1)
X0 <- runif(nrow(gSRDF)*3)
X1 <- matrix(sample(X0), ncol=3)
X1 <- prcomp(X1, center=FALSE, scale.=FALSE)$x
gSRDF$X1 <- X1[,1]
gSRDF$X2 <- X1[,2]
gSRDF$X3 <- X1[,3]
yrn <- rnorm(nrow(gSRDF))
gSRDF$yrn <- sample(yrn)
bw <- gwr.sel(yrn ~ X1 + X2 + X3, data=gSRDF, method="cv", adapt=FALSE, verbose=FALSE)
bw
bw <- gwr.sel(yrn ~ X1 + X2 + X3, data=gSRDF, method="aic", adapt=FALSE, verbose=FALSE)
bw
bw <- gwr.sel(yrn ~ X1 + X2 + X3, data=gSRDF, method="cv", adapt=TRUE, verbose=FALSE)
bw
bw <- gwr.sel(yrn ~ X1 + X2 + X3, data=gSRDF, method="aic", adapt=TRUE, verbose=FALSE)
bw
}
}
\keyword{spatial}

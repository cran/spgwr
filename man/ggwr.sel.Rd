\name{ggwr.sel}
\alias{ggwr.sel}
\alias{ggwr.cv.adapt.f}
\alias{ggwr.cv.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Crossvalidation of bandwidth for generalised GWR}
\description{
 The function finds a bandwidth for a given generalised geographically weighted regression
by optimzing a selected function. For cross-validation, this scores the root mean square
prediction error for the generalised geographically weighted regressions, choosing the
bandwidth minimizing this quantity.
}
\usage{
ggwr.sel(formula, data = list(), coords, adapt = FALSE, gweight = gwr.Gauss,
 family = gaussian, verbose = TRUE, longlat = NULL, RMSE=FALSE,
 tol=.Machine$double.eps^0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{regression model formula as in \code{glm}}
  \item{data}{model data frame as in \code{glm}, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame object as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations}
  \item{adapt}{either TRUE: find the proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours), or FALSE --- find global bandwidth}
  \item{gweight}{geographical weighting function, at present 
  \code{gwr.Gauss()} default, or \code{gwr.gauss()}, the previous default or \code{gwr.bisquare()}}
  \item{family}{a description of the error distribution and link function to
          be used in the model, see \code{glm}}
  \item{verbose}{if TRUE (default), reports the progress of search for bandwidth}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself}
  \item{RMSE}{default FALSE to correspond with CV scores in newer references (sum of squared CV errors), if TRUE the previous behaviour of scoring by LOO CV RMSE}
  \item{tol}{the desired accuracy to be passed to \code{optimize}}
}

\value{
  returns the cross-validation bandwidth.
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley;
\url{http://gwr.nuim.ie/}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{
  The use of GWR on GLM is only at the initial proof of concept stage, nothing should be treated as an accepted method at this stage.
}
\seealso{\code{\link{gwr.sel}}, \code{\link{ggwr}}}
\examples{
\donttest{
if (require(rgdal)) {
xx <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
bw <- ggwr.sel(SID74 ~ I(NWBIR74/BIR74) + offset(log(BIR74)), data=xx,
  family=poisson(), longlat=TRUE)
bw
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}


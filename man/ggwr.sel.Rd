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
ggwr.sel(formula, data = list(), coords, adapt = FALSE, gweight = gwr.gauss, family = gaussian, verbose = TRUE, longlat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{regression model formula as in \code{glm}}
  \item{data}{model data frame as in \code{glm}, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame object as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations}
  \item{adapt}{either TRUE: find the proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours), or FALSE --- find global bandwidth}
  \item{gweight}{geographical weighting function, at present only
  \code{gwr.gauss()} or \code{gwr.bisquare()}}
  \item{family}{a description of the error distribution and link function to
          be used in the model, see \code{glm}}
  \item{verbose}{if TRUE (default), reports the progress of search for bandwidth}
  \item{longlat}{if TRUE, use distances on an ellipse with WGS84 parameters}
}

\value{
  returns the cross-validation bandwidth.
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley;
\url{http://www.nuim.ie/ncg/GWR/index.htm}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{
  The use of GWR on GLM is only at the initial proof of concept stage, nothing should be treated as an accepted method at this stage.
}
\seealso{\code{\link{gwr.sel}}, \code{\link{ggwr}}}
\examples{
library(maptools)
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], 
  IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
bw <- ggwr.sel(SID74 ~ I(NWBIR74/BIR74) + offset(log(BIR74)), data=xx,
  family=poisson(), longlat=TRUE)
bw
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}


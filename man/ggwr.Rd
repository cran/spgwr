\name{ggwr}
\alias{ggwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalised geographically weighted regression}
\description{
  The function implements generalised geographically weighted regression approach
to exploring spatial non-stationarity for given global bandwidth and chosen
weighting scheme.
}
\usage{
ggwr(formula, data = list(), coords, bandwidth, gweight = gwr.Gauss,
 adapt = NULL, fit.points, family = gaussian, longlat = NULL, type = 
c("working", "deviance", "pearson", "response"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{regression model formula as in \code{glm}}
  \item{data}{model data frame as in \code{glm}, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame object as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations}
  \item{bandwidth}{bandwidth used in the weighting function, possibly
  calculated by \code{ggwr.sel}}
  \item{gweight}{geographical weighting function, at present 
  \code{gwr.Gauss()} default, or \code{gwr.gauss()}, the previous default or \code{gwr.bisquare()}}
  \item{adapt}{either NULL (default) or a proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours)}
  \item{fit.points}{an object containing the coordinates of fit points; often an object from package \pkg{sp}; if missing, the coordinates given through the data argument object, or the coords argument are used}
  \item{family}{a description of the error distribution and link function to
          be used in the model, see \code{glm}}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself}
  \item{type}{the type of residuals which should be returned. The alternatives are: "working" (default), "pearson", "deviance" and "response"}
}

\value{
A list of class \dQuote{gwr}:
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package "sp") with fit.points, weights, GWR coefficient estimates, dispersion if a "quasi"-family is used, and the residuals of type "type" in its "data" slot.}
  \item{lhat}{Leung et al. L matrix, here set to NA}
  \item{lm}{GLM global regression on the same model formula.}  \item{bandwidth}{the bandwidth used.}
  \item{this.call}{the function call used.}

}

\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley;
\url{http://www.nuim.ie/ncg/GWR/index.htm}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\note{
  The use of GWR on GLM is only at the initial proof of concept stage, nothing should be treated as an accepted method at this stage.
}
\seealso{\code{\link{ggwr.sel}}, \code{\link{gwr}}}
\examples{
library(maptools)
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], 
  IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
bw <- 144.4813
\dontrun{
bw <- ggwr.sel(SID74 ~ I(NWBIR74/BIR74) + offset(log(BIR74)), data=xx,
  family=poisson(), longlat=TRUE)
}
nc <- ggwr(SID74 ~ I(NWBIR74/BIR74) + offset(log(BIR74)), data=xx,
  family=poisson(), longlat=TRUE, bandwidth=bw)
nc
\dontrun{
nc <- ggwr(SID74 ~ I(NWBIR74/10000) + offset(log(BIR74)), data=xx,
  family=poisson(), longlat=TRUE, bandwidth=bw)
nc
nc <- ggwr(SID74 ~ I(NWBIR74/10000) + offset(log(BIR74)), data=xx,
  family=quasipoisson(), longlat=TRUE, bandwidth=bw)
nc
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}


\name{gwr.morantest}
\alias{gwr.morantest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moran's I for gwr objects}
\description{
  The function returns Leung et al. (2000) three moment approximation for Moran's I, for a gwr object calculated with argument hatmatrix=TRUE. This implementation should not be regarded as authoritative, as it involves assumptions about implied methods and about estimated degrees of freedom.
}
\usage{
gwr.morantest(x, lw, zero.policy = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{gwr} object returned by \code{gwr()} with argument hatmatrix=TRUE}
  \item{lw}{a \code{listw} object created for example by \code{nb2listw} in the spdep package}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA}
}

\value{
a \dQuote{htest} object with the results of testing the GWR residuals
}
\references{Leung Y, Mei C-L, Zhang W-X 2000 Testing for spatial autocorrelation among the residuals of the geographically weighted regression, Environment and Planning A, 32, 871-890.}
\author{Roger Bivand}

\examples{
if (suppressWarnings(require(spData)) && suppressWarnings(require(spdep))) {
  data(columbus, package="spData")
  bw <- gwr.sel(CRIME ~ INC + HOVAL, data=columbus, coords=coords)
  col0 <- gwr(CRIME ~ INC + HOVAL, data=columbus, coords=coords,
    bandwidth=bw, hatmatrix=TRUE)
  gwr.morantest(col0, nb2listw(col.gal.nb))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}


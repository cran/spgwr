\name{LMZ.F3GWR.test}
\alias{LMZ.F3GWR.test}
\alias{LMZ.F2GWR.test}
\alias{LMZ.F1GWR.test}
\alias{BFC99.gwr.test}
\alias{BFC02.gwr.test}
\alias{anova.gwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Global tests of geographical weighted regressions}
\description{
 Four related test statistics for comparing OLS and GWR models based on bapers by Brunsdon, Fotheringham and Charlton (1999) and Leung et al (2000), and a development from the GWR book (2002).
}
\usage{
LMZ.F3GWR.test(go)
LMZ.F2GWR.test(x)
LMZ.F1GWR.test(x)
BFC99.gwr.test(x)
BFC02.gwr.test(x, approx=FALSE)
\method{anova}{gwr}(object, \dots, approx=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go, x, object}{a \code{gwr} object returned by \code{gwr()}}
  \item{\dots}{arguments passed through (unused)}
  \item{approx}{default FALSE, if TRUE, use only (n - tr(S)) instead of (n - 2*tr(S) - tr(S'S)) as the GWR degrees of freedom}
}
\details{
  The papers in the references give the background for the analyses of variance presented.
}
\value{
BFC99.GWR.test, BFC02.gwr.test, LMZ.F1GWR.test and LMZ.F2GWR.test return "htest" objects, 
LMZ.F3GWR.test a matrix of test results.
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley;
\url{http://gwr.nuim.ie/}}

\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Danlin Yu}

\seealso{\code{\link{gwr}}}

\examples{
data(columbus)
col.bw <- gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y))
col.gauss <- gwr(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y), bandwidth=col.bw, hatmatrix=TRUE)
BFC99.gwr.test(col.gauss)
BFC02.gwr.test(col.gauss)
BFC02.gwr.test(col.gauss, approx=TRUE)
anova(col.gauss)
anova(col.gauss, approx=TRUE)
\dontrun{
BFC99.gwr.test(col.bisq)
}
}

\keyword{spatial}


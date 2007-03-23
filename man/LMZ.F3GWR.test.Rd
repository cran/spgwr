\name{LMZ.F3GWR.test}
\alias{LMZ.F3GWR.test}
\alias{LMZ.F2GWR.test}
\alias{LMZ.F1GWR.test}
\alias{BFC99.gwr.test}
\alias{BFC02.gwr.test}
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
BFC02.gwr.test(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go, x}{a \code{gwr} object returned by \code{gwr()}}
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
\url{http://www.nuim.ie/ncg/GWR/index.htm}}

\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Danlin Yu}

\seealso{\code{\link{gwr}}}

\examples{
example(gwr)
BFC99.gwr.test(col.gauss)
BFC02.gwr.test(col.gauss)
BFC99.gwr.test(col.bisq)
}

\keyword{spatial}


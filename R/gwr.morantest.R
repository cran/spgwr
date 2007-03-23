gwr.morantest <- function(x, listw, zero.policy = FALSE) {
	if(class(x) != "gwr") stop(paste(deparse(substitute(x)),
		"not a gwr object"))
	if (is.null(x$lhat)) stop("hatmatrix=TRUE needed in gwr fit")
	if (!inherits(listw, "listw")) 
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
	n <- ncol(x$lhat)
	if (n != length(listw$neighbours)) stop("objects of different length")
	if (listw$style != "W") warning(deparse(substitute(listw)),
		"not row standardised")
	W <- listw2mat(listw2U(listw))
	N <- diag(n) - x$lhat
	e.gwr <- N %*% x$lm$y
	I0.gwr <- c((t(e.gwr) %*% W %*% e.gwr) / (t(e.gwr) %*% e.gwr))
	A <- t(N) %*% (W - I0.gwr*diag(n)) %*% N
	EQ.gwr <- sum(diag(A))
	tr2 <- sum(diag(A %*% A))
	tr3 <- sum(diag(A %*% A %*% A))
	varQ.gwr <- 2*tr2
	EQ.EQ3.gwr <- 8*tr3
	h.gwr <- tr2^3/tr3^2
	chi2.gwr <- 0
	p.gwr <- 0
	if (EQ.EQ3.gwr > 0) {
		chi2.gwr <- h.gwr - ((sqrt(2*h.gwr)*EQ.gwr)/(sqrt(varQ.gwr)))
		p.gwr <- 1 - pchisq(chi2.gwr, h.gwr)
		}
	else if (EQ.EQ3.gwr < 0) {
		chi2.gwr <- (h.gwr +
			((sqrt(2*h.gwr)*EQ.gwr)/(sqrt(varQ.gwr))))
		p.gwr <- pchisq(chi2.gwr, h.gwr)
		}
	GWRtest <- list(estimate=c(I = I0.gwr),
		statistic=c(statistic = chi2.gwr), parameter=c(df = h.gwr),
		p.value=p.gwr, data.name="GWR residuals",
		method="Leung et al. 2000 three moment approximation for Moran's I")
	class(GWRtest) <- "htest"
	GWRtest
}

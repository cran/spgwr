gwr.gauss <- function(dist2, bandwidth) {
	w <- exp((-dist2)/(bandwidth^2))
	w
}


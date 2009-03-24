# Copyright 2001-2008 Roger Bivand and Danlin Yu
# 

gwr <- function(formula, data = list(), coords, bandwidth, 
	gweight=gwr.Gauss, adapt=NULL, hatmatrix=FALSE, fit.points, 
	longlat=FALSE, se.fit=FALSE, weights, cl=NULL) {
	this.call <- match.call()
	p4s <- as.character(NA)
	Polys <- NULL
	if (is(data, "SpatialPolygonsDataFrame")) 
		Polys <- as(data, "SpatialPolygons")
	if (is(data, "Spatial")) {
		if (!missing(coords))
		    warning("data is Spatial* object, ignoring coords argument")
		coords <- coordinates(data)
		p4s <- proj4string(data)
		data <- as(data, "data.frame")
	}

	if (missing(coords))
		stop("Observation coordinates have to be given")
	if (is.null(colnames(coords))) 
		colnames(coords) <- c("coord.x", "coord.y")

    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data", "weights"), names(mf), 0)
    	mf <- mf[c(1, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())
    	mt <- attr(mf, "terms")
	dp.n <- length(model.extract(mf, "response"))

#	mt <- terms(formula, data = data)
#	mf <- lm(formula, data, method="model.frame", na.action=na.fail)

    	weights <- as.vector(model.extract(mf, "weights"))
# set up default weights
    	if (!is.null(weights) && !is.numeric(weights)) 
        	stop("'weights' must be a numeric vector")
    	if (is.null(weights)) weights <- rep(as.numeric(1), dp.n)
    	if (any(is.na(weights))) stop("NAs in weights")
    	if (any(weights < 0)) stop("negative weights")
	y <- model.extract(mf, "response")
	x <- model.matrix(mt, mf)

	lm <- lm.wfit(x, y, w=weights)
	lm$x <- x
	lm$y <- y
	if (missing(fit.points)) {
		fp.given <- FALSE
		fit.points <- coords
		colnames(fit.points) <- colnames(coords)
	} else fp.given <- TRUE
	griddedObj <- FALSE
	if (is(fit.points, "Spatial")) {
		Polys <- NULL
		if (is(fit.points, "SpatialPolygonsDataFrame")) {
			Polys <- Polygons(fit.points)
			fit.points <- coordinates(fit.points)
		} else {
			griddedObj <- gridded(fit.points)
			fit.points <- coordinates(fit.points)
		}
	}

	n <- NROW(fit.points)
	rownames(fit.points) <- NULL
	if (is.null(colnames(fit.points))) colnames(fit.points) <- c("x", "y")
#	if (is.null(fit.points)) fit.points <- coords
	yiybar <- (y - mean(y))
	m <- NCOL(x)
	if (NROW(x) != NROW(coords))
		stop("Input data and coordinates have different dimensions")
	if (missing(bandwidth) && is.null(adapt))
	    stop("Bandwidth must be given for non-adaptive weights")
	if (missing(bandwidth)) bandwidth <- NULL
	if (hatmatrix) se.fit <- TRUE
	lhat <- NA


	GWR_args <- list(fp.given=fp.given, hatmatrix=hatmatrix, 
	    longlat=longlat, bandwidth=bandwidth, adapt=adapt, se.fit=se.fit)

	if (!is.null(cl) && length(cl) > 1 && fp.given && !hatmatrix) {
	    if (length(grep("cluster", class(cl))) > 0 && 
		.Platform$OS.type == "unix" && exists("splitIndices")) {
#		&& require(snow)) 
		l_fp <- lapply(splitIndices(nrow(fit.points), length(cl)), 
		    function(i) fit.points[i,])
		clusterEvalQ(cl, library(spgwr))

		clusterExport_l <- function(cl, list) {
                    gets <- function(n, v) {
                        assign(n, v, env = .GlobalEnv)
                        NULL
                    }
                    for (name in list) {
                        clusterCall(cl, gets, name, get(name))
                    }
		}

		clusterExport_l(cl, list("GWR_args", "coords", "gweight", "y",
		    "x", "yiybar", "weights"))

		res <- parLapply(cl, l_fp, function(fp) .GWR_int(fit.points=fp,
		    coords=coords, gweight=gweight, y=y, x=x, yiybar=yiybar,
		    weights=weights, GWR_args=GWR_args))
		clusterEvalQ(cl, rm(list=c("GWR_args", "coords", "gweight", "y",
		    "x", "yiybar", "weights")))
		df <- as.data.frame(do.call("rbind", 
		    lapply(res, function(x) x$df)))
		bw <- do.call("c", lapply(res, function(x) x$bw))
	        results <- NULL
	    } else stop("snow package not loaded")
	} else { # cl

	    df <- .GWR_int(fit.points=fit.points, coords=coords, 
		gweight=gweight, y=y, x=x, yiybar=yiybar, weights=weights,
		GWR_args=GWR_args)
	    if (!fp.given && hatmatrix) lhat <- df$lhat
	    bw <- df$bw
	    df <- as.data.frame(df$df)
	    results <- NULL

	} # cl
	if (!fp.given && hatmatrix) {

	#This section calculates the effective degree of freedoms, edf;
	#the normalized residual sum of square, sigma2; the model
	#residual of squares, rss; and various version of AICs
	
	#As long as the hat matrix is obtained, many statistics can be 
	#calculated, such as the effective degree of freedom, AIC, etc.
	#Now calculate the effective degree of freedom of the residual
	#Reference: GWR book, page 55
	#obtain v1 and v2:
	
		v1 <- sum(diag(lhat))
		B2 <- t(lhat)%*%lhat
		v2 <- sum(diag(B2))
	
	#effective d.f. is n - 2*v1 + v2
	
		edf <- n - 2*v1 + v2
	
	#Follow Leung et al. EPA 2000 page 15, the estimate of sigma square
	#can be obtained through rss and delta1 (which is actually edf)
	#Calculate rss:
	
		B1 <- t(diag(n)-lhat)%*%(diag(n)-lhat)
		rss <- c(t(y)%*%B1%*%y)
		delta1 <- sum (diag (B1))
		sigma2 <- rss/delta1 #line 77
	
	#Now the problem is, there are several version of AIC's calculation
	#formula, the GWR book's (page 61,96), Brunsdon's handouts, 
	#and the one from Hurvich, Simonoff and Tsai (1998, page 276)
	#I will implement all of them, and called them
	#AICb, from the book, AICh, from Brunsdon, and AICc, from Hurvich et al.
	
	#AICb <- n*log(sigma2) + n*log(2*3.14) + (n * (n + v1) / (n - 2 - v1))
	#AICh <- n*log(sigma2) + ((n + v1) / (n + 2 - v1))
	
	#To calculate AICc, there are several interal parameters
	#delta1 has already been calculated, detailed formula see
	#Hurvich et al. 1998, page 275, 276
	#B1 (above) and B2, delta2, nu1, nu2:
	
		delta2 <- sum(diag(B1)^2)
		nu1 <- sum(diag(B2))
	#nu2 <- sum(diag(B2^)2)
	
	#AICc is from the formula in Hurvich et al. 1998 page 276
	#AICc1:
	
	#AICc <- n*log(sigma2) + n * ((delta1/delta2)*(n + nu1))/((delta1^2/delta2)-2)
	
	#One thing that I did not notice is that the sigma2 here I used is not
	#the same sigma2 used in the GWR book (detailed reference in page 96).
	#The sigma2 I used is from Leung et al (2000, p 15), calculated in line 77
	#The sigma2 in the GWR book is a maximum likelihood estimate
	#It should be: sigma2 <- rss/n instead of sigma2 <- rss/delta1
	#For this reason, a corrected AICb.b, AICh.b, AICc.b are therefore provided
	#followed by creating the sigma sqare used in the book, termed here sigma2.b
	#All the above unncessary calculation is then commented out.
	
		sigma2.b <- rss / n
		AICb.b <- 2*n*log(sqrt(sigma2.b)) + n*log(2*pi) + 
			(n * ((n + v1) / (n - 2 - v1)))
# NOTE 2* and sqrt() inserted for legibility
		AICh.b <- 2*n*log(sqrt(sigma2.b)) + n*log(2*pi) + n + v1
		AICc.b <- 2*n*log(sqrt(sigma2.b)) + n * 
			((delta1/delta2)*(n + nu1))/((delta1^2/delta2)-2)
		results <- list(v1=v1, v2=v2, delta1=delta1, delta2=delta2, 
			sigma2=sigma2, sigma2.b=sigma2.b, AICb=AICb.b, 
			AICh=AICh.b, AICc=AICc.b, edf=edf, rss=rss, nu1=nu1)
	}
#	df <- data.frame(sum.w=sum.w, gwr.b, gwr.R2, gwr.se, gwr.e)
		
	SDF <- SpatialPointsDataFrame(coords=fit.points, 
		data=df, proj4string=CRS(p4s))

	if (griddedObj) {
		gridded(SDF) <- TRUE
	} else {
		if (!is.null(Polys)) {
			df <- data.frame(SDF@data)
			rownames(df) <- sapply(slot(Polys, "polygons"),
                            function(i) slot(i, "ID"))
			SDF <- SpatialPolygonsDataFrame(Sr=Polys, data=df)
		}
	}
	z <- list(SDF=SDF, lhat=lhat, lm=lm, results=results, 
		bandwidth=bw, adapt=adapt, hatmatrix=hatmatrix, 
		gweight=deparse(substitute(gweight)), this.call=this.call)
	class(z) <- "gwr"
	invisible(z)
}


print.gwr <- function(x, ...) {
	if(class(x) != "gwr") stop("not a gwr object")
	cat("Call:\n")
	print(x$this.call)
	cat("Kernel function:", x$gweight, "\n")
	n <- nrow(coordinates(x$SDF))
	if (is.null(x$adapt)) cat("Fixed bandwidth:", x$bandwidth, "\n")
	else cat("Adaptive quantile: ", x$adapt, " (about ", 
		floor(x$adapt*n), " of ", n, ")\n", sep="")
	m <- length(x$lm$coefficients)
	cat("Summary of GWR coefficient estimates:\n")
	CM <- t(apply(as(x$SDF, "data.frame")[,(1+(1:m)), drop=FALSE], 2, summary))[,c(1:3,5,6)]
	if (is.null(dim(CM))) CM <- t(as.matrix(CM))
	CM <- cbind(CM, coefficients(x$lm))
	colnames(CM) <- c(colnames(CM)[1:5], "Global")
	printCoefmat(CM)
	if (x$hatmatrix) {
		cat("Number of data points:", n, "\n")
		cat("Effective number of parameters:", 2*x$results$v1 -
			x$results$v2, "\n")
		cat("Effective degrees of freedom:", x$results$edf, "\n")
		cat("Sigma squared (ML):", x$results$sigma2.b, "\n")
		cat("AICc (GWR p. 61, eq 2.33; p. 96, eq. 4.21):", x$results$AICb, "\n")
		cat("AIC (GWR p. 96, eq. 4.22):", x$results$AICh, "\n")
		cat("Residual sum of squares:", x$results$rss, "\n")
	}
	invisible(x)
}

.GWR_int <- function(fit.points, coords, gweight, y, x, yiybar, weights, 
	GWR_args) {
	    n <- nrow(fit.points)
 	    m <- NCOL(x)
	    if (GWR_args$se.fit) {
		df <- matrix(nrow=n, ncol=(2*m + 3))
	        colnames(df) <- c("sum.w", colnames(x), "R2", "gwr.e", 
		    paste(colnames(x), "se", sep="_"))
	    } else {
		df <- matrix(nrow=n, ncol=(m + 3))
	    	colnames(df) <- c("sum.w", colnames(x), "R2", "gwr.e")
	    }
	    if (!GWR_args$fp.given && GWR_args$hatmatrix) 
	        lhat <- matrix(nrow=n, ncol=n)
	    if (is.null(GWR_args$adapt)) {
		    bw <- GWR_args$bandwidth
		    bandwidth <- rep(GWR_args$bandwidth, n)
	    } else {
		bandwidth <- gw.adapt(dp=coords, fp=fit.points, 
		    quant=GWR_args$adapt, longlat=GWR_args$longlat)
		bw <- bandwidth
	    }
	    if (any(bandwidth < 0)) stop("Invalid bandwidth")
	    for (i in 1:n) {
		dxs <- spDistsN1(coords, fit.points[i,], 
		    longlat=GWR_args$longlat)
		if (any(!is.finite(dxs)))
			dxs[which(!is.finite(dxs))] <- 0
#		if (!is.finite(dxs[i])) dxs[i] <- 0
		w.i <- gweight(dxs^2, bandwidth[i])
		w.i <- w.i * weights
		if (any(w.i < 0 | is.na(w.i)))
        		stop(paste("Invalid weights for i:", i))
		lm.i <- lm.wfit(x, y, w.i)
		df[i, 1] <- sum(w.i)
		df[i, 2:(m+1)] <- coefficients(lm.i)
		ei <- residuals(lm.i)
# use of diag(w.i) dropped to avoid forming n by n matrix
# bug report: Ilka Afonso Reis, July 2005
		rss <- sum(ei * w.i * ei)
#		if (!GWR_args$fp.given && GWR_args$hatmatrix) {
		df[i, (m+3)] <- ei[i]
		df[i, (m+2)] <- 1 - (rss / sum(yiybar * w.i * yiybar))
#                } else is.na(df[i, (m+(2:3))]) <- TRUE
	        if (GWR_args$se.fit) {
		    p <- lm.i$rank
		    p1 <- 1:p
		    inv.Z <- chol2inv(lm.i$qr$qr[p1, p1, drop=FALSE])
		    df[i,(m+4):(2*m + 3)] <- sqrt(diag(inv.Z) * (rss/(n-p)))
		}
		if (!GWR_args$fp.given && GWR_args$hatmatrix) 
			lhat[i,] <- t(x[i,]) %*% inv.Z %*% t(x) %*% diag(w.i)
	    }
	    if (!GWR_args$fp.given && GWR_args$hatmatrix) 
		return(list(df=df, lhat=lhat, bw=bw))
	    else return(list(df=df, bw=bw))
} # GWR_int



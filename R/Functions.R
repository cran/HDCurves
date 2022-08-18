# Functions needed to fit Derivative Curve Clustering Methodology

tpower <- function(x, t, p) {
	# Function for truncated p-th power function
	return((x - t) ^ p * (x > t))
}
# Function for constructing a B-spline basis
# x: covariate
# ndx: number of segments (related to the number of knots)
# bdeg: degree of the polynomial (e.g. 3: Cubic B-splines)
bbase <- function(x, ndx, bdeg = 3, eps = 1e-5) {
	xl = min(x)
	xr = max(x)
	dx <- (xr - xl)/ndx
	knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by=dx)
	P <- outer(x, knots, tpower, bdeg)
	n <- dim(P)[2]
	D <- diff(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
	B <- (-1) ^ (bdeg + 1) * P %*% t(D)
	B[B < eps] <- 0
	attr(B,"knots") <- knots
	attr(B,"bdeg") <- bdeg
	attr(B,"eps") <- eps
	class(B) <- c("bbase")
	B
}
# Prediction function
predict.bbase <- function(object, newx) {
	knots <- attr(object,"knots")
	bdeg <- attr(object,"bdeg")
	eps <- attr(object,"eps")

	dx <- diff(knots)[1]
	P <- outer(newx, knots, tpower, bdeg)
	n <- dim(P)[2]
	D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1)*dx^bdeg)
	B <- (-1) ^ (bdeg + 1) * P %*% t(D)
	B[B < eps] = 0
	B
}

emp_deriv <- function(y, t, D, u, ztype=0){
  #	ztype - 0 use transformation using D and u,
  #		    - 1 use transformation using u only
  #		    - 2 transformation from de Brabantar

  nobs <- length(y)
  z <- rep(0, nobs)

  C.out <- .C("empderiv",
              as.integer(nobs),
              as.double(y), as.double(t),
              as.integer(ztype),
              as.integer(D),as.double(u),
              fprime = as.double(z))

  C.out$fprime
}

HDCurves <- function(y, t, trt, ids, balanced = TRUE, tpred = NULL,
					    ztype = 0, p = 0.5, au = 1, bu = 1, A = 2, as = 1, bs = 1,
					    atau = 1, btau = 1/0.005, ndx = 10, q = 3, verbose=FALSE,
					    draws = 1100, burn = 100, thin = 1) {

  # ztype: 0 - mine (employing both D and u)
	#		     1 - cote (employing only u)

	df <- data.frame(t = t, y = y, trt = trt, ids = ids)


	usubject <- unique(df$ids)		# Vector with ids
	nsubject <- length(usubject)	# Number of individuals

	g <- sapply(usubject, function(x, trt, subject) trt[subject == x][1], trt = df$trt, subject = df$ids)
	nobs <- sapply(usubject, function(x, subject) sum(subject == x), subject = df$ids) # Observations per individual

	ngroups <- length(unique(df$trt))

  tm <- df$t
	p_lag <- p

	# Best to feed in each subjects basis created from Cote's code rather than
	# creating code in C. Not the Hmat matrix depends on the balanced argument
	#
	# balanced 1 - All subjects have same number of measurements and time at which
	#              they were measure and so have same design matrix
	#		   0 - Subjects do not have same number of measurements and they are
	#		       measured at different time points so each subject needs their
	#              own design matrix

	# Unique time points
	tt <- sort(unique(df$t))

	# Create the basis
	G <- bbase(tt, ndx = ndx, bdeg = q)

	# For prediction
	if(is.null(tpred)) {
		tpred <- seq(min(df$t), max(df$t), length = 50)
	}
	Gp <- predict.bbase(G, tpred)

	D <- diff(diag(ncol(G)), differences = 2)
	K <- crossprod(D)

  	nb <- ncol(K)
	Hmat <- G

	if(!balanced) {
		Hmat <- NULL
		cnobs <- 0
		for(j in 1:nsubject){
			ttmp <- df$t[(cnobs+1):(cnobs + nobs[j])]
			bsb <- predict.bbase(G, ttmp)
			Hmat <- c(Hmat, c(t(bsb)))
			cnobs <- cnobs + nobs[j]

		}
	}

	Htheta <- G
	nrHt <- nrow(Htheta)

	nout <- (draws - burn)/thin
	modelPriors <- c(as, bs, atau, btau)

	beta <- matrix(0, nrow = nout, ncol = nb*nsubject)
	zout <- fprime <- matrix(0, nrow = nout, ncol = nsubject*max(nobs))
	u <- sig2 <- D <- matrix(1, nrow = nout, ncol = nsubject)

	theta <- matrix(0, nrow = nout, ncol = nb*ngroups)
	fgprime <- matrix(0, nrow=nout, ncol=ngroups*nrHt)
	lam <- tau2 <- matrix(1, nrow = nout, ncol = ngroups)

	out <- NULL

	C.out <- .C("mcmcloop",
              	as.integer(draws), as.integer(burn), as.integer(thin),
              	as.integer(nsubject), as.integer(nobs),
              	as.double(df$y), as.double(t(df$t)),
              	as.integer(nb),as.double(t(K)),
              	as.integer(g), as.integer(ngroups),
              	as.double(A), as.double(p_lag), as.double(au), as.double(bu),
              	as.double(modelPriors),as.double(Hmat), as.double(t(Htheta)),
              	as.integer(nrHt), as.integer(balanced),as.integer(ztype),
	              as.integer(verbose),
              	beta.draws = as.double(beta), theta.draws = as.double(theta),
              	lam.draws = as.double(lam), tau2.draws = as.double(tau2),
              	u.draws = as.double(u),
              	D.draws = as.integer(D), z.draws = as.double(zout),
              	sig2.draws = as.double(sig2), fprime.draws = as.double(fprime),
              	fgprime.draws = as.double(fgprime))


	out$z <- array(C.out$z.draws, c(max(nobs), nsubject, nout))
	out$beta <- array(C.out$beta.draws, c(nsubject, nb, nout))

	out$tpred <- tpred

	# Obtain derivative curves for the desired times
	if(!is.null(tpred)){
		out$fprime <- apply(out$beta, c(1,3), function(x, G) G%*%x, G = Gp)
	} else {
		out$fprime <- array(C.out$fprime.draws, c(max(nobs),nsubject,nout))
	}

	out$u <- matrix(C.out$u.draws, nrow = nout, byrow = TRUE)
	out$d <- matrix(C.out$D.draws, nrow = nout, byrow = TRUE)
	out$sig2 <- matrix(C.out$sig2.draws, nrow = nout, byrow = TRUE)

  	out$theta <- array(C.out$theta.draws, c(ngroups, nb, nout))

  	# Obtain derivative curves for the desired times
  	out$fgprime <- apply(out$theta, c(1,3), function(x, G) G%*%x, G = Gp)

	out$lam <- matrix(C.out$lam.draws, nrow = nout, byrow = TRUE)
	out$tau2 <- matrix(C.out$tau2.draws, nrow = nout, byrow = TRUE)

 	out$HmatBySubject <- Hmat
	out$Htheta <- G
	out$tpred <- tpred
  out

}




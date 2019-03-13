ERA <- function(y=y, X=X, nvar=nvar, dist=1, const = 0, lambda = 0, it = 0, ceps = 0.0001){
	
	if(dist==1){ Yfam = "gaussian"
	} else if (dist==2) { Yfam = "binomial" #(link = "logit")
	} else if (dist==3) { Yfam = "poisson" #(link = "log")
	} else if (dist==4) { Yfam = "Gamma" #(link = "inverse")
	} else if (dist==5) { Yfam = "inverse.gaussian" #(link = "1/mu^2")
	}
	
	data = cbind(y,X)
	ndset = length(nvar)
	ncase = nrow(y)
	ny = ncol(y)
	sum_nvar = sum(nvar)
	# ------------------
	# data normalization
	# ------------------
	dataN <- scale(data)/sqrt(ncase-1)
	y <- dataN[,1]
	X <- dataN[,-1]
	
	# ------------------
	# initialization
	# ------------------
	W0 <- matrix(0, nrow=sum_nvar, ncol=ndset)

	kk = 0
	for (j in 1:ndset){
		k = kk + 1
		kk = kk + nvar[,j]
		if ( nvar[,j]==1 ) {
			W0[k:kk,j] = 1
		} else {
			W0[k:kk,j] = 99
		}
	}
	windex <- which(W0 == 99)
	num_windex = length(windex)
	W = W0
	W[windex] <- runif(num_windex)
	F = X%*%W
	if( const==1 ) {
		demF <- chol(t(F)%*%F)
		F <- t(solve(t(demF),t(F)))
	} else {
		for (j in 1:ndset) { F[,j] <- F[,j] / norm(F[,j,drop=FALSE], type = "2") }
	}
	A <- solve(t(F)%*%F,t(F)%*%y)
	lp <- F%*%A
	adj.DV <- linkfun(X, lp, y, family = Yfam)
	mu <- adj.DV$mu #getmu(lp,dist)
	v <- adj.DV$mvar #getvmu(mu,dist)
	V <- diag(v)
	z <- adj.DV$z # lp + (y-mu)/v
	
	# ------------------
	# IRLS Algorithm
	# ------------------
	vecW <- W[windex]
	est_new <- c(vecW,A)
	est_old <- rep(0, length(est_new))

	while ( sum(abs(est_new - est_old)) > ceps ) {
		it = it + 1
		if ( it > 1000 ) {
			print("Message: not converged within 1000 iterations")
			break
		}
		est_old <- est_new
		
		# Step1: Update weights
		M0 <- kronecker(t(A),X)
		M <- M0[,windex]
		temp <- t(M)%*%V%*%M + lambda*diag(num_windex)
		vecW <- solve(temp, t(M)%*%V%*%z)
		W[windex] <- vecW
		F = X%*%W
		if( const==1 ) {
			demF <- chol(t(F)%*%F)
			F <- t(solve(t(demF),t(F)))
		} else {
			for (j in 1:ndset) { F[,j] <- F[,j] / norm(F[,j,drop=FALSE], type = "2") }
		}
		
		# Step2: Update A
		A <- solve(t(F)%*%V%*%F, t(F)%*%V%*%z)
		
		# Step3: Update V and z
		lp <- F%*%A
		adj.DV <- linkfun(X, lp, y, family = Yfam)
		mu <- adj.DV$mu #getmu(lp,dist)
		v <- adj.DV$mvar #getvmu(mu,dist)
		V <- diag(v)
		z <- adj.DV$z # lp + (y-mu)/v
		est_new <- c(vecW,A)
	}
	vecW <- vecW
	A <- A

	if (dist == 1) {
		sigma2 <- sum((y-mu)^2)/(ncase-ndset) + 1E-5
	} else if (dist == 4) {
		sigma2 <- sum(((y - mu)/mu)^2)/(ncase-ndset) + 1E-5
	} else {
		sigma2 <- 1
	}
	
	# ------------------
	# Calculate the standard errors
	# ------------------
	G <- diag(num_windex)
	SS <- M
	H11 <- -t(F)%*%V%*%F
	H22 <- -t(SS)%*%V%*%SS - lambda*G
	DBDA <- matrix(0,nrow=ndset,ncol=num_windex)

	require(MASS)
	for (j in 1:ndset) {
		e <- matrix(0, nrow=ndset, ncol=1)
		e[j] <- 1
		KM <- kronecker(t(e),X)
		KMM <- KM[,windex]
		dBdA <- ginv(t(SS)%*%V%*%SS + lambda*G) %*% (t(KMM)%*%V%*%z - t(KMM)%*%V%*%SS%*%(ginv(t(SS)%*%V%*%SS + lambda*G)%*%(t(SS)%*%V%*%z)))
		DBDA[j,] <- -t(dBdA)
	}
	H12 <- DBDA%*%H22
	H21 <- t(H12)
	INF <- -rbind(cbind(H11,H12),cbind(H21,H22))
	cov_est <- solve(INF)
	se_est <- sqrt(sigma2)*sqrt(diag(cov_est))
	stdA <- se_est[1:ndset]
	stdW <- se_est[-(1:ndset)]
	DS <- getdeviance(y,mu,dist)
	LL = -DS/2
	DF = num_windex + ndset
	AIC = -2*LL + 2*DF;
	BIC = -2*LL + log(ncase)*DF
	CAIC = -2*LL + log(ncase+1)*DF

	# ------------------
	# Calculate the "estfun" for MOB
	# ------------------
	pm <- X%*%W%*%A
	dev <- y-pm
	scoreX <- X # estfun() for W
	for (i in 1:ncol(scoreX)) {	scoreX[,i] <- dev*scoreX[,i] }
	scoreF <- F # estfun() for A
	for (i in 1:ncol(scoreF)) { scoreF[,i] <- dev*scoreF[,i] }
	
	output.ERA <- list(it=it, lambda=lambda, sigma2=sigma2, X=X, scoreX=scoreX, scoreF=scoreF,
					lp=lp, mu=mu, V=V, adj.DV=z, W=W, F=F, vecW=vecW, A=A,
					INF=INF, cov_est=cov_est, stdA=stdA, stdW=stdW,
					DS=DS, LL=LL, DF=DF, AIC=AIC, BIC=BIC, CAIC=CAIC)
	res <- output.ERA
	res$it
	res$lambda
	res$sigma2
	res$X
	res$scoreX
	res$scoreF
	res$lp
	res$mu
	res$V
	res$adj.DV
	res$W
	res$F
	res$vecW
	res$A
	res$INF
	res$cov_est
	res$stdA
	res$stdW
	res$DS
	res$LL
	res$DF
	res$AIC
	res$BIC
	res$CAIC
	
	class(res) <- c("ERA", class(res))
	invisible(res)
}

getdeviance <- function(y, mu, datatype) {
	if (datatype == 1) {
		DS = sum((y - mu)^2)

	} else if (datatype == 2) {
	
	}
	
	return(DS)
}


linkfun <- function(x, lp, y, weights = NULL, offset = NULL, family = family,
					mustart = NULL, etastart = NULL, start = NULL)
{
	# Ref: stats::glm() in R
	# In R, family objects: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/glm.R
	# IWSL with adj DVs: https://www.statistics.ma.tum.de/fileadmin/w00bdb/www/czado/lec2.pdf
	# From linkfun(): adj.DV$w (weights in IWLS) is all 1 when gaussian
	
	x <- as.matrix(x)
	nvars <- ncol(x)
	nobs <- nrow(x)
	if (is.null(weights)) {	weights <- rep.int(1, nobs) }
    if (is.null(offset)) { offset <- rep.int(0, nobs) }
	
	## Define family:	
	if(is.character(family)) {
		family <- get(family, mode = "function", envir = parent.frame())
	}
	
	if(is.function(family)) family <- family()
	
	if(is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}
	
	## Get family functions:
	variance <- family$variance
	linkinv  <- family$linkinv
		 if (!is.function(variance) || !is.function(linkinv) ) {
			stop("'family' argument seems not to be a valid family object", call. = FALSE)
		}
    mu.eta <- family$mu.eta
	
	unless.null <- function(x, if.null) if(is.null(x)) if.null else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu  <- unless.null(family$validmu,  function(mu) TRUE)
	
	if(is.null(mustart)) {
        ## calculates mustart and may change y and weights and set n (!)
        eval(family$initialize)
    } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
	
	## Get (initial) eta
	coefold <- NULL
	eta <- lp
	mu <- linkinv(eta)
	mvar <- variance(mu) # the variance as a function of the mean
		
	if (!(validmu(mu) && valideta(eta))) {
		stop("cannot find valid starting values: please specify some", call. = FALSE)
	}

	## Calculate z and w
	good <- weights > 0
	mu.eta.val <- mu.eta(eta)
		if (any(is.na(mu.eta.val[good]))) {stop("NAs in d(mu)/d(eta)")}
	
	# drop observations for which w will be zero
	good <- (weights > 0) & (mu.eta.val != 0)
	z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
		# z: adjusted dependent variable
	w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
		# w: weights

	# BETA in GLM: Fisher Scoring
	# beta <- solve(t(X) %*% diag(w) %*% X) %*% (t(X) %*% diag(w) %*% z)
	
	# Outs
	return( list(mu=mu, mvar=mvar, z=z, w=w) )
}


coef.ERA <- function(object){
	structure(t(rbind(object$vecW,object$A)), class="coef")
}

logLik.ERA <- function(object){
	structure(object$LL, df=object$DF, class="logLik")
}

estfunW.ERA <- function(object, ...){
	require(sandwich)
    structure(object$scoreX, class="estfun")
}

estfunA.ERA <- function(object, ...){
	require(sandwich)
    structure(object$scoreF, class="estfun")
}
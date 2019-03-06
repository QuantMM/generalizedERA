Zgen <- function(N0=N, delt=delt, GN=GN){

	if (GN == 1) { 
		gratio <- c(1/3, 1/3, 1/3) # all equal group size
		} else if (GN == 2) {
		gratio <- c(1/2, 1/4, 1/4) # G1 = 50%
		} else {
		gratio <- c(2/3, 1/6, 1/6) # G1 = 66.67%
	}

	n_g1 <- ceiling(N0*gratio[1])
	n_g2 <- ceiling(N0*gratio[2])
	n_g3 <- ceiling(N0*gratio[3])
	N <- sum(n_g1, n_g2, n_g3)

	# True prob. for groups:
	P1 = n_g1/N				# P(Z1=T)=P1 & P(Z1=F)=q1
	P2 = n_g2/(n_g2+n_g3)	# P(Z2=T|Z1)=P2 & P(Z2=F|Z1)=q2
	p1 = P1 + (1-P1)*delt
	p2 = P2 + (1-P2)*delt
	
	# Generate Z1~Bin(p1), Z2~Bin(p2), Z3(random)
	Z1p <- rbinom(N,1,p1); Z1p <- Z1p[order(-Z1p)]
	Z2p <- rbinom(sum(n_g2+n_g3),1,p2) ; Z2p <- Z2p[order(-Z2p)]

	Z1 <- factor(Z1p) # Define "factor" for using mob() function
	Z2 <- factor( c( rep(1,n_g1) , Z2p ) )
	#Z2 <- factor( c( rbinom(n_g1, 1, 0.5) , Z2p ) )
	#Z3 <- rbinom(N, 1, 0.5)
	Z3 <- runif(N, min = 0, max = 1)
	
	if(delt==1) {
		G_id <- rep(1,N)
	} else {
		G_id <- c( rep(1,n_g1), rep(2,n_g2), rep(3,n_g3) )
	}
	
	Zdata <- data.frame(G_id=G_id,Z1=Z1,Z2=Z2,Z3=Z3)
	
	return(Zdata)
}

ERAgen <- function(corX=corX, b1=b1, b2=b2, npred=4, w1=c(0.7, 0.6, 0.5, 0.4), w2=c(0.6, 0.5, 0.4, 0.3), r=0){

	b_vec <- matrix(c(b1,b2),nrow=2) # b1=0.3 fixed
	
	# cov of predictors
	cx <- matrix(1,nrow=npred,ncol=npred)*corX - diag(1,npred)*corX + diag(1,npred)
	cx[1:2,3:4] <- cx[3:4,1:2] <- 0
	cz <- cx
	
	# cross-corr among predictors and Y
	W1 <- cbind(w1,0)
	W2 <- cbind(0,w2)
	W <- rbind(W1,W2)
	k <- 2 # ncomp
	CX <- kronecker(diag(k), cx)
	SdCX <- sqrt(t(W)%*%CX%*%W)
	Ws <- W *(1/SdCX[1])
	W_vec <- Ws[Ws!=0]
	
	crossCor <- CX%*%Ws%*%b_vec

	# pop. cov matrix
	COVX <- cbind(CX, crossCor)
	COVY <- cbind(t(crossCor),1)
	COV <- rbind(COVX,COVY)
	
	out.ERA_gen <- list(b_vec=b_vec,W_vec=W_vec,COV=COV)
	res <- out.ERA_gen
	res$b_vec
	res$W_vec
	res$COV
	invisible(res)
}

memPred <- function(Tree=Tree){

	nodeT <- predict(Tree, type="node") # predicted memberships
	MemT <- length(unique(nodeT))
	
	Tempo <- data.frame(MobG=row.names(coef(Tree)), sign(coef(Tree)))
	Tempo <- transform(Tempo, TrueG = ifelse(xF1==xF2, 1, ifelse(xF1 > xF2, 2, 3)))
	Tempo <- Tempo[ order(Tempo$TrueG), ]
	Tempo <- transform(Tempo, G_N = NA)
	for (i in 1: nrow(Tempo)) Tempo$G_N[i] <- sum(Tempo$MobG[i] == nodeT)
	memPred <- rep(Tempo[ ,"TrueG"],Tempo[ ,"G_N"])

	return(memPred)
}


ERA_simul <- function(y=y, X=X, nvar=nvar, dist = 1, const = 0, lambda = 0, it = 0, ceps = 0.0001){
	
	RawD = cbind(y,X)
	ndset = length(nvar)
	ncase = nrow(RawD)
	sum_nvar = sum(nvar)
	# data normalization:
	dataN <- scale(RawD)/sqrt(ncase-1)
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
	mu <- getmu(lp,dist)
	v <- getvmu(mu,dist)
	V <- diag(v)
	z <- lp + (y-mu)/v
	
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
		mu <- getmu(lp,dist)
		v <- getvmu(mu,dist)
		V <- diag(v)
		z <- lp + (y-mu)/v
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

	# Data output: F
	Fnames <- paste("F", c(1:length(A)) , sep="")
	F <- data.frame(F)
	colnames(F) <- Fnames
	ERAD <- cbind(y,X,F)

	output.ERA <- list(F=F, A=A, adj.DV=z)
	res <- output.ERA
	res$F
	res$A
	res$adj.DV
	invisible(res)
}

coef.ERA <- function(object){
	totalcoef <- rbind(object$vecW,object$A)
	Xv <- colnames(object$X)
	Fv <- paste("F", c(1:length(object$A)) , sep="")
	cvec <- as.vector(totalcoef)
	names(cvec) <- c(Xv,Fv)
	structure(cvec, class="coef")
}

logLik.ERA <- function(object){
	structure(object$LL, df=object$DF, class="logLik")
}

getmu <- function(lp, datatype){
	if (datatype == 1) {
		mu = lp
		
	} else if (datatype == 2) {
		mu = exp(lp) / (1 + exp(lp))
		nan = is.na(mu) | mu>1000
		if (sum(nan)!=0) {mu[which(nan == TRUE)] <- 1000}
	
	} else if (datatype == 3) {
		mu = exp(lp)
	
	} #else if (datatype == 4) {
	#	mu = 1/lp
	#	mu = 
	#}

	return(mu)
}

getvmu <- function(mu, datatype) {
	if (datatype == 1) {
		vmu = rep(1, nrow(mu))

	} else if (datatype == 2) {
	
	}
	return(vmu)
}

getdeviance <- function(y, mu, datatype) {
	if (datatype == 1) {
		DS = sum((y - mu)^2)

	} else if (datatype == 2) {
	
	}
	
	return(DS)
}

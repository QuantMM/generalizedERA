library(partykit);library(strucchange);library(sandwich);library(mvtnorm)
source("ERAfunction_simulation.r")

#---------------------
# Random ERA data generate: Becker, Rai, & Rigdon (2013)
# "Predictive validity and formative measurement in structural equation modeling: Embracing practical relevance."
#---------------------
ncomp <- 2 								# num. of components
r <- 0									# corr(F1,F2)
B1 <- 0.3								# reg coeff for F1
R2 <- 0.4								# R2 of endo- composite
B2 <- 0.5568							# obtained reg coeff for F2 based on (R2, b1, r)
corstr <- 1								# corr structure of predictors (1= nearby pair)
npred <- 4 								# num. of predictors per comp
corX <- 0								# corr among predictors
w1 <- c(0.7, 0.6, 0.5, 0.4)				# initial (i.e., unstandardized) values for F1 weights
w2 <- c(0.6, 0.5, 0.4, 0.3)				# initial values for F2 weights
Nvar <- npred*ncomp						# total num. of predictors
nvar_t <- Nvar+1						# total num. of observed variables

nrep <- 1000 							# num. of replications

#---------------------
# Design factors
#---------------------
corX <- 0.4								# corr among predictors; standardized weights will vary
gsize <- c(1, 2, 3) 					# 1=equal/balanced group size, 2=G1,50%, 3=G1,66.67%
delta <- 0				 				# (0=heterogeniety, 1=homogeneity, by 0.2)
startN <- c(60,90,120,180,300,500)		# sample sizes

#---------------------
# Combinations? for each senario
#---------------------
ERA_senario <- expand.grid(F1coef = B1, F2coef = B2, weightCoef = corX)
Tree_senario <- expand.grid(N = startN, GroupN = gsize)

lmfun <- function(y, x, start=NULL, weights=NULL, offset=NULL, ...){lm(y~0+x,...)}

Result <- vector("list", nrow(Tree_senario))
names(Result) <- paste("Tree_senario", 1:nrow(Tree_senario), sep="")

for (sn in 1:nrow(Tree_senario)) {

	startN <- Tree_senario[sn,"N"]
	gsize <- Tree_senario[sn,"GroupN"]
	
	# Create covariates (Z1, Z2, Z3) with true group memberships in "G_id"
	Zdata <- Zgen(N0=startN, delt=delta, GN=gsize)
	
	G1 <- Zdata[Zdata$G_id==1,]
	G2 <- Zdata[Zdata$G_id==2,]
	G3 <- Zdata[Zdata$G_id==3,]
	
	b10 <- ERA_senario[,"F1coef"]
	b20 <- ERA_senario[,"F2coef"]
	corrX <- ERA_senario[,"weightCoef"]
	
	# Save space for the outcome measures
	Result[[sn]]$Zsplit1 <- matrix(,ncol=nrep)
	Result[[sn]]$Power <- matrix(,ncol=nrep)
	Result[[sn]]$Accuracy <- matrix(,ncol=nrep)


	for (np in 1:nrep) { # 1000 replications
			# Group1
			ERAG1 <- ERAgen(corX=corrX, b1=b10, b2=b20)
			G1d <- MASS::mvrnorm(n = nrow(G1), mu = rep(0, nrow(ERAG1$COV)), Sigma = ERAG1$COV)
			G1y <- G1d[,9] # in COV or ERAgen(), the last column is for Y
			G1X <- G1d[,1:8]
			nvar <- matrix(c(4,4), ncol=2)
			G1ERA <- ERA_simul(y=G1y, X=G1X, nvar=nvar)
			G1Y <- G1ERA$adj.DV
			G1F <- G1ERA$F
			G1D <- cbind(G1Y,G1F) # test: lm(G1Y~0+F1+F2, data=G1D) & G1ERA$A
			colnames(G1D)[1] <- "Y" 

			# Group2
			ERAG2 <- ERAgen(corX=corrX, b1=b10, b2=-b20)
			G2d <- MASS::mvrnorm(n = nrow(G2), mu = rep(0, nrow(ERAG2$COV)), Sigma = ERAG2$COV)
			G2y <- G2d[,9] # in COV or ERAgen(), the last column is for Y
			G2X <- G2d[,1:8]
			G2ERA <- ERA_simul(y=G2y, X=G2X, nvar=nvar)
			G2Y <- G2ERA$adj.DV
			G2F <- G2ERA$F
			G2D <- cbind(G2Y,G2F) # test: lm(G2Y~0+F1+F2, data=G2D) & G2ERA$A
			colnames(G2D)[1] <- "Y" 
			
			# Group3
			ERAG3 <- ERAgen(corX=corrX, b1=-b10, b2=b20)
			G3d <- MASS::mvrnorm(n = nrow(G3), mu = rep(0, nrow(ERAG3$COV)), Sigma = ERAG3$COV)
			G3y <- G3d[,9] # in COV or ERAgen(), the last column is for Y
			G3X <- G3d[,1:8]
			G3ERA <- ERA_simul(y=G3y, X=G3X, nvar=nvar)
			G3Y <- G3ERA$adj.DV
			G3F <- G3ERA$F
			G3D <- cbind(G3Y,G3F) # test: lm(G3Y~0+F1+F2, data=G3D) & G3ERA$A
			colnames(G3D)[1] <- "Y" 
			
			# Final data and ERA
			ERAd <- data.frame(rbind(G1D, G2D, G3D))
			Data <- cbind(Zdata, ERAd)
			
			# MOB with maxdept=3
			MOBf <- mob(Y ~ 0+F1+F2 | Z1+Z2+Z3, data=Data, fit=lmfun,
					control = mob_control(ordinal = "L2", alpha=0.1, prune="AIC")
					)
			
			nodeP <- predict(MOBf, type="node")
			MemP <- length(unique(nodeP))		# predicted memberships
			MemT <- length(unique(Data$G_id))	# true memberships
			
			if (MemP == 1) { # no split, 1 group
				divided <- FALSE
				firstsplit <- NA
				CramV <- NA # the "q" below would be 1, leading to the 0 value for denominator
			} else {
				divided <- TRUE
			
				sct <- sctest(MOBf)$"1"
				sct <- sct[,complete.cases(sct["p.value",]),drop=FALSE]
				firstsplit <- colnames(sct)[ min(sct["p.value",]) == sct["p.value",]]
				if( length(firstsplit) != 1) firstsplit <- firstsplit[1]
						
				# Cramer's V
				for (i in 1:length(unique(nodeP))) nodeP[nodeP == unique(nodeP)[i]] <- (i+99)
				for (i in 1:length(unique(nodeP))) nodeP[nodeP == unique(nodeP)[i]] <- i
				crossT <- table(Data$G_id, nodeP)
				n <- sum(crossT)
				q <- min( nrow(crossT), ncol(crossT) )
				chis <- unname(chisq.test(crossT, correct=FALSE)$statistic)
				CramV <- sqrt( chis / (n*(q-1)))
			}
	
		# Outcomes
		Result[[sn]]$Zsplit1[1,np] <- firstsplit
		Result[[sn]]$Power[1,np] <- divided
		Result[[sn]]$Accuracy[1,np] <- CramV
		
		} # replications end
	
}


comp <- Result

saveM <- matrix(,nrow=nrow(Tree_senario),ncol=4,dimnames=list(NULL,c("whichZ","howmany","POWER","CramerV")))
saveM <- cbind(Tree_senario, saveM)

for (i in 1:nrow(saveM)){

	ttt <- table(comp[[i]]$Zsplit1)
	if ( length(ttt) == 0 ) {saveM[i,"whichZ"] <- NA
		} else {
		tt <- as.matrix(ttt[ttt==max(ttt)])
		saveM[i,"whichZ"] <- rownames(tt)
		}
	saveM[i,"howmany"] <- as.numeric(tt)
	saveM[i,"POWER"] <- mean(as.vector(comp[[i]]$Power))
	saveM[i,"CramerV"] <- mean(comp[[i]]$Accuracy, na.rm=T)

}

saveM



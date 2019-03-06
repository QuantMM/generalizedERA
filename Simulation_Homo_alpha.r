#---------------------
# Final Simulation Setups
#---------------------
library(partykit);library(strucchange);library(sandwich);library(mvtnorm)
source("ERAfunction_simulation.r")
lmfun <- function(y, x, start=NULL, weights=NULL, offset=NULL, ...){lm(y~0+x,...)} # for MOB

#---------------------
# Random ERA data generate: Becker, Rai, & Rigdon (2013)
# "Predictive validity and formative measurement in structural equation modeling: Embracing practical relevance."
#---------------------
ncomp <- 2 							# num. of components
r <- 0								# corr(F1,F2)
B1 <- 0.3							# reg coeff for F1
R2 <- 0.4							# R2 of endo- composite
B2 <- 0.5568							# obtained reg coeff for F2 based on (R2, b1, r)
corstr <- 1							# corr structure of predictors (1= nearby pair)
npred <- 4 							# num. of predictors per comp
corX <- 0							# corr among predictors
w1 <- c(0.7, 0.6, 0.5, 0.4)					# initial (i.e., unstandardized) values for F1 weights
w2 <- c(0.6, 0.5, 0.4, 0.3)					# initial values for F2 weights
Nvar <- npred*ncomp						# total num. of predictors
nvar_t <- Nvar+1						# total num. of observed variables

homo <- 1			 				# 1=homogeneity

#---------------------
# Simulation Design factors
#---------------------
gsize <- c(1, 2, 3) 						# 1=equal/balanced group size, 2=G1,50%, 3=G1,66.67%
startN <- c(60,90,120,180,300,500)				# sample sizes
# all possible senarios
Tree_senario <- expand.grid(N = startN, GroupN = gsize)

nrep <- 1000							# num. of replications

#---------------------
# Save the results
#---------------------
Result <- vector("list", nrow(Tree_senario))
names(Result) <- paste("Tree_senario", 1:nrow(Tree_senario), sep="")

#---------------------
# Simulation starts here
#---------------------
for (sn in 1:nrow(Tree_senario)) {

	startN <- Tree_senario[sn,"N"]
	gsize <- Tree_senario[sn,"GroupN"]
	
	# Create covariates (Z1, Z2, Z3) with true group memberships in "G_id"
	Zdata <- Zgen(N0=startN, delt=homo, GN=gsize)
	
	# Save space for the outcome measures
	Result[[sn]]$parval <- paste(startN, gsize, sep="/")
	Result[[sn]]$Zsplit1 <- matrix(,ncol=nrep)
	Result[[sn]]$Alpha <- matrix(,ncol=nrep)

	for (np in 1:nrep) { # 1,000 replications
	
		# Group1, 2, and 3: (B1, B2)
		ERA <- ERAgen(corX=corX, b1=B1, b2=B2)
		ERAd <- MASS::mvrnorm(n = nrow(Zdata), mu = rep(0, nrow(ERA$COV)), Sigma = ERA$COV)
		ERAy <- ERAd[,9] # in COV or ERAgen(), the last column is for Y
		ERAX <- ERAd[,1:8]
		nvar <- matrix(c(4,4), ncol=2)
		ERAfit <- ERA_simul(y=ERAy, X=ERAX, nvar=nvar)
		Data <- data.frame(Y=ERAfit$adj.DV,ERAfit$F,ERAX) # test: lm(Y~0+F1+F2, data=Data) & ERAfit$A

		Data <- data.frame(Zdata, Data)
		Data <- data.frame(Data, Pred_Node1 = NA)

		# Parameter Instability Tests: (1)Split or not? (2)Group memberships 
		
		# MOB with maxdept=2
		Node1 <- mob(Y ~ 0+F1+F2 | Z1+Z2+Z3, data=Data, fit=lmfun,
			control = mob_control(bonferroni = TRUE, ordinal = "chisq", maxdept=2)
		)

		Data$Pred_Node1 <- nodeP <- predict(Node1, type="node") # predicted memberships
		MemP <- length(unique(nodeP))
		
		if (MemP==1) { # no split

			firstsplit <- NA
			divided <- FALSE

		} else { # mistakenly split

			sct <- data.frame(sctest(Node1)[[1]])
			sct <- sct[,!is.na(colSums(sct)),drop=FALSE] # discard any Zs having NA
			Zselected <- which(min(sct["p.value",]) == sct["p.value",])
			# when ties in p-values
			if( length(Zselected) !=1 ){ firstsplit <- colnames(sct)[Zselected[1]] 
			} else { firstsplit <- colnames(sct)[Zselected] }
				
			divided <- TRUE
		}
		
	# Save results
	Result[[sn]]$Zsplit1[1,np] <- firstsplit
	Result[[sn]]$Alpha[1,np] <- divided

	}

}

#---------------------
# If you'd like to save the results in .RData
#---------------------
save(list = ls(all.names = TRUE), file = "AlphaSimulation.RData", envir = .GlobalEnv)

#---------------------
# Compute the performance measures
#---------------------

comp <- Result

saveM <- matrix(,nrow=nrow(Tree_senario), ncol=3, dimnames=list(NULL,c("whichZ","howmany","Alpha")))
saveM <- cbind(Tree_senario,saveM)

for (i in 1:nrow(saveM)){

	ttt <- table(comp[[i]]$Zsplit1)
	if ( length(ttt) == 0 ) {saveM[i,"whichZ"] <- NA
		} else {
		tt <- as.matrix(ttt[ttt==max(ttt)])
		saveM[i,"whichZ"] <- rownames(tt)
		}
	saveM[i,"howmany"] <- as.numeric(tt)
	saveM[i,"Alpha"] <- mean(as.vector(comp[[i]]$Alpha))

}




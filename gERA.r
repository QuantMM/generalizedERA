####################
## Basic Settings ##
####################

# (1) Load packages:
library(partykit)
library(strucchange)
library(sandwich)
library(mvtnorm)

# (2) Get the fitERA() function in your current R session
source("fitERA.r")
source("ERAfunction_simulation.r")

#####################
## Example Dataset ## Let's generate an example dataset for ERA
##################### Recall the figure on the front page!
B1 <- 0.3								# reg coeff for F1
B2 <- 0.5568							# obtained reg coeff for F2 based on (R2, b1, r)
corX <- 0								# corr among predictors
COVgenerate <- ERAgen(corX=corX, b1=B1, b2=B2)

startN <- 500							# sample sizes
ERAdata <- MASS::mvrnorm(n = startN, mu = rep(0, nrow(COVgenerate$COV)), Sigma = COVgenerate$COV)
ERAdata <- data.frame(ERAdata)
colnames(ERAdata)[9] <- "Y"				# in ERAgen(), the last column is for Y
head(ERAdata)

########################
## Data setup for ERA ##
########################
y <- ERAdata[,"Y",drop=FALSE]
x1 <- ERAdata[,c(paste0("X", 1:4)),drop=FALSE]
x2 <- ERAdata[,c(paste0("X", 5:8)),drop=FALSE]
X <- cbind(x1,x2)
nvar <- matrix( c(ncol(x1),ncol(x2)), ncol=2 )

#############
## Fit ERA ##
#############
fit <- ERA(y=y, X=X, nvar=nvar, dist=1) 

# To define family objective, refer to the below:
#	if(dist==1){ Yfam = "gaussian"
#	} else if (dist==2) { Yfam = "binomial" #(link = "logit")
#	} else if (dist==3) { Yfam = "poisson" #(link = "log")
#	} else if (dist==4) { Yfam = "Gamma" #(link = "inverse")
#	} else if (dist==5) { Yfam = "inverse.gaussian" #(link = "1/mu^2")
#	}

str(fit)
fit$A
fit$W	
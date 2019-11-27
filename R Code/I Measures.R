#-------------------------------------------------------------------
# I Measures
#	This contains code for computing I_PO and I_YP.
#
#	This code assumes that your data is in the the following form:
#		- Column 1 = survival time (time)
#		- Column 2 = censoring indicator (censor)
#		- Columns 3+ = genes
#
#	Output: 
#		I_PO, outPO (I, I test statistic, I p-value)
#		I_YP, outYP (I, I test statistic, I p-value)
#
#	Note:  There is separate code that adjusts for age
#		and stage in I_PO, if desired. 
#
# Required R Packages: survival, timereg, YPmodel
#-------------------------------------------------------------------

library(survival)
library(timereg)
library(YPmodel)

# Create the following data set if you'd like to run an example
# 	gg <- matrix(rnorm(2000, 0, 1), 200, 10)	# Generate gene expression
# 	stime <- rexp(200,1)	# Generate survival time
# 	cens <- sample(c(0,1), replace = T, 200)	# Generate censoring indicator
# 	dat <- as.data.frame(cbind(stime, cens, gg))
# 	colnames(dat)[1:2] <- c("time", "censor")

#------------------------------------------------------------------------
# I_PO 
#------------------------------------------------------------------------

# Fit the PO model and obtain the beta coefficients (code similar to what
#	is seen in 'Model Fits & GOF')

doPropOdds <- function(x){
	propOdds <- prop.odds(Event(time, censor) ~ x, data = dat)
	coeff <- coef(propOdds)[1]
	coeff_se <- coef(propOdds)[2]
	p <- coef(propOdds)[6]
	p_gof <- propOdds$pval.Prop

	return(c(coeff, coeff_se, p, p_gof))
	}

fitPO <- apply(dat[,3:ncol(dat)], 2, doPropOdds)
beta <- fitPO[1,]

#Function for obtaining the variance under the null

getVarBetaHo <- function(x){
	po0 <- prop.odds(Event(time, censor) ~ x, beta = 0, Nit = 1, data = dat)	
	sig2 <- po0$var.gamma [1]
	return(sig2)
	}

varBetaHo <- apply(dat[,3:ncol(dat)], 2, getVarBetaHo)  

allGenes <- cbind(varBetaHo,t(dat[,3:ncol(dat)]))
allGenes[allGenes == 0] <- .0000001

# Calculate I1_PO

I1_step1 <- apply(allGenes[,2:ncol(allGenes)], 2, function(x) x*beta-2+2*x*beta/(exp(x*beta)-1))
I1 <- apply(I1_step1, 1, function(x) sum(x, na.rm = T))
I1[I1 < 0] <- 0
I1[which(is.na(I1))] <- 0

# Calculate variance of I1_PO

varI1 <- apply(allGenes, 1, 
	function(x) sum(outer(2:length(x), 2:length(x),
      	function(i,j) 1/18 * x[i]^2 * x[j]^2 * x[1]^2   
	 ), na.rm = T))

#Calculate Test Statistic and p-value

testI1 <- I1^2/varI1
getP <- pchisq(testI1, 1, lower.tail = FALSE)

#Combine output into one summary data frame

outPO <- as.data.frame(cbind(I1_PO = I1, I1_PO_Test = testI1, I1_PO_pval = getP))

#------------------------------------------------------------------------
# I_PO - Adjusted for age/stage
#	Last 2 columns of data set must be age & stage
#------------------------------------------------------------------------

#Fit PO model adjusted for age & stage and retain beta coefficients

doPropOdds_AgeStage <- function(x){
	propOdds <- prop.odds(Event(time, censor) ~ x + Age + Stage, data = dat)
	coeff <- coef(propOdds)[1]
	coeff_se <- coef(propOdds)[6]
	p <- coef(propOdds)[26]
	p_gof <- propOdds$pval.Prop[1]

	return(c(coeff, coeff_se, p, p_gof))
	}

fitPO <- apply(dat[,3:(ncol(dat)-2)], 2, doPropOdds)
beta <- fitPO[1,]

#Function for obtaining the variance under the null

getVarBetaHo <- function(x){
	po0 <- prop.odds(Event(time, censor) ~ x + Age + Stage, beta = 0, Nit = 1, data = dat)	
	sig2 <- po0$var.gamma [1]
	return(sig2)
	}

varBetaHo <- apply(dat[,3:(ncol(dat)-2)], 2, getVarBetaHo)  

allGenes <- cbind(varBetaHo,t(dat[,3:(ncol(dat)-2)]))
allGenes[allGenes == 0] <- .0000001

# Calculate I1_PO

I1_step1 <- apply(allGenes[,2:ncol(allGenes)], 2, function(x) x*beta-2+2*x*beta/(exp(x*beta)-1))
I1 <- apply(I1_step1, 1, function(x) sum(x, na.rm = T))
I1[I1 < 0] <- 0
I1[which(is.na(I1))] <- 0

# Calculate variance of I1_PO

varI1 <- apply(allGenes, 1, 
	function(x) sum(outer(2:length(x), 2:length(x),
      	function(i,j) 1/18 * x[i]^2 * x[j]^2 * x[1]^2   
	 ), na.rm = T))

#Calculate Test Statistic and p-value

testI1 <- I1^2/varI1
getP <- pchisq(testI1, 1, lower.tail = FALSE)

#Combine output into one summary data frame

outPO <- as.data.frame(cbind(I1_PO = I1, I1_PO_Test = testI1, I1_PO_pval = getP))


#------------------------------------------------------------------------
# I_YP
#------------------------------------------------------------------------

# Dichotomize data 

split <- function(x){  
	y <- as.matrix(x) 
	cut(  
		y,
		breaks <- c(-Inf, median(y, na.rm = T), Inf),
		include.lowest = TRUE,
		labels = c("0", "1")
		)
	}

splitdat <- apply(dat[,3:ncol(dat)], 2, split)

# Functions for calculating I2_YP

I2func <- function(x) {
	(exp(beta*x)*(beta*x*exp(gamma*x)-gamma*x-2) - exp(gamma*x)*(gamma*x*exp(beta*x)-beta*x-2))/(exp(beta*x)-exp(gamma*x))
	}

#Obtain YP model estimates

beta <- c()
gamma <- c()

for(i in 1:ncol(splitdat)){	
	thisdat <- as.data.frame(cbind(V1 = as.numeric(dat[,1]), V2 = as.numeric(dat[,2]), V3 = as.numeric(splitdat[,i])))
	yp <- YPmodel(thisdat)
	yp2 <- yp$Estimate
	beta <- c(beta, yp2$beta[1])
	gamma <- c(gamma, yp2$beta[2])
}

# Calculate I2_YP 
#	Note: change the 0/1 split to 1/2 for I_YP calculation

newsplit <- ifelse(splitdat == 0, 1, 2)
allGenes <- t(newsplit)
I2_step1 <- apply(allGenes[,1:ncol(allGenes)], 2, I2func)
I2_YP <- apply(I2_step1, 1, function(x) sum(x, na.rm = T))

# Calculate variance of I2_YP

varbeta <- c()
vargamma <- c()
for(i in 1:ncol(splitdat)){
	pickthese <- which(!is.na(splitdat[,i]))
	thisdat <- as.data.frame(cbind(V1 = as.numeric(dat[pickthese,1]), V2 = as.numeric(dat[pickthese,2]), V3 = as.numeric(splitdat[pickthese,i])))
	yp <- YPmodel.estimate(thisdat, startPoint = c(0,0), nm = c(0,0), maxIter1 = 0, maxIter2 = 0 )
	varbeta <- c(varbeta, yp$variance.beta1^2)
	vargamma <- c(vargamma, yp$variance.beta2^2)
}

var_step1 <- apply(allGenes, 1, function(x) sum(outer(1:length(x), 1:length(x),
      	function(i,j) x[i]^2 * x[j]^2 ), na.rm = T))

varI2 <- 1/36*var_step1*(2*varbeta^2+varbeta*vargamma+2*vargamma^2)

# Calculate test statistic & p-value

testI2_YP <- I2_YP^2/varI2
pval_I2_YP <- pchisq(testI2_YP, 1, lower.tail = FALSE)

#Combine output into one summary data frame

outYP <- as.data.frame(cbind(YPgamma = beta, YPbeta = gamma, I2_YP = I2_YP, I2_YP_Test = testI2_YP, I2_YP_pval = pval_I2_YP))


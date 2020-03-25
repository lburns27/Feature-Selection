#----------------------------------------------------------------------------
# EXAMPLE - This example does the following: 
#	1.  Creates a simulated data set (Scheme 1 - W, 33% censoring)
#	2.  Obtains PH, PO & YP model fits and GOF
#	3.  Computes all proposed measures for feature selection:
#		I measures - I_PO, I_YP
#		R2_measures - R2_LR, R2_IPO, R2_IPH, R2_PO, R2_CO, R2_ModCH
#		R2_measures (by Rouam et al. 2010, 2011) - R2_CH, R2_PH
#	4.  Computes other existing measures (concreg, Uno's C, R2_G & R2_SH)
#	5.  Computes Sensitivity, Specificity, Youden & AUC for each measure.
#	6.  Creates venn diagrams showing overlaps between measures.
#
#	Note:  Before running this example, some functions need to be run
#		from the other R code in this repository.  All required functions
#		are noted thoughout below.
#
#	Required packages: survival, timereg, YPmodel, concreg, survAUC, pec, 
#				MESS, gplots, VennDiagram, latex2exp
#-----------------------------------------------------------------------------

library(survival)
library(timereg)
library(YPmodel)
library(concreg)
library(survAUC)
library(pec)
library(MESS)
library(gplots)
library(VennDiagram)
library(latex2exp)

#-----------------------------------------------------------------------------
# Run simulation function of interested.  Example: makeDataW function in Scheme 1
# Create data set with 200 significant genes (100 small effect, 100 large effect)
#	and 500 insignificant genes.  Use 33% censoring.
#-----------------------------------------------------------------------------

numSmall <- 100
numLarge <- 100
notsig <- 500
N <- 100

dat <- makeDataW(N, numSmall, numLarge, 5.3, notsig, 1,.5)

#-----------------------------------------------------------------------------
# Run functions for model fitting:  doCox, doPropOdds, split 
#	(from 'Model Fits & GOF' code)
#-----------------------------------------------------------------------------

#Fit Cox model
fitCox <- t(apply(dat[,3:ncol(dat)], 2, doCox))   

#Fit PO model
fitPO <- t(apply(dat[,3:ncol(dat)], 2, doPropOdds))
colnames(fitPO)[4] <- "PO_GOF_p"

# Split data
	splitdat <- apply(dat[,3:ncol(dat)], 2, split)

# Fit YP model via loop (code from 'Model Fits & GOF')
	yp_beta <- c()
	yp_gamma <- c()
	yp_beta_se <- c()
	yp_gamma_se <- c()
	yp_p_sig <- c()
	yp_p_fit <- c()

	for(i in 1:ncol(splitdat)){	
		thisdat <- as.data.frame(cbind(V1 = as.numeric(dat[,1]), V2 = as.numeric(dat[,2]), V3 = as.numeric(splitdat[,i])))
		yp <- YPmodel(thisdat)
		yp2 <- yp$Estimate
		yp_beta <- c(yp_beta, yp2$beta[1])
		yp_gamma <- c(yp_gamma, yp2$beta[2])
		yp_beta_se <- c(yp_beta_se, sqrt(abs(yp2$variance.beta1)))
		yp_gamma_se <- c(yp_gamma_se, sqrt(abs(yp2$variance.beta2)))
		yp_p_sig <- c(yp_p_sig, yp$Adlgrk$pval)
		yp_p_fit <- c(yp_p_fit, yp$LackFitTest$pvalu1)

	}

outYP <- as.data.frame(cbind(YP_gamma = yp_gamma, YP_beta = yp_beta, YPbeta_se = yp_beta_se, 
		YPgamma_se = yp_gamma_se, YP_p_sig = yp_p_sig, YP_p_fit = yp_p_fit))

#-----------------------------------------------------------------------------
# Compute I measures (I_PO & I_YP).  Code copied from 'I measures'
#-----------------------------------------------------------------------------

# I_PO -------------------------------------------------------------------------

	beta <- fitPO[,1]

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

	outIPO <- as.data.frame(cbind(I_PO = I1, I_PO_Test = testI1, I_PO_p = getP))

# I_YP -------------------------------------------------------------------------

	# Functions for calculating I2_YP

	I2func <- function(x) {
		(exp(yp_beta*x)*(yp_beta*x*exp(yp_gamma*x)-yp_gamma*x-2) - exp(yp_gamma*x)*(yp_gamma*x*exp(yp_beta*x)-yp_beta*x-2))/(exp(yp_beta*x)-exp(yp_gamma*x))
		}

	# Calculate I2_YP.  Note: change the 0/1 split to 1/2 for I_YP calculation.

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

	outIYP <- as.data.frame(cbind(I_YP = I2_YP, I_YP_Test = testI2_YP, I_YP_p = pval_I2_YP))

#-----------------------------------------------------------------------------
# Compute R2_Measures
#-----------------------------------------------------------------------------

# R2_LR & R2_I Measures: Run functions r2LR and r2_I from 'R2_LR & R2_I Measures'

	R2_LR <- apply(dat[,3:ncol(dat)], 2, function(x) r2LR(x, dat$time, dat$censor))
	R2_I <- t(apply(dat[,3:ncol(dat)], 2, function(x) r2_I(x, dat$time, dat$censor)))

# Pseudo R2 Measures: Run functions Index.PO, Index.CO & Index.ModCH from 'Pseudo R2 Measures'

	R2_PO <- apply(dat[,3:ncol(dat)], 2, function(x) Index.PO(x, dat$time, dat$censor))
	R2_CO <- apply(dat[,3:ncol(dat)], 2, function(x) Index.CO(x, dat$time, dat$censor))
	R2_ModCH <- apply(dat[,3:ncol(dat)], 2, function(x) Index.modCH(x, dat$time, dat$censor))

# Existing Pseudo R2 Measures: Run functions Index.PH & Index.CH from 'Pseudo R2 Measures'

	R2_CH <- apply(dat[,3:ncol(dat)], 2, function(x) Index.CH(x, dat$time, dat$censor))
	R2_PH <- apply(dat[,3:ncol(dat)], 2, function(x) Index.PH(x, dat$time, dat$censor))

outR2 <- cbind(R2_LR, R2_I, R2_PO, R2_CO, R2_ModCH, R2_CH, R2_PH)

#------------------------------------------------------------------------------
# Compute Other Existing Measures
#	Note - First run the following functions from 'Other Existing Measures':
#		predictSurvProb.propodds, doSurev, doExisting
#------------------------------------------------------------------------------

# R2_SH

	R2_SH <- apply(dat[,3:ncol(dat)], 2, function(x) doSurev(dat$time, dat$censor, x))

# Concreg, Uno's C & R2_G (PH & PO)

	other <- t(apply(dat[,3:ncol(dat)], 2, function(x) doExisting(dat$time, dat$censor, x)))
	colnames(other)[1] <- "ConcEffect"

outOther <- cbind(R2_SH, other)

#------------------------------------------------------------------------------
# Combine all results & measures into one data frame
#------------------------------------------------------------------------------

allOut <- as.data.frame(cbind(fitCox, fitPO, outYP, outIPO, outIYP, outR2, outOther))
head(allOut)


#------------------------------------------------------------------------------
# Compute Sensitivity, Specificity, Youden, & AUC for each measure
#	Note - Run getResults function in 'Youden & AUC' first.
#------------------------------------------------------------------------------

# Specify measures of interest

	measures <- c("I_YP", "I_PO", "ConcEffect", "R2_LR", "R2_IPO", "R2_IPH", "R2_PO",
		"R2_CO", "R2_ModCH", "R2_CH", "R2_PH", "R2G_PH", "R2G_PO", "R2_SH")              

# Calcaulte Sensitivity, Specificity, Youden, & AUC for numbers of interest
	
	allResults <- c()
	roww <- 0
	for(i in which(colnames(allOut) %in% measures)){
		roww <- roww + 1
		allResults <- rbind(allResults, getResults(allOut[,i], 200))
		rownames(allResults)[roww] <- colnames(allOut)[i]
	}

	rownames(allResults) <- measures


#------------------------------------------------------------------------------
# Create Venn Diagrams (code copied from 'Venn Diagrams'
#------------------------------------------------------------------------------

# Rank genes by each measure ------------------------------------------------

	out <- allOut

	# Choose number of top selected genes to examine
	n <- 50

	r2po <- rownames(out[order(out$R2_PO, decreasing = T),])[1:n]
	r2lr <- rownames(out[order(out$R2_LR, decreasing = T),])[1:n]
	r2ipo <- rownames(out[order(out$R2_IPO, decreasing = T),])[1:n]
	r2co <- rownames(out[order(out$R2_CO, decreasing = T),])[1:n]
	r2modch <- rownames(out[order(out$R2_ModCH, decreasing = T),])[1:n]
	i1po <- rownames(out[order(out$I_PO, decreasing = T),])[1:n]
	concreg <- rownames(out[order(out$ConcEffect, decreasing = T),])[1:n]
	out$Conc2 <- .5 + abs(out$ConcEffect -.5)	#Ensure that the absolute effect size is being used
	concreg2 <- rownames(out[order(out$Conc2, decreasing = T),])[1:n]
	i2yp <- rownames(out[order(out$I_YP, decreasing = T),])[1:n]

	# Create Venn Diagrams (continuous expression) -------------------------------------

	listR2PO <- list(r2po, r2ipo, r2lr)
	venn.plot <- venn.diagram(listR2PO , NULL, category.names=c(TeX('$\\textit{R^2_{PO}}$'), TeX('\\textit{$R^2_{\\tilde{I}_{PO}}}$'), TeX('$\\textit{R^2_{LR}}$')), 
			rotation.degree =60, cat.cex = 2.2, cex = 2.2, lwd = 3, cat.pos = c(-130,0,140), cat.dist = c(.015, .05, .045), scaled = F, euler.d = F)
	grid.draw(venn.plot)

	listR2 <- list(r2po, r2modch, r2co)
	venn.plot <- venn.diagram(listR2 , NULL, category.names=c(TeX('$\\textit{R^2_{PO}}$'), TeX('$\\textit{R^2_{ModCH}}$'), TeX('$\\textit{R^2_{CO}}$')), 
			rotation.degree =60, cat.cex = 2.2, cex = 2.2, lwd = 3, cat.pos = c(-130,0,140), cat.dist = c(.015, .05, .045), scaled = F, euler.d = F)
	grid.draw(venn.plot)

	listOther <- list(i1po, i2yp, concreg2)
	venn.plot <- venn.diagram(listOther , NULL, category.names=c(TeX('$\\textit{I_{PO}}$'), TeX('$\\textit{I_{YP}}$'), TeX('$\\hat{c}^\'_{+}$')), 
			rotation.degree =60, cat.cex = 2.2, cex = 2.2, lwd = 3, cat.pos = c(-130,0,140), cat.dist = c(.015, .04, .045), scaled = F, euler.d = F)
	grid.draw(venn.plot)


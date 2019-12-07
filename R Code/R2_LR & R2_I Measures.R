#-------------------------------------------------------------------
# R2_LR & R2_I Measures
#	This code contains functions for R2_LR, R2_IPO & R2_IPH.
#		r2LR: Computes R2_LR
#		r2_I: Computes R2_IPO & R2_IPH
#	Input: gene, survival time, censoring indicator
#	Output: R2 values 
#
#	Note:  There are separate functions that adjust for age
#		and stage, if desired. 
#
#	*Examples at end of code
#
# Required R Packages: survival, timereg
#-------------------------------------------------------------------

library(survival)
library(timereg)

#-------------------------------------------------------------------
# R2_LR (O'Quigley)
#-------------------------------------------------------------------

# NOT adjusted for age and stage 

r2LR <- function(gene, tt, cens){

	df <- as.data.frame(cbind(time = tt, censor = cens, x = gene))

	poNull <- prop.odds(Event(time, censor) ~ x, beta = 0, Nit = 1, data = df)
	po <- prop.odds(Event(time, censor) ~ x, data = df)
	k <- sum(cens)

	logNull <- poNull$loglik
	logBeta <- po$loglik

	quigley <- 1-exp(-2/k*(logBeta - logNull))
	
	return(R2_LR = quigley)
}

# Adjusted for age and stage

r2LR_AgeStage <- function(gene, tt, cens, age, stage){

	df <- as.data.frame(cbind(time = tt, censor = yetat, x = gene))

	poNull <- prop.odds(Event(time, censor) ~ x + age + stage, beta = 0, Nit = 1, data = df)
	po <- prop.odds(Event(time, censor) ~ x + age + stage, data = df)
	k <- sum(cens)

	logNull <- poNull$loglik
	logBeta <- po$loglik

	quigley <- 1-exp(-2/k*(logBeta - logNull))
	
	return(R2_LR = quigley)
}

#-------------------------------------------------------------------
# R2_LR (O'Quigley)
#-------------------------------------------------------------------

# NOT adjusted for age and stage

r2_I <- function(g, tt, cens){

	df <- as.data.frame(cbind(time = tt, censor = cens, g = g))

	po <- prop.odds(Event(tt, cens) ~ scale(g), data = df)
	cox <- coxph(Surv(tt, cens) ~ scale(g))	
	coxSum <- summary(cox)
	po_coeff <- coef(po)[1]
	cox_coeff <- coxSum$coeff[1]

	ipo <- 1-exp(-1/3*po_coeff^2)
	iph <- 1-exp(-2*(exp(cox_coeff^2/2)-1))
	
	return(c(R2_IPO = ipo, R2_IPH = iph))
}

# Adjusted for age and stage

r2_I_AgeStage <- function(g, tt, cens, age, stage){

	df <- as.data.frame(cbind(time = tt, censor = cens, g = g))

	po <- prop.odds(Event(tt, cens) ~ scale(g) + age + stage, data = df)
	cox <- coxph(Surv(tt, cens) ~ scale(g) + age + stage)	
	coxSum <- summary(cox)
	po_coeff <- coef(po)[1]
	cox_coeff <- coxSum$coeff[1]

	ipo <- 1-exp(-1/3*po_coeff^2)
	iph <- 1-exp(-2*(exp(cox_coeff^2/2)-1))
	
	return(c(R2_IPO = ipo, R2_IPH = iph))
}


#-------------------------------------------------------------------
# Examples
#-------------------------------------------------------------------

gg <- rnorm(200,0,1)	# Generate gene expression
stime <- rexp(200,1)	# Generate survival time
cens <- sample(c(0,1), replace = T, 200)	# Generate censoring indicator

#R2_LR
r2LR(gg, stime, cens)

#R2_IPO & R2_IPH
r2_I(gg, stime, cens)



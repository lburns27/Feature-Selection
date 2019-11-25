#----------------------------------------------------------------------
# Includes functions for fitting the following models:
#	Cox (R package: survival)
#	PO (R package: timereg)
#	YP (R package: YPmodel)
#
#	Input: predictor, x (i.e. gene expression)
#	Output: 
#		coeff - model coefficient (beta)
#		se - coefficient standard error (
#		p - significance p-value
#		p_gof - GOF p-value
#
#	Note: The Cox & PO models have separate functions adjusting for
#		age & stage, if desired. 
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Cox & PO adjusted for age & stage
#----------------------------------------------------------------------

doPropOdds_AgeStage <- function(x){
	propOdds <- prop.odds(Event(time, censor) ~ x + Age + Stage, data = dat)
	coeff <- coef(propOdds)[1]
	coeff_se <- coef(propOdds)[6]
	p <- coef(propOdds)[26]
	p_gof <- propOdds$pval.Prop[1]

	return(c(coeff, coeff_se, p, p_gof))
	}

doCox_AgeStage <- function(x){
	cox <- coxph(Surv(time, censor) ~ x + Age + Stage, data = dat)	
	gof <- cox.zph(cox)
	coxSum <- summary(cox)
	coeff <- coxSum$coefficients[1]
	coeff_se <- coxSum$coefficients[11]
	p <- coxSum$coefficients[21]
	p_gof <- gof$table[13]

	return(c(Coeff = coeff, coeff_se,  p, p_gof))
	}

# Example: dat (column 1 = time, column 2 = censor, columns 3+ = gene expression)
#	po <- apply(dat[,3:ncol(dat), 2, doPropOdds_AgeStage)
# 	ph <- apply(dat[,3:ncol(dat), 2, doCox_AgeStage)

#----------------------------------------------------------------------
# Cox & PO (NOT adjusted for age & stage)
#----------------------------------------------------------------------


doCox <- function(x){
	cox <- coxph(Surv(time, censor) ~ x, data = dat)	
	gof <- cox.zph(cox)
	coxSum <- summary(cox)
	coeff <- coxSum$coefficients[1]
	coeff_se <- coxSum$coefficients[3]
	p <- coxSum$coefficients[5]
	p_gof <- gof$table[3]

	return(c(Coeff = coeff, coeff_se,  p, p_gof))
	}

doPropOdds <- function(x){
	propOdds <- prop.odds(Event(time, censor) ~ x, data = dat)
	coeff <- coef(propOdds)[1]
	coeff_se <- coef(propOdds)[2]
	p <- coef(propOdds)[6]
	p_gof <- propOdds$pval.Prop

	return(c(coeff, coeff_se, p, p_gof))
	}

# Example: dat (column 1 = time, column 2 = censor, columns 3+ = gene expression)
#	po <- apply(dat[,3:ncol(dat), 2, doPropOdds)
# 	ph <- apply(dat[,3:ncol(dat), 2, doCox)

#----------------------------------------------------------------------
# YP model
#	1. Data must be dichotomized. We create a function 'split' to split
#		data by median.
#	2. apply() function will not work here, so we create a loop
#----------------------------------------------------------------------

# Split function

	split <- function(x){  
		y <- as.matrix(x) 
		cut(  
			y,
			breaks <- c(-Inf, median(y, na.rm = T), Inf),
			include.lowest = TRUE,
			labels = c("0", "1")
			)
		}

# Split data
	splitdat <- apply(dat[,3:ncol(dat)], 2, split)

# Fit YP model via loop
	beta <- c()
	gamma <- c()
	beta_se <- c()
	gamma_se <- c()
	p_sig <- c()
	p_fit <- c()

	for(i in 1:ncol(splitdat)){	
		thisdat <- as.data.frame(cbind(V1 = as.numeric(dat[,1]), V2 = as.numeric(dat[,2]), V3 = as.numeric(splitdat[,i])))
		yp <- YPmodel(thisdat)
		yp2 <- yp$Estimate
		beta <- c(beta, yp2$beta[1])
		gamma <- c(gamma, yp2$beta[2])
		beta_se <- c(beta_se, sqrt(abs(yp2$variance.beta1)))
		gamma_se <- c(gamma_se, sqrt(abs(yp2$variance.beta2)))
		p_sig <- c(p_sig, yp$Adlgrk$pval)
		p_fit <- c(p_fit, yp$LackFitTest$pvalu1)

	}


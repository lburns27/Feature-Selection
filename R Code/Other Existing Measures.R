#--------------------------------------------------------------------------
# This code contains functions for computing the following existing
#	measures in the literature:
#		Concreg (Dunkler et al. 2010)
#		Uno's C (Uno et al. 2011)
#		R2_G - PH & PO cases (Graf et al. 1999; Gerds & Schumacher 2006)
#		R2_SH (Schemper & Henderson 2000)

#	Input: gene, survival time, censoring indicator
#	Output: Concreg (absolute effect size), Uno's C, R2_G (PH), R2_G (PO)
#			and R2_SH.
#
#	*Examples at end of code
#
# Required R Packages: concreg, survAUC, survival, pec, timereg
#--------------------------------------------------------------------------
#Load packages.................
library(survival)
library(concreg)
library(survAUC)
library(timereg)
library(pec)

#------------------------------------------------------------------------------
# Predict function based on PO model (used in pec function later)
#------------------------------------------------------------------------------

predictSurvProb.propodds <- function(object,newdata,times,...){
 	#  require(timereg)
  	##  The time-constant effects first
  	const <- c(object$gamma)
  	#  names(const) <- substr(dimnames(object$gamma)[[1]],6,nchar(dimnames(object$gamma)[[1]])-1)
  	names(const) <-dimnames(object$gamma)[[1]]
  	constant.part <- t(newdata[,names(const)])*const
  	constant.part <- exp(colSums(constant.part))
  	##  Then extract the time-varying effects
  	time.coef <- data.frame(object$cum)
  	ntime <- nrow(time.coef)
  	objecttime <- time.coef[,1,drop=TRUE]
  	ntimevars <- ncol(time.coef)-2
  	time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
  	nobs <- nrow(newdata)
  	time.part <- .C("survest_cox_aalen",timehazard=double(ntime*nobs),as.double(unlist(time.coef[,-1])),as.double(unlist(time.vars)),as.integer(ntimevars+1),as.integer(nobs),as.integer(ntime),PACKAGE="pec")$timehazard
  	time.part <- matrix(time.part,ncol=ntime,nrow=nobs)
  	## dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
  	surv <- pmin(exp(-time.part*constant.part),1)
  	if (missing(times)) times <- sort(unique(objecttime))
 	p <- surv[,prodlim::sindex(objecttime,times)]
  	if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    	stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
 	p
}

#------------------------------------------------------------------------------
# Function for calculation R2_SH
#------------------------------------------------------------------------------

doSurev <- function(tt, cens, x, type = "kaplan-meier")

{
	coxfit <- coxph(Surv(tt, cens) ~ x) 

	status <- coxfit$y[, 2]
	times <- coxfit$y[, 1]
	km <- survfit(Surv(times, status)~1, type = type)
	km.inv <- survfit(Surv(time = coxfit$y[, 1], event = abs(coxfit$y[, 2] - 1))~1)
	orderedtimes <- km$time[km$n.event != 0]
	km.curve <- km$surv[km$n.event != 0]
	basefit <- survfit(coxfit, type = type)$surv[status == 1]
	G <- km.inv$surv
	G.inv <- 1/G[km$n.event != 0]
	dj <- km$n.event[km$n.event != 0]
	w <- sum(G.inv * dj)	

	## Next line calculates predicted survival for every subject at every failure time.
	## The results are stored in a matrix, each line representing survival for one subject.

	Stx <- matrix(rep(basefit, each = length(coxfit$linear))^exp(coxfit$linear), byrow = F, nrow = length(coxfit$linear))
	oneszeros <- matrix(1 * (rep(times, each = length(basefit)) > orderedtimes), byrow = T, ncol = length(basefit))
	ties <- matrix(1 * (rep(times, each = length(basefit)) == orderedtimes), byrow = T, ncol = length(basefit)) * abs(status - 1)
	oneszeros <- oneszeros + ties
	Sti <- Stx * oneszeros
	Sti <- (Sti == 0) + Sti
	Stix <- apply(Sti, 1, min)
	extrapolation <- Stx + (Stx/Stix) - 2 * (Stx^2/Stix)
	extrapolation <- extrapolation * abs(status - 1)
	extrapolation <- extrapolation * abs(oneszeros - 1)
	selection <- 1 * (extrapolation == 0)	

	## We now calculate M_hat at every time point. First we need its components (of
	## which we later take the average). Of these we select those that need not to be
	## extrapolated and we than add the extrapolated parts

	M.comp <- (abs(oneszeros - Stx)) * selection + extrapolation	

	## M_hat is obtained by averaging columns of M.comp

	M.hat <- apply(M.comp, 2, mean)	

	## We also have to calculate M_hat under the null model (Kaplan-Meier)

	St0 <- matrix(rep(km.curve, length(coxfit$linear)), byrow = T, ncol = length(basefit))
	St0i <- St0 * oneszeros
	St0i <- (St0i == 0) + St0i
	St0ix <- apply(St0i, 1, min)
	extrapolation0 <- St0 + (St0/St0ix) - 2 * (St0^2/St0ix)
	extrapolation0 <- extrapolation0 * abs(status - 1)
	extrapolation0 <- extrapolation0 * abs(oneszeros - 1)
	selection0 <- 1 * (extrapolation0 == 0)
	M0.comp <- (abs(oneszeros - St0)) * selection0 + extrapolation0
	M0.hat <- apply(M0.comp, 2, mean)

	D <- sum(G.inv * dj * M0.hat)/w
	Dx <- sum(G.inv * dj * M.hat)/w

 	return(c("V" = (D - Dx)/D))
}

#------------------------------------------------------------------------------
# Function for calculating concreg, Uno's C & R2_G
#------------------------------------------------------------------------------


doExisting <- function(tt, cens, gene){

	dd2 <- as.data.frame(cbind("time" = tt, "status" = cens, "g" = gene))
  	dd <- dd2[order(dd2$time),]

	#Concreg 
	concMod <- concreg(data = dd, Surv(time, status) ~ g)
	conc1 <- exp(concMod$coef[1])/(1+exp(concMod$coef[1]))
	conc2 <- .5 + abs(conc1-.5)
	
	#Uno's C
	uu <- UnoC(Surv(tt, cens), Surv(tt, cens), gene)
	UC <- .5 + abs(uu - .5)

	#R2_G
	tobs <- dd$time[which(dd$status == 1)]
  	fitCox <- coxph(Surv(time, status) ~ g, data=dd, x = TRUE)
  	fitPO <- prop.odds(Event(time, status) ~ g, data=dd)
  	class(fitPO) <- "propodds"

  	grafCox <- ibs(pec(list(fitCox), Surv(time, status) ~ 1, times = tobs, data = dd, cens.model = "marginal"), times = max(tobs))
  	grafPO <- ibs(pec(list(fitPO), Surv(time, status) ~ 1, times = tobs, data = dd, cens.model = "marginal"), times = max(tobs))
  
  	r2grafC <- 1-grafCox[2]/grafCox[1]
  	r2grafP <- 1-grafPO[2]/grafPO[1]

	allMeasures <- c("ConcEffect" = conc2, "UnoC" = UC, "R2G_PH" = r2grafC, "R2G_PO" = r2grafP)

	return(allMeasures)
}

#-------------------------------------------------------------------
# Examples
#-------------------------------------------------------------------

gg <- rnorm(200,0,1)	# Generate gene expression
stime <- rexp(200,1)	# Generate survival time
cens <- sample(c(0,1), replace = T, 200)	# Generate censoring indicator

#R2_SH
doSurev(stime, cens, gg)

#Concreg, Uno's C, R2_G (PH & PO)
doExisting(stime, cens, gg)


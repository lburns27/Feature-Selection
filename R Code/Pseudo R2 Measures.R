#-------------------------------------------------------------------
# Pseudo-R2 Measures
#	Cases: PO, CO, CH, ModCH, PH
#	Input: gene, survival time, censoring indicator
#
#	Note: Examples at end of code
#
# Required R Packages: survival
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# Index.PO(): R2_PO
#-------------------------------------------------------------------

Index.PO <- function(gene, tt, cens){

	n <- length(tt)
	o <- order(tt)

	s.reg <- summary(survfit(Surv(tt, cens)~1))
	ni <- s.reg$n.risk

	tps.deces <- s.reg$time
	delti <- cens

	mat <- as.data.frame(matrix(nrow=n, ncol=7))
	names(mat) <- c( "tps" , "n.risk", "delta", "x.event", 
			"x.deces", "x.risk", "poids" )
	mat$tps <- round(tt[o], 2)
	mat$n.risk <- n:1
	mat$delta <- delti[o]
	geneo <- gene[o]
	etato <- cens[o]
	mat$x.deces <- as.numeric(geneo*etato)

	mat$x.event <- as.numeric(geneo)
	mat$x.risk <- rev(cumsum(rev(as.numeric(geneo))))

	St0 <- cumprod((mat$n.risk-mat$delta)/mat$n.risk)
	omega <- c(1, St0[-length(St0)])
	mat$poids <- omega*mat$delta

	Y <- matrix(1, ncol=n, nrow=n)
	Y[upper.tri(Y)] <- 0

	Ui.tmp <- mat$delta*(mat$x.deces -(mat$x.risk/mat$n.risk))
	Ui <- Ui.tmp*mat$poids
	
	v1 <- mat$delta/mat$n.risk
	m1 <- matrix(v1, nrow=n, ncol=n, byrow=T)*Y
	n1 <- m1*mat$x.event
	v2 <- mat$delta/(mat$n.risk*mat$n.risk)
	m2 <- matrix(v2, nrow=n, ncol=n, byrow=T)*Y
	n2 <- matrix(mat$x.risk, ncol=length(mat$tps) ,nrow=length(mat$tps), byrow=T)*m2
	EUij.tmp <- n1-n2
	EUij <- EUij.tmp*matrix(mat$poids, ncol=length(mat$tps), nrow=length(mat$tps), byrow=T)
	EUi <- apply(EUij, 1, sum)
	Wi <- Ui - EUi
	Cox <- sum(Wi)^2/sum(Wi^2)
	Ind <- Cox/length(s.reg$time)
	return(Ind)
	}

#-------------------------------------------------------------------
# Index.CO(): R2_CO
#-------------------------------------------------------------------

Index.CO <- function(gene, tt, cens){

	n <- length(tt)
	o <- order(tt)

	s.reg <- summary(survfit(Surv(tt, cens)~1))
	ni <- s.reg$n.risk

	tps.deces <- s.reg$time
	delti <- cens

	mat <- as.data.frame(matrix(nrow=n, ncol=7))
	names(mat) <- c( "tps" , "n.risk", "delta", "x.event", 
			"x.deces", "x.risk", "poids" )
	mat$tps <- round(tt[o], 2)
	mat$n.risk <- n:1
	mat$delta <- delti[o]
	geneo <- gene[o]
	etato <- cens[o]
	mat$x.deces <- as.numeric(geneo*etato)

	mat$x.event <- as.numeric(geneo)
	mat$x.risk <- rev(cumsum(rev(as.numeric(geneo))))

	St0 <- cumprod((mat$n.risk-mat$delta)/mat$n.risk)
	omega <- c(1, St0[-length(St0)])
	w <- ifelse(omega == 1, -7.2, 1 + omega - omega*log(omega/(1-omega)))
	mat$poids <- w*mat$delta

	Y <- matrix(1, ncol=n, nrow=n)
	Y[upper.tri(Y)] <- 0

	Ui.tmp <- mat$delta*(mat$x.deces -(mat$x.risk/mat$n.risk))
	Ui <- Ui.tmp*mat$poids
	
	v1 <- mat$delta/mat$n.risk
	m1 <- matrix(v1, nrow=n, ncol=n, byrow=T)*Y
	n1 <- m1*mat$x.event
	v2 <- mat$delta/(mat$n.risk*mat$n.risk)
	m2 <- matrix(v2, nrow=n, ncol=n, byrow=T)*Y
	n2 <- matrix(mat$x.risk, ncol=length(mat$tps) ,nrow=length(mat$tps), byrow=T)*m2
	EUij.tmp <- n1-n2
	EUij <- EUij.tmp*matrix(mat$poids, ncol=length(mat$tps), nrow=length(mat$tps), byrow=T)
	EUi <- apply(EUij, 1, sum)
	Wi <- Ui - EUi
	Cox <- sum(Wi)^2/sum(Wi^2)
	Ind <- Cox/length(s.reg$time)
	return(Ind)
	}

#-------------------------------------------------------------------
# Index.CH(): R2_CH
#-------------------------------------------------------------------

Index.CH <- function (gene, tt, cens) {

 	n <- length(tt)
 	o <- order(tt)

 	s.reg <- summary(survfit(Surv(tt , cens )~1))
 	ni <- s.reg$n.risk

 	tps.deces <- s.reg$time
 	delti <- cens

 	mat <- as.data.frame(matrix(nrow=n , ncol=7))

	 names(mat) <- c( "tps", "n.risk", "delta", "x.event", "x.deces", "x.risk", "poids")
 	mat$tps <- round(tt[o] ,2)
 	mat$n.risk <- n:1
 	mat$delta <- delti[o]
 	geneo <- gene[o]
 	etato <- cens[o]
 	mat$x.deces <- as.numeric(geneo*etato)
 	mat$x.event <- as.numeric(geneo)
 	mat$x.risk <- rev(cumsum(rev(as.numeric(geneo))))

	 Kt0 <- cumsum(mat$delta/mat$n.risk)
 	Kt2 <- c(0, Kt0[-length(Kt0)])
 	omega <- ifelse(Kt2==0, 1, 1+log(Kt2))
 	mat$poids <- omega*mat$delta

 	Y <- matrix(1 , ncol=n, nrow=n)
 	Y[upper.tri(Y)] <- 0

 	Ui.tmp <- mat$delta*(mat$x.deces -(mat$x.risk /mat$n.risk))
 	Ui <- Ui.tmp*mat$poids

 	v1 <- mat$delta/mat$n.risk
 	m1 <- matrix(v1, nrow=n, ncol=n, byrow=T)*Y
 	n1 <- m1*mat$x.event
 	v2 <- mat$delta/(mat$n.risk*mat$n.risk)
 	m2 <- matrix(v2, nrow=n, ncol=n, byrow=T)*Y
 	n2 <- matrix(mat$x.risk, ncol=length(mat$tps), nrow=length(mat$tps ), byrow=T)*m2

 	EUij.tmp <- n1-n2
 	EUij <- EUij.tmp*matrix(mat$poids, ncol=length(mat$tps), nrow=length(mat$tps), byrow=T)
 	EUi <- apply(EUij , 1 ,sum)
 	Wi <- Ui - EUi
 	Cox <- sum(Wi)^2/sum(Wi^2)
 	Ind <- Cox/length(s.reg$time)
 	return (Ind)
 }

#-------------------------------------------------------------------
# Index.ModCH(): R2_ModCH
#-------------------------------------------------------------------

Index.modCH <- function (gene , tt, cens) {

 	n <- length(tt)
 	o <- order(tt)

 	s.reg <- summary(survfit(Surv(tt , cens )~1))
 	ni <- s.reg$n.risk

 	tps.deces <- s.reg$time
 	delti <- cens

 	mat <- as.data.frame(matrix(nrow=n , ncol=7))

	 names(mat) <- c( "tps", "n.risk", "delta", "x.event", "x.deces", "x.risk", "poids")
 	mat$tps <- round(tt[o] ,2)
 	mat$n.risk <- n:1
 	mat$delta <- delti[o]
 	geneo <- gene[o]
 	etato <- cens[o]
 	mat$x.deces <- as.numeric(geneo*etato)
 	mat$x.event <- as.numeric(geneo)
 	mat$x.risk <- rev(cumsum(rev(as.numeric(geneo))))

	 Kt0 <- cumsum(mat$delta/mat$n.risk)
 	Kt2 <- c(0, Kt0[-length(Kt0)])
 	omega <- ifelse(Kt2==0, -8.2, 1+log(Kt2))

 	mat$poids <- omega*mat$delta

 	Y <- matrix(1 , ncol=n, nrow=n)
 	Y[upper.tri(Y)] <- 0

 	Ui.tmp <- mat$delta*(mat$x.deces -(mat$x.risk /mat$n.risk))
 	Ui <- Ui.tmp*mat$poids

 	v1 <- mat$delta/mat$n.risk
 	m1 <- matrix(v1, nrow=n, ncol=n, byrow=T)*Y
 	n1 <- m1*mat$x.event
 	v2 <- mat$delta/(mat$n.risk*mat$n.risk)
 	m2 <- matrix(v2, nrow=n, ncol=n, byrow=T)*Y
 	n2 <- matrix(mat$x.risk, ncol=length(mat$tps), nrow=length(mat$tps ), byrow=T)*m2

 	EUij.tmp <- n1-n2
 	EUij <- EUij.tmp*matrix(mat$poids, ncol=length(mat$tps), nrow=length(mat$tps), byrow=T)
 	EUi <- apply(EUij , 1 ,sum)
 	Wi <- Ui - EUi
 	Cox <- sum(Wi)^2/sum(Wi^2)
 	Ind <- Cox/length(s.reg$time)
 	return (Ind)
 }

#-------------------------------------------------------------------
# Index.PH(): R2_PH
#-------------------------------------------------------------------

Index.PH <- function(gene , tt, cens){

	cens.fake <- rep(1 ,length(tt))
 	n <- length(tt)
 	o <- order(tt)

 	s.reg <- summary(survfit(Surv(tt, cens)~1))
	s.fake <- summary(survfit(Surv(tt ,cens.fake)~1))
	di <- s.reg$n.event
	ni <- s.reg$n.risk[di > 0]

	ei <- s.fake$n.event
	tps.tot <- s.fake$time
	tps.deces <- s.reg$time
	delti <- as.numeric(is.element(tps.tot, tps.deces))

	mat<- as.data.frame(matrix(nrow=length(tps.tot), ncol=10))
	names(mat) <- c("tps", "n.risk", "n.event", "n.deces", "delta", "x.deces", "x.event", "x.risk")
	mat$tps <- tps.tot
	mat$n.risk <- s.fake$n.risk
	mat$n.event <- ei
	mat$delta <- delti
	deb <- c(1, cumsum(ei[-length(ei)]) + 1)
	finn <- deb + ei - 1
	geneo <- gene[o]
	etato <- cens[o]
	genetat <- geneo*etato
	for(i in 1:length(deb)){
		mat$n.deces[i] <- sum(etato[deb[i]:finn[i]])
		mat$x.deces[i] <- sum(genetat[deb[i]:finn[i]])
		mat$x.event[i] <- sum(geneo[deb[i]:finn[i]])
		}

	mat$x.risk <- rev(cumsum(rev(mat$x.event)))
	
	Y <- matrix (1, ncol=length(mat$tps), nrow=length(mat$tps))
	Y[upper.tri(Y)] <- 0
	
	Ui <- mat$delta*(mat$x.deces - (mat$n.deces/mat$n.risk)*mat$x.risk)
	
	v1.1 <- mat$delta*mat$n.deces/mat$n.risk
	m1.1 <- matrix(v1.1, nrow=length(mat$tps), ncol=length(mat$tps), byrow=T)*Y
	n1.1 <- m1.1*mat$x.event
	v2.1 <- mat$delta*mat$n.deces/(mat$n.risk*mat$n.risk)
	m2.1 <- matrix(v2.1, nrow=length(mat$tps), ncol=length(mat$tps), byrow=T)*Y

	tmp2.1 <- m2.1*matrix(mat$n.event, ncol=length(mat$tps), nrow=length(mat$tps))
	n2.1 <- matrix(mat$x.risk, ncol=length(mat$tps), nrow=length(mat$tps), byrow=T)*tmp2.1
	EUij <- n1.1 - n2.1
	EUi <- apply(EUij, 1, sum)
	Wi <- Ui - EUi
	
	Cox.PH <- sum(Wi)^2/sum(Wi^2)
	Ind.PH <- Cox.PH/length(s.reg$time)
	return(Ind.PH)
}

#-------------------------------------------------------------------
# Example
#-------------------------------------------------------------------

gg <- rnorm(200,0,1)	# Generate gene expression
stime <- rexp(200,1)	# Generate survival time
cens <- sample(c(0,1), replace = T, 200)	# Generate censoring indicator

#R2_PO
Index.PO(gg, stime, cens)

#R2_CO
Index.CO(gg, stime, cens)

#R2_CH
Index.CH(gg, stime, cens)

#R2_ModCH
Index.modCH(gg, stime, cens)

#R2_ModCH
Index.PH(gg, stime, cens)

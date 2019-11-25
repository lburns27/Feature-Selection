#----------------------------------------------------------------------------------------
# Data generating function: W (Weibull, same shape, different scale)
# Incorporates correlations between genes
#	N: Number of subjects
#	small: # of genes with small effect
#	large: # of genes with large effect
#	tau: censoring parameter (see examples below)
#	NoSig: number of genes with no effect
#	lam1, lam2: scale parameters
#----------------------------------------------------------------------------------------

makeDataWMV <- function(N, small, large, tau, NoSig, lam1, lam2){
	Bo_small <- .289
	Bo_large <- .577
	Bt_small <- function(time) {rep(Bo_small*log(lam1/lam2), length(time))} 
	Bt_large <- function(time) {rep(Bo_large*log(lam1/lam2), length(time))} 
	Bt_noSig <- function(time) {rep(0, length(time))}

	p <- small + large + NoSig

	v <- runif(N, min=0, max=1)
	y <- (-log(v))^(1/2)/0.5

	if(tau != "inf"){
		z <- runif(N, 0.00001, tau)
		t <- pmin(y, z)
		d <- (y <= z) + 0
		}
	if(tau == "inf"){ 
		t <- y
		d <- rep(1,N) 
		}

	survDat_unsorted <- cbind(t, d)
	survDat <- as.data.frame(survDat_unsorted[order(t),])
	dat <- survDat

	geneMatrix <- matrix(0, N, p)

	for(j in 1:N){
		u1 <- runif(1, 0, 1)
		u2 <- runif(1, 0, 1)
		u3 <- runif(1, 0, 1)
		for(g in 1:p){
			e <- rnorm(1,0,1)
			if(j <= 0.5*N & g <= 0.05*p){geneMatrix[j,g] <- -1 + e}
			if(j > 0.5*N & g <= 0.05*p){geneMatrix[j,g] <- 1 + e}
			if(0.05 < g & g <= 0.1*p){geneMatrix[j,g] <- 1.5*((u1 < 0.4) + 0) + e}
			if(0.1 < g & g <= 0.2*p){geneMatrix[j,g] <- 0.5*((u2 < 0.7) + 0) + e}
			if(0.2 < g & g <= 0.3*p){geneMatrix[j,g] <- 1.5*((u3 < 0.3) + 0) + e}
			else{geneMatrix[j,g] <- e}
		}
	}

	x <- matrix(0, N, p)
	R <- 1:N
	a <- geneMatrix
	for(i in 1:N){
		if(length(R) == 1){ x[i,] <- a[R,] }
		else if (survDat$d[i] == 0){
			j <- sample(R, 1)
			x[i,] <- a[j,]
			R <- R[R != j]
			}
		else{
			smallGenes <- apply(Bt_small(survDat$t[i])*a[R,1:small], 1, sum)
			largeGenes <- apply(Bt_large(survDat$t[i])*a[R,(small+1):(small + large)], 1, sum)
			noSigGenes <- apply(Bt_noSig(survDat$t[i])*a[R,(small+large+1):p], 1, sum)

			topSum <- smallGenes + largeGenes + noSigGenes
			#probs <- exp(topSum)/sum(exp(topSum))

			#If there are infinity's in exp(topSum) use:
			expSum <- exp(topSum)
			expSum[which(expSum == Inf)] <- 1E99
			expSum[which(expSum == -Inf)] <- -1E99
			probs <- expSum/sum(expSum)

			j <- sample(R, 1, prob = probs)
			x[i,] <- a[j,]
			R <- R[R != j]
			}
		}

	d <- cbind(dat, x)
	colnames(d) <- c("time", "censor", paste("S", 1:small, sep = ""), paste("L", 1:large, sep = ""), paste("NoSig", 1:NoSig, sep = ""))
	return(d)

}

# Example----------------------------------------------------------------------------------------------------------

numSmall <- 200
numLarge <- 200
notsig <- 4600
N <- 200

# 0% censoring
dat0 <- makeDataWMV(N, numSmall, numLarge, "inf", notsig, 1,.5)

# 33% censoring
dat33 <- makeDataWMV(N, numSmall, numLarge, 5.3, notsig, 1,.5)

# 67% censoring
dat67 <- makeDataWMV(N, numSmall, numLarge, 2.4, notsig, 1,.5)

# 80% censoring
dat80 <- makeDataWMV(N, numSmall, numLarge, 1.7, notsig, 1,.5)





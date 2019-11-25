#----------------------------------------------------------------------------------------
# Data generating function: LN (log-normal)
#	N: Number of subjects
#	small: # of genes with small effect
#	large: # of genes with large effect
#	tau: censoring parameter (see examples below)
#	NoSig: number of genes with no effect
#----------------------------------------------------------------------------------------


makeDataLN <- function(N, small, large, tau, NoSig){
	Bo_small <- .346
	Bo_large <- 2.098

	y <- rlnorm(N)

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

	Bo <- Bo_small
	Bt <- function(time) {Bo*(time^2-1)} 
	for(i in 1:small){
		a <- rnorm(N, 0, 1)
		x <- c()
		R <- 1:N

		for(i in 1:N){
			if(length(R) == 1){ x[i] <- a[R] }
			else if (survDat$d[i] == 0){
				j <- sample(R, 1)
				x[i] <- a[j]
				R <- R[R != j]
				}
			else{
				probs <- exp(Bt(survDat$t[i])*a[R])/sum(exp(Bt(survDat$t[i])*a[R]))
				probs[is.nan(probs)] <- 1E-200
				j <- sample(R, 1, prob = probs)
				x[i] <- a[j]
				R <- R[R != j]
				}
			}
		dat <- cbind(dat, x)
		}

	Bo <- Bo_large
	Bt <- function(time) {Bo*(time^2-1)} 
	for(i in 1:large){
		a <- rnorm(N, 0, 1)
		x <- c()
		R <- 1:N

		for(i in 1:N){
			if(length(R) == 1){ x[i] <- a[R] }
			else if (survDat$d[i] == 0){
				j <- sample(R, 1)
				x[i] <- a[j]
				R <- R[R != j]
				}
			else{
				probs <- exp(Bt(survDat$t[i])*a[R])/sum(exp(Bt(survDat$t[i])*a[R]))
				probs[is.nan(probs)] <- 1E-200
				j <- sample(R, 1, prob = probs)
				x[i] <- a[j]
				R <- R[R != j]
				}
			}
		dat <- cbind(dat, x)
		}

	Bo <- 0
	Bt <- function(time) {Bo*(time^2-1)} 
	for(i in 1:NoSig){
		a <- rnorm(N, 0, 1)
		x <- c()
		R <- 1:N

		for(i in 1:N){
			if(length(R) == 1){ x[i] <- a[R] }
			else if (survDat$d[i] == 0){
				j <- sample(R, 1)
				x[i] <- a[j]
				R <- R[R != j]
				}
			else{
				probs <- exp(Bt(survDat$t[i])*a[R])/sum(exp(Bt(survDat$t[i])*a[R]))
				j <- sample(R, 1, prob = probs)
				x[i] <- a[j]
				R <- R[R != j]
				}
			}
		dat <- cbind(dat, x)
		}


	colnames(dat) <- c("time", "censor", paste("CrossS", 1:small, sep = ""), 
				paste("CrossL", 1:large, sep = ""), paste("NoSig", 1:NoSig, sep = ""))
	return(dat)	
}

# Example----------------------------------------------------------------------------------------------------------

numSmall <- 100
numLarge <- 100
notsig <- 500
N <- 100

# 0% censoring
dat0 <- makeDataLN(N, numSmall, numLarge, "inf", notsig)

# 33% censoring
dat33 <- makeDataLN(N, numSmall, numLarge, 3.59, notsig)

# 67% censoring
dat67 <- makeDataLN(N, numSmall, numLarge, 1.36, notsig)

# 80% censoring
dat80 <- makeDataLN(N, numSmall, numLarge, 0.85, notsig)
#----------------------------------------------------------------------------
# getResults: Function to calculate Sens, Spec, Youden & AUC for measure
#	of interest
#	
#	out: Data frame with n rows (number of genes) and 1 column (measure)
#	effectGenes: number of significant genes 
#
#	Returns: Sensitivity, Specificity, Youden, & AUC 
#
#	Required package: MESS
#-----------------------------------------------------------------------------

library(MESS)

getResults <- function(out, effectGenes){

	# Create data frame with measure as 1st column
	# Second column = arbitrary place holder
	out <- as.data.frame(cbind(out,1))
	colnames(out) <- c("Measure", "NotNeeded")

	n <- effectGenes
	m <- nrow(out)

	goodgenes <- rownames(out)[1:n]
	badgenes <- rownames(out)[(n+1):nrow(out)]

	sortMeasure<- out[order(abs(out$Measure), decreasing = T),]

	#Sens & Spec

	sm <- rownames(sortMeasure)[1:n]

	goodpick <- length(sm[! sm %in% badgenes])
	badpick <- length(sm[sm %in% badgenes])
	sens <- goodpick/n
	spec <- 1-badpick/(m-n)
	youden <- sens + spec - 1	


	#AUC

	choosenum <- seq(0, m, 10)
	tpr <- rep(0, length(choosenum)) 
	fpr <- rep(0, length(choosenum)) 

	count <- 0
	for(i in choosenum){
		count <- count+1
		sm <- rownames(sortMeasure)[1:i]
		goodpick <- length(sm[! sm %in% badgenes])
		badpick <- length(sm[sm %in% badgenes])
		tpr[count] <- goodpick/n
		fpr[count] <- badpick/(m-n)
		}
	
	area <- auc(fpr, tpr, 0, 1) 
		
	changeArea <- function(area, lower, upper){
		area <- auc(fpr, tpr, lower, upper)
		return(area)
	}
	
	if(sum(is.na(area)) > 0){area <- changeArea(area, .001, 1)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .9999)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .999)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .99)}
	if(sum(is.na(area)) > 0){area <- changeArea(area, .01, .9)}
	
	spitOut <- c(Sens = sens, Spec = spec, Youden = youden, AUC = area)

	return(spitOut)
}


#----------------------------------------------------------------------------
# This code creates Venn Diagrams showing various interestions between
#	different feature selection measures
#	
#	Venn Diagrams created:
#		I_PO, IYP, & concreg (*dichotomized case is coded separately*)
#		R2_PO, R2_IPO & R2_LR 
#		R2_PO, R2_ModCH & R2_CO
#
#	Note:  
#		Obtain a data frame named "out" with n rows (number of genes) and 
#		1 column for each of the measures listed above (except I_YP).  All
# 		measures here are based on continuous gene expression. 
#
#		Obtain a data frame named "out2" with n rows (number of genes) and 
#		1 column for each of the measures based on dichotomized expression
#		(I_YP, I_PO, concreg). 
#	
#	Required packages: gpplots, VennDiagram, latex2exp
#-----------------------------------------------------------------------------

#Required packages
library(gplots)
library(VennDiagram)
library(latex2exp)


# Obtain data frame 'out' and 'out2' as described above prior to running this code.

# Rank genes by each measure ------------------------------------------------

n <- 500

# Continuous expression
r2po <- rownames(out[order(out$R2_PO, decreasing = T),])[1:n]
r2lr <- rownames(out[order(out$R2_LR, decreasing = T),])[1:n]
r2ipo <- rownames(out[order(out$R2_IPO, decreasing = T),])[1:n]
r2co <- rownames(out[order(out$R2_CO, decreasing = T),])[1:n]
r2modch <- rownames(out[order(out$R2_ModCH, decreasing = T),])[1:n]
i1po <- rownames(out[order(out$I_PO, decreasing = T),])[1:n]
concreg <- rownames(out[order(out$ConcEffect, decreasing = T),])[1:n]
out$Conc2 <- .5 + abs(out$ConcEffect -.5)	#Ensure that the absolute effect size is being used
concreg2 <- rownames(out[order(out$Conc2, decreasing = T),])[1:n]

# Dichotomized expression
i2yp <- rownames(out[order(out2$I_YP, decreasing = T),])[1:n]
i1po_dich <- rownames(out[order(out2$I_PO_dich, decreasing = T),])[1:n]
out2$Conc2 <- .5 + abs(out2$ConcEffect_dich -.5)
concreg2_dich <- rownames(out2[order(out2$Conc2, decreasing = T),])[1:n]

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

# Create Venn Diagrams (dichotomized expression) -------------------------------------

listOther <- list(i1po_dich, i2yp, concreg2_dich)
venn.plot <- venn.diagram(listOther , NULL, category.names=c(TeX('$\\textit{I_{PO}}$'), TeX('$\\textit{I_{YP}}$'), TeX('$\\hat{c}^\'_{+}$')), 
			rotation.degree =60, cat.cex = 2.2, cex = 2.2, lwd = 3, cat.pos = c(-130,0,140), cat.dist = c(.015, .04, .045), scaled = F, euler.d = F)
grid.draw(venn.plot)


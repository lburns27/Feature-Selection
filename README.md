# Feature-Selection
This repository contains R code for various feature selection methods. The code is organized as outlined below.

* **Youden & AUC**: Computes Youden & AUC values based on gene ranking by a specified feature selection method
  * Required: R packages (MESS)
  * Inputs: 
    * data frame with n rows (number of genes) and 1 column (feature selection measure)
    * effectGenes: number of significant genes
  * Output option: Specificity, Sensitivity, Youden, AUC & AUC standard deviation

* **Simulations**: Contains R code for simulating data
  * Scheme 1: Univariate aproach; genes linked to survival one at a time
  * Scheme 2: Multivariate approach; incorporates correlation between features
  * For both schemes, there are options to simulate from the following models: LN, LL1, LL2, W1, W2

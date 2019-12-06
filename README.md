# Feature-Selection
This repository contains R code for various feature selection methods. The code is organized as outlined below.

* **Model Fits & GOF**: Includes functions for fitting the Cox, PO & YP models, as well as testing for goodness-of-fit.  For the Cox and PO cases, there are options for adjusting for age and stage, if desired. 
  * Required: R packages (survival, timereg, YPmodel)
  * Inputs: 
    * predictor, x (gene expression)
    * data frame with n rows (number of subjects) and columns with survival information (time = survival time, censor = censoring information) and p genes   
  * Output: model coefficient (beta), coefficient standard error, significance p-value, GOF p-value
  
* **Pseudo-R<sup>2</sup> Measures**: Includes functions for calculating each pseudo-R<sup>2</sup> measures (PO, CO, CH, ModCH and PH).  Note, each measure has a separate function.  In the future, we plan to combine the code into one R<sup>2</sup> function with option `type = c("PO", "CO", "CH", "ModCH", "PH")`. 
  * Required: R package (survival)
  * Inputs: 
    * predictor, x (gene expression)
    * survival time
    * censoring indicator
  * Output: R<sup>2</sup> measure
  
* **R<sup>2</sup><sub>LR</sub> & R<sup>2</sup><sub>I</sub> Measures**: This code contains functions for R<sup>2</sup><sub>LR</sub>, R<sup>2</sup><sub>I<sub>PO</sub></sub> & R<sup>2</sup><sub>I<sub>PH</sub></sub>. There are options for adjusting for age and stage, if desired.
  * Required: R packages (survival, timereg)
  * Inputs: 
    * predictor, x (gene expression)
    * survival time
    * censoring indicator
  * Output: R<sup>2</sup> measures
  
* **I Measures**: This code contains functions for I<sub>PO</sub> and I<sub>YP</sub>. There are options for adjusting for age and stage for I<sub>PO</sub>, if desired.
  * Required: R packages (survival, timereg, YPmodel)
  * This code assumes that your data is in the the following form:
    * Column 1 = survival time (time)
    * Column 2 = censoring indicator (censor)
    * Columns 3+ = genes
  * Output: 
    * I<sub>PO</sub>, outPO (I, I test statistic, I p-value)
    * I<sub>YP</sub>, outYP (I, I test statistic, I p-value)
  
* **Youden & AUC**: Computes Youden & AUC values based on gene ranking by a specified feature selection method
  * Required: R packages (MESS)
  * Inputs: 
    * data frame with n rows (number of genes) and 1 column (feature selection measure)
    * effectGenes: number of significant genes
  * Output option: Specificity, Sensitivity, Youden & AUC 

* **Venn Diagrams**: This code creates Venn Diagrams showing various interestions between different feature selection measures.
  * Required: R packages (gpplots, VennDiagram, latex2exp)
  * Inputs: 
    * Obtain a data frame named "out" with n rows (number of genes) and	1 column for each of the measures listed above (except I<sub>YP</sub>). All measures here are based on continuous gene expression. 
    * Obtain a data frame named "out2" with n rows (number of genes) and 1 column for each of the measures based on dichotomized expression (I<sub>YP</sub>, I<sub>PO</sub>, concreg). 
  * Venn Diagrams created: 
    * I<sub>PO</sub>, I<sub>YP</sub>, & concreg (*dichotomized case is coded separately*)
    * R<sup>2</sup><sub>PO</sub>, R<sup>2</sup><sub>I<sub>PO</sub></sub> & R<sup>2</sup><sub>LR</sub> 
    * R<sup>2</sup><sub>PO</sub>, R<sup>2</sup><sub>ModCH</sub> & R<sup>2</sup><sub>CO</sub>
    
* **Other Existing Measures**: This code contains functions for computing some existing measures in the literature.
  * Required: R packages (concreg, survAUC, survival, pec, timereg)
  * Measures computed:
    * Concreg (Dunkler et al. 2010)
    * Uno's C (Uno et al. 2011)
    * R<sup>2</sup><sub>G</sub> - PH & PO cases (Graf et al. 1999; Gerds & Schumacher 2006)
    * R<sup>2</sup><sub>SH</sub> (Schemper & Henderson 2000)
  * Inputs: 
    * predictor, x (gene expression)
    * survival time
    * censoring indicator
  * Output: Concreg (absolute effect size), Uno's C, R<sup>2</sup><sub>G</sub> (PH), R<sup>2</sup><sub>G</sub> (PO) and R<sup>2</sup><sub>SH</sub>.
    
* **Simulations**: Contains R code for simulating data
  * Scheme 1: Univariate aproach; genes linked to survival one at a time
  * Scheme 2: Multivariate approach; incorporates correlation between features
  * For both schemes, there are options to simulate from the following models: LN, LL1, LL2, W1, W2

# Copyright & Citations
Copyright © Lauren Spirko-Burns and Karthik Devarajan 

Spirko, L.N., Devarajan, K. Unified methods for variable selection in large-scale genomic studies with censored survival outcomes. Under review. COBRA pre-print series, Article 120 (June 2019). http://biostats.bepress.com/cobra/art120. 

Spirko, L. (2017). Variable Selection and Supervised Dimension Reduction for Large-Scale Genomic Data with Censored Survival Outcomes. Ph.D. Dissertation. Department of Statistical Science, Temple University, Philadelphia.

# License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Lauren Spirko Burns and Karthik Devarajan</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

# 

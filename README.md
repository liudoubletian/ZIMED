# ZIMED (The effect decomposition of mediation analysis for zero-inflated models in microbiome data)

# Introduction
The characteristics of microbiome data are complicated: sparsity and over-dispersion. However, there are few existing causal mediation methods specifically designed to handle sparse and high dimensional microbiome data. We develop a novel zero-inflated mediation effect decomposition model (ZIMED) to estimate and test the mediation effects of the microbiome utilizing the zero-inflated negative-binomial (ZINB) regression model.

# Installation
```r
install.packages("devtools")  
devtools::install_github("liudoubletian/ZIMED")  
library(ZIMED)  
```
# Basic Usage
## mediation analysis for microbiome data
```r
zimed(M_mat,Treat,Outcome,method="joint",ci.method="delta)
```
* `M_mat` : an OTU table with n rows (samples) and m columns (taxa)
* `Treat` : a n-vector of group indicators
* `method` : joint, HDMT or DACT
* `ci.method` : bootstrap or delta method

it returns a list of results:  
* `NIE` : the estimated natural indirect effect
* `NIE.p` : the calculated p value for NIE
* `NIE.ci` : the calculated confidence interval for NIE
* `NIEA` : the estimated natural indirect effect through changes of abundance
* `NIEA.p` : the calculated p value for NIEA
* `NIEA.ci` : the calculated confidence interval for NIEA
* `NIEP` : the estimated natural indirect effect through changes of absence/presence
* `NIEP.p` : the calculated p value for NIEP
* `NIEP.ci` : the calculated confidence interval for NIEP





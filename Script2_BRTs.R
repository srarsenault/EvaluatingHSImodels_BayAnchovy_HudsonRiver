#=============================================================================
#
# Evaluating effects of data quality and variable weighting on habitat suitability modelling 
# Script 2) Boosted Regression Trees (BRTs)
#
# Authors: Stephanie Arsenault, Robyn Linner, and Yong Chen
# Date: 02.12.2025
#
## LOAD REQUIRED PACKAGES =====================================================
#
l_packages <- c("dismo", "gbm3", "dplyr")
for (p in l_packages){
  if(! p %in% installed.packages()){
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = T)
}
#

## BOOSTED REGRESSION TREES ===================================================
#
# Load data
df=read.csv() #add path

# Create matrix where the results from the BRT runs will be save, ncol = 10 (or however many simulations), nrow = number of habitat variables
BRT_res=matrix(data=0, ncol=10, nrow=4) 

# Names for the environmental variables in df
rownames(BRT_res)=c("Depth", "Temp", "DO","Salinity")

# Boosted Regression Trees
## test different combinations of learning rates and tree complexities
## select the combination that minimizes the cross-validated deviance
for (i in 1:10){          
  LRS_brt = gbm.step(data = LRS_1982, gbm.x = c(1:4), gbm.y = 5, family = "gaussian", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.75)
  q = summary(LRS_brt)
  reorder = c("RIVER.DEPTH", "WATER.TEMPERATURE","DISSOLVED.OXYGEN","SALINITY")
  q = q[match(reorder, q$var),]
  BRT_res[,i] = q$rel.inf
}
#
#Take the mean of the relative influences from the BRTs
LRS_inf=as.data.frame(rowMeans(BRT_res)/100)
summary(LRS_inf)


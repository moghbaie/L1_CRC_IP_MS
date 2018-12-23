## Mehrnoosh Oghbaie
## 12/11/2018
## Imputate missing values

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","plyr","magrittr","igraph","sqldf","stringr","corrplot","ggplot2","ggridges","gridExtra")
bioconductor.packages <- c("biomaRt")
source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

## Imputation consists of three stages:
##  1. Remove proteins not identified in both cases and controls
##  2. Impute values for proteins with zero intensity (cases and controls)
##  3. Impute values for proteins that have at least one non-zero intensity

load(file="../Result/2. After_preparation/log_transformed.RData")
set.seed(123)
#####################################################################################################
##  1. Remove proteins not identified in both cases and controls
# Ovary

list_experiment$Ovary <- list_experiment$Ovary[rowSums(list_experiment$Ovary[,-(1:2)], na.rm =T)!=0,]
dim(list_experiment$Ovary)

# Liver
list_experiment$Liver <- list_experiment$Liver[rowSums(list_experiment$Liver[,-(1:2)], na.rm =T)!=0,]
dim(list_experiment$Liver)

# Colon
list_experiment$Colon <- list_experiment$Colon[rowSums(list_experiment$Colon[,-(1:2)], na.rm =T)!=0,]
dim(list_experiment$Colon)

list_experiment[["PreImputation"]] <- list_experiment
#####################################################################################################
### Imputation method 6_2

#Method_2:
#  For each protein, which has zeros in each replica, replica with the least number of zeros is chosen.
#  For it mean and standard deviation of non-zero proteins intensities are counted. 
#  New intensity for this protein in this replica is sampled from uniform distribution with parameters: 
#  start = mu - 3*sd, end = mu - 2*sd. Then other replicas for this protein are imputed as in 6th step
#Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_6:
# Impute values for proteins, which have zero intensities only in some replicates:
#  Build distribution of deltas for all non zero proteins, where 
# delta =(Intrep1-Intrep2) mean(Intrep1,Intrep2)
# Calculate mudelta , sddelta
# Calculate new delta and new Intensity:
# deltanew=rnorm(mu=mudelta, sd=sddelta*sqrt(2)*mean(correlations))
# Inew=mean(Intother)*abs(1+deltanew)

run_order <- rbind(cbind(rep("Ovary",2),c("Igg","Tumor")),
                   cbind(rep("Liver",3),c("Igg","Normal","Tumor")),
                   cbind(rep("Colon",3),c("Igg","Normal","Tumor")))

##  2. Impute values for proteins with zero intensity (cases and controls)
#for(i in 1:dim(run_order)[1]){
#  list_experiment[[run_order[i,1]]] <- mutate_all_zeros_2(list_experiment[[run_order[i,1]]],run_order[i,2])
#}

##  3. Impute values for proteins that have at least one non-zero intensity

#for(i in 1:dim(run_order)[1]){
#  list_experiment[[run_order[i,1]]] <- imutate_partial_zero_6(list_experiment[[run_order[i,1]]],run_order[i,2])
#}

#####################################################################################################
### Imputation method 7_1 (* We implement this method)

#Method_1:
# For each replica mean and standard deviation of non-zero proteins intensities are counted. 
# New intensity for each missing value in each replica is sampled from uniform distribution with parameters:
#  start = mu - 3*sd, end = mu - 2*sd
# Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_7:
# For outliers (proteins, which have zero values in all but one replica) this method implies one of the methods for imputation of all zero replicas.
# For other proteins uses method6

for(i in 1:dim(run_order)[1]){
  list_experiment[[run_order[i,1]]] <- imutate_partial_zero_7_1(list_experiment[[run_order[i,1]]],run_order[i,2])
}

###############################################################################################
## Compare distribution of data before and after imputation.

## Colon tissue
png(filename="../Image/Imputation/comparison_before_after_imputation_Colon_7_1.png",width = 1100, height = 600)
my.env <- plot_density(data=list_experiment, tissue ="Colon")
p <- grid.arrange(grobs=list(my.env$p1,my.env$p2,my.env$p3,
                        my.env$p4,my.env$p5,my.env$p6,
                        my.env$p7,my.env$p8,my.env$p9),top = "Colon tissue")
print(p)
dev.off()

## Liver tissue 
png(filename="../Image/Imputation/comparison_before_after_imputation_Liver_7_1.png",width = 1100, height = 600)
my.env <-  plot_density(data=list_experiment, tissue ="Liver")
p2<-arrangeGrob(my.env$p1,my.env$p2,my.env$p3,
                 my.env$p4,my.env$p5,my.env$p6,
                 my.env$p7,my.env$p8,my.env$p9)
grid.arrange(p2,top = "Liver tissue")
print(p2)
dev.off()

## Ovary tissue
png(filename="../Image/Imputation/comparison_before_after_imputation_Ovary_7_1.png",width = 733, height = 600)
my.env <- plot_density(data=list_experiment, tissue ="Ovary")
p3<-arrangeGrob(my.env$p1,my.env$p2,my.env$p3,
                my.env$p4,my.env$p5,my.env$p6)
grid.arrange(p3,top = "Ovary tissue")
print(p3)
dev.off()
##########################################################################
## Plot average value

## Average intensities before imputation
png(filename="../Image/Imputation/average_intensities_before_imputation.png",width = 1000, height = 600)
r1 <- draw_mean_impute(data=list_experiment$PreImputation ,tissue="Colon")
r2 <- draw_mean_impute(data=list_experiment$PreImputation ,tissue="Liver")
r3 <- draw_mean_impute(data=list_experiment$PreImputation , tissue="Ovary")
r <- grid.arrange(grobs=list(r1,r2,r3),ncol=3, top="Before imputation")
print(r)
dev.off()

## Average intensities after imputation
png(filename="../Image/Imputation/average_intensities_after_imputation.png",width = 1000, height = 600)
q1 <- draw_mean_impute(data=list_experiment,tissue="Colon")
q2 <- draw_mean_impute(data=list_experiment,tissue="Liver")
q3 <- draw_mean_impute(data=list_experiment,tissue="Ovary")
q <- grid.arrange(grobs=list(q1,q2,q3),ncol=3,top="After imputation")
print(q)
dev.off()


####################################################################################################################
## Save the output in "../Result/3. After_imputation/"
####################################################################################################################


save(KM050217_merged,list_experiment,file="../Result/3. After_imputation/after_imputation.RData")

## Mehrnoosh Oghbaie
## 12/11/2018
## Imputate missing values

## Imputation consists of three stages:
##  1. Remove proteins not identified in both cases and controls
##  2. Impute values for proteins with zero intensity (cases and controls)
##  3. Impute values for proteins that have at least one non-zero intensity

#load(file="../Result/2. After_preparation/log_transformed.RData")
set.seed(123)
#####################################################################################################
##  1. Remove proteins not identified in both cases and controls

Template$set("public","PreImputation", list())

Template$set("public","removeAllZeros", function(run_order){
  self$PreImputation <- self$experiment
  for ( i in names(self$PreImputation)){
    self$PreImputation[["Ovary"]] <- self$PreImputation[["Ovary"]][rowSums(self$PreImputation[["Ovary"]][,-(1:2)], na.rm =T)!=0,]
  }
}
)

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

Template$set("public","experimentImputed", list())

Template$set("public","imputeAll", function(run_order){
  self$experimentImputed <- self$PreImputation
  for(i in 1:dim(run_order)[1]){
    self$experimentImputed[[run_order[i,1]]] <- imutate_partial_zero_7_1(self$experimentImputed[[run_order[i,1]]],run_order[i,2])
  }
}
)



###############################################################################################
## Compare distribution of data before and after imputation.

Template$set("public","drawComparisonAfterImputation", function(){
  ## Colon tissue
  png(filename="../Image/Imputation/comparison_before_after_imputation_Colon_7_1.png",width = 1100, height = 600)
  my.env <- plot_density(data=self, tissue ="Colon")
  p <- grid.arrange(grobs=list(my.env$p1,my.env$p2,my.env$p3,
                          my.env$p4,my.env$p5,my.env$p6,
                          my.env$p7,my.env$p8,my.env$p9),top = "Colon tissue")
  print(p)
  dev.off()
  
  ## Liver tissue 
  png(filename="../Image/Imputation/comparison_before_after_imputation_Liver_7_1.png",width = 1100, height = 600)
  my.env <-  plot_density(data=self,  tissue ="Liver")
  p2<-arrangeGrob(my.env$p1,my.env$p2,my.env$p3,
                   my.env$p4,my.env$p5,my.env$p6,
                   my.env$p7,my.env$p8,my.env$p9)
  grid.arrange(p2,top = "Liver tissue")
  print(p2)
  dev.off()
  
  ## Ovary tissue
  png(filename="../Image/Imputation/comparison_before_after_imputation_Ovary_7_1.png",width = 733, height = 600)
  my.env <- plot_density(data=self,  tissue ="Ovary")
  p3<-arrangeGrob(my.env$p1,my.env$p2,my.env$p3,
                  my.env$p4,my.env$p5,my.env$p6)
  grid.arrange(p3,top = "Ovary tissue")
  print(p3)
  dev.off()
})
##########################################################################
## Plot average value
Template$set("public","drawAverageImputated", function(){
  ## Average intensities before imputation
  png(filename="../Image/Imputation/average_intensities_before_imputation.png",width = 1000, height = 600)
  r1 <- draw_mean_impute(data=self$PreImputation ,tissue="Colon")
  r2 <- draw_mean_impute(data=self$PreImputation ,tissue="Liver")
  r3 <- draw_mean_impute(data=self$PreImputation , tissue="Ovary")
  r <- grid.arrange(grobs=list(r1,r2,r3),ncol=3, top="Before imputation")
  print(r)
  dev.off()
  
  ## Average intensities after imputation
  png(filename="../Image/Imputation/average_intensities_after_imputation.png",width = 1000, height = 600)
  q1 <- draw_mean_impute(data=self$experimentImputed,tissue="Colon")
  q2 <- draw_mean_impute(data=self$experimentImputed,tissue="Liver")
  q3 <- draw_mean_impute(data=self$experimentImputed,tissue="Ovary")
  q <- grid.arrange(grobs=list(q1,q2,q3),ncol=3,top="After imputation")
  print(q)
  dev.off()
})


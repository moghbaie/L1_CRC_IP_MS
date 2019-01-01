# Mehrnoosh Oghbaie
# 12/10/2018
# Repository for all the functions
# This file includes most of the general functions that are going to be used in this project

######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  #if (length(bioconductor.packages) > 0) {
  #  source("http://bioconductor.org/biocLite.R")
  #}
  
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p,version = "3.8") 
    }
    library(p,lib.loc="~/R/win-library/3.5",character.only=T)
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
    }	
    #library(p, lib.loc="~/R/win-library/3.5",character.only=T) 
    library(p, lib.loc="~/R/win-library/3.5",character.only=T)
  }
}


###################################################################################
### Imputation functions
##################################################################################

### function : finding column with minimum zeros
count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_mutate <- function(na_zeros , mu,sd){
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}

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



### Method2
mutate_all_zeros_2 <- function(data=list_experiment$Liver,condition="Normal"){
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  na_zeros <- rowSums(data[,colnames],na.rm=T)==0
  min_zero_col <- min_col(apply(data[,colnames],2,count_zeros))
  sample <- unname(unlist(data[,min_zero_col]))
  mu <- mean(sample, na.rm=T)
  sd <- sd(sample, na.rm=T)
  data[na_zeros,min_zero_col] <- na_zeros_mutate(na_zeros , mu,sd)
  return(data)
}

### Funtion that Imputate partial rows

### Method6

imutate_partial_zero_6 <- function(data=list_experiment$Colon,condition="Tumor"){
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  count_na <- apply(data[,colnames],1, count_zeros)
  y <- data[count_na>0,colnames]
  for(i in 1:dim(y)[1]){
    col_NAs <- names(y[i,is.na(unlist(y[i,]))])
    for(j in col_NAs){
      col_NA <- j
      col_select <- names(y[i,colnames(y)!=col_NA])
      sample <- y[complete.cases(y[,col_select]),col_select]
      delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
      mu <- mean(unlist(unname(delta)))
      std <- sd(unlist(unname(delta)))
      cor <- cor(data[count_na==0,grepl(condition,colnames(data))])
      mean_cor <- mean(cor[col_NA,col_select])
      deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
      y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
    }
  }
  data[count_na>0,colnames] <- y
  return(data)
} 

#Method_1:
# For each replica mean and standard deviation of non-zero proteins intensities are counted. 
# New intensity for each missing value in each replica is sampled from uniform distribution with parameters:
#  start = mu - 3*sd, end = mu - 2*sd
# Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_7:
# For outliers (proteins, which have zero values in all but one replica) this method implies one of the methods for imputation of all zero replicas.
# For other proteins uses method6

### function: mutate all NA values for specific condition( control or case)
### Method1
mutate_all_zeros_1 <- function(data,condition){
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  na_zeros <- rowSums(data[,colnames],na.rm=T)==0
  for(i in 1:length(colnames)){
    sample <- data[,colnames[i]]
    sample <- sample[complete.cases(sample),][[colnames[i]]]
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[na_zeros,colnames[i]] <- na_zeros_mutate(na_zeros , mu,sd)
  }
  
  return(data)
}

### Method7_1
imutate_partial_zero_7_1 <- function(data=list_experiment$Colon,condition="Tumor"){
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  count_na_initial <- apply(data[,colnames],1, count_zeros)
  OutVals = boxplot(data[[colnames[1]]])$out
  data[count_na_initial>1&(data[[colnames[1]]] %in% OutVals), colnames][!is.na(data[count_na_initial>1&(data[[colnames[1]]] %in% OutVals), colnames])] <- NA 
  data <- mutate_all_zeros_1(data,condition)
  count_na <- apply(data[,colnames],1, count_zeros)
  y <- data[count_na>0,colnames]
  for(i in 1:dim(y)[1]){
    col_NAs <- names(y[i,is.na(unlist(y[i,]))])
    for(j in col_NAs){
      col_NA <- j
      col_select <- names(y[i,colnames(y)!=col_NA])
      sample <- y[complete.cases(y[,col_select]),col_select]
      delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
      mu <- mean(unlist(unname(delta)))
      std <- sd(unlist(unname(delta)))
      cor <- cor(data[count_na==0,colnames])
      mean_cor <- mean(cor[col_NA,col_select])
      deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
      y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
    }
  }
  data[count_na>0,colnames] <- y
  return(data)
} 




### Function that plot intensity comparison before and after imputation

plot_density<- function(data , tissue="Ovary"){
  g<-dim(data[["PreImputation"]][[tissue]][,-(1:2)])
  my.env <- new.env()
  for(i in 1:g[2]){
    df <- data.frame(rbind(
      cbind(unlist(unname(data[["PreImputation"]][[tissue]][,-(1:2)][,i])),rep("Before_imputation", g[1])),
      cbind(unlist(unname(data[["experimentImputed"]][[tissue]][,-(1:2)][,i])),rep("After_imputation", g[1]))
    ))
    colnames(df) <- c("value","group")
    df$value<- as.numeric(as.character(df$value))
    mu <- plyr::ddply(df, "group", summarise, grp.mean=mean(value,na.rm=TRUE))
    assign(paste0("p",i),ggplot(df, aes(x=value, y=group, fill=factor(..quantile..))) +
             stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.25, 0.75)) +
             scale_fill_manual(
               name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
               labels = c("(0, 0.25]", "(0.25, 0.75]", "(0.75, 1]")
             )+labs(title = paste("Log intensity",names(data[[tissue]][,-(1:2)][,i]))),envir = my.env )
    
  }
  return(my.env)
}


### Function that plot average intensity  after imputation
draw_mean_impute <- function(data, tissue="Colon"){
  df <- data[[tissue]][,-(1:2)]
  df_average<- data.frame(matrix(NA, ncol=2, nrow=0))
  for (k in run_order[run_order[,1]==tissue,2]){
    df_average<- rbind(df_average,
                       cbind(apply(df[,grepl(k,colnames(df))],1,mean),rep(k,dim(df)[1])))
  }
  df_average$V1 <- as.numeric(as.character(df_average$V1))
  colnames(df_average) <- c("Log_intensity","Condition")
  
  p <- ggplot(df_average, aes(x = Log_intensity, y = Condition)) +
    geom_density_ridges(
      jittered_points = TRUE, quantile_lines = TRUE, scale = 0.6, alpha = 0.9,
      vline_size = 1.2, vline_color = "red", vline_alpha = 1,
      point_size = 0.8, point_alpha = 1,
      position = position_raincloud(adjust_vlines = TRUE)
    )+labs(title = paste("Average Log intensity",tissue))
  return(p)
}


###################################################################################
### Anova functions
##################################################################################
#### Volcano plot

draw_volcanoplot <- function(data, condition="Liver_Tumor_Igg"){
  ds <- data[[condition]]
  fold_cutoff = 1
  pvalue_cutoff = 0.05
  col <- ifelse(ds$Significant == "Yes"|ds$ID=="LORF1", "red", "grey50")
  lab <- ifelse(ds$Significant=="Yes"|ds$ID=="LORF1",ds$ID,"")
  png(paste0("../Image/Volcano_plot/Volcano_plot_",condition,".png"),width = 800, height = 600)
  
  p <- ggplot(ds) +
    geom_point(aes(logfold, -log10(p.adj)),color = col) +
    geom_text_repel(aes(logfold, -log10(p.adj),label = lab))+
    geom_vline(xintercept = fold_cutoff, col = "blue")+
    geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
    ggtitle(condition)
  
  
  print(p)
  dev.off()
}


#######################################################################################
#### Normalize functions
######################################################################################

numerize<- function(x) ifelse(x=="Yes",1,0)

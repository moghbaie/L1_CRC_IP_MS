## Mehrnoosh Oghbaie
## 12/11/2018
## Integrate normalized Values

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

## install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","magrittr","igraph","sqldf","stringr","corrplot","ggplot2","ggrepel","limma","rgl")
bioconductor.packages <- c("biomaRt")
source("Functions/All_Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
load(file="../Result/4. ANOVA_result/Significant.RData")

## Normalization consists of two stages:
##  1. Normalize intensities by LORF1 intensity
##  2. Take average of normalized intensity on all cases

run_order <- rbind(cbind(rep("Ovary",2),c("Igg","Tumor")),
                   cbind(rep("Liver",3),c("Igg","Normal","Tumor")),
                   cbind(rep("Colon",3),c("Igg","Normal","Tumor")))


average_list <- data.frame(matrix(NA, ncol=0,nrow=dim(KM050217_merged)[1]))
average_list$id <- KM050217_merged$id
average_list$uniprotID <- KM050217_merged$uniprotID
  
for(i in unique(run_order[,1])){
  df <- list_experiment[[i]]
  df <- df[grepl("id|uniprotID|Tumor",colnames(df))]
  df[,-(1:2)] <- df[grepl("Tumor",colnames(df))]/t(df[df[["id"]]=="LORF1",grepl("Tumor",colnames(df))])
  df[["mean"]] <- apply(df[,-(1:2)],1,mean)
  average_list[[paste0("Average.",i)]]<- df$mean[match(average_list$id,df$id)]
}

average_list[is.na(average_list)] <- 0
average_list[-(1:2)][average_list[-(1:2)]>1] <- 1
list_experiment[["Averaged.normalized"]]<- average_list

####################################################################################################################
### Draw venn diagram
venn_list <- list_experiment[["Averaged.normalized"]][-(1:2)]
venn_list[venn_list!=0] <- 1
png(filename="../Image/Integrated_plot/VennDiagram.png",width = 600, height = 600)
a <- vennCounts(venn_list)
vennDiagram(a,circle.col = c("darkmagenta", "darkblue",  "orange"), main="Comparison between expressed proteins")
dev.off()


###################################################################################################################
### Draw heatmap of average value for significant ones
heatmap_average <- list_experiment[["Averaged.normalized"]][(1:2)]

numerize<- function(x) ifelse(x=="Yes",1,0)

for(i in unique(run_order[,1])){
  coln <- names(list_experiment)[grepl(i,names(list_experiment))&grepl("_",names(list_experiment))]
  orlist <- list_experiment[[coln[1]]][,(1:2)]
  for(j in coln){
    orlist[[j]] <- list_experiment[[j]]$Significant
  } 
  heatmap_average[[paste0("Average.",i)]] <- apply(orlist[-(1:2)],1,function(x) ifelse("Yes" %in% x,1,0))[match( heatmap_average$id,orlist$ID)]
}
heatmap_average[is.na(heatmap_average)]<-0

dl <- cbind(melt(heatmap_average%>%filter(uniprotID%in% Significant_list$uniprotID)), 
            melt(average_list%>%filter(uniprotID%in% Significant_list$uniprotID))[-(1:3)]%>%dplyr::rename("value.2"=value))
colnames(dl) <- c("Gene_Name","uniprotID","Tissues", "Significance","Expression")

colon_significant <- unique(unlist(c(list_experiment[["Colon_Tumor_Igg"]] %>% filter(Significant=="Yes") %>% dplyr::select(ID),
                       list_experiment[["Colon_Tumor_Normal"]] %>% filter(Significant=="Yes") %>% dplyr::select(ID)
                       )))
dl2 <- dl %>% filter(Gene_Name %in% colon_significant)

sequence_length = length(unique(dl2$Gene_Name))
first_sequence = c(1:(sequence_length%/%2)) 
second_sequence = c((sequence_length%/%2+1):sequence_length) 
first_angles = c(90 - 180/length(first_sequence) * first_sequence)
second_angles = c(-90 - 180/length(second_sequence) * second_sequence)

## polar plot with only colon significant
png(filename="../Image/Integrated_plot/heatmap_colon_significant.png",width = 1200, height = 1200)
q <- ggplot(dl2, aes(Tissues, Gene_Name)) + 
  geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
  scale_alpha(range = c(0.5, 1))+
  geom_tile(alpha=0.6-dl2$Significance*(3/5), fill="gray")+
  #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
  scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
  scale_y_discrete(position = "right")+
  coord_polar(theta = "y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle= c(first_angles,second_angles),size=12),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=16,face="bold"))
 
print(q)
dev.off()

## heatmap with all significant
png(filename="../Image/Integrated_plot/heatmap_all significant.png",width = 400, height = 3000)
u <- ggplot(dl, aes(Tissues, Gene_Name)) + 
  geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
  scale_alpha(range = c(0.5, 1))+
  geom_tile(alpha=0.6-dl$Significance*(3/5), fill="gray")+
  #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
  scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
  scale_y_discrete(position = "right")+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=8,face="bold"),
        axis.title=element_text(size=14,face="bold"))

print(u)
dev.off()
###################################################################################################################
### MDS plot
dz <- average_list%>%filter(uniprotID%in% Significant_list$uniprotID)
rownames(dz) <- dz$id
d4 = dist(dz[-(1:2)])

fit4 = cmdscale(d4, eig=TRUE, k=3) #k is number of dimensions
fit4
dim1 = fit4$points[,1]
dim2 = fit4$points[,2]
dim3 = fit4$points[,3]

colfunc <- colorRampPalette(c("green", "red3"))
color <- data.frame(seq(10)/10,colfunc(10))
color <- color$colfunc.10.[match(round(dz$Average.Colon,1),color$seq.10..10)]

plot3d(dim1, dim2, dim3, type="s", size =1, lwd=4,col=color)
text3d(dim1, dim2, dim3+0.08, rownames(dz), fontweight="bold", cex= 0.7, col="black")
#create a spinning object
s = spin3d(axis = c(0,0,1), rpm=3)
#play 3d
play3d(s, duration=15)
### You need to have Imagemagick installed and give the path to convert.exe to it
imconvertstring<-"\"c:\\Program Files\\ImageMagick-7.0.8-Q16\\convert.exe\" -delay 1x%d %s*.png %s.%s"
movie3d(s, duration=15, dir="../Image/Integrated_plot/", clean=T, convert = imconvertstring, type = "gif")




####################################################################################################################
## Save the output in "../Result/5. Integrated_normalized_data/"
####################################################################################################################


save(KM050217_merged,Significant_proteins,Significant_list,list_experiment, file="../Result/5. Integrated_normalized_data/Normalized.RData")

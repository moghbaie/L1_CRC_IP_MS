## Mehrnoosh Oghbaie
## 12/11/2018
## Integrate normalized Values

## Normalization consists of two stages:
##  1. Normalize intensities by LORF1 intensity
##  2. Take average of normalized intensity on all cases

Template$set("public","experimentAveragedNormalized", list())
Template$set("public","experimentAveraged", list())
Template$set("public","calculateAverageNormalized", function(run_order2){
  average_list <- data.frame(matrix(NA, ncol=0,nrow=dim(self[["input_merged"]])[1]))
  average_list$id <- self[["input_merged"]]$id
  average_list$uniprotID <- unlist(lapply(self[["input_merged"]]$uniprotID,function(x) strsplit(x, ";")[[1]][1]))
  average_list$uniprotID <- unlist(lapply(self[["input_merged"]]$uniprotID,function(x) strsplit(x, "-")[[1]][1]))
  average_unnormalized_list <- average_list
  for(i in unique(run_order2[,1])){
    print(i)
    iBAQ_col <- colnames(self[["input_merged"]])[grepl("iBAQ.",colnames(self[["input_merged"]]))]
    
    #df <- self[["experimentImputed"]][[i]]
    con_list <- run_order2[run_order2[,1]==i,2]
    
    for(j in con_list[grepl("T_ORF1",con_list)]){
      dt <- self$input_merged[,c("id","uniprotID",iBAQ_col[grepl(j,iBAQ_col)])]
      
      for(s in colnames(dt[,-c(1:2)])){
        dt[s] <- log(unlist(unname(dt[s])))
        dt[is.infinite(unname(unlist(dt[s]))),s] <- 0
      }
      dy <- dt[grepl("T_ORF1",colnames(dt))]/t(apply(dt[grepl("Q9UN81|L1RE1",dt$uniprotID),grepl("T_ORF1",colnames(dt))],2,sum,na.rm="T"))
      dw <-  dt[grepl("T_ORF1",colnames(dt))]
      dt[["mean"]] <- apply(dy,1, function(x) mean(x[x!=0],na.rm=T))
      dt[["mean_unnormalized"]] <- apply(dw,1,function(x) mean(x[x!=0],na.rm=T))
      
      
      average_list[[paste0("Average.",j)]]<- dt$mean
      average_unnormalized_list[[paste0("Average.unnormalized.",j)]]<- dt$mean_unnormalized
    }
  }
  average_list[is.na(average_list)] <- 0
  average_list <- average_list[self$input_merged$Potential.contaminant!="+"|self$input_merged$Reverse!="+",]
  
  average_unnormalized_list[is.na(average_unnormalized_list)] <- 0
  average_unnormalized_list <- average_unnormalized_list[self$input_merged$Potential.contaminant!="+",]
  
  
  names <- names(average_list)[c(-1,-2)]
  names <- gsub("144T","Tumor-A",names)
  names <- gsub("159T","Tumor-B",names)
  names <- gsub("163T","Tumor-C",names)
  names(average_list)[c(-1,-2)] <- names
  names(average_unnormalized_list)[c(-1,-2)] <- gsub("Average","Average.unnormalized",names)
  self[["experimentAveragedNormalized"]]<- average_list
  self[["experimentAveraged"]] <-  average_unnormalized_list
}
)

####################################################################################################################
### Draw venn diagram
Template$set("public","drawCommonVenndiagram", function(){
  venn_list <- self[["experimentAveragedNormalized"]][-(1:2)]
  venn_list[venn_list!=0] <- 1
  png(filename="../Image/Integrated_plot/VennDiagram.png",width = 1000, height = 800)
  a <- vennCounts(venn_list)
  vennDiagram(a,circle.col = c("darkmagenta", "darkblue",  "orange"), main="Comparison between expressed proteins")
  dev.off()
})

###################################################################################################################
### Draw heatmap of average value for significant ones
### Draw heatmap of average value for significant ones
Template$set("public","drawHeatmap", function(){
  
  heatmap_average <- self[["experimentAveragedNormalized"]][,c(1:2)]
  #heatmap_average <- heatmap_average[heatmap_average$uniprotID!="0",]
  for(i in unique(run_order2[,1])){
    coln <- names(self[["experimentImputateSignificant"]])[grepl(i,names(self[["experimentImputateSignificant"]]))&grepl("_",names(self[["experimentImputateSignificant"]]))]
    
    for(j in coln){
      orlist <- self[["experimentImputateSignificant"]][[j]][,c(1:2)]
      orlist$uniprotID <- unlist(lapply(orlist$uniprotID, function(x) strsplit(x,";")[[1]][1]))
      orlist$uniprotID <- unlist(lapply(orlist$uniprotID, function(x) strsplit(x,"-")[[1]][1]))
      orlist[[j]] <- self[["experimentImputateSignificant"]][[j]]$Significant
      heatmap_average[[paste0("Average.",j)]] <- apply(orlist[-(1:2)], 1, function(x) ifelse("Yes" %in% x,1,0))[match(heatmap_average$uniprotID, orlist$uniprotID)]
    } 
  }
  heatmap_average[is.na(heatmap_average)]<-0
  
  heatmap_average[["Average.159T_ORF1_3"]] <- ifelse(heatmap_average$Average.Liver_Tumor_159T_IgG_4+heatmap_average$Average.Liver_Tumor_159N_ORF1_5 >0, 1, 0)
  heatmap_average[["Average.144T_ORF1_1"]] <- heatmap_average[["Average.Ovary_Tumor_144T_IgG_2"]]
  heatmap_average[["Average.144T_ORF1_9"]] <- heatmap_average[["Average.Ovary_Tumor_144T_IgG_10"]]
  heatmap_average[["Average.163T_ORF1_6"]] <- ifelse(heatmap_average$Average.Colon_Tumor_163T_IgG_7+heatmap_average$Average.Colon_Tumor_163N_ORF1_8 >0, 1, 0)
  
  sel_cols <- c("id", "uniprotID",colnames(heatmap_average)[grepl("T_ORF1_",colnames(heatmap_average))])
  
  heatmap_average <- heatmap_average[,grepl(paste0(sel_cols,collapse="|"),colnames(heatmap_average))]
  uniprotID_selected <- heatmap_average[apply(heatmap_average[,(3:6)],1, max)>0,"uniprotID"]
  
  uniprotID_sel2 <- heatmap_average[apply(heatmap_average[,(3:6)],1, sum)>1,"uniprotID"]
  uniprotID_sel2 <- unique(c(uniprotID_sel2 , eLife_list$uniprotID))
  
  
  significants <- unlist(lapply(self[["Significant_list"]]$uniprotID, function(x) strsplit(x, ";")[[1]][1]))
  significants <- unlist(lapply(significants, function(x) strsplit(x, "-")[[1]][1]))
  
  colnames(heatmap_average) <- gsub("144T","Tumor-A",colnames(heatmap_average))
  colnames(heatmap_average) <- gsub("159T","Tumor-B",colnames(heatmap_average))
  colnames(heatmap_average) <- gsub("163T","Tumor-C",colnames(heatmap_average))
  
  dl <- cbind(melt(heatmap_average[,-1]%>%dplyr::filter(uniprotID%in% significants)), 
              melt(self[["experimentAveragedNormalized"]][,-1]%>%dplyr::filter(uniprotID%in% significants))[-(1:2)]%>%dplyr::rename("value.2"=value))
  colnames(dl) <- c("uniprotID","condition", "Significance","Expression")
  dl<- dl[dl$uniprotID %in%  uniprotID_selected,]
  
  gene_protein <- read.delim("protein.tab")
  gene_protein <- gene_protein %>% dplyr::mutate(Gene.names = as.character(Gene.names),
                                                 Entry = as.character(Entry))
  gene_protein$Gene.names  <- apply(gene_protein,1, function(x) strsplit(x[["Gene.names"]], " ")[[1]][1])
  #ds$ID <- gene_protein$Gene.names[match(ds$uniprotID,gene_protein$Entry)]
  #ensembl=useMart("ensembl")
  #ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  #gene_protein <- getBM(attributes=c('uniprot_gn', 'external_gene_name'), 
  #                      filters = 'uniprot_gn', 
  #                      values = dl[["uniprotID"]], 
  #                      mart = ensembl)
  geneName_list <- gene_protein$Gene.names[match(dl$uniprotID, gene_protein$Entry)]
  
  dl[["geneName"]] <- ifelse(is.na(geneName_list),dl$uniprotID, geneName_list)
  dl$geneName[dl$geneName=="Q9UN81"] <- "L1RE1"
  
  colon_significant <- unique(unlist(c(self[["experimentImputateSignificant"]][["Colon_Tumor_163T_IgG_7"]] %>% dplyr::filter(Significant=="Yes") %>% dplyr::select(uniprotID),
                                       self[["experimentImputateSignificant"]][["Colon_Tumor_163N_ORF1_8"]] %>% dplyr::filter(Significant=="Yes") %>% dplyr::select(uniprotID))))
  colon_significant <- unlist(lapply(colon_significant, function(x) strsplit(x, ";")[[1]][1]))
  colon_significant <- unlist(lapply(colon_significant, function(x) strsplit(x, "-")[[1]][1]))
  
  dl2 <- dl %>% filter(uniprotID %in% uniprotID_sel2)
  a <-ifelse(unique(dl2[["geneName"]])[order(unique(dl2[["geneName"]]))] %in% unique(eLife_list$geneName),"red","black")
  sequence_length = length(unique(dl2$uniprotID))
  first_sequence = c(1:(sequence_length%/%2)) 
  second_sequence = c((sequence_length%/%2+1):sequence_length) 
  first_angles = c(90 - 180/length(first_sequence) * first_sequence)
  second_angles = c(-90 - 180/length(second_sequence) * second_sequence)
  
  
  ## polar plot with only colon significant
  png(filename="../Image/Integrated_plot/heatmap_significant_eLife.png",width = 1200, height = 1200)
  q <- ggplot(dl2, aes(condition, geneName)) + 
    geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
    scale_alpha(range = c(0.5, 1))+
    geom_tile(alpha=0.6-dl2$Significance*(3/5), fill="gray")+
    #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
    scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
    scale_y_discrete(position = "right")+
    coord_polar(theta = "y")+
    theme_bw()+
    theme(axis.text.x = element_text(angle= c(first_angles,second_angles),size=10, colour = a),
          axis.text.y = element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  print(q)
  dev.off()
  
  ## heatmap with all significan
  png(filename="../Image/Integrated_plot/heatmap_all_significant.png",width = 600, height = 3700)
  u <- ggplot(dl, aes(condition, geneName)) + 
    geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
    scale_alpha(range = c(0.5, 1))+
    geom_tile(alpha=0.6-dl$Significance*(3/5), fill="gray")+
    #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
    scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
    scale_y_discrete(position = "right")+
    theme_bw()+
    theme(axis.text.x = element_text(size=12,angle = 90, hjust = 1),
          axis.text.y = element_text(size=12,face="bold"),
          axis.title=element_text(size=14,face="bold"))
  print(u)
  dev.off()
})

###################################################################################################################
### Draw barplot of value for significant ones
Template$set("public","drawHeatmapUnnormalized", function(){
  
  heatmap_average <- self[["experimentAveraged"]][,(1:2)]
  for(i in unique(run_order2[,1])){
    coln <- names(self[["experimentImputateSignificant"]])[grepl(i,names(self[["experimentImputateSignificant"]]))&grepl("_",names(self[["experimentImputateSignificant"]]))]
    
    for(j in coln){
      orlist <- self[["experimentImputateSignificant"]][[j]][,c(1:2)]
      orlist$uniprotID <- unlist(lapply(orlist$uniprotID, function(x) strsplit(x,";")[[1]][1]))
      orlist$uniprotID <- unlist(lapply(orlist$uniprotID, function(x) strsplit(x,"-")[[1]][1]))
      orlist[[j]] <- self[["experimentImputateSignificant"]][[j]]$Significant
      heatmap_average[[paste0("Average.",j)]] <- apply(orlist[-(1:2)],1,function(x) ifelse("Yes" %in% x,1,0))[match(heatmap_average$uniprotID,orlist$uniprotID)]
    } 
  }
  heatmap_average[is.na(heatmap_average)]<-0
  
  heatmap_average[["Average.unnormalized.159T_ORF1_3"]] <- ifelse(heatmap_average$Average.Liver_Tumor_159T_IgG_4 + heatmap_average$Average.Liver_Tumor_159N_ORF1_5>0,1,0)
  heatmap_average[["Average.unnormalized.144T_ORF1_1"]] <- heatmap_average[["Average.Ovary_Tumor_144T_IgG_2"]]
  heatmap_average[["Average.unnormalized.144T_ORF1_9"]] <- heatmap_average[["Average.Ovary_Tumor_144T_IgG_10"]]
  heatmap_average[["Average.unnormalized.163T_ORF1_6"]] <- ifelse(heatmap_average$Average.Colon_Tumor_163T_IgG_7 + heatmap_average$Average.Colon_Tumor_163N_ORF1_8>0,1,0)
  
  
  sel_cols <- c("id", "uniprotID",colnames(heatmap_average)[grepl("T_ORF1_",colnames(heatmap_average))])
  heatmap_average <- heatmap_average[,grepl(paste0(sel_cols,collapse="|"),colnames(heatmap_average))]
  uniprotID_selected <- heatmap_average[apply(heatmap_average[,(3:6)],1, max)>0,"uniprotID"]
  
  significants <- unlist(lapply(self[["Significant_list"]]$uniprotID, function(x) strsplit(x, ";")[[1]][1]))
  significants <- unlist(lapply(significants, function(x) strsplit(x, "-")[[1]][1]))
  
  colnames(heatmap_average) <- gsub("144T","Tumor-A",colnames(heatmap_average))
  colnames(heatmap_average) <- gsub("159T","Tumor-B",colnames(heatmap_average))
  colnames(heatmap_average) <- gsub("163T","Tumor-C",colnames(heatmap_average)) 
  
  dl <- cbind(melt(heatmap_average[,-1]%>%dplyr::filter(uniprotID%in% significants)), 
              melt(self[["experimentAveraged"]][,-1]%>%dplyr::filter(uniprotID%in% significants))[-(1:2)]%>%dplyr::rename("value.2"=value))
  colnames(dl) <- c("uniprotID","condition", "Significance","Expression")
  dl<- dl[dl$uniprotID %in%  uniprotID_selected,]
  
  dl <- dl[!dl$uniprotID %in% self$input_merged$uniprotID[self$input_merged$Potential.contaminant=="+"],]
  
  gene_protein <- read.delim("protein.tab")
  gene_protein <- gene_protein %>% dplyr::mutate(Gene.names = as.character(Gene.names),
                                                 Entry = as.character(Entry))
  gene_protein$Gene.names  <- apply(gene_protein,1, function(x) strsplit(x[["Gene.names"]], " ")[[1]][1])
  
  #ensembl=useMart("ensembl")
  #ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  #gene_protein <- getBM(attributes=c('uniprot_gn', 'external_gene_name'), 
  #                      filters = 'uniprot_gn', 
  #                      values = dl[["uniprotID"]], 
  #                      mart = ensembl)
  
  geneName_list <- gene_protein$Gene.names[match(dl$uniprotID, gene_protein$Entry)]
  #geneName_list <- gene_protein$external_gene_name[match(dl$uniprotID, gene_protein$uniprot_gn)]
  
  dl[["geneName"]] <- ifelse(is.na(geneName_list),dl$uniprotID, geneName_list)
  dl$geneName[dl$geneName=="Q9UN81"] <- "L1RE1"
  colon_significant <- unique(unlist(c(self[["experimentImputateSignificant"]][["Colon_Tumor_163T_IgG_7"]] %>% dplyr::filter(Significant=="Yes") %>% dplyr::select(uniprotID),
                                       self[["experimentImputateSignificant"]][["Colon_Tumor_163N_ORF1_8"]] %>% dplyr::filter(Significant=="Yes") %>% dplyr::select(uniprotID))))
  
  colon_significant <- unlist(lapply(colon_significant, function(x) strsplit(x, ";")[[1]][1]))
  colon_significant <- unlist(lapply(colon_significant, function(x) strsplit(x, "-")[[1]][1]))
  
  uniprotID_sel2 <- heatmap_average[apply(heatmap_average[,(3:6)],1, sum)>1,"uniprotID"]
  uniprotID_sel2 <- unique(c(uniprotID_sel2 , eLife_list$uniprotID))
  
  dl2 <- dl %>% filter(uniprotID %in%
                         unlist(lapply(colon_significant,function(x) strsplit(x, "-")[[1]][1])))
  a <-ifelse(unique(dl2[["geneName"]])[order(unique(dl2[["geneName"]]))] %in% unique(eLife_list$geneName),"red","black")
  sequence_length = length(unique(dl2$uniprotID))
  first_sequence = c(1:(sequence_length%/%2)) 
  second_sequence = c((sequence_length%/%2+1):sequence_length) 
  first_angles = c(90 - 180/length(first_sequence) * first_sequence)
  second_angles = c(-90 - 180/length(second_sequence) * second_sequence)
  
  ## polar plot with only colon significant
  png(filename="../Image/Integrated_plot/heatmap_significant_eLife_unnormalized.png",width = 1200, height = 1200)
  q <- ggplot(dl2, aes(condition, geneName)) + 
    geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
    scale_alpha(range = c(0.5, 1))+
    geom_tile(alpha=0.6-dl2$Significance*(3/5), fill="gray")+
    #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
    scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
    scale_y_discrete(position = "right")+
    coord_polar(theta = "y")+
    theme_bw()+
    theme(axis.text.x = element_text(angle= c(first_angles,second_angles),size=10, colour = a),
          axis.text.y = element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  print(q)
  dev.off()
  
  
  ## heatmap with all significan
  png(filename="../Image/Integrated_plot/heatmap_all_significant_unnormalized.png",width = 600, height = 3700)
  u <- ggplot(dl, aes(condition, geneName)) + 
    geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
    scale_alpha(range = c(0.5, 1))+
    geom_tile(alpha=0.6-dl$Significance*(3/5), fill="gray")+
    #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
    scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
    scale_y_discrete(position = "right")+
    theme_bw()+
    theme(axis.text.x = element_text(size=12,angle = 90, hjust = 1),
          axis.text.y = element_text(size=12,face="bold"),
          axis.title=element_text(size=14,face="bold"))
  print(u)
  dev.off()
})


###################################################################################################################
### MDS plot
Template$set("public","drawMDSplot", function(){
  average_list <- self$experimentAveragedNormalized
  dz <- average_list%>%filter(uniprotID%in% self$Significant_list$uniprotID)
  rownames(dz) <- dz$id
  d4 = dist(dz[-(1:2)])
  
  fit4 = cmdscale(d4, eig=TRUE, k=3) #k is number of dimensions
  dim1 = fit4$points[,1]
  dim2 = fit4$points[,2]
  dim3 = fit4$points[,3]
  
  colfunc <- colorRampPalette(c("green", "red3"))
  color <- data.frame(seq(10)/10,colfunc(10))
  color <- color$colfunc.10.[match(round(dz$Average.163T_ORF1_6,1),color$seq.10..10)]
  
  gene_protein <- read.delim("protein.tab")
  gene_protein <- gene_protein %>% dplyr::mutate(Gene.names = as.character(Gene.names),
                                                 Entry = as.character(Entry))
  gene_protein$Gene.names  <- apply(gene_protein,1, function(x) strsplit(x[["Gene.names"]], " ")[[1]][1])
  
  name <- gene_protein$Gene.names[match(dz$uniprotID, gene_protein$Entry)]
  
  dz[["GeneName"]] <- ifelse(!is.na(name), name , dz$uniprotID)
  dz[["GeneName"]][dz[["GeneName"]]=="Q9UN81"] <- "L1RE1"
  plot3d(dim1, dim2, dim3, type="s", size =1, lwd=4,col=color)
  text3d(dim1, dim2, dim3+0.08, ifelse(!is.na(color),dz[["GeneName"]],""), fontweight="bold", cex= 0.7, col="black")
  #create a spinning object
  s = spin3d(axis = c(0,0,1), rpm=2)
  #play 3d
  play3d(s, duration=33)
  ### You need to have Imagemagick installed and give the path to convert.exe to it
  imconvertstring<-"\"c:\\Program Files\\ImageMagick-7.0.8-Q16\\convert.exe\" -delay 1x%d %s*.png %s.%s"
  movie3d(s, duration=33, dir="../Image/Integrated_plot/", clean=T, convert = imconvertstring, type = "gif")
  
})
#######################################################################


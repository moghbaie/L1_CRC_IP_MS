## Mehrnoosh Oghbaie
## 12/11/2018
## ANOVA test between cases and controls



##  1. Filter proteins with less than two peptides (I'll skip this one for now)
##  2. Perform t-test between cases and controls
##  3. Adjust pvalueswith Benjamin Hochberg correction test
##  4. Select significant proteins with (p.adj < 0.05 & log2fold > 1)

Template$set("public","anovaAnalysis", function(run_order){
for( i in unique(run_order[,1])){
  y <- self[["experimentImputed"]][[i]]
  z <- run_order[run_order[,1]==i,2]
  for(j in z[z!="Tumor"]){
    
      dt <- data.frame(matrix(NA,ncol=0,nrow=dim(y)[1]))
      dt[["ID"]] <- y[["id"]]
      dt[["uniprotID"]] <- y[["uniprotID"]]
      dt[["p.value"]] <- NA
      dt[["logfold"]] <- NA
      dt[["Significant"]] <- NA
      for(k in 1:dim(dt)[1]){
        dt[[k,"p.value"]] <- t.test(unname(unlist(y[k,grepl("Tumor",colnames(y))])),unname(unlist(y[k,grepl(j,colnames(y))])))$p.value
        dt[[k,"logfold"]] <- mean(unname(unlist(y[k,grepl("Tumor",colnames(y))])))-mean(unname(unlist(y[k,grepl(j,colnames(y))])))
      }
      dt[["p.adj"]] <- p.adjust(dt[["p.value"]], method = "BH", n = length(dt[["p.value"]]))
      dt[["Significant"]] <- ifelse(dt[["p.adj"]]<0.05&dt[["logfold"]]>1,"Yes","No")
      self[["experimentImputed"]][[paste0(i,"_","Tumor","_",j)]] <-dt
  }
}
})

####################################################################################################################
## Venn diagram
Template$set("public","Significant_list", NA)
Template$set("public","drawVenDiagramSignificant", function(){
Significant_proteins <- data.frame(matrix(NA, nrow=length(self[["input_merged"]]$uniprotID), ncol= 2+length( names(self[["experimentImputed"]])[grep("_",names(self[["experimentImputed"]]))])))
colnames(Significant_proteins) <- c("id","uniprotID", names(self[["experimentImputed"]])[grep("_",names(self[["experimentImputed"]]))])
Significant_proteins$uniprotID <- self[["input_merged"]]$uniprotID
Significant_proteins$id <- self[["input_merged"]]$id


for(i in names(self[["experimentImputed"]])[grep("_",names(self[["experimentImputed"]]))]){
  Significant_proteins[[i]] <- self[["experimentImputed"]][[i]][,"Significant"][match(Significant_proteins$uniprotID,self[["experimentImputed"]][[i]][,"uniprotID"])]
  Significant_proteins[[i]][is.na(Significant_proteins[[i]])] <- "NO"
  Significant_proteins[[i]] <- ifelse(Significant_proteins[[i]]=="Yes",1,0)
}
png(paste0("../Image/Volcano_plot/Significant_Venndiagram.png"),width = 800, height = 600)
a <- vennCounts(Significant_proteins[,-(1:2)])
vennDiagram(a,circle.col = c("darkmagenta", "darkblue", "pink1", "skyblue", "orange"), main="Comparison significant proteins")
print((a))
dev.off()

Significant_list <- Significant_proteins%>%
  dplyr::filter(.[[3]]|.[[4]]|.[[5]]|.[[6]]|.[[7]]==1)%>%
  dplyr::select(id,uniprotID)
self[["Significant_list"]] <- Significant_list 

})


####################################################################################################################
## 5. Volcano plot "../Image/Volcano_plot/"
####################################################################################################################
Template$set("public","drawVolcanoPlot", function(){
  for(condition in names(self[["experimentImputed"]])[grep("_",names(self[["experimentImputed"]]))]){
    draw_volcanoplot(data=self[["experimentImputed"]],condition=condition)
    Sys.sleep(2)
  }
})




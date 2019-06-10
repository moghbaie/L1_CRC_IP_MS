## Mehrnoosh Oghbaie
## 05/21/2019
## Draw mutation from consensus in ORF1s sequences


Template$set("public","drawMutation", function(){
  detected_orf1 <- CRC$Significant_list$uniprotID[grepl("ORF1P",CRC$Significant_list$uniprotID)]
  detected_orf1 <- c(detected_orf1)
  
  peptides <- read.delim(paste0("../Input_data/",strsplit(as.character(CRC$input$V2[1]),"/")[[1]][1],"/txt/peptides.txt"))
  
  unlist(lapply(pep$Proteins, function(x) lapply(strsplit(as.character(x),";")[[1]], function(y) any(grepl(paste0(c(detected_orf1),collapse = "|"),y, fixed =TRUE)))))
 
  pep <- peptides[(peptides$Leading.razor.protein %in% detected_orf1)&(peptides$Missed.cleavages==0),
                c(colnames(peptides)[c(1,2,3,34,35,36,37,38)],colnames(peptides)[grepl("Intensity.", colnames(peptides))])]

  orf1_dir <- "../Input_data/orf1.fasta"
  mySequences <- readAAStringSet(orf1_dir)
  name_selected <- names(mySequences)[unlist(lapply(names(mySequences), function(x) strsplit(x,"\\|")[[1]][2])) %in% detected_orf1]
  myseq <- mySequences[name_selected]
  
  myseq1 <- myseq
  names(myseq1) <- paste(lapply(names(myseq1), function(x) strsplit(x,"\\|")[[1]][2]),lapply(names(myseq1), function(x) strsplit(x,"\\.")[[1]][5]), sep="-")
  
  orf1_align <- msa(myseq1)
  orf1_align2 <- msaConvert(orf1_align, type="seqinr::alignment")
  
  d <- dist.alignment(orf1_align2, "identity")
  orf1Tree <- nj(d)
  png(filename="../Image/Phylogenetic_tree/tree.png",width = 1000, height = 900)
  plot(orf1Tree, main="ORF1_mutations")
  dev.off()
  
  names(myseq) <- paste(lapply(names(myseq), function(x) strsplit(x,"\\|")[[1]][2]),lapply(names(myseq), function(x) strsplit(x,"\\.")[[1]][5]), sep="-")
  
  pep_observe <- list()
  for(name in names(CRC$experimentImputateSignificant)){
    print(name)
   
    case1more <- CRC[["experimentPreImputateSignificant"]][[name]][["Case"]][match(CRC$experimentImputateSignificant[[name]][["uniprotID"]],CRC[["experimentPreImputateSignificant"]][[name]][,"uniprotID"])]
    significant_name <- CRC$experimentImputateSignificant[[name]][["uniprotID"]][CRC$experimentImputateSignificant[[name]][["Significant"]]=="Yes"&case1more>1]

    significant_name_ORF <- significant_name[grepl("ORF1",significant_name)]
   # if("ORF1P.76" %in% significant_name){significant_name <- c(significant_name, "ORF1P.2")}
   # if("ORF1P.46" %in% significant_name){significant_name <- c(significant_name, "ORF1P.81","ORF1.P1")}
    
    myseq_Tumor_name <- myseq[lapply(names(myseq), function(x) strsplit(x,"-")[[1]][1]) %in% significant_name_ORF]
    if(length(myseq_Tumor_name)>1){
      orf1_align_name <- msa(myseq_Tumor_name)
      orf1_align_name2 <- msaConvert(orf1_align_name, type="seqinr::alignment")
      
      aln = data.frame(
        letter=strsplit(paste0(orf1_align_name2$seq, collapse=""), "")[[1]], 
        mutation = rep(orf1_align_name2$nam, each=338),
        x       = rep(1:338, orf1_align_name2$nb)
      )
      
      peptide_consensus <- as.character(peptides[peptides$Leading.razor.protein=="Q9UN81","Sequence"])
      
      pep[is.na(pep)] <- 0
      pep_orf1 <- pep[as.character(pep$Leading.razor.protein) %in% lapply(orf1_align_name2$nam, function(x) strsplit(x,"-")[[1]][1]) & (apply(pep[,colnames(pep)[grepl(paste0(strsplit(name, "_|N|T")[[1]][4],"T_ORF1_"), colnames(pep))]],1,sum)!=0), c("Sequence","Leading.razor.protein", "Proteins","Start.position","End.position" , colnames(pep)[grepl(paste0(strsplit(name, "_|N|T")[[1]][4],"T_ORF1_"), colnames(pep))])]
      pep_orf2 <- cbind(pep_orf1[,c(1,2,3,4,5)], apply(pep_orf1[,c(-1,-2,-3,-4,-5)],1, mean))
      pep_orf3 <- pep_orf2[pep_orf2$Sequence %in% peptide_consensus,]
      
      pep_observe[[name]] <- pep_orf2
      
      aln[["mut"]] = 'no'
      for(i in as.character(pep_orf2[["Sequence"]])){
        aln[(aln[["x"]] %in% c(as.numeric(as.character(pep_orf2[as.character(pep_orf2[["Sequence"]])==i, "Start.position"])):as.numeric(as.character(pep_orf2[as.character(pep_orf2[["Sequence"]])==i, "End.position"]))))& lapply(as.character(aln[["mutation"]]), function(x) strsplit(x,"-")[[1]][1])== as.character(pep_orf2[as.character(pep_orf2[["Sequence"]])==i,  "Leading.razor.protein" ]),"mut"] <- "yes"
      }
      
      aln$consensus <- "yes"
      Lorf1 <- strsplit("MGKKQNRKTGNSKTQSASPPPKERSSSPATEQSWMENDFDELREEGFRRSNYSELREDIQTKGKEVENFEKNLEECITRITNTEKCLKELMELKTKARELREECRSLRSRCDQLEERVSAMEDEMNEMKREGKFREKRIKRNEQSLQEIWDYVKRPNLRLIGVPESDVENGTKLENTLQDIIQENFPNLARQANVQIQEIQRTPQRYSSRRATPRHIIVRFTKVEMKEKMLRAAREKGRVTLKGKPIRLTADLSAETLQARREWGPIFNILKEKNFQPRISYPAKLSFISEGEIKYFIDKQMLRDFVTTRPALKELLKEALNMERNNRYQPLQNHAKM","")[[1]]
      aln$LORF1 <- NA
      for( i in as.character(aln$mutation)){
        aln$LORF1[aln$mutation==i] <- Lorf1[aln$x[aln$mutation==i]]
      }
      aln$consensus <- ifelse(as.character(aln$letter)==as.character(aln$LORF1),"yes","no")
      
      dt <- data.frame(cbind(Lorf1,rep("L1RE1",338),c(1:338),rep("yes",338),rep("yes",338),Lorf1))
      colnames(dt) <- colnames(aln)
      dt$x <- as.integer(as.character(dt$x))
      dt$mut <- as.character(dt$mut)
      dt$consensus <- as.character(dt$consensus)
      dt$LORF1 <- as.character(dt$LORF1)
      aln <- rbind(aln,dt)
      aln[["Change"]] <- ifelse(aln[["consensus"]]=="no", paste0(aln[["LORF1"]],aln[["x"]],aln[["letter"]]),"")
      aln[["Change2"]] <- ""
      for(z in as.numeric(rownames(aln)[aln[["mut"]]=="yes"])){

        if(z-1>1){
          aln[z -1 , "Change2"] <- ifelse(aln[z-1, "Change"]!="",aln[z-1, "Change"],"")
        }
        if(z+1< dim(aln)[1]){
          aln[z +1 , "Change2"] <- ifelse(aln[z+1, "Change"]!="",aln[z+1, "Change"],"")        
          }
        aln[z , "Change2"] <- ifelse(aln[z, "Change"]!="",aln[z, "Change"],"")        
        
      }
      sequence_length = length(unique(aln$x))
      first_sequence = c(1:(sequence_length%/%2)) 
      second_sequence = c((sequence_length%/%2+1):sequence_length) 
      first_angles = c(90 - 180/length(first_sequence) * first_sequence)-90
      second_angles = c(-90 - 180/length(second_sequence) * second_sequence)+90
      png(filename=paste0("../Image/Phylogenetic_tree/",name,".png"),height = (orf1_align_name2$nb+2)*25, width =2200)
      #png(filename=paste0("../Image/Phylogenetic_tree/",name,".png"),width = 1900, height =1900)
      
      # Generate the sequence alignment
      p2 = ggplot(aln, aes(x, mutation)) +
        #geom_text(aes(label=letter, color=consensus, size=mut, angle= rep(c(first_angles,second_angles),length(unique(aln$mutation))))) + 
        geom_text(aes(label=letter, color=consensus, size=mut)) + 
        scale_x_continuous(breaks=1:10, expand = c(0.0105, 0)) + xlab('') + 
       # expand_limits(y=c(-15-(length(unique(aln$mutation))-2)/2, 10+(length(unique(aln$mutation))-2)/4))+
        scale_color_manual(values=c('red', "black")) + 
        scale_size_manual(values=c(2.5, 5)) + 
        theme_logo() + 
        #coord_polar(theta = "x")+
        theme(legend.position = 'none')
      print(p2)
      dev.off()
      
      dy <- data.frame(table(aln[aln$consensus=="no"& aln$mut=="yes",c("Change2")]))
      dy <- dy[order(-dy$Freq),]
      dy[["loc"]] <- as.numeric(as.character(gsub('[A-Z]+|-', '', dy$Var1)))
      dy$loc <- as.factor(dy$loc)
      png(filename=paste0("../Image/Phylogenetic_tree/mutation_fre_",name,".png"),width = 200+500* length(unique(dy$Var1))/50, height =200+800* length(unique(dy$Var1))/80)
      p3 <- ggplot(dy, aes(fill=Var1, y=Freq, x=loc,label=Var1)) + 
        geom_bar( stat="identity",colour="black")+
        geom_text(size = 3, position = position_stack(vjust = 0.5))+
        scale_y_reverse()+
        coord_flip()+ggtitle(paste0("Mutations in ",name))
      
      print(p3)
      dev.off()
      
    }
  }
  
  write_xlsx(pep_observe,"../Image/Phylogenetic_tree/mutated_peptides.xlsx")
  
})
             




## Mehrnoosh Oghbaie
## 12/11/2018
## Preparing data from MaxQuant or ...

## Data preparation consists of four stages:
##  1. Calculate iBAQ intensity
##  2. Remove contaminants and reverse proteins
##  3. Log transformation
##  4. Separate different experiment


Template <- R6Class("Template",
                   private = list(
                     input.info.dir = "../Input_data/Input.info",
                     df = NA,
                     df_log=NA
                   ),
                   active = list(
                     input_dir = function() {
                       dirname(private$input.info.dir)
                     }, 
                     input = function(){
                       read.delim(private$input.info.dir, header=FALSE)
                     },
                     data_type = function(){
                       self$input %>%
                         dplyr::filter(V1=="data_type") %>%
                         dplyr::mutate(V2=as.character(V2))%>% 
                         dplyr::select(V2) %>% 
                         .$V2
                     }, 
                    
                     input_merged = function(){
                       file_name <- self$input %>%
                         dplyr::filter(V1=="MaxQuant_Output") %>%
                         dplyr::mutate(V2=as.character(V2))%>% 
                         dplyr::select(V2) %>% 
                         .$V2

                      tbl <- read.delim(paste(self$input_dir,file_name, sep="/"))
                      tbl$uniprotID <- apply(tbl, 1, function(x) strsplit(as.character(x[["Protein.IDs"]]), ";")[[1]][1])
                         #tbl$uniprotID <- apply(tbl, 1, function(x) strsplit(as.character(x[["uniprotID"]]), "-")[[1]][1])
                         #tbl$uniprotID <- ifelse(grepl("ORF1|ORF2",toupper(as.character(tbl$Protein.IDs))),as.character(tbl$Protein.IDs),tbl$uniprotID)
 
                       return(tbl)
                     }
                     ,
                     contaminant_list =function(){
                       contaminants_tbl <- read.table(paste(self$input_dir,
                                                            self$input %>%
                                                              dplyr::filter(V1=="Contaminants") %>%
                                                              dplyr::mutate(V2=as.character(V2))%>% 
                                                              dplyr::select(V2) %>% 
                                                              .$V2, sep="/"),
                                                      sep=";", quote="\"")
                       contaminant_list <- unlist(lapply(as.character(contaminants_tbl[,1])[grepl(">", as.character(contaminants_tbl[,1]))], function(x) strsplit(x," |>")[[1]][2]))
                       return(contaminant_list)
                       } 
                      )
                   )

Template$set("public","non_unique_ORF", list())

Template$set("public","removeContaminant", function(){
  cols <- colnames(self$input_merged)
  dl <-  self[["input_merged"]]
  dl <- dl %>% dplyr::mutate(Peptide.counts..unique. = as.character(Peptide.counts..unique.),
                               Protein.IDs = as.character(Protein.IDs),
                               Potential.contaminan = as.character(Potential.contaminant))
  dl$Peptide.counts..unique. <-  apply(dl,1,function(x) strsplit(x[["Peptide.counts..unique."]],";")[[1]][1])
  self[["non_unique_ORF"]] <- dl  %>% dplyr::filter(Peptide.counts..unique.=="0" & grepl("ORF1",Protein.IDs)) %>% .$uniprotID
  private$df <- dl[!(dl$Potential.contaminant=="+"|dl$Reverse=="+"|grepl( "P01857",dl$uniprotID)) ,c("id","uniprotID",cols[grepl(self$data_type,cols)])]    
  #private$df <- dl[!(dl$Potential.contaminant=="+"||dl$Reverse=="+"|grepl( "P01857|P02671|P02675|P02679|B9A064|A2NJV5|P01871",dl$uniprotID)) ,c("id","uniprotID",cols[grepl(self$data_type,cols)])]    

  })


Template$set("public","logTransformation", function(){
  private$df <- private$df
  tp <- lapply(private$df, class)
  df_log <- private$df
  for(i in colnames(private$df)[tp =="numeric"]){
    df_log[i] <- log(unlist(unname(df_log[i])))
    df_log[is.infinite(unname(unlist(df_log[i]))),i] <- NA
  }
  colnames(df_log) <- c("id","uniprotID",unlist(lapply(colnames(df_log)[-(1:2)], function(x) strsplit(x,"\\.")[[1]][3])))
  
  private$df_log <- df_log
 
}
)

Template$set("public","experiment", list())
                      
Template$set("public","separatedList", function(){
  selected_col <- colnames(private$df_log)
  #selected_col <- c("uniprotID", "id", colnames(private$df_log)[grepl(self$data_type,colnames(private$df_log))])
  self$experiment[["Ovary"]] <- private$df_log[,selected_col[grepl("144|uniprotID|id",selected_col)]]
  self$experiment[["Liver"]] <- private$df_log[,selected_col[grepl("159|uniprotID|id",selected_col)]]
  self$experiment[["Colon"]] <- private$df_log[,selected_col[grepl("163|uniprotID|id",selected_col)]]

  }
)





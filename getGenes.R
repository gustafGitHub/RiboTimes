# Class getGenes

getGenes <- setClass("getGenes",
         representation(annotationFile = "character",
                        startPos = "integer",
                        endPos = "integer",
                        genes = "character",
                        direction = "logical",
                        nrGenes = "integer"),
         prototype(annotationFile = "/Users/gustafullman/Documents/Data/annotation/all_gene_sequences.txt",
                   startPos = integer(0),
                   endPos = integer(0),
                   genes = "ATG",
                   direction = TRUE,
                   nrGenes = integer(0))
         )

setGeneric(name = "loadGenes",
           def = function(theObject){
             standardGeneric("loadGenes")
           }
)
setMethod(f = "loadGenes",
          signature = "getGenes",
          definition = function(theObject){
            exeptionFlagPKset = FALSE
            geneTable = read.table(theObject@annotationFile, 
                                   stringsAsFactors = F,header=F,sep="\t")
            
            if(exeptionFlagPKset){
              nrGenesTemp = length(geneTable[1]$V1)
              genesTemp <- geneTable[1]$V1
              theObject@nrGenes = as.integer(nrGenesTemp/2)
              theObject@genes = genesTemp[1:nrGenesTemp %% 2 == 0]
              characterList = strsplit(geneTable[[1]][1:nrGenesTemp %% 2 != 0],split = ":|-|\\(|\\)")
              
              theObject@startPos = strtoi(do.call("rbind", 
                                                  lapply(characterList, "[[", 5)))
              theObject@endPos = strtoi(do.call("rbind", 
                                                lapply(characterList, "[[", 6)))
              theObject@direction = grepl("\\+", 
                                          do.call("rbind", lapply(characterList, "[[", 7)))
              #browser()
              
            }
            else{
              #browser()
              theObject@nrGenes = length(geneTable[2]$V2)
              theObject@genes <- geneTable[2]$V2
              characterList = strsplit(geneTable[[1]],split = ":|-|\\(|\\)")
              
              theObject@startPos = strtoi(do.call("rbind", 
                                                  lapply(characterList, "[[", 2)))
              theObject@endPos = strtoi(do.call("rbind", 
                                                lapply(characterList, "[[", 3)))
              theObject@direction = grepl("\\+", 
                                          do.call("rbind", lapply(characterList, "[[", 4)))
            }
            

            return(theObject)
          }
)
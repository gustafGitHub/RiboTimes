# loadCounts
source("getGenes.r")
source("codonToIndex.r")

loadCountsClassInstance <- setClass("loadCountsClass",
                             representation("geneName" = "character",
                                            "geneNameCounts" = "character",
                                            "geneSequence" = "character",
                                            "codons" = "list",
                                            "codonIndex" = "list",
                                            "codonIndexCount" = "list",
                                            "fileNameGenes" = "character",
                                            "fileNamesCounts" = "list",
                                            "startPos" = "integer",
                                            "startPosTest" = "integer",
                                            "endPos" = "integer",
                                            "counts" = "list",
                                            "nucleotideIndex" = "list",
                                            "nucleotides" = "list",
                                            "SequencesLoaded" = "logical",
                                            "geneDataExists" = "logical",
                                            "FA_concentrations" = "numeric",
                                            "countsPerCodon" = "list",
                                            "countsPerNucleotide" = "list",
                                            "timePerCodonTot" = "list",
                                            "timePerCodonPeptTrans" = "list",
                                            "timePerCodonTransLoc" = "list",
                                            "highestExpressedGenes" = "integer",
                                            "elongationTimes" = "numeric",
                                            "numberCodons" = "numeric",
                                            "numberGenes" = "numeric",
                                            "codonOffset" = "numeric"),

                             prototype(
                               geneName = "Initiation",
                               geneNameCounts = "Initiation",
                               geneSequence = "Initiation",
                               codons = list(),
                               codonIndex = list(),
                               codonIndexCount = list(),
                               fileNameGenes = paste0(getwd(),"/data/1_Genes_CDS.bed"),
                               fileNamesCounts = list(),
                               startPos = integer(0),
                               startPosTest = integer(0),
                               endPos = integer(0),
                               counts = list(),
                               nucleotideIndex = list(),
                               nucleotides = list(),
                               SequencesLoaded = FALSE,
                               geneDataExists = FALSE,
                               FA_concentrations = c(0,0.2,0.5,1),
                               countsPerCodon = list(),
                               countsPerNucleotide = list(),
                               timePerCodonTot = list(),
                               timePerCodonPeptTrans = list(),
                               timePerCodonTransLoc = list(),
                               highestExpressedGenes = integer(0),
                               numberCodons = 64,
                               elongationTimes = numeric(64),
                               numberGenes = 0,
                               codonOffset = 0)
                             )

setGeneric(name = "getSequenceAll",
           def = function(theObject){
             standardGeneric("getSequenceAll")
           }
)

setMethod(f = "getSequenceAll",
          signature = "loadCountsClass",
          definition = function(theObject){
            genesTable = read.table(theObject@fileNameGenes,stringsAsFactors = F,header=F,sep="\t")
            geneName = genesTable[[4]]
            theObject@geneName = geneName
            getGenes1 <- getGenes()
            loadGenes1 <- loadGenes(getGenes1)
            genes = loadGenes1@genes
            lengthGenes = length(genes)
            theObject@geneSequence = genes
            theObject@startPos = genesTable[[2]]
            theObject@endPos = genesTable[[3]]
            theObject@geneDataExists = rep(FALSE, lengthGenes)
            lengthGenes2 = length(genesTable[[2]])
            codons = list()
            nucleotides = list()
            codonIndex = list()
            codonIndexCount = list()
            for (i in 1:lengthGenes2){
              # if(i == 267){browser()}
              ii = which((loadGenes1@startPos == theObject@startPos[i]) & (loadGenes1@endPos == theObject@endPos[i]))
              startNucleotide = nchar(genes[ii]) %% 3
              gene_i = substring(genes[ii], seq(1 + startNucleotide, nchar(genes[ii]), 3), seq(1 + startNucleotide, nchar(genes[ii]), 3) + 2)
              nucl_i = substring(genes[ii], seq(1 + startNucleotide, nchar(genes[ii]), 1), seq(1 + startNucleotide, nchar(genes[ii]), 1))
              codons[[i]] = gene_i
              nucleotides[[i]] = nucl_i
              # browser()
              codonIndex_i = codonToIndex(gene_i)
              codonIndex[[i]] = codonIndex_i
              codonIndexCount[[i]] = codonsIndexCount(codonIndex_i)
            }
            theObject@codons = codons
            theObject@codonIndex = codonIndex
            theObject@codonIndexCount = codonIndexCount
            theObject@SequencesLoaded = TRUE
            theObject@nucleotides = nucleotides
            return(theObject)
            
          }
)

setGeneric(name = "setupLoadCounts",
           def = function(theObject){
             standardGeneric("setupLoadCounts")
           }
)

setMethod(f = "setupLoadCounts",
          signature = "loadCountsClass",
          definition = function(theObject){
            dirRoot = getwd()

            fname = 'calibrated_counts.bed'
            
            FA_conc = 0
            fileNamesCounts = list()
            
            fileNamesCounts[[1]] = paste0(dirRoot, "/data/", fname)
            
            theObject@fileNamesCounts = fileNamesCounts
            theObject@FA_concentrations = FA_conc
            
            return(theObject)
          }
)

setGeneric(name = "loadFiles",
           def = function(theObject){
             standardGeneric("loadFiles")
           }
)

setMethod(f = "loadFiles",
          signature = "loadCountsClass",
          definition = function(theObject){
            if (!theObject@SequencesLoaded) {stop("You must load Sequence first")}
            
            countsAll = list()
            countsAllNuc = list()
            countsPerCodon = list()
            countsPerNucleotide = list()
            
            for(i in 1:length(theObject@fileNamesCounts)){
            #  browser()
            cat(theObject@fileNamesCounts[[i]], '\n')
            reads=read.table(theObject@fileNamesCounts[[i]],stringsAsFactors=F,header=F,sep="\t")
            reads=reads[order(reads[,1],reads[,2]),] # sort again to be sure
            maxCol=ncol(reads)
            genes_u=unique(reads[,4])
            for(k in 1:length(genes_u)){
              #cat(k,'\n')
              idx=which(reads[,4]==genes_u[k])
              # reverse if on minus strand
              if(reads[idx[1],6]=='-'){
                idx=rev(idx)
              }
              counts = reads[idx,maxCol]
              lengthCounts = length(counts)
              countBias = lengthCounts%%3
              
              countsPerCodon[[k]] = counts[seq(1+countBias,lengthCounts,3)] + 
                counts[seq(2+countBias,lengthCounts,3)] + counts[seq(3+countBias,lengthCounts,3)]
              
              countsPerNucleotide[[k]] = counts
            }
            countsAll[[i]] = countsPerCodon
            countsAllNuc[[i]] = countsPerNucleotide
            
            }
            
            theObject@countsPerCodon = countsAll
            theObject@countsPerNucleotide = countsAllNuc
            theObject@numberGenes = length(genes_u)
            
            return(theObject)
          }
)

setGeneric(name = "loadCounts",
           def = function(theObject){
             standardGeneric("loadCounts")
           })


setMethod(f = "loadCounts",
          signature = "loadCountsClass",
          definition = function(theObject){
            theObject<-getSequenceAll(theObject)
            theObject<-loadFiles(theObject)
            return(theObject)
            
          })

setGeneric(name = "getHighestExpressedGenes",
           def = function(theObject){
             standardGeneric("getHighestExpressedGenes")
           })


setMethod(f = "getHighestExpressedGenes",
          signature = "loadCountsClass",
          definition = function(theObject){
            C = theObject
            numberGenes = C@numberGenes
            numberCounts = numeric(numberGenes)
            geneIndex = 1:numberGenes
            for(j in 1:numberGenes){
              numberCounts[j] = sum(C@countsPerCodon[[1]][[j]])/length(C@countsPerCodon[[1]][[j]])
            }
            highestExpressed <- sort.int(numberCounts, index.return=TRUE, decreasing = TRUE)
            theObject@highestExpressedGenes = highestExpressed$ix
            return(theObject)
            
          })


loadAll<-function(){
  
  C<-loadCountsClassInstance()
  C<-setupLoadCounts(C)
  C<-loadCounts(C)
  C<-getHighestExpressedGenes(C)
  return(C)
  
}

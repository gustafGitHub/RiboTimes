# loadCountsCsv
source("codonToIndex.R")

loadCountsTableinstance <- setClass("loadCountsTable", representation("codonName" = "list", "codonIndex" = "list", "geneName" = "character", "countsPerCodon" = "list"),
         prototype(codonName = list(), codonIndex = list(), geneName = "lpp", countsPerCodon = list()))


loadCountsTable <- function(fileName){
  
countsTable <- read.table(file = fileName, sep = "\t", stringsAsFactors = F)
CodonIndexes <- read.table("./data/CodonIndexes.txt")
codonNameOrder = substr(CodonIndexes[[2]], 1,3)
codonName <- countsTable[,3]
codonIndex <- codonToIndexVar(countsTable[,3], codonNameOrder)
geneName <- countsTable[,2]
countsPerCodon <- countsTable[,4]

codonName_list = list()
codonIndex_list = list()
countsPerCodon_list = list()

geneUnique <- unique(geneName)

for(i in 1:length(geneUnique)){
  
  gene_i <- which(geneUnique[i] == geneName)
  codonName_list[[i]] = codonName[gene_i]
  codonIndex_list[[i]] = codonIndex[gene_i]
  countsPerCodon_list[[i]] = countsPerCodon[gene_i]
}

CT <- new("loadCountsTable")

CT@codonName = codonName_list
CT@codonIndex = codonIndex_list
CT@geneName = geneUnique
CT@countsPerCodon = countsPerCodon_list

return(CT)

}


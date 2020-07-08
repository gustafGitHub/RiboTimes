#library("xlsx")

nGenes = C@numberGenes

geneName = character(0)
nucleotideNames = character(0)
counts22 = numeric(0)
counts23 = numeric(0)
counts24 = numeric(0)
counts25 = numeric(0)
counts26 = numeric(0)
counts27 = numeric(0)
counts28 = numeric(0)


saveNgeneFragments <- function(N, C){
  for(ii in 1:N){
    
    i = C@highestExpressedGenes[ii]
      
    geneName1<-gsub(".*\\.","",C@geneName[i])
    nucleotideNames_t = C@nucleotides[[i]]
    
    counts22_t = C@countsPerNucleotide[[1]][[i]]
    counts23_t = C@countsPerNucleotide[[2]][[i]]
    counts24_t = C@countsPerNucleotide[[3]][[i]]
    counts25_t = C@countsPerNucleotide[[4]][[i]]
    counts26_t = C@countsPerNucleotide[[5]][[i]]
    counts27_t = C@countsPerNucleotide[[6]][[i]]
    counts28_t = C@countsPerNucleotide[[7]][[i]]
    
    geneName_t = rep(geneName1, length(counts22_t))
    
    if(length(counts22_t) == length(nucleotideNames_t)){
      
      geneName = c(geneName, geneName_t)
      nucleotideNames = c(nucleotideNames, nucleotideNames_t)
      
      counts22 = c(counts22, counts22_t)
      counts23 = c(counts23, counts23_t)
      counts24 = c(counts24, counts24_t)
      counts25 = c(counts25, counts25_t)
      counts26 = c(counts26, counts26_t)
      counts27 = c(counts27, counts27_t)
      counts28 = c(counts28, counts28_t)
      
      }
    }
  
  dataFrame = data.frame(geneName, nucleotideNames, counts22, 
                         counts23, counts24, counts25,
                         counts26, counts27, counts28)
  
  write.csv(dataFrame, file = paste0("data/",N ,"_geneFragments_FA0"))
  
}

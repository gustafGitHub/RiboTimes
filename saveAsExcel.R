# library("xlsx")

nGenes = C@numberGenes

geneName = character(0)
codonNames = character(0)
counts = numeric(0)

saveNgenes <- function(N, C){
  for(ii in 1:N){
    
    i = C@highestExpressedGenes[ii]
      
      geneName1<-gsub(".*\\.","",C@geneName[i])
      codonNames_t = C@codons[[i]]
      
      counts_t = C@countsPerCodon[[1]][[i]]
      
      geneName_t = rep(geneName1, length(counts_t))
      geneName = c(geneName, geneName_t)
      codonNames = c(codonNames, codonNames_t)
      
      counts = c(counts, counts_t)

  }
  
  dataFrame = data.frame(geneName, codonNames, counts)
  write.csv(dataFrame, file = paste0("data/",N ,"_genes_FA0"))
  
}

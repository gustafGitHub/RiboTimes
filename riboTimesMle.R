#RiboTimesMleClassInstance <- setClass("RiboTimesMleClass",
#                                    representation(),
#                                    prototype())
source("loadCountsTable.R")
fileName = "./data/Old_FA0_Set1.txt"
CT<-loadCountsTable(fileName)

Mle <- function(CT){
  
  
  
  
  U_exp = 0
  z_k = 0
  g_mod = 0
  
  ngenes = length(CT@geneName)
  
  codon_context = 15
  posAsite = 8
  number_codons = 64
  
  Omega_exp = matrix(nrow = number_codons, ncol = codon_context)
  
  
  for(i in 1:ngenes){
    
    # Calculate U_exp
    U_exp[i] = sum(CT@countsPerCodon[[i]])
    length_gene = length(CT@countsPerCodon[[i]])
    counts_gene = CT@countsPerCodon[[i]]
    
    for(j in posAsite:(length_gene - posAsite)){
      # Calculate Omega_exp
      j_relAsite = j - posAsite + 1 # See figure 2 in manuscript
      counts_j = counts_gene[j]
      Omega_exp[CT@codonIndex[kA], j_relAsite] = Omega_exp[CT@codonIndex[kA], j_relAsite] + counts_j
    }
    
    
  }
  
  
  
  Mle_out <- new("RiboTimesMleClass")
  
  return(Mle_out)
  
}
library(Rcpp)
sourceCpp("./C/RiboTimesMLE.cpp")
source("loadCountsTable.R")

InputFile = "/Users/gustafullman/Documents/src/RiboTimes/data/ZI_30000_FA0.txt"
PathToOutput = "/Users/gustafullman/Documents/src/RiboTimes/data/"
fileName = "/Users/gustafullman/Documents/src/RiboTimes/data/ZI_30000_FA0_table.txt"

CountsObj <- function(){
  CountsObj<-loadCountsTable(fileName)
  return(CountsObj)
}

MleInitiate <- function(CountsObj){
  
  U_exp = 0
  z_k = 0
  g_mod = 0
  
  ngenes = length(CountsObj$geneName)
  
  pNumber = 15
  posAsite = 8
  cNumber = 64
  rowLenght = pNumber * cNumber
  
  Omega_exp = matrix(nrow = pNumber, ncol = cNumber, data = 0)
  n_ipc = array(dim = c(ngenes, pNumber, cNumber), data = 0)
  #n_ipc2 = array(dim = c(ngenes, pNumber, cNumber), data = 0)
  delta_c_seq_i = list()
  
  #tGeneElong = matrix(nrow = pNumber, ncol = ngenes, data = 0)
  muCoff = matrix(nrow = pNumber, ncol = cNumber, data = 0)
  zHat = matrix(nrow = pNumber, ncol = cNumber, data = 0)
  mu_initial = 0
  mu_new = 0
  
  
  for(i in 1:ngenes){
    
    # Calculate U_exp
    U_exp[i] = 0
    length_gene = length(CountsObj$countsPerCodon[[i]])
    counts_gene = CountsObj$countsPerCodon[[i]]
    delta_c_seq_i = matrix(nrow = length(counts_gene), ncol = pNumber * cNumber)
    codonIndex = CountsObj$codonIndex[[i]]
    tGeneElong1 = length_gene - pNumber
    
    n_pc = matrix(nrow = pNumber, ncol = cNumber, data = 0)
    
    for(j in 1:(length_gene - pNumber)){
      # Calculate Omega_exp
      # Calculate delta_{c, seq_i(p + i - pA)
      j_relAsite = j + posAsite - 1
      
      counts_j = counts_gene[j_relAsite]
      U_exp[i] = U_exp[i] + counts_j
      
      for(k in 1:pNumber){
        
        p = j_relAsite + k - posAsite
        
        codonIndex_p = codonIndex[p]
        
        Omega_exp[k, codonIndex_p] = Omega_exp[k, codonIndex_p] + counts_j
        n_pc[k, codonIndex_p] = n_pc[k, codonIndex_p] + 1

      }
      
    }
    
    mu_initial = mu_initial + U_exp[i]/tGeneElong1
    n_ipc[i,,] = n_pc
    
  }

  
  maxIter = 10
  
  deltaZ = 1e-22
  
  z_next = matrix(nrow = pNumber, ncol = cNumber, data = 1)
  
  for(n_iter in 1:maxIter){
    
    z_previous = z_next
    
    for(p in 1:pNumber){
      for(cc in 1:cNumber){
        
        sum_i = 0
        for(i in 1:ngenes){
          
          sum_c = 0
          n_ipc_temp = n_ipc[i, p,]
          for(c in 1:cNumber){
            
            sum_c = sum_c + z_previous[p, c] * n_ipc_temp[c]
            
          }
          
          sum_i = sum_i + U_exp[i]*n_ipc[i, p, cc] / sum_c
        }
        z_next[p, cc] = Omega_exp[p, cc]/sum_i
      }
      
    }
    
    # Normalize
    
    
    
    for(p in 1:pNumber){
      mu_new[p] = 0
      for(i in 1:ngenes){
        length_gene = length(CountsObj$countsPerCodon[[i]])
        codonIndex = CountsObj$codonIndex[[i]]
        tGeneElong = 0
        for(j in 1:(length_gene - pNumber)){
          j_relAsite = j + posAsite - 1
          codonIndex_p = codonIndex[j_relAsite + p - posAsite]
          tGeneElong = tGeneElong + z_next[p, codonIndex_p]
        }
        mu_new[p] = mu_new[p] + U_exp[i] /  tGeneElong
      }
      for(c in 1:cNumber){
        
        coffB = mu_new[p] / mu_initial
        z_next[p, c] = z_next[p, c] * coffB
        if(is.nan(z_next[p, c])){z_next[p, c] = 0}
      }
    }
    
    diff = sum((z_next - z_previous)^2)/rowLenght
    if(diff < deltaZ){break}
    
  }
  
  mRow = as.integer(0)
  for(i in 1:length(CountsObj$countsPerCodon)){mRow = mRow + length(CountsObj$countsPerCodon[[i]])}
  
  #MleObj <- new("MleClass")
  MleObj <- list()
  
  MleObj$Omega = Omega_exp
  MleObj$U = U_exp
  MleObj$cNumber = as.integer(cNumber)
  MleObj$pNumber = as.integer(pNumber)
  MleObj$ngenes = as.integer(ngenes)
  MleObj$z_initial = z_next
  MleObj$codonIndex = C$codonIndex
  MleObj$countsPerCodon = C$countsPerCodon
  MleObj$mRow = mRow
  
  return(MleObj)
  
}

MleAlgorithm <- function(){
  
  outputList1 = runMLE(InputFile, PathToOutput)
  
  outputList$geneCodeTable = outputList1$geneCodeTable[2:length(outputList1$geneCodeTable)]
  outputList$geneCodeTableOrdered = outputList1$geneCodeTableOrdered[2:length(outputList1$geneCodeTableOrdered)]
  outputList$rCodonSeq = outputList1$rCodonSeq[2:length(outputList1$rCodonSeq)]
  outputList$geneStart = outputList1$geneStart[2:length(outputList1$geneStart)]
  outputList$geneEnd = outputList1$geneEnd[2:length(outputList1$geneEnd)]
  outputList$geneCodons = outputList1$geneCodons[2:length(outputList1$geneCodons)]
  outputList$geneRPFtotal = outputList1$geneRPFtotal[2:length(outputList1$geneRPFtotal)]
  outputList$geneRPFdensity = outputList1$geneRPFdensity[2:length(outputList1$geneRPFdensity)]
  outputList$gene_Elong_AA = outputList1$gene_Elong_AA[2:length(outputList1$gene_Elong_AA)]
  outputList$geneRPF_Elong = outputList1$geneRPF_Elong[2:length(outputList1$geneRPF_Elong)]
  
  return(outputList1)
  
}



UnitTests <- function(MleObj){
  
  
}

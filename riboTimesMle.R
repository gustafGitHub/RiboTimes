

CountsObj <- function(){
  source("loadCountsTable.R")
  fileName = "./data/Counts.txt"
  CountsObj<-loadCountsTable(fileName)
  return(CountsObj)
}

MleInitiate <- function(CountsObj){
  
  MleClassInstance <- setClass("MleClass", 
                   
                   representation("Omega" = "matrix", "U" = "numeric"), 
                   
                   prototype("Omega" = matrix(0), "U" = 0
                             ), 
                   
                   contains = "loadCountsTable") 
  
  U_exp = 0
  z_k = 0
  g_mod = 0
  
  ngenes = length(CountsObj@geneName)
  
  pNumber = 15
  posAsite = 8
  ncodons = 64
  
  Omega_exp = matrix(nrow = pNumber, ncol = ncodons, data = 0)
  delta_c_seq_i = list()
  
  
  for(i in 1:ngenes){
    
    # Calculate U_exp
    U_exp[i] = 0
    length_gene = length(CountsObj@countsPerCodon[[i]])
    counts_gene = CountsObj@countsPerCodon[[i]]
    delta_c_seq_i = matrix(nrow = length(counts_gene), ncol = pNumber*ncodons)
    codonIndex = CountsObj@codonIndex[[i]]
    
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
        
      }
    }
    
  }
  
  MleObj <- new("MleClass")
  
  MleObj@Omega = Omega_exp
  MleObj@U = U_exp
  
  return(MleObj)
  
}

initialGuess <- function(CountsObj, U_exp, Omega_exp, ncodons, codon_context){
 
  
  # Initial guess 
  z_previous = matrix(nrow = ncodons, ncol = codon_context, data = 1)

  for(p in 1:codon_context){
    for(c in 1:ncodons){
      
      sum_i_expr = 0
      
      nipc = CountsObj@CountsPerCodon
      
      for(i in 1:ngenes){
        
        nipc = CountsObj@CountsPerCodon[[i]]
        
        c_sum = 0
        
        for(c2 in 1:codon_context){
          
          
          
        }
        
        sum_i_expr_temp = U_exp[i]/(c_sum*nipc)
        sum_i_expr = sum_i_expr + sum_i_expr_temp
        
      }
      
      z_next = 1
      
    }
  }
  
}


calculateHessian <- function(){
  
  
  
}


MleAlgorithm <- function(){
  
  
  
}


UnitTests <- function(MleObj){
  
  Omega_ref <- read.table("./data/Omega.txt")
  Omega = MleObj@Omega
  
  test1 = sum(Omega - Omega_ref) == 0
  
  U_ref = read.table("./data/U.txt")
  U = MleObj@U
  
  test2 = sum(U - U_ref) == 0
  
  testaAll = c(test1, test2)
  
  return(testaAll)
  
}
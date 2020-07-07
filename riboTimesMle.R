#RiboTimesMleClassInstance <- setClass("RiboTimesMleClass",
#                                    representation(),
#                                    prototype())
source("loadCountsTable.R")
fileName = "./data/Old_FA0_Set1.txt"
CT<-loadCountsTable(fileName)

Mle <- function(CT){
  
  Omega_exp = 0
  
  # Calculate U_exp
  U_exp = 0
  z_k = 0
  g_mod = 0
  
  
  
  Mle_out <- new("RiboTimesMleClass")
  
  return(Mle_out)
  
}
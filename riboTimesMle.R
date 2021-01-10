library(Rcpp)
sourceCpp("./C/RiboTimesMLE.cpp")
#source("loadCountsTable.R")

InputFile = "/Users/gustafullman/Documents/src/RiboTimes/data/ZI_30000_FA0.txt"
PathToOutput = "/Users/gustafullman/Documents/src/RiboTimes/data/"
fileName = "/Users/gustafullman/Documents/src/RiboTimes/data/ZI_30000_FA0_table.txt"

CountsObj <- function(){
  CountsObj<-loadCountsTable(fileName)
  return(CountsObj)
}

MleAlgorithm <- function(InputFile, PathToOutput){
  
  outputList = list()
  
  outputList1 = runMLE(InputFile, PathToOutput)
  
  outputList$dataSetName = outputList1$dataSetName
  outputList$dataSubSetName = outputList1$dataSubSetName
  outputList$dataSubSet_RPF_Total = outputList1$dataSubSet_RPF_Total
  outputList$doublingTime = outputList1$doublingTime
  outputList$pAsite = outputList1$pAsite
  outputList$pNumber = outputList1$pNumber
  outputList$cNumber = outputList1$cNumber
  outputList$jTot = outputList1$jTot
  outputList$global_Time_Factor = outputList1$global_Time_Factor
  
  outputList$geneCodeTable = outputList1$geneCodeTable[2:length(outputList1$geneCodeTable)]
  outputList$geneCodeTableOrdered = outputList1$geneCodeTableOrdered[2:length(outputList1$geneCodeTableOrdered)]
  outputList$rCodonSeq = outputList1$rCodonSeq[2:length(outputList1$rCodonSeq)]
  outputList$geneName = outputList1$geneName[2:length(outputList1$geneName)]
  outputList$geneStart = outputList1$geneStart[2:length(outputList1$geneStart)]
  outputList$geneEnd = outputList1$geneEnd[2:length(outputList1$geneEnd)]
  outputList$geneCodons = outputList1$geneCodons[2:length(outputList1$geneCodons)]
  outputList$geneRPFtotal = outputList1$geneRPFtotal[2:length(outputList1$geneRPFtotal)]
  outputList$geneRPFdensity = outputList1$geneRPFdensity[2:length(outputList1$geneRPFdensity)]
  outputList$gene_Elong_AA = outputList1$gene_Elong_AA[2:length(outputList1$gene_Elong_AA)]
  outputList$gene_Ci_Exper = outputList1$gene_Ci_Exper[[1]][2:length(outputList1$gene_Ci_Exper[[1]])]
  
  outputList$gene_di_Exper = outputList1$gene_di_Exper[[1]][2:length(outputList1$gene_di_Exper[[1]])]
  outputList$gene_fi_Model = outputList1$gene_fi_Model[[1]][2:length(outputList1$gene_fi_Model[[1]])]
  outputList$gene_Gi_Model = outputList1$gene_Gi_Model[[1]][2:length(outputList1$gene_Gi_Model[[1]])]
  outputList$gene_Gi_Model_Time = outputList1$gene_Gi_Model_Time[[1]][2:length(outputList1$gene_Gi_Model_Time[[1]])]
  outputList$gene_Ti_Model_Time_Abs = outputList1$gene_Ti_Model_Time_Abs[[1]][2:length(outputList1$gene_Ti_Model_Time_Abs[[1]])]
  
  outputList$rpfOmegaExper = outputList1$rpfOmegaExper[2:length(outputList1$rpfOmegaExper)]

  outputList$iORFcodons = outputList1$iORFcodons[2:length(outputList1$iORFcodons)]
  outputList$indCodonOrder = outputList1$indCodonOrder[2:length(outputList1$indCodonOrder)]
  outputList$nRPF = outputList1$nRPF[[1]][2:length(outputList1$nRPF[[1]])]
  
  outputList$sij_Exper = outputList1$sij_Exper[[1]][2:length(outputList1$sij_Exper[[1]])]
  outputList$sij_Model = outputList1$sij_Model[[1]][2:length(outputList1$sij_Model[[1]])]
  outputList$sij_Model_Sigma = outputList1$sij_Model_Sigma[[1]][2:length(outputList1$sij_Model_Sigma[[1]])]
  outputList$sij_Model_Time = outputList1$sij_Model_Time[[1]][2:length(outputList1$sij_Model_Time[[1]])]
  outputList$sij_Model_Time_Sigma = outputList1$sij_Model_Time_Sigma[[1]][2:length(outputList1$sij_Model_Time_Sigma[[1]])]
  
  outputList$gij_Model = outputList1$gij_Model[[1]][2:length(outputList1$gij_Model[[1]])]
  outputList$gij_Model_Sigma = outputList1$gij_Model_Sigma[[1]][2:length(outputList1$gij_Model_Sigma[[1]])]
  
  outputList$gij_Model_Time = outputList1$gij_Model_Time[[1]][2:length(outputList1$gij_Model_Time[[1]])]
  outputList$gij_Model_Time_Sigma = outputList1$gij_Model_Time_Sigma[[1]][2:length(outputList1$gij_Model_Time_Sigma[[1]])]
  
  outputList$tij_Model_Time_Abs = outputList1$tij_Model_Time_Abs[[1]][2:length(outputList1$tij_Model_Time_Abs[[1]])]
  outputList$tij_Model_Time_Abs_Sigma = outputList1$tij_Model_Time_Abs_Sigma[[1]][2:length(outputList1$tij_Model_Time_Abs_Sigma[[1]])]
  
  #outputList$strModel_Info = outputList1$strModel_Info
  
  #outputList$rDwellTime = outputList1$rDwellTime
  
  outputList$nCodCount = outputList1$nCodCount[2:length(outputList1$nCodCount)]
  outputList$codUsage = outputList1$codUsage[2:length(outputList1$codUsage)]
  outputList$codRPFsAver = outputList1$codRPFsAver[2:length(outputList1$codRPFsAver)]
  outputList$codRPFsSigma = outputList1$codRPFsSigma[2:length(outputList1$codRPFsSigma)]
  
  vFP_Full = outputList1$zFP[2:length(outputList1$zFP)]
  vFP_Sigma_Full = outputList1$zFP_Sigma[2:length(outputList1$zFP_Sigma)]
  
  vFP_Full2 = outputList1$tML_long[2:length(outputList1$zFP)]
  vFP_Sigma_Full2 = outputList1$tML_Sigma_long[2:length(outputList1$zFP_Sigma)]
  
  zFP = matrix(0, nrow = 15, ncol = 64)
  zFP_Sigma = matrix(0, nrow = 15, ncol = 64)
  tML_long = matrix(0, nrow = 15, ncol = 64)
  tML_Sigma_long = matrix(0, nrow = 15, ncol = 64)
  
  k = 0
  for(iPos in 1:15){
    for (jCod in 1:64) {
      k = k + 1
      #print(k)
      zFP[iPos, jCod] = vFP_Full[k]
      tML_long[iPos, jCod] = vFP_Full2[k]
      zFP_Sigma[iPos, jCod] = vFP_Sigma_Full[k]
      tML_Sigma_long[iPos, jCod] = vFP_Sigma_Full2[k]
    }
  }
  
  outputList$zFP = zFP
  outputList$zFP_Sigma = zFP_Sigma
  outputList$tML_long = tML_long
  outputList$tML_Sigma_long = tML_Sigma_long
  
  return(outputList)
  
}

UnitTests <- function(MleObj){
  
  
}

library(Rcpp)
sourceCpp("./C/RiboTimesMLE.cpp")

InputFile = "./data/Demo_input_file.txt"
PathToOutput = "./output/"

if(!file.exists(PathToOutput)){
  dir.create(PathToOutput)
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
  
  outputList$nCodCount = outputList1$nCodCount[2:length(outputList1$nCodCount)]
  outputList$codUsage = outputList1$codUsage[2:length(outputList1$codUsage)]
  outputList$codRPFsAver = outputList1$codRPFsAver[2:length(outputList1$codRPFsAver)]
  outputList$codRPFsSigma = outputList1$codRPFsSigma[2:length(outputList1$codRPFsSigma)]
  
  zFP_temp = outputList1$zFP
  zFP_Sigma_temp = outputList1$zFP_Sigma
  zFP_HAT_temp = outputList1$zFP_HAT
  zFP_HAT_Sigma_temp = outputList1$zFP_HAT_Sigma
  
  pNumber = outputList$pNumber
  cNumber = outputList$cNumber
  
  zFP = matrix(0, nrow = pNumber, ncol = cNumber)
  zFP_Sigma = matrix(0, nrow = pNumber, ncol = cNumber)
  zFP_HAT = matrix(0, nrow = pNumber, ncol = cNumber)
  zFP_HAT_Sigma = matrix(0, nrow = pNumber, ncol = cNumber)
  
  for(i in 1:pNumber){
    zFP[i,] = zFP_temp[[i + 1]][2:length(zFP_temp[[i + 1]])]
    zFP_Sigma[i,] = zFP_Sigma_temp[[i + 1]][2:length(zFP_Sigma_temp[[i + 1]])]
    zFP_HAT[i,] = zFP_HAT_temp[[i + 1]][2:length(zFP_HAT_temp[[i + 1]])]
    zFP_HAT_Sigma[i,] = zFP_HAT_Sigma_temp[[i + 1]][2:length(zFP_HAT_Sigma_temp[[i + 1]])]
  }
  
  colnames(zFP) <- outputList$geneCodeTableOrdered
  colnames(zFP_Sigma) <- outputList$geneCodeTableOrdered
  colnames(zFP_HAT) <- outputList$geneCodeTableOrdered
  colnames(zFP_HAT_Sigma) <- outputList$geneCodeTableOrdered

  outputList$zFP = zFP
  outputList$zFP_Sigma = zFP_Sigma
  outputList$zFP_HAT = zFP_HAT
  outputList$zFP_HAT_Sigma = zFP_HAT_Sigma
  
  return(outputList)
  
}

UnitTests <- function(MleObj){
  
  
}

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List Hes_Pos_ML_Refine_zFactors(List MleObject) {
  
  List output;
  
  int pNumber = MleObject["pNumber"]; // Number of positions
  
  int ngenes = MleObject["ngenes"];
  int jTot = ngenes; // Number of genes
  
  int pAsite = 8; // Position of A-site
  int cNumber = 64; //Number of codons
  int mRow = MleObject["mRow"];
  // double rowLength = pNumber * cNumber;
  // int jGeneStartShift = 0;
  int jEndAcc, jEnd;
  int nCodon_Elong, nRPF_kA;
  double geneRPF_Elong, t_Gene_Elong;
  
  NumericMatrix zFP = MleObject["z_initial"];
  List geneCounts = MleObject["countsPerCodon"];
  List codonIndex = MleObject["codonIndex"];
  
  int kA = 0;
  double t_pAsite;
  int indCodon;
  
  NumericMatrix zFP_Old(pNumber, cNumber);
  NumericMatrix ijH_A(pNumber, cNumber);
  NumericVector codonIndexGene;
  NumericVector counts;
  NumericVector tModel(mRow);
  
  NumericVector OutTemp1(200);
  NumericVector OutTemp2(200);
  NumericVector OutTemp3(200);
  
    // Refines z-Factors using Position-Diagonal Hessian Approximation  by Marquardt-Like approach

    //Convert z-factors to v-factors by normalyzing each position by z(p,1)
    
  int i = 0;
  for(int iPos = 0; iPos < pNumber; iPos++){
    double z1 = zFP(iPos, 0);
    for(int indCodon = 0; indCodon < cNumber; indCodon++){
    zFP(iPos, indCodon) = zFP(iPos, indCodon) / z1;
    zFP_Old(iPos, indCodon) = zFP(iPos, indCodon);
    ijH_A(iPos, indCodon) = 0;
      if(zFP(iPos, indCodon) > 0){ //the position/codon  is active
        i = i + 1;
        ijH_A(iPos, indCodon) = i; //Transforms Pos/Codon into active Vector Position
      }
    }
  }
  
  //output["zFP"] = zFP;
  //output["ijH_A"] = ijH_A;
  
  //return output;
  
  /*  
  int nActive = i;
  //Prepare the arrays
  int nActiveShort = nActive - pNumber; //The dimention of active subspace
  NumericVector vDir(nActive); */
  
    //Active movement spase
    
  
  
  //Calculate model gene times from v-factors
  //Get Initial  model times from zFactors
  jEndAcc = 0;
  jEnd = 0;
  
  for(int j = 0; j < jTot; j++){
    counts = geneCounts[j];
    codonIndexGene = codonIndex[j];
    jEndAcc = jEndAcc + jEnd;
    jEnd = counts.size();
    for(int k = 0; k < (jEnd - pNumber); k++){
      t_pAsite = 1; //time with a particular codon in the A-site
      for(int i = 0; i < pNumber; i++){
        indCodon = codonIndexGene(i + k);
        t_pAsite = t_pAsite * zFP(i, indCodon - 1);
      }
      kA = pAsite + k + jEndAcc; //A-site in the original data set
      tModel(kA) = t_pAsite;
    }
  }
      
  NumericMatrix rpfOmegaExper(pNumber, cNumber);
  NumericMatrix rpfOmegaModel(pNumber, cNumber);
  NumericMatrix nCodPos(pNumber, cNumber);
  NumericMatrix tFetha(pNumber, cNumber);
  NumericVector uGeneElong(ngenes);
  NumericVector nCodGeneElong(ngenes);
      
      //Prepare the main Omega(p,c) array of RPF fingerprint and additional arrays
      
  double muInitial = 0;
  double muStandard = 0;
  jEndAcc = 0;
  jEnd = 0;
  
  for(int j = 0; j < jTot; j++){
    counts = geneCounts[j];
    codonIndexGene = codonIndex[j];
    jEndAcc = jEndAcc + jEnd;
    jEnd = counts.size();
    nCodon_Elong = 0;
    geneRPF_Elong = 0;
    t_Gene_Elong = 0;
    for(int k = 0; k < (jEnd - pNumber); k++){
      kA = pAsite + k + jEndAcc; //Count from the A-site
      nCodon_Elong = nCodon_Elong + 1;
      nRPF_kA = counts[k + pAsite - 1];
      geneRPF_Elong = geneRPF_Elong + nRPF_kA;
      t_Gene_Elong = t_Gene_Elong + tModel(kA);
      for(int iPos = 0; iPos < pNumber; iPos++){
        indCodon = indCodon = codonIndexGene(iPos + k);
        rpfOmegaExper(iPos, indCodon - 1) = rpfOmegaExper(iPos, indCodon - 1) + nRPF_kA;
        nCodPos(iPos, indCodon - 1) = nCodPos(iPos, indCodon - 1) + 1;
      }
    }
    uGeneElong(j) = geneRPF_Elong;  //RPFs in the inner gene region
    nCodGeneElong(j) = nCodon_Elong; //Codon Number in the inner gene region
    muInitial = muInitial + geneRPF_Elong / t_Gene_Elong;
    muStandard = muStandard + geneRPF_Elong / nCodon_Elong;
  }
  
  output["U"] = uGeneElong;
  output["Omega"] = rpfOmegaExper;
  
  //-----------------------------
  //main iterations starts =========================================================
  
  double alfa = 3;  //Marquard parameter
  double beta = 0;  //Marquard parameter
  double stepSize_Coff_Iter = 0.8;
  int iTerNumber = 19;
  
  // Iterations Start
    for(int Iter = 0; Iter < iTerNumber; Iter++){
      
      
    }
  
  
  return output;

}


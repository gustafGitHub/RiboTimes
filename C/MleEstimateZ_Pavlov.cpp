#include <Rcpp.h>
//#include "RcppArmadillo.h"
#include <vector>
#include <math.h>
#include <cmath>
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
  double rowLength = pNumber * cNumber;
  // int jGeneStartShift = 0;
  int jEndAcc, jEnd;
  int nCodon_Elong, nRPF_kA;
  double geneRPF_Elong, t_Gene_Elong;
  
  NumericMatrix zFP = MleObject["z_initial"];
  List geneCounts = MleObject["countsPerCodon"];
  List codonIndex = MleObject["codonIndex"];
  //List LogLikelihoodOut, LogLikelihoodList;
  
  int kA = 0;
  double t_pAsite;
  int indCodon, indCod;
  
  NumericVector zFP_Refined_long(rowLength);
  NumericVector zFP_Sigma_Refined_long(rowLength);
  NumericVector wFP_Refined_long(rowLength);
  
  NumericMatrix zFP_Old(pNumber, cNumber);
  NumericMatrix ijH_A(pNumber, cNumber);
  NumericVector codonIndexGene;
  NumericVector counts;
  NumericVector tModel(mRow);
  NumericVector rpfModel(mRow);

  
  NumericMatrix rpfOmegaExper(pNumber, cNumber);
  NumericMatrix rpfOmegaModel(pNumber, cNumber);
  NumericMatrix nCodPos(pNumber, cNumber);
  NumericMatrix tFetha(pNumber, cNumber);
  NumericMatrix gradFull(pNumber, cNumber);
  NumericVector uGeneElong(ngenes);
  NumericVector tGeneElong(ngenes);
  NumericVector nCodGeneElong(ngenes);
  //arma::Cube<double> fyHessian(pNumber, cNumber, cNumber);
  double fyHessian[pNumber][cNumber][cNumber];
  double mH_Full[pNumber][cNumber][cNumber];
  double mH_Full_Modified[pNumber][cNumber][cNumber];
  
  NumericMatrix mH_Full_Modified_iPos1(cNumber, cNumber);
  
  NumericMatrix fyHessian_temp(cNumber,cNumber);
  NumericMatrix zShift(pNumber,cNumber);
  NumericVector indCod1(pNumber);
  NumericVector zShift_Long(rowLength);
  
  NumericMatrix diag_mH_Full(pNumber,cNumber);
  
  NumericMatrix muCoff(pNumber,cNumber);
  NumericMatrix nPosCodInGene(pNumber,cNumber);
  NumericMatrix diagHM1_Full(pNumber,cNumber);
  
  NumericMatrix mU;
  NumericMatrix mL;
  NumericVector vP;
  
  Function FactorLUPA("FactorLUPA");
  Function SolveLUPA("SolveLUPA");
  Function InvertLUPA("InvertLUPA");
  
  Environment base = Environment("package:base");
  Function readline = base["readline"];
  
  int iA, jA, nPos, j;

  
  double zFP_pc;
  
  double wtFetha;
  double coff_Fy;
  double tModel_kA;
  double muNew;
  
  int k;
  int i;
  int imax = 1000;
  double gradActiveNorm2;
  NumericVector iGradActive(imax);
  NumericVector ijH_Active(imax);
  NumericVector gradActive(imax);
  NumericVector zPosAverage(pNumber);
  
  double gradActiveNorm;
  double vDir_Norm, gradActive_Norm, vDir_LogLkh_Derivative_Min, cos_vDirOld_GradNew;
  double vDir_LogLkh_Derivative, cos_vDir_Grad, stepSize;
  double tLhd_Total_Old, zDiff2, zPosCodDiff, muWeight, zAverTot;
  
  double vDirGradProjection;
  double vDir_Norm2;
  double gradActive_Norm2;
  double tLhd_delta = 0;
  double tLhd_Total = 0;
  double zPosAver, wCodPos, coffB, coffZ;
  
  
  // Refines z-Factors using Position-Diagonal Hessian Approximation  by Marquardt-Like approach

  //Convert z-factors to v-factors by normalyzing each position by z(p,1)
    
  i = 0;
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
  
  
  int nActive = i;
  Rcout << "nActive: " << nActive << "\n";
  //Prepare the arrays
  int nActiveShort = nActive - pNumber; //The dimention of active subspace
  NumericVector vDir(nActive);
  NumericVector grad_PS(nActive);
  NumericMatrix mH_PS(nActive, nActive);
  NumericMatrix mH_PS_Pos1(nActive, nActive);
  NumericVector vDir_Short(nActiveShort);
  NumericVector vDir_PS;
  NumericVector mHM1;
  NumericVector gradActive_Short(nActiveShort);
  
  List FactorLupaOut;
  List SolveLupaOut;
  
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
      kA = pAsite + k + jEndAcc - 1; //A-site in the original data set
      tModel(kA) = t_pAsite;
    }
  }
      

      
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
      kA = pAsite + k + jEndAcc - 1; //Count from the A-site
      nCodon_Elong = nCodon_Elong + 1;
      nRPF_kA = counts[k + pAsite - 1];
      geneRPF_Elong = geneRPF_Elong + nRPF_kA;
      t_Gene_Elong = t_Gene_Elong + tModel(kA);
      for(int iPos = 0; iPos < pNumber; iPos++){
        indCodon = codonIndexGene(iPos + k);
        rpfOmegaExper(iPos, indCodon - 1) = rpfOmegaExper(iPos, indCodon - 1) + nRPF_kA;
        nCodPos(iPos, indCodon - 1) = nCodPos(iPos, indCodon - 1) + 1;
      }
    }
    uGeneElong(j) = geneRPF_Elong;  //RPFs in the inner gene region
    tGeneElong(j) = t_Gene_Elong;
    nCodGeneElong(j) = nCodon_Elong; //Codon Number in the inner gene region
    muInitial = muInitial + geneRPF_Elong / t_Gene_Elong;
    muStandard = muStandard + geneRPF_Elong / nCodon_Elong;
  }
  
  output["U"] = uGeneElong;
  output["Omega"] = rpfOmegaExper;
  
  //Function Get_Log_Likelihood("Get_Log_Likelihood");
  //LogLikelihoodOut = Get_Log_Likelihood(LogLikelihoodList);
  
  //-----------------------------
  //main iterations starts =========================================================
  
  double alfa = 3;  //Marquard parameter
  double beta = 0;  //Marquard parameter
  //double stepSize_Coff_Iter = 0.8;
  int iTerNumber = 1; //19;
  
  // Iterations Start
  for(int Iter = 0; Iter < iTerNumber; Iter++){
    
    //Initite the above arrays
    for(int iPos = 0; iPos < pNumber; iPos++){
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
        tFetha(iPos, indCodon) = 0;
        rpfOmegaModel(iPos, indCodon) = 0;
        gradFull(iPos, indCodon) = 0;
        for(int iCod = 0; iCod < cNumber; iCod++){
          fyHessian[iPos][indCodon][iCod] = 0;
        }
      }
    }
    
    jEndAcc = 0;
    jEnd = 0;
    
    //Get "tFetha(iPos, indCodon)", "rpfOmegaModel(iPos, indCodon)"
    // and "fyHessian(iPos, indCodon, iCod)"
    muNew = 0;
    for(int jGene = 0; jGene < jTot; jGene++){
    //jStart = .geneStart(jGene) + jGeneStartShift:
    counts = geneCounts[jGene];
    codonIndexGene = codonIndex[jGene];
    jEndAcc = jEndAcc + jEnd;
    jEnd = counts.size();
    wtFetha = uGeneElong(jGene) / tGeneElong(jGene); //U(j)/T(j)
    coff_Fy = wtFetha / tGeneElong(jGene); //U(j)/(T(j)*T(j))
    muNew = muNew + wtFetha;
    
    // Get "tFetha" and model counts "rpfModel" for gene j
    for(int k = 0; k < (jEnd - pNumber); k++){
      kA = pAsite + k + jEndAcc - 1; //Count from the A-site
      //Rcout << "kA" << kA << "\n";
      tModel_kA = tModel(kA);
      for(int iPos = 0; iPos < pNumber; iPos++){
        indCodon = codonIndexGene(iPos + k); //.iORFcodons(iPos + k - 1):
        tFetha(iPos, indCodon - 1) = tFetha(iPos, indCodon - 1) + tModel_kA;
      }
      rpfModel(kA) = wtFetha * tModel_kA; //Expected model RPFs
    }
        
    // get model RPF  count fingerprint rpfOmega and fyHessian
    for(int iPos = 0; iPos < pNumber; iPos++){
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
        rpfOmegaModel(iPos, indCodon) = rpfOmegaModel(iPos, indCodon) + 
          wtFetha * tFetha(iPos, indCodon);
        //Rcout << "wtFetha: " << wtFetha << " tFetha: " << tFetha(iPos, indCodon) << "\n";
        
        for(int iCod = 0; iCod < cNumber; iCod++){
          fyHessian[iPos][indCodon][iCod] = fyHessian[iPos][indCodon][iCod] + 
            coff_Fy * tFetha(iPos, indCodon) * tFetha(iPos, iCod);
          if(iPos == 0){
            fyHessian_temp(indCodon, iCod) = fyHessian_temp(indCodon, iCod) + 
              coff_Fy * tFetha(iPos, indCodon) * tFetha(iPos, iCod);
          }
        }
      }
    }
      
    //re-initiate tFetha for the next gene
    for(int iPos = 0; iPos < pNumber; iPos++){
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
        tFetha(iPos, indCodon) = 0;
      }
    }
  }
  
  //Get Full and Active Gradient Vectors
  k = 0;
  i = 0;
  gradActiveNorm2 = 0;
  
  for(int iPos = 0; iPos < pNumber; iPos++){
  indCod1(iPos) = i + 1;
    for(int indCodon = 0; indCodon < cNumber; indCodon++){
      zShift_Long(k) = 0;
      k++; 
      zShift(iPos, indCodon) = 0;
      gradFull(iPos, indCodon) = rpfOmegaExper(iPos, indCodon) - rpfOmegaModel(iPos, indCodon);
      if(zFP(iPos, indCodon) > 0){
        gradFull(iPos, indCodon) = gradFull(iPos, indCodon) / zFP(iPos, indCodon);
      }
      i = ijH_A(iPos, indCodon);
      if(i > 0){ //it is active vector position
        iGradActive(i - 1) = k;
        ijH_Active(i - 1) = i;
        gradActive(i - 1) = 0;
        if(zFP(iPos, indCodon) > 0){
          gradActive(i - 1) = gradFull(iPos, indCodon);
          gradActiveNorm2 = gradActiveNorm2 + gradActive(i) * gradActive(i);
        }
      }
    }
  }
  
  /*output["gradFull"] = gradFull;
  output["iGradActive"] = iGradActive;
  output["ijH_Active"] = ijH_Active;
  output["gradActive"] = gradActive;
  output["gradActiveNorm2"] = gradActiveNorm2;
  output["rowLength"] = rowLength;*/
  
    
  gradActiveNorm2 = gradActiveNorm2 / nActive;
  gradActiveNorm = sqrt(gradActiveNorm2);
      
  //Check the Projection of the New Gradient on the Old seach directin vDir:
  //It Should be close to Zero
  vDirGradProjection = 0;
  vDir_Norm2 = 0;
  gradActive_Norm2 = 0;
  
  for(int iA = 0; iA < nActive; iA++){
    vDirGradProjection = vDirGradProjection + vDir(iA) * gradActive(iA);
    vDir_Norm2 = vDir_Norm2 + vDir(iA) * vDir(iA);
    gradActive_Norm2 = gradActive_Norm2 + gradActive(iA) * gradActive(iA);
  }
  vDir_Norm = sqrt(vDir_Norm2);
  gradActive_Norm = sqrt(gradActive_Norm2);
    
  if(vDir_Norm > 0){
    vDir_LogLkh_Derivative_Min = vDirGradProjection / vDir_Norm;
  }
  else{
    vDir_LogLkh_Derivative_Min = 0;
  }
    
  cos_vDirOld_GradNew = vDir_LogLkh_Derivative_Min / gradActive_Norm;
  
  //Construct  Full Block-Hessian Approximation
  for(int iPos = 0; iPos < pNumber; iPos++){
    for(int iC = 0; iC < cNumber; iC++){
      zFP_pc = zFP(iPos, iC);
      if(zFP_pc > 0){
        for(int iCod = (iC + 1); iCod < cNumber; iCod++){
          mH_Full[iPos][iC][iCod] = 0;
          if(zFP(iPos, iCod) > 0){
              mH_Full[iPos][iC][iCod] = 
            fyHessian[iPos][iC][iCod] / (zFP_pc * zFP(iPos, iCod));
          }
          mH_Full_Modified[iPos][iC][iCod] = mH_Full[iPos][iC][iCod];
          mH_Full[iPos][iCod][iC] = mH_Full[iPos][iC][iCod];
          mH_Full_Modified[iPos][iCod][iC] = mH_Full[iPos][iC][iCod];
        }
        mH_Full[iPos][iC][iC] = 
        (fyHessian[iPos][iC][iC] - rpfOmegaExper(iPos, iC)) / (zFP_pc * zFP_pc);
        diag_mH_Full(iPos, iC) = mH_Full[iPos][iC][iC];
        mH_Full_Modified[iPos][iC][iC] = mH_Full[iPos][iC][iC];
      }
    }
  }
    
      
    //Increase the Block-Hessian diagonal entries to tilt the shift direction towards the gradient
    for(int iPos = 0; iPos < pNumber; iPos++){
      for(int iC = 0; iC < cNumber; iC++){
    mH_Full_Modified[iPos][iC][iC] = mH_Full[iPos][iC][iC] 
    + alfa * (diag_mH_Full(iPos, iC) - 1) - beta * gradActiveNorm;
      }
    }
    
    //Make most frequent codon (indexed as  1) inactive at  all positions of the local context
    for(int iPos = 0; iPos < pNumber; iPos++){
      i = ijH_A(iPos, 1);
      ijH_Active(i - 1) = 0;
    }
      
    //Accelerated solution using the block-diagonal structure of Hessian
    k = 0;
    
    for(int iPos = 0; iPos < pNumber; iPos++){
      //Get Hessian and Gradient for each p-position
        //Get Gradient as one-dimension vector from pos/codon matrix
        iA = 0;
        for(int iCod = 1; iCod < cNumber; iCod++){
          i = ijH_A(iPos, iCod);
          if(i > 0){
            //Rcout << "grad_PS, iA: " << iA << "\n";
            grad_PS(iA) = gradFull(iPos, iCod);
            iA = iA + 1;
          }
        }
        nPos = iA;
        
        //Re-number Hessian accrding to grad_PS
        iA = 0;
        for(int iC = 1; iC < cNumber; iC++){
          i = ijH_A(iPos, iC);
          if(i > 0){
            jA = 0;
            for(int iCod = 1; iCod < cNumber; iCod++){
              j = ijH_A(iPos, iCod);
              if(j > 0){
                mH_PS(iA, jA) = mH_Full_Modified[iPos][iC][iCod];
                if(iPos == 0){
                  mH_PS_Pos1(iA, jA) = mH_PS(iA, jA);
                  if(jA == 0){
                    Rcout << "mH_Full_Modified: " << mH_Full_Modified[iPos][iC][iCod] << "\n";
                  }
                }
                jA = jA + 1;
              }
            }
            iA = iA + 1;
          }
        }
        
    //Solve Equations for the current context position iPos:
    //Solve  (mH_PS+(alfa-1)*diag_mH_PS)*vDir=gradActive

    //FactorLUPA(1, mH_PS, vP, mL, mU); 
    
    FactorLupaOut = FactorLUPA(1, mH_PS); //LUP factorise mH_PS
    //Call Print_LUP_Factors(mH, vP, mL, mU)
    // SolveLUPA(mL, mU, vP, vDir_PS, grad_PS); //Solve to find the direction
    
    vDir_PS = FactorLupaOut["vDir_PS"];
    NumericMatrix mL = FactorLupaOut["mL"];
    NumericMatrix mU = FactorLupaOut["mU"];
    vP = FactorLupaOut["vP"];
    
    grad_PS = SolveLUPA(mL, mU, vP, vDir_PS); //Solve to find the direction
    
    // restore gradActive short and vDir_Short
      //Change  vDir_Short into the oposite, assent direction in the full active space
    for(int i = 1; i < nPos; i++){
      k++;
      vDir_Short(k) = -vDir_PS(i);
      gradActive_Short(k) = grad_PS(i);
    }
    }
      
      //Restore vDir and zShift_Long in the full total space
      iA = 0;
    for(int i = 0; i < nActive; i++){
      vDir(i) = 0;
      if(ijH_Active(i) > 0){
        iA = iA + 1;
        vDir(i) = vDir_Short(iA);
      }
      k = iGradActive(i);
      zShift_Long(k) = vDir(i);
    }
      
      // Check vDir direction legibility
      vDirGradProjection = 0;
      vDir_Norm2 = 0;
    gradActive_Norm2 = 0;
    for(int iA = 0; iA < nActiveShort; iA++){
      vDirGradProjection = vDirGradProjection + vDir_Short(iA) * gradActive_Short(iA);
      vDir_Norm2 = vDir_Norm2 + vDir_Short(iA) * vDir_Short(iA);
      gradActive_Norm2 = gradActive_Norm2 + gradActive_Short(iA) * gradActive_Short(iA);
    }
      vDir_Norm = sqrt(vDir_Norm2);
      gradActive_Norm = sqrt(gradActive_Norm2);
      
      vDir_LogLkh_Derivative = vDirGradProjection / vDir_Norm;
      cos_vDir_Grad = vDir_LogLkh_Derivative / gradActive_Norm;
      
    if(cos_vDir_Grad >= 0){
      stepSize = 0.5 * tLhd_delta / vDir_LogLkh_Derivative;
      if(std::abs(stepSize) > 0.8){stepSize = 0.8;}
      
      if(std::abs(vDirGradProjection) < 1000){
        alfa = alfa * 0.1;
        beta = beta * 0.5;
        stepSize = 0.8;
      }
    }
    else{ //vDir direction is bad, use gradActive direction instead
      vDirGradProjection = 0;
      vDir_Norm2 = 0;
    for(int iA = 0; iA < nActiveShort; iA++){
      vDir_Short(iA) = gradActive_Short(iA);
      vDirGradProjection = vDirGradProjection + vDir_Short(iA) * gradActive_Short(iA);
      vDir_Norm2 = vDir_Norm2 + vDir_Short(iA) * vDir_Short(iA);
    }
    vDir_Norm = sqrt(vDir_Norm2);
    vDir_LogLkh_Derivative = vDirGradProjection / vDir_Norm;
    stepSize = 0.05 * tLhd_Total / vDirGradProjection; //Try to increase LL by 5%
      //stepSize = stepSize;
      
      //Restore vDir and zShift_Long
      iA = 0;
    for(int i = 0; i < nActive; i++){
      vDir(i) = 0;
      if(ijH_Active(i) > 0){
      iA = iA + 1;
      vDir(i) = vDir_Short(iA);
      }
      k = iGradActive(i);
      zShift_Long(k) = vDir(i);
    }
    }
      
      //record zFP_Old and transform  vector zShift_long into matrix zShift
      k = 0;
      for(int iPos = 0; iPos < pNumber; iPos++){
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
      zFP_Old(iPos, indCodon) = zFP(iPos, indCodon);
      k = k + 1;
    zShift(iPos, indCodon) = zShift_Long(k);
      }
      }
      
      //Refine stepSize by minimizing log-likelihood along vDir
      //Call RefineStepSize(jSet, newRPFdata, zFP, tModel, tGeneElong, _
      //                      tLhd_Total, stepSize, zShift_Long)
      
      //Do zFP shift
      for(int iPos = 0; iPos < pNumber; iPos++){
        for(int indCodon = 0; indCodon < cNumber; indCodon++){
        zFP(iPos, indCodon) = zFP_Old(iPos, indCodon) + stepSize * zShift(iPos, indCodon);
        }
      }
      
      //Get new Log-Likelihood function, new model times and gene times
      //Call Get_Log_Likelihood_Only(jSet, newRPFdata, zFP, tModel, tGeneElong, tLhd_Total)
      
      tLhd_Total_Old = tLhd_Total;
      
      zDiff2 = 0;
    for(int iPos = 0; iPos < pNumber; iPos++){
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
        zPosCodDiff = (zFP(iPos, indCodon) - zFP_Old(iPos, indCodon));
        zDiff2 = zDiff2 + zPosCodDiff * zPosCodDiff;
      }
    }
    zDiff2 = zDiff2 / rowLength;
      
      
    
    if(stepSize == 0){break;}
  }

    //
    //Main iteration ends ==================================================
    //-----------------------------
    //
    
    //Finalize z-factors
      
    //Zero muCoff=Effective weights of Codons
      for(int iPos = 0; iPos < pNumber; iPos++){
        for(int indCod = 0; indCod < cNumber; indCod++){
          muCoff(iPos, indCod) = 0;
        }
      }
      
      muWeight = 0;
      muNew = 0;
      muStandard = 0;
      for(int jGene = 0; jGene < jTot; jGene++){
        wtFetha = uGeneElong(jGene) / tGeneElong(jGene); //U(j)/T(j)
        muWeight = muWeight + nCodGeneElong(jGene) * wtFetha;
        muNew = muNew + wtFetha;
        muStandard = muStandard + uGeneElong(jGene) / nCodGeneElong(jGene);
        //jStart = .geneStart(jGene) + jGeneStartShift;
        jEnd = 0; //.geneEnd(jGene);
        
        //Zero nPosCodInGene
        for(int iPos = 0; iPos < pNumber; iPos++){
          for(int indCod = 0; indCod < cNumber; indCod++){
            nPosCodInGene(iPos, indCod) = 0;
          }
        }
        
        //Obtain nPosCodInGene
        for(int k = 0; k < jEnd - pNumber; k++){
          for(int iPos = 0; iPos < pNumber; iPos++){
          indCod = 0; //.iORFcodons(iPos + k - 1);
          nPosCodInGene(iPos, indCod) = nPosCodInGene(iPos, indCod) + 1;
          }
        }
        
        //Use nPosCodInGene to get muCoff=Effective weights of Codons
        for(int iPos = 0; iPos < pNumber; iPos++){
          for(int indCod = 0; indCod < cNumber; indCod++){
          muCoff(iPos, indCod) = muCoff(iPos, indCod) +
          nPosCodInGene(iPos, indCod) * wtFetha;
          }
        }
      }
      
      //Finalize  muCoff=Effective weights of nucleotides
      for(int iPos = 0; iPos < pNumber; iPos++){
        for(int indCod = 0; indCod < cNumber; indCod++){
          muCoff(iPos, indCod) = 1000 * muCoff(iPos, indCod) / muWeight;
        }
      }
      
      //Get position Hessian and Gradient
      k = 0;
    zAverTot = 1;
    //get the dimentions of the position Hessian
      for(int iPos = 0; iPos < pNumber; iPos++){
      iA = 0;
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
        i = ijH_A(iPos, indCodon);
        if(i > 0){
          iA = iA + 1;
        }
      }
      nPos = iA - 1; //Check
      
      
      //Get Possition Hessian with one position excluded
      
      
      //prepare the diagonal
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
        diagHM1_Full(iPos, indCodon) = 0;
      }
      
      
      for(int kTest = 0; kTest < 1; kTest++){
        iA = 0;
        for(int indCodon = 0; indCodon < cNumber; indCodon++){
          if(indCodon != kTest){
            i = ijH_A(iPos, indCodon);
            if(i > 0){
              iA = iA + 1;
              jA = 0;
              for(int iCod = 0; iCod < cNumber; iCod++){
                if(iCod != kTest){
                j = ijH_A(iPos, iCod);
                if(j > 0){
                  jA = jA + 1;
                  mH_PS(iA, jA) = mH_Full[iPos][indCodon][iCod];
                }
              }
            }
          }
        }
      }
      
      
      //invert the curtailed position Hessian
      FactorLupaOut = FactorLUPA(1, mH_PS);
      NumericMatrix mL = FactorLupaOut["mL"];
      //Call Print_LUP_Factors(mH_PS, vP, mL, mU)
      //InvertLUPA(mL, mU, vP, mHM1);
      mHM1 = InvertLUPA(mL, mU, vP);
      
      
      iA = 0;
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
      if(indCodon != kTest){
      i = ijH_A(iPos, indCodon);
      if(i > 0){
      iA = iA + 1;
      diagHM1_Full(iPos, indCodon) = diagHM1_Full(iPos, indCodon) + mHM1(iA, iA);
      //iA = iA;
      }
      }
      }
      }
      
      for(int kTest = 0; kTest < 1; kTest++){
        diagHM1_Full(iPos, kTest) = 2 * diagHM1_Full(iPos, kTest);
      }
      
      
      //Get position averaged z-factors
      zPosAver = 0;
      wCodPos = 0;
      for(int indCod = 0; indCod < cNumber; indCod++){
        wCodPos = wCodPos + muCoff(iPos, indCod);
        zPosAver = zPosAver + zFP(iPos, indCod) * muCoff(iPos, indCod);
      }
      zPosAverage(iPos) = zPosAver / wCodPos;
      zAverTot = zAverTot * zPosAverage(iPos);
      }
      
      // Normalization Coffs
      coffB = muNew / muStandard;
      coffB = zAverTot * coffB;
      coffZ = exp(log(coffB) / pNumber);
      
      //Finalize the refined z-factors
      k = 0;
      for(int iPos = 0; iPos < pNumber; iPos++){
        zPosAver = zPosAverage(iPos);
        //coffZ = 1
        for(int indCod = 0; indCod < cNumber; indCod++){
          k = k + 1;
          zFP_Refined_long(k) = coffZ * zFP(iPos, indCod) / zPosAver;
          zFP_Sigma_Refined_long(k) = coffZ * sqrt(abs(diagHM1_Full(iPos, indCod))) / zPosAver;
          wFP_Refined_long(k) = muCoff(iPos, indCod);
        }
      }
    
     
    
  output["rpfModel"] = rpfModel;
  output["tGeneElong"] = tGeneElong;
  output["tModel"] = tModel;
  output["rpfOmegaModel"] = rpfOmegaModel;
  output["fyHessian"] = fyHessian_temp;
  output["tFetha"] = tFetha;
  output["diag_mH_Full"] = diag_mH_Full;
  output["mH_PS_Pos1"] = mH_PS_Pos1;
  output["ijH_A"] = ijH_A;
  //output["mH_Full_Modified_iPos1"] = mH_Full_Modified_iPos1;
  
  return output;

}

// [[Rcpp::export]]
NumericVector SolveLUPA(NumericMatrix mL, NumericMatrix mU, NumericVector vP, NumericVector vB){
  
  //===========================================================================
  //
  //Public Sub SolveLUPA(mL() As Double, mU() As Double, vP() As Double, _
  //                       vX() As Double, vB() As Double)
  //
  //Given equation A*x=b and P*A=L*U factorization
  //Solve the equation L*U*x=PM1*b
  //first  solve Lv=PM1*b by forward substitutions
  //
  //
  
  //List LupaOut;
  int iDim, jDim, j;
  double jP, vj, mUjj, zj;
  
  double zTol = 1E-09;
  
  iDim = mL.nrow();
  jDim = mL.ncol();

  NumericVector v(jDim);
  NumericVector vX(jDim);
  NumericVector z(jDim);
  NumericVector W(jDim);
  
  //
  // prepare intermediate vectors
  //
  for(int j = 0; j < jDim; j++){
    v(j) = 0;
    vX(j) = 0;
    z(j) = 0;
    jP = vP(j); //permutate vb components
    W(j) = vB(jP);
  }
  // get L*v=w solution by forward substitutions
  for(int j = 0; j < jDim; j++){
    vj = 0;
    if(abs(W(j)) > 0){vj = W(j) / mL(j, j);}
    v(j) = vj;
    if(abs(vj) > 0){
      for(int i = j; i < jDim; i++){
        W(i) = W(i) - vj * mL(i, j);
      }
    }
  }
  //then solve mU*z=v by backsubstitutions
  for(int jj = 0; jj < jDim; jj++){
    j = jDim - jj - 1;
    mUjj = mU(j, j);
    zj = 8;
    if(abs(mUjj) > zTol){zj = v(j) / mUjj;}
    z(j) = zj;
    vX(j) = zj;
    for(int i = 0; i < j; i++){
      v(i) = v(i) - zj * mU(i, j);
    }
  }

  return vX;
  
}

// [[Rcpp::export]]
List FactorLUPA(int iP, NumericMatrix mA){
  
  //Sub FactorLUPA(iP As Long, mA() As Double, vP() As Double, mL() As Double, mU() As Double)
  //
  //Factorise the matrix as A*P=L*U where L is low and U is upper triangular
  //And keep the row rearrangements in the A matrix as a permutation matrix
  // or rather a permiutation vector P (vP)
  // Note that the algorithm with permutations is described in Sprang and
  //
  // This is a row oriented algorithm which first permutate rows of A to find the pivot
  // It first calculates a new column of mL then a new row of R and Gauss eliminate  rows of A
  // It has an advantage that it does not stop if a pivot is zero but just skips
  // the elimination step
  //

  List LupaOut;
  
  int mRow = mA.nrow(); 
  int nCol = mA.ncol();
  int iRowMax, jR, j1,iFlag;
  double rRowMax, rRowMaxAbs, rRowAbs, rRow, mUjj, mLij;
  
  NumericMatrix mU(mRow, nCol);
  NumericMatrix mL(mRow, nCol);
  NumericVector vP(mRow);
  NumericVector mUjj_temp(nCol);

  // prepare mL and mU matrices
  for(int i = 0; i < mRow; i++){
    for(int j = 0; j < nCol; j++){
      mU(i, j) = mA(i, j);
    }
    vP(i) = i;
  }
  // main factorization cycle

  double epsTol = 1E-31;
  
  for(int j = 0; j < nCol; j++){
    // at step j find the largest pivot in column j below j
    iRowMax = j;
    rRowMax = mU(j, j);
    rRowMaxAbs = std::abs(rRowMax);
    for(int i = j; i < mRow; i++){
    rRowAbs = std::abs(mU(i, j));
      if(rRowMaxAbs < rRowAbs) { //a larger pivot is found
        rRowMaxAbs = rRowAbs;
        iRowMax = i;
      }
    }

    if(iP == 0){iRowMax = j;} //Run without pivoting
    if(iRowMax > j){ //Swap  rows *iRowMax* and *j* in mL+mR
    //record the row swaping in A in the vector vP
    jR = vP(j);
    vP(j) = vP(iRowMax);
    vP(iRowMax) = jR;
    //actually swap the rows j and iRowMax (both in L and in the remaining of A)
      for(int k = 0; k < nCol; k++){
        rRow = mU(j, k);
        mU(j, k) = mU(iRowMax, k);
        mU(iRowMax, k) = rRow;
        Rcout << "Here!" << "\n";
      }
    }
    //Find new column j of a low triangular mL matrix
    //The k element of column j mLkj contains elimination coefficients
    //for subtracting row j from from row k creating a zero subcolumn in
      //modified A

    mUjj = mU(j, j); //get the current pivot
    mUjj_temp(j) = mUjj;
    j1 = j + 1;
    if (std::abs(mUjj) > epsTol) {//run the elimination
      for(int i = j1; i < nCol; i++){
          mLij = mU(i, j) / mUjj; //the new component of the column mL
          mU(i, j) = mLij; //save elimination coefficients in a low part of mU
          //Rcout << "mU: " << mU(i, j) << " i: " << i << " j: " << j << "\n";
          //subtract row j multiplied by Lij from row i of modified A
          for(int k = j1; k < mRow; k++){
            mU(i, k) = mU(i, k) - mLij * mU(j, k);
          }
        }
      }
      else{ //put real zeroes to stress that mUjj=0 and skip the elimination step
        iFlag = 2;
        for(int i = j; i < nCol; i++){
          mU(i, j) = 0;
        }
      }
    }
    // separate mU into mL and mR matrices
    if(nCol <= mRow){
      for(int j = 0; j < nCol; j++){
        for(int i = j + 1; i < mRow; i++){
          mL(i, j) = mU(i, j);
          mU(i, j) = 0;
        }
        mL(j, j) = 1;
      }
    }
    else{
    mL(0, 0) = 0;
    for(int i = 1; i < mRow; i++){
      for(int j = 0; j < i - 1; j++){
      mL(i, j) = mU(i, j);
      mU(i, j) = 0;
      }
      for(int j = i; j < nCol; j++){
        mL(i, j) = 0;
      }
      mL(i, i) = 1;
    }
  }
  
  LupaOut["vP"] = vP;
  LupaOut["mL"] = mL; 
  LupaOut["mU"] = mU;
  

  return LupaOut;
  
}

List Get_Log_Likelihood(int jTot){

  List output;
  
/*Sub Get_Log_Likelihood(jSet As Long, newRPFdata As RPFdataSet, zFP() As Double, tModel() As Double, _
                         tGeneElong() As Double, tLhd_Total As Double, tLhd_Gene() As Double, _
                         tLhd_Gene_UB() As Double, tLhd_Total_UB As Double) */
  
  //Common Block
  
  /*Dim rRPFgeneAver As Double
  Dim jGene As Long, jAsiteCodon As Long, iPos As Long
  Dim tModel_kA As Double, t_pAsite As Double
  
  Dim cNumber As Long, pNumber As Long, pAsite As Long
  Dim jStart As Long, jEnd As Long, jGeneStartShift As Long
  Dim jTot As Long, mRow As Long, indCodon As Long
  Dim i As Long, j As Long, k As Long, kA As Long*/
  
  
  //Likelyhood block
  //Dim gene_RPF_ln_tModel As Double, gene_RPF_ln_RPF As Double
  
  //With newRPFdata
  //jTot = UBound(.geneEnd): 'number of genes in data set
  //mRow = .geneEnd(jTot): 'number of codons in fata set
  //ReDim tModel(mRow), rpfModel(mRow)
  
  //jGeneStartShift = .jGeneStartShift
  //Calculate z-Factors from time ML Fingerprints

  int pAsite = 15;
  //pNumber = .pNumber:
  int cNumber = 64; //Number of codons
  double t_pAsite;
  int kA, indCodon;
  NumericVector tModel;

  //Calculate model gene times from v-factors
  // calculate Initial  model times from zFactors
  for(int jGene = 0; jGene < jTot; jGene++){
    //jStart = .geneStart(jGene) + jGeneStartShift:
    //jEnd = .geneEnd(jGene)
    for(int k = 0; k < jEnd){
      t_pAsite = 1; //time barrie with a particular codon in the A-site
      for(int i = 0; i < pNumber; i++){
        indCodon = //.iORFcodons(i + k - 1)
        t_pAsite = t_pAsite * zFP(i, indCodon);
      }
      kA = pAsite + k - 1; //A-site in the original data set
      tModel(kA) = t_pAsite;
    }
  }
  
/*Dim nCodon_Elong As Long
  Dim geneRPF_Elong As Double
  Dim t_Gene_Elong As Double, gene_RPF_per_Codon As Double, gene_Time_per_Codon As Double
  Dim dRPF_kA As Double, nRPF_kA As Long*/
  
  
  tLhd_Total = 0: //total likelyhood function
  tLhd_Total_UB = 0: //Upper Bound of the total likelyhood function
  For jGene = 1 To jTot
  jStart = .geneStart(jGene) + jGeneStartShift:
  jEnd = .geneEnd(jGene)
  nCodon_Elong = 0:
  geneRPF_Elong = 0:
  t_Gene_Elong = 0:
  gene_RPF_ln_tModel = 0: 'Likelyhood calculation
  gene_RPF_ln_RPF = 0: 'Likelyhood Upper Bound Calculation
  
  For k = jStart To jEnd - pNumber
  kA = pAsite + k - 1: 'Count from the A-site
  nCodon_Elong = nCodon_Elong + 1:
  nRPF_kA = .nRPF(kA, jSet)
  
  dRPF_kA = nRPF_kA
  geneRPF_Elong = geneRPF_Elong + nRPF_kA:
  t_Gene_Elong = t_Gene_Elong + tModel(kA):
  If (tModel(kA) > 0) Then
  gene_RPF_ln_tModel = gene_RPF_ln_tModel + dRPF_kA * Log(tModel(kA)):
  Else
  i = i
  End If
  If (nRPF_kA > 0) Then
  gene_RPF_ln_RPF = gene_RPF_ln_RPF + dRPF_kA * Log(dRPF_kA):
  Else
  i = i
  End If
  Next k
  
  gene_Time_per_Codon = t_Gene_Elong / nCodon_Elong
  tLhd_Gene(jGene) = gene_RPF_ln_tModel - geneRPF_Elong * Log(gene_Time_per_Codon):
  tLhd_Total = tLhd_Total + tLhd_Gene(jGene):
  tLhd_Gene_UB(jGene) = 0
If (geneRPF_Elong > 0) Then
  gene_RPF_per_Codon = geneRPF_Elong / nCodon_Elong:
  tLhd_Gene_UB(jGene) = gene_RPF_ln_RPF - geneRPF_Elong * Log(gene_RPF_per_Codon):
  End If
  tLhd_Total_UB = tLhd_Total_UB + tLhd_Gene_UB(jGene):
  tGeneElong(jGene) = t_Gene_Elong:
  Next jGene
  End With
  End Sub*/

  return output;

}
  
// [[Rcpp::export]]
List RefineStepSize(NumericMatrix zFP, double tLhd_Total, double stepSize, NumericVector zShift_Long){
  
  List output;
  
//  '
//'=============================
  //Sub RefineStepSize(jSet As Long, newRPFdata As RPFdataSet, zFP() As Double, tModel() As Double, _
  //                     tGeneElong() As Double, tLhd_Total As Double, stepSize As Double, zShift_Long() As Double)
  
  
  //The Sub refines the stepSize by minimizing likelihood along the vDir_Short
  
  //Standard block

  
  int nStepMax = 19;
  double stepSizeMin = abs(stepSize) / 10;
  double tLhd_Total_Old = tLhd_Total;
  double yTol = stepSizeMin * 0.001;

  NumericVector stepSizeSequence(nStepMax); 
  NumericVector fValueSequence(nStepMax);
  
  stepSizeSequence(0) = 0;
  fValueSequence(0) = tLhd_Total;
  
  //do the shift to Get new z-factor and estimate the convergence
  

  
  int pNumber = 15;
  int cNumber = 64; //Number of codons
  double y1, y2,y3,v1,v2,v3;

  double d32, d21, d31, aCoff, bCoff, yI;
  double vI;
  double y31_Length, y32_Length, y21_Length;
  double v31_Diff, v23_Diff, v21_Diff;
  
  NumericMatrix zFP_New(pNumber, cNumber);
  NumericMatrix zShift(pNumber, cNumber);
  
  //Get zShift
  int k = 0;
  for(int iPos = 0; iPos < pNumber; iPos++){
    for(int indCodon = 0; indCodon < cNumber; indCodon++){
      k++;
      zShift(iPos, indCodon) = zShift_Long(k);
    }
  }
  
  //Do stepping
  double stepSizeAdd = stepSizeMin;
  int iStepLast;
  bool iFlag_Zero, iFlag;
  
  stepSize = 0;
  for(int iStepSize = 1; iStepSize < nStepMax; iStepSize++){
    iStepLast = iStepSize;
    //safeguard from negatives in the new zFP
    stepSize = stepSize + stepSizeAdd;
    iFlag_Zero = 0;
  for(int iA34 = 0; iA34 < 5; iA34++){
    iFlag = 0;
    for(int iPos = 0; iPos < pNumber; iPos++){
      for(int indCodon = 0; indCodon < cNumber; indCodon++){
      if((zFP(iPos, indCodon) + stepSize * zShift(iPos, indCodon)) < 0){
          iFlag = 1;
          iFlag_Zero = 1;
          break;
        }
      }
      
      if(iFlag == 1){break;}
    }
    
    
    if(iFlag == 1){ //reduce the step size
      stepSize = stepSize / 2;
    }
    if(iFlag_Zero == 1 && iFlag == 0){
    //stepSize = stepSize
      return output;
    }
    }
    
    //do a tryal shift
      
      for(int iPos = 0; iPos < pNumber; iPos++){
        for(int indCodon = 0; indCodon < cNumber; indCodon++){
          zFP_New(iPos, indCodon) = zFP(iPos, indCodon) + stepSize * zShift(iPos, indCodon);
        }
      }
      //k = k
      
      // Get new  value of the L-function
      // Call Get_Log_Likelihood_Only(jSet, newRPFdata, zFP_New, tModel, tGeneElong, tLhd_Total)
      double tLhd_Total_Old;
      iFlag = 0;
      if(tLhd_Total > tLhd_Total_Old){ //The step is legitimate; Record
        stepSizeSequence(iStepSize) = stepSize;
        fValueSequence(iStepSize) = tLhd_Total;
        tLhd_Total_Old = tLhd_Total;
      }
      else{
      stepSizeSequence(iStepSize + 1) = stepSize;
      fValueSequence(iStepSize + 1) = tLhd_Total;
      stepSize = stepSize - stepSizeAdd / 2;
      //try an intermediate step
        for(int iPos = 0; iPos < pNumber; iPos++){
          for(int indCodon = 0; indCodon  < cNumber; indCodon++){
            zFP_New(iPos, indCodon) = zFP(iPos, indCodon) + stepSize * zShift(iPos, indCodon);
          }
        }
        //Call Get_Log_Likelihood_Only(jSet, newRPFdata, zFP_New, tModel, tGeneElong, tLhd_Total)
        stepSizeSequence(iStepSize) = stepSize;
        fValueSequence(iStepSize) = tLhd_Total;
        
        tLhd_Total = tLhd_Total;
        if(fValueSequence(iStepSize) > fValueSequence(iStepSize - 1)){
        stepSize = stepSizeSequence(iStepSize);
        tLhd_Total = fValueSequence(iStepSize);
        //Use 3-point bracket, v2 is max
        y1 = stepSizeSequence(iStepSize - 1);
        y2 = stepSizeSequence(iStepSize);
        y3 = stepSizeSequence(iStepSize + 1);
        v1 = fValueSequence(iStepSize - 1);
        v2 = fValueSequence(iStepSize);
        v3 = fValueSequence(iStepSize + 1);
        }
        else{
        stepSize = stepSizeSequence(iStepSize - 1);
        tLhd_Total = fValueSequence(iStepSize - 1);
        //Use 3-point bracket, v2 is max
        y1 = stepSizeSequence(iStepSize - 2);
        y2 = stepSizeSequence(iStepSize - 1);
        y3 = stepSizeSequence(iStepSize);
        v1 = fValueSequence(iStepSize - 2);
        v2 = fValueSequence(iStepSize - 1);
        v3 = fValueSequence(iStepSize);
        
        }
        iFlag = 1;
      break;
    }
        
        stepSizeAdd = 2 * stepSizeAdd;
  }
  if(iFlag == 0){stepSize = stepSizeSequence(iStepLast);}
  if(iStepLast < 3){
    stepSize = 0;
    return output;
  }
  
  // Refine using quadratic interpolation
  for(int iDummy = 0; iDummy < 6; iDummy++){
    d32 = (v3 - v2) / (y3 - y2);
    d21 = (v2 - v1) / (y2 - y1);
    d31 = (v3 - v1) / (y3 - y1);
    aCoff = (d32 - d21) / (y3 - y1);
    bCoff = d31 - aCoff * (y3 + y1);
    yI = -bCoff / (2 * aCoff);
  
  //Get Log-Likelihood (=vI) with yI step size
  for(int iPos = 0; iPos < pNumber; iPos++){
    for(int indCodon = 0; indCodon < cNumber; indCodon++){
      zFP_New(iPos, indCodon) = zFP(iPos, indCodon) + yI * zShift(iPos, indCodon);
    }
  }
  //Call Get_Log_Likelihood_Only(jSet, newRPFdata, zFP_New, tModel, tGeneElong, vI)
  
  //Four major cases
  iFlag = 0;

  if(yI < y2 && vI < v2){
    y1 = yI; v1 = vI; iFlag = 1;
  } //replace (1) with (I)
    
  if(yI < y2 && vI > v2){  
    y3 = y2; v3 = v2; y2 = yI; v2 = vI; iFlag = 1;
  } //replace (3) with (2) and (2) with (I); (1) the same
      
  if (yI > y2 && vI > v2){ //replace (1) with (2) and (2) with (I); (3) the same
    y1 = y2; v1 = v2; y2 = yI; v2 = vI; iFlag = 1;
  }
        
  if(yI > y2 && vI < v2){ y3 = yI; v3 = vI; iFlag = 1;}//replace (1) with (I)
        
        if (iFlag == 0){break;}
        y31_Length = y3 - y1;
        y32_Length = y3 - y2;
        y21_Length = y2 - y1;
        v31_Diff = v3 - v1;
        v23_Diff = v2 - v3;
        v21_Diff = v2 - v1;
        if(y31_Length < yTol || y32_Length < yTol || y21_Length < yTol) { //we are at the bracket boundary
        //i = i
        break;
        }
    }
    stepSize = y2;
    tLhd_Total = v2;
    //v2 = v2
    return output;
  
}


// [[Rcpp::export]]            
List InvertLUPA(){
  
  List output;
                        
  //===========================================================================
  //
  //Public Sub InvertLUPA(mL() As Double, mU() As Double, vP() As Double, mAM1() As Double)
  //
  //Invert matrix A that has been LUP factorised by solving n-equations
  //Given equation A*AM1(k)=e(k) where AM1(k) is a k-column of A inverse
  //Solve the equation L*U*AM1(k)=PM1*e(k)
  //first  solve Lv=PM1*e(k) by forward substitutions
  //
  //
  /*Dim i As Long, j As Long, k As Long, jP As Long, kP As Long
  Dim iDim As Long, jDim As Long
  Dim z() As Double, v() As Double, W() As Double, e() As Double
  Dim vj As Double, zj As Double, mUjj As Double
  Dim zTol As Double*/
  double zTol = 1E-09;
          
  /* iDim = UBound(mL, 1)
            jDim = UBound(mL, 2)
            ReDim v(jDim)
            ReDim z(jDim)
            ReDim W(jDim)
            ReDim mAM1(jDim, jDim)
            
            '
          ' prepare intermediate vectors
            '
          For k = 1 To jDim
            For j = 1 To jDim:
            v(j) = 0:
            mAM1(j, k) = 0:
            z(j) = 0:
            W(j) = 0:
            jP = vP(j): 'permutate e(k) components
            If (jP = k) Then kP = k
            Next j
            W(kP) = 1:
            ' Solve L*v=W for v by forward substitutions
              For j = kP To jDim:
              vj = 0:
              If (Abs(W(j)) > 0) Then vj = W(j) / mL(j, j):
              v(j) = vj:
              If (Abs(vj) > 0) Then
              For i = j To jDim:
              W(i) = W(i) - vj * mL(i, j):
              Next i
              End If
              Next j
              'then solve mU*z=v by backsubstitutions
              For j = jDim To 1 Step -1
            mUjj = mU(j, j):
              zj = 8:
              If (Abs(mUjj) > zTol) Then zj = v(j) / mUjj:
              z(j) = zj:
              mAM1(j, k) = zj:
              For i = 1 To j:
              v(i) = v(i) - zj * mU(i, j):
              Next i
              Next j
              Next k
              End Sub */

  return output;

}

// [[Rcpp::export]]              
List Get_Log_Likelihood_Only(){
  
  List output;
                            
              /*Sub Get_Log_Likelihood_Only(jSet As Long, newRPFdata As RPFdataSet, zFP() As Double, tModel() As Double, _
                                            tGeneElong() As Double, tLhd_Total As Double)
              
              'Common Block
              
              Dim rRPFgeneAver As Double
              Dim jGene As Long, jAsiteCodon As Long, iPos As Long
              Dim tModel_kA As Double, t_pAsite As Double
              
              Dim cNumber As Long, pNumber As Long, pAsite As Long
              Dim jStart As Long, jEnd As Long, jGeneStartShift As Long
              Dim jTot As Long, mRow As Long, indCodon As Long
              Dim i As Long, j As Long, k As Long, kA As Long
              
              
              'Likelyhood block
              Dim gene_RPF_ln_tModel As Double, gene_RPF_ln_RPF As Double
              
              Dim nCodon_Elong As Long
              Dim geneRPF_Elong As Double
              Dim t_Gene_Elong As Double, gene_RPF_per_Codon As Double, gene_Time_per_Codon As Double
              Dim dRPF_kA As Double, nRPF_kA As Long
              
              With newRPFdata
              jTot = UBound(.geneEnd): 'number of genes in data set
              mRow = .geneEnd(jTot): 'number of codons in fata set
              ReDim tModel(mRow), rpfModel(mRow)
              
              jGeneStartShift = .jGeneStartShift
              '
            ' Calculate z-Factors from time ML Fingerprints
              '
            pAsite = .pAsite:
              pNumber = .pNumber:
              cNumber = 64: 'Number of codons
              '
            'Calculate model gene times from v-factors
              ' calculate Initial  model times from zFactors
              For jGene = 1 To jTot
              jStart = .geneStart(jGene) + jGeneStartShift:
              jEnd = .geneEnd(jGene)
              For k = jStart To jEnd - pNumber
              t_pAsite = 1: 'time barrie with a particular codon in the A-site
              For i = 1 To pNumber
              indCodon = .iORFcodons(i + k - 1):
              t_pAsite = t_pAsite * zFP(i, indCodon):
              Next i
              kA = pAsite + k - 1: 'A-site in the original data set
              tModel(kA) = t_pAsite:
              Next k
              Next jGene
              '
            tLhd_Total = 0: 'total likelyhood function
              For jGene = 1 To jTot
              jStart = .geneStart(jGene) + jGeneStartShift:
              jEnd = .geneEnd(jGene)
              nCodon_Elong = 0:
              geneRPF_Elong = 0:
              t_Gene_Elong = 0:
              gene_RPF_ln_tModel = 0: 'Likelyhood calculation
              gene_RPF_ln_RPF = 0: 'Likelyhood Upper Bound Calculation
              
              For k = jStart To jEnd - pNumber
              kA = pAsite + k - 1: 'Count from the A-site
              nCodon_Elong = nCodon_Elong + 1:
              nRPF_kA = .nRPF(kA, jSet)
              
              dRPF_kA = nRPF_kA
              geneRPF_Elong = geneRPF_Elong + nRPF_kA:
              t_Gene_Elong = t_Gene_Elong + tModel(kA):
              If (tModel(kA) > 0) Then
              gene_RPF_ln_tModel = gene_RPF_ln_tModel + dRPF_kA * Log(tModel(kA)):
              Else
              i = i
              End If
              If (nRPF_kA > 0) Then
              gene_RPF_ln_RPF = gene_RPF_ln_RPF + dRPF_kA * Log(dRPF_kA):
              Else
              i = i
              End If
              Next k
              
              gene_Time_per_Codon = t_Gene_Elong / nCodon_Elong
              tLhd_Total = tLhd_Total + gene_RPF_ln_tModel - geneRPF_Elong * Log(gene_Time_per_Codon):
              tGeneElong(jGene) = t_Gene_Elong:
              Next jGene
              End With
              End Sub*/
              
  return output;
}
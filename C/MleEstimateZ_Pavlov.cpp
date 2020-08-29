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
  int indCodon;
  
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
  
  Function FactorLUPA("FactorLUPA");
  Function SolveLUPA("SolveLUPA");
  
  Environment base = Environment("package:base");
  Function readline = base["readline"];
  
  int iA, jA, nPos, j;
  NumericVector grad_PS(60);
  NumericMatrix mH_PS(60, 60);
  NumericMatrix mH_PS_Pos1(60, 60);
  
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
  
  double gradActiveNorm;
  double vDir_Norm, gradActive_Norm, vDir_LogLkh_Derivative_Min, cos_vDirOld_GradNew;
  
  double vDirGradProjection;
  double vDir_Norm2;
  double gradActive_Norm2;
  
  
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
  //Prepare the arrays
  //int nActiveShort = nActive - pNumber; //The dimention of active subspace
  NumericVector vDir(nActive);
  
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
    //List FactorLupaOut, SolveLupaOut;
    
    //FactorLupaOut = FactorLUPA(mH_PS); //LUP factorise mH_PS
    //Call Print_LUP_Factors(mH, vP, mL, mU)
    // SolveLUPA(mL, mU, vP, vDir_PS, grad_PS); //Solve to find the direction
    
    //SolveLupaOut = SolveLUPA(grad_PS); //Solve to find the direction
      
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
List SolveLUPA(){
  
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
  
  List LupaOut;
  
  //double zTol = 1E-09;
  
  /*iDim = UBound(mL, 1)
  jDim = UBound(mL, 2)

  '
  ' prepare intermediate vectors
  '
  For j = 1 To jDim:
  v(j) = 0:
  vX(j) = 0:
  z(j) = 0:
  jP = vP(j): 'permutate vb components
  W(j) = vB(jP):
  Next j
  ' get L*v=w solution by forward substitutions
  For j = 1 To jDim:
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
  vX(j) = zj:
  For i = 1 To j:
  v(i) = v(i) - zj * mU(i, j):
  Next i
  Next j
  End Sub */

  return LupaOut;
  
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
  NumericMatrix mU_temp0(mRow, nCol);
  NumericMatrix mU_temp1(mRow, nCol);
  NumericMatrix mU_temp2(mRow, nCol);
  NumericMatrix mL(mRow, nCol);
  NumericVector vP(mRow);
  NumericVector mUjj_temp(nCol);

  // prepare mL and mU matrices
  for(int i = 0; i < mRow; i++){
    for(int j = 0; j < nCol; j++){
      mU(i, j) = mA(i, j);
      mU_temp0(i, j) = mU(i, j);
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
          mU_temp1(i, j) = mU(i, j);
          mU(i, j) = mLij; //save elimination coefficients in a low part of mU
          mU_temp2(i, j) = mLij;
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
  LupaOut["mU_temp0"] = mU_temp0;
  LupaOut["mU_temp1"] = mU_temp1;
  LupaOut["mU_temp2"] = mU_temp2;
  LupaOut["mUjj_temp"] = mUjj_temp;
  

  return LupaOut;
  
}
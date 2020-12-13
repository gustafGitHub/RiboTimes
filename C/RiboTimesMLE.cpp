#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <Rcpp.h>

using std::iostream;
using std::string;
using std::vector;
using std::fstream;
//[[Rcpp::plugins(cpp11)]]

using MatrixDouble=vector<vector <double> >;
using MatrixLong=vector<vector <long> >;
using D3_VectorDouble=vector<vector <vector<double> > >;
using D3_VectorLong=vector<vector <vector<long> > >;
using namespace Rcpp;

struct DS_Corr {
		string Name;
        double setA_Aver;
		double setA_Sigma;
        double setB_Aver;
		double setB_Sigma;
        double setAB_Corr_Raw;
		double setAB_Corr_Pearson;
        double setAB_R2;
	};

struct RPFdataSet {
    string dataSetName;
    vector<string> dataSubSetName;
    vector<double> dataSubSet_RPF_Total;
    vector<double> doublingTime;        	// Time factor extracted from doubling time
    long nFus;                           	// Nubmer of data sub-sets; FA conc used, for instance
    vector<double> concFA;              	// Conc of FA for a data subset
    long pAsite;                     		// Position of the A-site in the cocal context
    long pNumber;                     		// The length of the local contaxt in codons
    long jGeneStartShift;              		// If = 1, indicates thet we want to start from codon 2 of ORF
    long cNumber;                     		// Number of codons = 64
    vector<string> geneCodeTable;           // Standard Genetic Code Table
    vector<string> geneCodeTableOrdered;    // Re-Ordered Standard Genetic Code Table
    vector<long> jCodon_to_Nucleotide_Table;
    string rSeqName;
    vector<string> rSequence;
    vector<string> rCodonSeq;          		// Codon sequence of Data Set
    string rSeqMutatedName;            		// Name of mutated ORF
    vector<string> rSeqMutated;             // Mutated Sequence  of an ORF
    string rSeqOpt;
    vector<string> geneName;
    vector<long> geneStart;                 // Gene Start in the DataSet in codons
    vector<long> geneEnd;                   // Gene End in the DataSet in codons
    long jTot;                    			// Number of Genes in the Data set
    vector<long> geneCodons;               	// Total Number of codons in a gene
    MatrixDouble geneRPFtotal;           	// Total Number of RPFs in a gene
    MatrixDouble geneRPFdensity;
    vector<long> gene_Elong_AA;             // Number of Codons in the inner Gene par
    MatrixDouble geneRPF_Elong;            	// Number of RPFs in the inner Gene part
    MatrixDouble geneU_Elong;
    MatrixDouble geneRPF_Elong_Dens;        //RPFs/codon=Ui/nCodons in the inner Gene part
    vector<double> gene_Model_Density;      // Ui/Ti for the inner Gene part
    vector<double> gene_Elong_N2_Time;
    vector<double> growth_Rate_Standard;
    MatrixDouble rpfOmegaExper;             //Experimental RPF fingerprint for this DataSubSet
    vector<long> indexRowsIncluded;

    vector<double> rGeneWeight;
    vector<double> relTimeGene;
    vector<string> rGeneSequence;
    vector<string> geneProtSeq;
    vector<string> strCodonGeneName;
    vector<int> iORFcodons;                 //Index (from 1 to 64) that correspond to codon name
    vector<long> indCodonOrder;
    MatrixDouble nRPF;                      //RPFs ascribed to the codon in the dataset (A-site)
    long nRPFadded;                          // RPFs added per codom; Normally zero
    MatrixDouble pauseScore;
    MatrixDouble pauseScoreModel;
    MatrixDouble pauseScoreModelSigma;
    long p_NarrowModel_First;
    long p_NarrowModel_Last;
    vector<long> p_NarrowModel_Included;
    string strModel_Narrow_Info;
    vector<double> timeTransLoc;
    vector<double> timeTransLocR2;
    vector<double> timeTransLocSigma;
    vector<double> timeTransLoc_Model;
    vector<double> timeTransLoc_Model_Sigma;
    vector<double> timeTransLoc_Model_R2;
    vector<double> rT3ptcTime;
     MatrixDouble rDwellTime;
     MatrixDouble rDwellTimeModel;
     MatrixDouble rDwellTimeModelSigma;
    vector<long> nCodCount;
    vector<double> codUsage;
    vector<double> codRPFsAver;
    vector<double> codRPFsSigma;

    double dK_FA ;  //Inhibitory Constant for FA
    double dK_FA_Sigma;  //Its Sigma
};

// Rcpp interface

List outputList;

List runMLE();

RCPP_MODULE(mod) {
  function( "runMLE", &runMLE );
}

// SUB Declarations

void Get_Sequence_RPFs(int , const string, RPFdataSet&);

void FindArrayOrder(const vector<double>& , vector<long>& );

void Print_DataSet_Statistics(const string, RPFdataSet&);

void Get_ML_zHAT(const string, long, long, long, long, RPFdataSet& ,
        vector<double>& , vector<double>& , vector<double>& );

void Get_Log_Likelihood(const long& , RPFdataSet& , MatrixDouble& ,
    vector<double>&, vector<double>& , double& ,
    vector<double>&, vector<double>& , double& );

void Get_Log_Likelihood_Only(const long& , const RPFdataSet&,
                             const MatrixDouble& , double& );

void RefineStepSize_New(const int&, const RPFdataSet& , const MatrixDouble& ,
                                                const MatrixDouble& , double& );

void Hes_Pos_ML_Refine_zFactors(const string, long, long, long, long,
		RPFdataSet& , vector<double>& , vector<double>& ,
		vector<double>& , vector<double>& ,
		vector<double>& ,vector<double>& );

void Print_vFP_Full_As_MatrixB(const string, string, const string, const int,
        int, int, int, RPFdataSet&, const long ,
		vector<double>& , vector<double>& , vector<double>& );

void VB_Gauss_Solve(int , const MatrixDouble& ,
								vector<double>& , const vector<double>& );

void VB_LUPA_Invert(int, const MatrixDouble& , MatrixDouble& );

void Get_PCorr_Rows_mA_New(long& , long& , const MatrixDouble& ,
								const MatrixDouble& , vector<double>& ,
									vector<double>& , MatrixDouble& );

void Get_PCorr_Cols_mA_New(long , long , const MatrixDouble& ,
								const MatrixDouble& , vector<double>& ,
									vector<double>& , MatrixDouble& );
void Report_R2_zFP_Statistics(string , long , string , string ,
			string , string , double , int& , int& ,
            RPFdataSet& , vector<double>& , vector<double>& );

void Get_X_Y_Correlation(long , long , vector<double>& , vector<double>& ,
	vector<double>& , double& , double& , double& , double& ,double& );

void VB_printVector(string infoVector, vector <double>& vB){
	long j=0;
    long nCol = vB.size()-1;
        std::cout<<"Vector: " <<infoVector<<std::endl;
        std::cout<< "nCol="<<nCol<<std::endl;
    for (j = 1; j<= nCol; j++){
    std::cout<<vB.at(j)<<" ,";
    }
		std::cout<<std::endl;
}

List runMLE(){

	string strDataSetFileName=
	"/Users/gustafullman/Documents/src/RiboTimes/data/ZI_30000_FA0.txt";
	int nRPFadd=0;
	RPFdataSet DS;
  Get_Sequence_RPFs(nRPFadd,strDataSetFileName,DS);

	string outPutFileStatistics=
	"/Users/gustafullman/Documents/src/RiboTimes/data/ZI_OutPut.txt";
	Print_DataSet_Statistics(outPutFileStatistics,DS);

	long jSet=0, pAsite=8, pNumber=15, jGeneStartShift=0;
	vector<double> tML_long, tML_Sigma_long, wML_long;
	vector<double> zFP_long, zFP_Sigma_long, wFP_long;

	string reportOut_zHAT=
	"/Users/gustafullman/Documents/src/RiboTimes/data/Z_HAT_Report.txt";

	Get_ML_zHAT(reportOut_zHAT,jSet,pAsite,pNumber,jGeneStartShift,DS,
        tML_long,tML_Sigma_long,wML_long);

    // print zHAT results
    string zHAT_OutputFile=
    "/Users/gustafullman/Documents/src/RiboTimes/data/zFP_HAT_OutPut.txt";
        string strText="First guess zHAT coefficients";
        string strNorm="NATIVE";
        int iPrint_rpfOmega=1;
        int iPrint_Col_Corr=0;
        int pC_First=1;
        int pC_Last=15;
    Print_vFP_Full_As_MatrixB(zHAT_OutputFile,strText, strNorm, iPrint_rpfOmega,
        iPrint_Col_Corr, pC_First, pC_Last, DS, jSet,
		tML_long,tML_Sigma_long,wML_long);

    string zFP_Refinementr_Log=
    "/Users/gustafullman/Documents/src/RiboTimes/data/Refinement_Log.txt";
    Hes_Pos_ML_Refine_zFactors(zFP_Refinementr_Log, jSet, pAsite, pNumber, jGeneStartShift,DS,
        tML_long,tML_Sigma_long,wML_long, zFP_long, zFP_Sigma_long, wFP_long);


    // more refinement?
    int kIter=0, kEnd=0;
    do{
        std::cout<<"GradNorn<0.5? Stop? To stop enter 1; to continue enter 0" << std::endl;
        //std::cin >>  kEnd;
        kEnd = 0;
        if(kEnd==1) {break;}

    Hes_Pos_ML_Refine_zFactors(zFP_Refinementr_Log, jSet, pAsite, pNumber, jGeneStartShift,DS,
         zFP_long, zFP_Sigma_long, wFP_long, zFP_long, zFP_Sigma_long, wFP_long);
        kIter++;
    } while(kIter<7||kEnd==1);

    // print zFP refined results
    string zFP_OutputFile=
    "/Users/gustafullman/Documents/src/RiboTimes/data/zFP_Refined_OutPut.txt";
       strText="Refined zFP factors";
        strNorm="NATIVE";
        iPrint_rpfOmega=0;
        iPrint_Col_Corr=0;
        pC_First=1;
        pC_Last=15;
    Print_vFP_Full_As_MatrixB(zFP_OutputFile,strText, strNorm, iPrint_rpfOmega,
        iPrint_Col_Corr, pC_First, pC_Last, DS, jSet,
		zFP_long, zFP_Sigma_long, wFP_long);

	// print R2 gene statistics
    string R2_OutPutFile=
    "/Users/gustafullman/Documents/src/RiboTimes/data/R2_OutPut.txt";
        string strMode="";
        string strMark="MARK";
        string strWeight="WEIGHTED";
        string strInfo="R2 Output using zHAT";
        double tol_R2=0.1;
    Report_R2_zFP_Statistics(R2_OutPutFile, jSet, strMode,  strMark,strWeight, strInfo,
                             tol_R2, pC_First, pC_Last, DS,  zFP_long, zFP_Sigma_long);
    
    
    // Add some elements of DS to outputList
    
    outputList["dataSetName"] = DS.dataSetName;
    outputList["DS.dataSubSet_RPF_Total"] = DS.dataSubSet_RPF_Total;
    outputList["doublingTime"] = DS.doublingTime;
    outputList["pAsite"] = DS.pAsite;
    outputList["pNumber"] = DS.pNumber;
    outputList["jGeneStartShift"] = DS.jGeneStartShift;
    outputList["cNumber"] = DS.cNumber;
    outputList["geneCodeTable"] = DS.geneCodeTable;
    outputList["geneCodeTableOrdered"] = DS.geneCodeTableOrdered;
    outputList["rCodonSeq"] = DS.rCodonSeq;
    outputList["geneStart"] = DS.geneStart;
    outputList["geneEnd"] = DS.geneEnd;
    outputList["jTot"] = DS.jTot;
    outputList["geneCodons"] = DS.geneCodons;
    outputList["geneRPFtotal"] = DS.geneRPFtotal;
    outputList["rpfOmegaExper"] = DS.rpfOmegaExper;
    outputList["iORFcodons"] = DS.iORFcodons;
    outputList["nRPF"] = DS.nRPF;

    // Add some other variables of interest to outputList
    
    outputList["zFP"] = zFP_long;
    outputList["zFP_Sigma"] = zFP_Sigma_long;
    
    return outputList;
}

// ================================
void Get_Sequence_RPFs(int nRPFadd, const string strDataSetFileName, RPFdataSet& DS) {
	//
	// transforms information from dataset file into
	// dataset structure DS
	// nRPFadd: normaly=0;
	// strDataSetFileName: the name of a text file containing RPF dataset
	//
	// NOTE: we use positive indexes in long vectors: index "0" is ignored
	// NOTE: to this end the size of vectors is increased by one
	//       the matrix dimentions are also increased by one
	// Implemented by Michael Pavlov

	long mRow;
	string strCodon, strCodonORF, strAA;
	long i=0, j=0, k=0, iCodon;
	long cNumber, indCodon, jSet=0, nFus, kRow, jFus;

	// Gene Code Table
	cNumber = 64; //Number of codons
	int cNumber1=cNumber+1; //dimention of codon arrays

	vector<string> geneCodeTableOrdered(cNumber1);
	// A standard Table of Genetic code
	vector<string> geneCodeTable{"UUU_Phe_F_tF1_GAA" ,
	 "UUU_Phe_F_tF1_GAA", "UUC_Phe_F_tF1_GAA", "UUA_Leu_L_tL4_UAA", "UUG_Leu_L_tL5_CAA",
	 "CUU_Leu_L_tL2_GAG", "CUC_Leu_L_tL2_GAG", "CUA_Leu_L_tL3_UAG", "CUG_Leu_L_tL1_CAG",
	 "AUU_Ile_I_tI1_GAU", "AUC_Ile_I_tI1_GAU", "AUA_Ile_I_tI2_LAU", "AUG_Met_M_tM1_CAU",
	 "GUU_Val_V_tV2_GAC", "GUC_Val_V_tV2_GAC", "GUA_Val_V_tV1_VAC", "GUG_Val_V_tV1_VAC",
	 "UCU_Ser_S_tS5_GGA", "UCC_Ser_S_tS5_GGA", "UCA_Ser_S_tS1_VGA", "UCG_Ser_S_tS2_CGA",
	 "CCU_Pro_P_tP2_GGG", "CCC_Pro_P_tP2_GGG", "CCA_Pro_P_tP3_VGG", "CCG_Pro_P_tP1_CGG",
	 "ACU_Thr_T_tT3_GGU", "ACC_Thr_T_tT3_GGU", "ACA_Thr_T_tT4_UGU", "ACG_Thr_T_tT2_CGU",
	 "GCU_Ala_A_tA2_GGC", "GCC_Ala_A_tA2_GGC", "GCA_Ala_A_tA1_VGC", "GCG_Ala_A_tA1_VGC",
	 "UAU_Tyr_Y_tY1_QUA", "UAC_Tyr_Y_tY1_QUA", "UAA_Ohr_O_rF1", "UAG_Amb_A_rF1",
	 "CAU_His_H_tH1_QUG", "CAC_His_H_tH1_QUG", "CAA_Gln_Q_tQ1_SUG", "CAG_Gln_Q_tQ2_CUG",
	 "AAU_Asn_N_tN1_GUU", "AAC_Asn_N_tN1_GUU", "AAA_Lys_K_tK1_SUU", "AAG_Lys_K_tK1_SUU",
	 "GAU_Asp_D_tD1_QUC", "GAC_Asp_D_tD1_QUC", "GAA_Glu_E_tE2_SUC", "GAG_Glu_E_tE2_SUC",
	 "UGU_Cys_C_tC1_GCA", "UGC_Cys_C_tC1_GCA", "UGA_Opa_O_rf2", "UGG_Trp_W_tW1_CCA",
	 "CGU_Arg_R_tR2_ICG", "CGC_Arg_R_tR2_ICG", "CGA_Arg_R_tR2_ICG", "CGG_Arg_R_tR3_CCG",
	 "AGU_Ser_S_tS3_GCU", "AGC_Ser_S_tS3_GCU", "AGA_Arg_R_tR4_MCU", "AGG_Arg_R_tR5_CCU",
	 "GGU_Gly_G_tG3_GCC", "GGC_Gly_G_tG3_GCC", "GGA_Gly_G_tG2_UCC", "GGG_Gly_G_tG1_CCC"
	 };

	 // Read from a Text File
	string strLine,  strProtein;
	vector<string> strLinesInFile(1,"");
	mRow=1;
	nFus=1;
  std::ifstream dataInput(strDataSetFileName); // define the input object
    //Note line "0" in the input file corresponds to strLinesInFile vector index "1"
    if(dataInput.is_open()){ //Fill lines from the file into string vector
        while(getline(dataInput,strLine)){
           strLinesInFile.push_back(strLine); // construct the vector of file-lines
           mRow++;
        }
    }
    dataInput.close();
  // The dataset is in, and the number of its lines is determined
        std::cout<<"kLast="<<mRow<<std::endl;
		int mRow1=mRow+1;

  //Parse the header lines
  string strField, strField1, strField2;
	long dataSetRPF_Total;
	double setDoublingTime;
	std::stringstream ssLine; // assigne the streaming/parcing buffer
	ssLine.clear();
    ssLine.str("");
	ssLine<<strLinesInFile.at(1);
	ssLine>>strField>>strField1; // Parse the first header line
  if(strField == "dataSet") {// data set Field is found
		DS.dataSetName =strField1;

		// second line
		ssLine.clear();ssLine.str("");
		ssLine<<strLinesInFile.at(2);
		ssLine>>strField>>dataSetRPF_Total;
		DS.dataSubSet_RPF_Total.push_back(dataSetRPF_Total);

		// third line
		ssLine.clear();ssLine.str("");
		ssLine<<strLinesInFile.at(3);
		ssLine>>strField>>setDoublingTime;
		DS.doublingTime.push_back(setDoublingTime);

		// fouth line
		ssLine.clear();ssLine.str("");
		ssLine<<strLinesInFile.at(4);
		ssLine>>strField>>strField1>>strField2;
		DS.dataSubSetName.push_back(strField2);

	// Initialize jGene Arrays putting  0 at position zero
	vector<long> geneStart(1,0);
	vector<long> geneEnd(1,0);
	vector<string> geneName(1, "   ");
	vector<double> geneRPFtotal(1,0.0);

	// Prepare major Array
	vector<string> rCodonSeq(mRow1, "   ");
	vector<int> indCodonSt(mRow1, 0);
	MatrixDouble nRPF(nFus, vector<double>(mRow1,0.0));

	// Parse the data lines from the dataset file
	// convert DNA sequence of the codons  in the first column of the Range into the RNA sequence
	// convert the codons into their indexes in the standard gene-code Table
	long jGene;
	long jStart=0, jEnd=0, jTot=0, iFlagStop;
	int numberRPFs=0;
	string strORF, strSuffix;
			ssLine.clear();ssLine.str("");
			j = 1;
			long iA=0;
			geneStart.push_back(j); //the start of the first codon in the first gene is 1!
		for(i = 5; i<=  mRow; i++) {    //mRow is  the total number of lines=rows in the DataSet
				ssLine<<strLinesInFile.at(i);
				ssLine>>strORF>>strCodonORF>>numberRPFs; // Parse the common line
				ssLine.clear();
				ssLine.str("");
			if(i==5){geneName.push_back(strORF);}  // the name of the first gene
				 //cout<<strLinesInFile[k]<<endl;
			// convert DNA to RNA
			if (strCodonORF[0]== 'T') {strCodonORF[0] = 'U';}
			if (strCodonORF[1]== 'T') {strCodonORF[1] = 'U';}
			if (strCodonORF[2]== 'T') {strCodonORF[2] = 'U';}

			// recod mRNA codon sequence and RPFs per codon
				iA=iA+1;
				rCodonSeq.at(iA) = strCodonORF;
				nRPF[0][iA]=numberRPFs;

			// find and recod codon index in Genetic Code Table
				iCodon = 0;
			for (k = 1; k <= cNumber; k++) {
					strCodon = geneCodeTable.at(k).substr(0,3);
			   if(strCodon == strCodonORF) {
					iCodon = k;
					break;
			   }
			}
				indCodonSt.at(iA) = iCodon; // record the codon index

			if (!(strORF == geneName.at(j))) {// next gene is found
					geneEnd.push_back(iA - 1);
				 j = j + 1; // increase gene counter
				if (strORF == "STOP") {// the last dummy record is found
					break;
				}
					geneStart.push_back (iA);     //Record Gene start
					geneName.push_back (strORF);  //Record Gene Name
			}
		}
				jTot = j - 1; //Total Gener number in the data set
				mRow = geneEnd.at(jTot); // get the number of real data rows
				mRow1=mRow+1; // refined the number of core data lines

	// ============ The main Input is done
	//Check ORFs for consistency (the absence of internal stop codons + proper last stop)
    for (jGene = 1; jGene <= jTot; jGene++) {
            jStart = geneStart.at(jGene);
            jEnd = geneEnd.at(jGene);
            strORF = geneName.at(jGene);
            strSuffix = "";
            iFlagStop = 1;
        for (i = jStart;  i <= jEnd; i++) {
                strCodonORF =rCodonSeq.at(i);
            if ((strCodonORF == "UAA") || (strCodonORF == "UGA") || (strCodonORF == "UAG")) {
				 if (i == jEnd) {
					iFlagStop = 0;
					break;
				 }
				strSuffix = "_ORF Stop"; // a stop inside ORF is found
			}
        }
        if (iFlagStop ==1) {// a gene with no ORF end is found
            strSuffix = strSuffix + "_No End Stop";
        }
        if (!(strSuffix == "")) { // mark non-standard ORFs in the data set
            geneName.at(jGene) = "_" + strORF + strSuffix;
            k = k;
        }
	}
	//
	// get statistics about codon frequences in the data set
	vector<int> nCodStat(cNumber1,0);

		for (k = 1; k <= mRow; k++) {
				iCodon = indCodonSt[k]; // codon index
			nCodStat[iCodon] = nCodStat[iCodon] + 1;
		}
	//
	// find the decending order of codon frequences and create the corresponding
	// Re-Ordere genetic code Table so that the codons are nubered according to their frequences
	// => the codon with new index 1 is the most frequent codon in the dataset
	//
	long nMax, jMax;
	vector<int> indCodonOrder(cNumber1,0);
	vector<int> indCodonReOrder(cNumber1,0);
	vector<int> indORFcodon(mRow1,0);

    for(i = 1; i <= cNumber;  i++) { //step recording the ascending order of codon frequences
          nMax = -1;
		  jMax = 0;
        for (j = 1; j <= cNumber; j++) {
			if (nCodStat[j] > nMax) {
                nMax = nCodStat[j];
                jMax = j;
			}
        }
		indCodonOrder.at(i) = jMax;
		indCodonReOrder.at(jMax) = i;
		geneCodeTableOrdered.at(i) = geneCodeTable.at(jMax);
		nCodStat.at(jMax) = -2; //remove current nMax
    }
	//
	// Change the codon indexes to match the ordered geneTable
	// and record the basic dataeet information in DS structure
	 DS.iORFcodons.resize(mRow1,0);
	 DS.rCodonSeq.resize(mRow1);

	for(i = 1; i <= mRow; i++) {
		strCodonORF = rCodonSeq.at(i);
			DS.rCodonSeq.at(i) = rCodonSeq.at(i);  //Record the codon sequence
			indCodon = indCodonSt.at(i);
				k = indCodonReOrder.at(indCodon);
			indORFcodon.at(i) = k;
			DS.iORFcodons.at(i) = indORFcodon.at(i); //record the corresponding indexes
	}

	// NOTE the exception in the jFus data sub-set indexing in nRPF matrix: it starts from "0"
	DS.nRPF.resize(nFus, vector<double>(mRow1,0.0));
		for( jFus =0; jFus < nFus; jFus++) {
			for( i = 1; i <= mRow; i++) {
					DS.nRPF[jFus][i] = nRPF[jFus][i] + nRPFadd;
			}
		}

	// Do thourough data set statistics
	  double diffCounts=0.0;
	  vector<long> nCodCount(cNumber1,0);
	  vector<double> codRPFsAver(cNumber1,0.0);
	  vector<double> codRPFsSigma(cNumber1,0.0);

    for( i = 1; i <= mRow; i++) {
        iCodon = indORFcodon[i];
        nCodCount[iCodon] = nCodCount[iCodon] + 1;
        codRPFsAver[iCodon] = codRPFsAver[iCodon] + nRPF[0][i];
    }

	// averaged RPFs per codon in the data sub-set jFus=0
    for (k = 1; k <= cNumber; k++) {
        if (nCodCount[k] > 0) {codRPFsAver[k] = codRPFsAver[k] / nCodCount[k];}
    }

    // calculate RPF sigma per codon
    for( i = 1; i <= mRow; i++) {
            iCodon = indORFcodon[i];
            diffCounts =  nRPF[jSet][i] - codRPFsAver[iCodon];
            codRPFsSigma[iCodon] = codRPFsSigma[iCodon] + diffCounts * diffCounts;
    }
    for (k = 1; k <= cNumber; k++) {
        if (nCodCount[k] > 0) { codRPFsSigma[k] = std::sqrt(codRPFsSigma[k] / nCodCount[k]);}
	}

	//Store  the basic and statistic information on the Data Set in DS
	double rRPFgeneAver, pauseScoreTotal;

		DS.nRPFadded = nRPFadd;
		DS.cNumber = cNumber;
		DS.jTot = jTot;
	  DS.geneCodeTable.resize(cNumber1);
	  DS.geneCodeTableOrdered.resize(cNumber1);
	  DS.indCodonOrder.resize(cNumber1); DS.nCodCount.resize(cNumber1);
	  DS.codRPFsAver.resize(cNumber1); DS.codRPFsSigma.resize(cNumber1);
		for (k = 1; k <= cNumber; k++) {
			DS.indCodonOrder.at(k) = indCodonOrder.at(k);
			DS.geneCodeTableOrdered.at(k) = geneCodeTableOrdered.at(k);
			DS.geneCodeTable.at(k) = geneCodeTable.at(k);
			DS.nCodCount.at(k) = nCodCount.at(k);
			DS.codRPFsAver.at(k) = codRPFsAver.at(k);
			DS.codRPFsSigma.at(k) = codRPFsSigma.at(k);
		}
		int jTot1=jTot+1;
	  DS.geneName.resize(jTot1);
	  DS.geneStart.resize(jTot1);
	  DS.geneEnd.resize(jTot1);
    for(jGene = 1; jGene<= jTot; jGene++) {
        DS.geneName.at(jGene) = geneName.at(jGene);
        DS.geneStart.at(jGene) = geneStart.at(jGene);
        DS.geneEnd.at(jGene) = geneEnd.at(jGene);
    }

  // Initialize Vector/Matrices of Data Structure DS
	DS.geneRPFtotal.resize(nFus, vector<double>(jTot1, 0.0));
	DS.geneRPFdensity.resize(nFus, vector<double>(jTot1, 0.0));
	DS.geneCodons.resize(jTot1,0);
	DS.geneRPF_Elong.resize(nFus, vector<double>(jTot1, 0.0));
    DS.pauseScore.resize(nFus, vector<double>(mRow1, 0.0));
	DS.rDwellTime.resize(nFus, vector<double>(mRow1, 0.0));
    DS.rDwellTimeModel.resize(nFus, vector<double>(mRow1, 0.0));
	DS.rDwellTimeModelSigma.resize(nFus, vector<double>(mRow1, 0.0));
    DS.pauseScoreModel.resize(nFus, vector<double>(mRow1, 0.0));
	DS.pauseScoreModelSigma.resize(nFus, vector<double>(mRow1, 0.0));

	// Transform RPFs into rude Pausing  Scores and record them
	// Determine also the codon Usage frequency by the ribosomes
   DS.codUsage.resize(cNumber1,0.0);
    double growthRateStandard =0.0;
    double geneRPF_Total=0;
    double weightTotal =0.0;
    double rWeightNorm=0.0;
        long kORF=0;
	for(jFus =0; jFus < nFus; jFus++) {
			kRow = 0;
			growthRateStandard =0.0;
            geneRPF_Total=0;
			weightTotal =0.0;
		for (jGene = 1; jGene <=jTot; jGene++) {
				jStart = DS.geneStart.at(jGene);
				jEnd = DS.geneEnd.at(jGene);
				DS.geneCodons.at(jGene) = jEnd - jStart + 1;
				kORF = 0;
				geneRPF_Total=0;
			for (k = jStart; k <= jEnd; k++) {
					kORF = kORF + 1;
					kRow = kRow + 1;
				geneRPF_Total =geneRPF_Total+ DS.nRPF[jFus][k];
			}
				DS.geneRPFtotal[jFus][jGene] = geneRPF_Total;
				rRPFgeneAver = 0;
			  if (kORF > 0) {
				rRPFgeneAver = geneRPF_Total / kORF;
			  }
				growthRateStandard = growthRateStandard + rRPFgeneAver;
				DS.geneRPFdensity[jFus][jGene] = rRPFgeneAver; //relative gene weight as average  RPF/codon

				pauseScoreTotal = 0;
			if (rRPFgeneAver > 0) {
				for (k = jStart; k <= jEnd; k++) {
					DS.pauseScore[jFus][k] = DS.nRPF[jFus][k] / rRPFgeneAver;
					pauseScoreTotal = pauseScoreTotal + DS.pauseScore[jFus][k];
					weightTotal = weightTotal +rRPFgeneAver;
					j = DS.iORFcodons.at(k);
					DS.codUsage.at(j) = DS.codUsage.at(j)+rRPFgeneAver;
				}
			}
		}
	}

	// normalize codon usage frequency
		rWeightNorm = weightTotal / kRow;
    for (k = 1; k <= cNumber; k++) {
        DS.codUsage.at(k) = DS.codUsage.at(k)/rWeightNorm;
	}
 }
}


// =================================
 void FindArrayOrder(const vector<double>& vArray, vector<long>& iOrderArray) {
	// It finds a descending order in array vArray without modifying it
	// vArray: input array
	// iOrderArray: output/ j= iOrderArray(i)  j-rearanges the i-order of vArray
	// NOTE: index "0" is ignored here!
	// Implemented by Michael Pavlov

	long i, k;
	long mL1 = vArray.size();
    iOrderArray.resize(mL1,0);
	vector<double> bArray(mL1,0.0);
	long mL=mL1-1;

	//	find vArray boundaries and copy original vArray in bArray
	double maxValue, minValue, bAmaxi;
	long iMax, iMin, iMaxi;

	// Find max and min in vArray
        maxValue = vArray.at(1);
        minValue = vArray.at(1);
        iMax = 1;
        iMin = 1;
	for (i = 1;  i <=  mL; i++) { //step through the values
            bArray.at(i) = vArray.at(i);
            iOrderArray.at(i) = i;
        if (vArray[i] > maxValue) {maxValue = vArray[i]; iMax = i;}
        if (vArray[i] < minValue) {minValue = vArray[i]; iMin = i;}
	}

	// find the decending order of bArray values
	for (k = 1; k <= mL; k++) { // step throug vArray to order it
			bAmaxi = minValue - 1.0; //Set initial vAmaxi well below minValue
			iMaxi = 0;
        for (i = 1;  i <=  mL; i++) {
           if (bArray[i] > bAmaxi) {
                bAmaxi = bArray[i];
                iMaxi = i;
            }
        }
		iOrderArray[k] = iMaxi;
		bArray.at(iMaxi) = minValue - 11.0; //remove the current bAmaxi from consideration
    }
}

// =============================================================
 void Print_DataSet_Statistics(const string strOutPutFile, RPFdataSet& DS){
//
// Print Codon and RPF statistics in the dataset recoded in DS structure
//
// strOutPutFile: 	path to the output file and its name
// DS: 				dataset structure
// Implemented by Michael Pavlov

long i , j , k , iCodon, nCodonTotal, geneLength;
int jGene;

	// Prepare to do the statistics
		int nFus = DS.geneRPFtotal.size(); 		// The number of subsets in the dataset
		int jSet=0;
        int cNumber = DS.geneCodeTableOrdered.size()-1;  // Number of codon in the Gene Code Table
		int cNumber1=cNumber+1;
		int jTot=DS.geneEnd.size()-1;     		//Total Gene number in the data set
		int jTot1=jTot+1;
		long mRow=DS.geneEnd.at(jTot);				// The number of core data rows

	// Do thourough data set statistics
	// get averaged RPFs per codon in the data set jSet
		double diffCounts=0.0;
		vector<long> nCodCount(cNumber1,0);
		vector<double> codRPFsAver(cNumber1,0.0);

    for( i = 1; i <= mRow; i++) {
        iCodon = DS.iORFcodons.at(i);
        nCodCount.at(iCodon) = nCodCount.at(iCodon) + 1;
        codRPFsAver.at(iCodon) = codRPFsAver.at(iCodon) + DS.nRPF[jSet][i];
    }
    for (k = 1; k <= cNumber; k++) {
        if (nCodCount.at(k) > 0) {codRPFsAver.at(k) = codRPFsAver.at(k) / nCodCount.at(k);}
    }

    // also calculate sigma per codon
		vector<double> codRPFsSigma(cNumber1,0.0);
    for( i = 1; i <= mRow; i++) {
            iCodon = DS.iORFcodons.at(i);
            diffCounts =  DS.nRPF[jSet][i] - codRPFsAver.at(iCodon);
            codRPFsSigma.at(iCodon) = codRPFsSigma.at(iCodon) + diffCounts * diffCounts;
    }
    for (k = 1; k <= cNumber; k++) {
        if (nCodCount.at(k) > 0) { codRPFsSigma.at(k) = std::sqrt(codRPFsSigma.at(k) / nCodCount.at(k));}
	}
		nCodonTotal = 0;
		double totalRPFcodon = 0;
		double totalRPF = 0;
    for (k = 1; k <= cNumber; k++) {
         nCodonTotal = nCodonTotal + nCodCount.at(k);
         totalRPFcodon = totalRPFcodon + codRPFsAver.at(k);
         totalRPF = totalRPF + nCodCount.at(k) * codRPFsAver.at(k);
    }

	// OUTPUT for DATASET CODONS
	std::ofstream dataOutput(strOutPutFile);
	 std::stringstream ssLine;
		ssLine.clear(); ssLine.str("");  //Clear Stream Buffer
	ssLine << "Data Set Basic Statistics: #Codons in DataSet = " << nCodonTotal <<
                ";  Total RPFs in Data Set= " << totalRPF<< std::endl;
    dataOutput<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
    ssLine <<"Codon ; Rank  ; #Codons ; Fraction ; codUsage ; Fraction ; RPF/Codon ; "
	<< "Sigm_RPFs ; Codon_Identity"<<std::endl;
	dataOutput<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    for (k = 1; k <= cNumber; k++) {
        ssLine <<DS.geneCodeTableOrdered.at(k).substr(0,6)<<" ; "
        <<k<<" ; "<<nCodCount.at(k)
        <<" ; " << (1000 * nCodCount.at(k) / nCodonTotal)<<" ; "
        <<DS.codUsage.at(k)<<" ; "<<(1000 * DS.codUsage.at(k))<<" ; "
        <<codRPFsAver.at(k)<<" ; "<< codRPFsSigma.at(k)<<" ;  "
        <<DS.geneCodeTableOrdered.at(k)<<std::endl;
      dataOutput<<ssLine.str();
      std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
    }
    dataOutput.close();

	// OUTPUT for DATASET GENES
	std::ofstream geneDataOutput("/Users/gustafullman/Documents/src/RiboTimes/data/ZI_GeneOutput.txt");

    vector<double> aVector(jTot1,0.0);
	vector<long> iVectorOrder(jTot1,0);

    ssLine << "Original Ordering  "<< std::endl;
    geneDataOutput<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

	ssLine << "Gene# ; Gene_Name ; Starts ; #AAs  ; "<<
	"Total_RPF ; RPF/ORF_Codon"<<std::endl;
	geneDataOutput<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    for (jGene = 1; jGene <= jTot; jGene++) {
            geneLength = DS.geneEnd.at(jGene) - DS.geneStart.at(jGene) + 1;
        ssLine <<jGene<<" ; "<< DS.geneName.at(jGene)<< " ; "
        <<DS.geneStart.at(jGene)<<" ; "
        <<geneLength<<" ;  " <<DS.geneRPFtotal[jSet][jGene]<< " ; "
		<<(DS.geneRPFtotal[jSet][jGene] / geneLength)<<std::endl;
      geneDataOutput<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str(""); //do output and clear buffer
        aVector.at(jGene) = DS.geneRPFtotal[jSet][jGene] / geneLength;
        iVectorOrder.at(jGene) = jGene;
    }

	FindArrayOrder(aVector, iVectorOrder);  //mind assending array order

    ssLine<< "Genes Ordered by RPF/ORF_Codon"<<std::endl;
    geneDataOutput<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine<< "Rank ; Gene ; #AAs ; Total_RPF ; RPF/ORF_Codon" << std::endl;
    geneDataOutput<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    for (j= 1; j <= jTot; j++) {
        jGene = iVectorOrder.at(j);
		geneLength = DS.geneEnd.at(jGene) - DS.geneStart.at(jGene) + 1;
        ssLine <<j<<"  "<< DS.geneName.at(jGene)<< "  "
        <<geneLength<<"   " <<DS.geneRPFtotal[jSet][jGene]<< "   "
		<<(DS.geneRPFtotal[jSet][jGene] / geneLength)<<std::endl;
      geneDataOutput<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str(""); //do output and clear buffer
    }
     geneDataOutput.close();
 }

//=============== Principle Simple FingerPrint Sub
void Get_ML_zHAT(const string strReportFile,long jSet, long pAsite , long pNumber , long jGeneStartShift,
    RPFdataSet& DS, vector<double>& tML_long, vector<double>& tML_Sigma_long, vector<double>& wML_long) {
	//
	// Calculates ML Fingerprints
	// Implemented by Michael Pavlov

	int iPrint = 1;

    int nFus = DS.nRPF.size();
    long mRow1 = DS.iORFcodons.size();
	long mRow=mRow1-1;
    int jTot1 = DS.geneEnd.size();
	int jTot=jTot1-1;   //number of genes in data set
    int cNumber = 64;  //Number of codons
	int cNumber1=cNumber+1;
    long rowLength = pNumber * cNumber;
	long rowLength1=rowLength+1;

	// resize main output arrays
	tML_long.resize(rowLength1,0.0);
	tML_Sigma_long.resize(rowLength1,0.0);
	wML_long.resize(rowLength1,0.0);

        DS.pAsite = pAsite;
        DS.pNumber = pNumber;
        DS.jGeneStartShift = jGeneStartShift;

	int pNumber1=pNumber+1;

	long i, j, k, kA, Iter, iA, jA;
	string strCodon;
	long  nActive, nCol =0, pShift =0;
	long jStart =0, jEnd =0;
	int jGene=0, iPos =0;
	double rpfOmega_PC =0.0, muCoff_PC =0.0, tModel_k=0.0, muCheck=0.0;
	double nRPF_kA =0.0, dRPF_kA =0.0;
	long  nPosCodActive;
	int  indCodon, iCod, indCod, nCodon_Elong, kRow;
	double wtFetha, coff_Fy, muCoff_P_Cod_Sum;

    MatrixLong nPosCodInGene(pNumber1, vector<long>(cNumber1,0));
    MatrixDouble rpfOmegaExper(pNumber1,vector<double>(cNumber1, 0.0));
    MatrixDouble zHAT(pNumber1, vector<double>(cNumber1, 0.0));
    MatrixDouble zHAT_Old(pNumber1, vector<double>(cNumber1, 0.0));
    MatrixDouble muCoff(pNumber1, vector<double>(cNumber1, 0.0));
    vector<double> muCoff_Cod_Sum(pNumber1);
    vector<double> gradActive(rowLength);

	double elongRPF_Total=0.0, geneRPF_Elong =0.0;
	double tPosCodon =0, tPosCodonDiff =0, tDiff2 =0;
	vector<long> lCodonGeneElong(jTot1);
	vector<double> uGeneElong(jTot1);
	MatrixDouble tGeneElong(pNumber1, vector<double>(jTot1,0.0));
	MatrixDouble tGeneElongNew(pNumber1, vector<double>(jTot1,0.0));

	DS.rpfOmegaExper.resize(pNumber1, vector<double>(cNumber1,0.0));

	//
	double  coffB, t_Gene_Elong;
	vector<double> muNew(pNumber1);
	MatrixDouble tLhd_Gene(pNumber1, vector<double>(jTot1,0.0));
	vector<double> tLhd_Total(pNumber1);
	vector<double> tLhd_Gene_UB(jTot1);

	double tLhd_Total_UB;
	double gene_RPF_ln_tModel, gene_RPF_per_Codon, gene_RPF_ln_RPF;
	double gradActiveNorm2, gradActiveNorm;

	//empty the arrays
    for (iPos = 1;  iPos <= pNumber; iPos++){
		for( indCod = 1;  indCod <= cNumber; indCod++) {
			rpfOmegaExper[iPos][indCod] = 0.0;
			zHAT_Old[iPos][indCod] = 1.0;
			zHAT[iPos][indCod] = 0.0;
		}
    }

	// Prepare rpfOmega(p,c) FingerPring, Standard Gene Times, log Likelihood
	//  and miscelanous arrays
		double muInitial = 0.0;
		kRow = 0;
	for(jGene = 1; jGene <= jTot; jGene++){
			   jStart = DS.geneStart.at(jGene) + jGeneStartShift;
			   jEnd = DS.geneEnd.at(jGene);
					nCodon_Elong = 0;
					geneRPF_Elong = 0.0;
				gene_RPF_ln_RPF = 0.0;

	   for (k = jStart; k<= jEnd - pNumber; k++) {
				kRow = kRow + 1; 		//Counts total number of used codons
				kA = pAsite + k - 1; 	//Count from the A-site
			nRPF_kA = DS.nRPF[jSet][kA];
			dRPF_kA = nRPF_kA;
				   nCodon_Elong = nCodon_Elong + 1;
				   geneRPF_Elong = geneRPF_Elong + nRPF_kA;

		   for (iPos = 1; iPos <= pNumber; iPos++) {
				indCodon = DS.iORFcodons.at(kA + iPos - pAsite);
				strCodon = DS.geneCodeTableOrdered.at(indCodon);
			rpfOmegaExper[iPos][indCodon] = rpfOmegaExper[iPos][indCodon] + nRPF_kA;
			}

			if (nRPF_kA > 1) {
				gene_RPF_ln_RPF = gene_RPF_ln_RPF + dRPF_kA * std::log(dRPF_kA);
			}
		}
			uGeneElong.at(jGene)= geneRPF_Elong;
			lCodonGeneElong.at(jGene) = nCodon_Elong;
			muInitial = muInitial + geneRPF_Elong / nCodon_Elong;
			elongRPF_Total = elongRPF_Total + geneRPF_Elong;
		for (iPos = 1; iPos <= pNumber; iPos++) {
			tGeneElong[iPos][jGene] = nCodon_Elong;
		}

			tLhd_Gene_UB.at(jGene) = 0;
		if (geneRPF_Elong > 0) {
			gene_RPF_per_Codon = geneRPF_Elong / nCodon_Elong;
		  tLhd_Gene_UB.at(jGene) = gene_RPF_ln_RPF - geneRPF_Elong * std::log(gene_RPF_per_Codon);
		}
		  tLhd_Total_UB = tLhd_Total_UB + tLhd_Gene_UB.at(jGene);
	}

	//record rpfOmega RPF FingerPrint
	for (iPos = 1; iPos <= pNumber; iPos++) {
		for (indCodon = 0; indCodon <= cNumber; indCodon++) {
			DS.rpfOmegaExper[iPos][indCodon] = rpfOmegaExper[iPos][indCodon];
		}
	}
	//
	// main iterations start =====================================
	iPrint = 1;
        std::ofstream reportOut(strReportFile);
        std::stringstream ssLine;
		ssLine.clear(); ssLine.str("");  //Clear Stream Buffer
    if (iPrint == 1) {
        ssLine<< "DataSet: "<< DS.dataSetName<<std::endl;
        reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "#_of_Genes= "<< jTot<<std::endl;
        reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "subSet_#= " << jSet<<std::endl;
        reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        //std::cout<< "SubSet: "<< DS.dataSubSetName[jSet]<<std::endl;
        ssLine<< "pAsite= "<< DS.pAsite << " pNumber= " << DS.pNumber<<std::endl;
        reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "Iter    Log_LiHo(pA)  Grag_Norm     Shift_Norm      Coff_B"<<std::endl;
        reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
	}

	//================================================================

	for (Iter = 1; Iter<=  21; Iter++) {

		//Get muCoff
		//Zero muCoff Array
		for (iPos = 1; iPos <= pNumber; iPos++) {
			for (indCodon = 1; indCodon <= cNumber; indCodon++) {
				muCoff[iPos][indCodon] = 0.0;
			}
		}

	   for(jGene = 1; jGene<= jTot; jGene++){
			 jStart = DS.geneStart.at(jGene) + jGeneStartShift;
			 jEnd = DS.geneEnd.at(jGene);

			//Zero nPosCodInGene Array
			for (iPos = 1; iPos <= pNumber; iPos++) {
				for (indCodon = 1; indCodon <= cNumber; indCodon++) {
					nPosCodInGene[iPos][indCodon] = 0;
				}
			}

			//Fill nPosCodInGene
			for (k = jStart; k<= jEnd - pNumber; k++) {
					kA = pAsite + k - 1; //Count from the A-site
				for (iPos = 1; iPos <= pNumber; iPos++) {
					indCodon = DS.iORFcodons.at(kA + iPos - pAsite);
				   nPosCodInGene[iPos][indCodon] = nPosCodInGene[iPos][indCodon] + 1;
				}
			}

			//Use nPosCodInGene to get muCoff=Effective weights of Codons
			for (iPos = 1; iPos<= pNumber; iPos++) {
				   wtFetha = uGeneElong.at(jGene) / tGeneElong[iPos][jGene]; //U(j)/T(j)
				for (indCodon = 1; indCodon <= cNumber; indCodon++) {
					muCoff[iPos][indCodon] = muCoff[iPos][indCodon] +
						 nPosCodInGene[iPos][indCodon] * wtFetha;
				}
			}
		}

			//
			for (iPos = 1; iPos <= pNumber; iPos++) {
					muCoff_P_Cod_Sum = 0;
				for (indCodon = 1; indCodon <= cNumber; indCodon++) {
					muCoff_P_Cod_Sum = muCoff_P_Cod_Sum + muCoff[iPos][indCodon];
				}
					muCoff_Cod_Sum.at(iPos) = muCoff_P_Cod_Sum;
			}
	//
	// Calculate the gradient
			nActive = 0;
			gradActiveNorm2 = 0.0;
	   for (iPos = 1; iPos <= pNumber; iPos++) {
		   for (indCodon = 1; indCodon <= cNumber; indCodon++) {
					tPosCodon = zHAT_Old[iPos][indCodon];
				if (tPosCodon > 0) {
						nActive = nActive + 1;
						muCoff_PC = muCoff[iPos][indCodon];
						rpfOmega_PC = rpfOmegaExper[iPos][indCodon];
					if (muCoff_PC > 0) {
						gradActive.at(nActive) = muCoff_PC - rpfOmega_PC / tPosCodon;
						gradActiveNorm2 = gradActiveNorm2 +
						            gradActive.at(nActive) * gradActive.at(nActive);
					}
				}
			}
		}
		gradActiveNorm = std::sqrt(gradActiveNorm2);

    //Get refined zHAT FingerPrints
			nPosCodActive = 0;
	   for (iPos = 1; iPos <= pNumber; iPos++) {
		   for (indCodon = 1; indCodon <= cNumber; indCodon++) {
					zHAT[iPos][indCodon] = 0;
				if (muCoff[iPos][indCodon] > 0) {
				   zHAT[iPos][indCodon] = rpfOmegaExper[iPos][indCodon]
														/ muCoff[iPos][indCodon];
				   nPosCodActive = nPosCodActive + 1;
				}
				else {
				  i = i;
				}
			}
		}

    //Calculate new Likelihood -function, recalculate gene elongation times and new growth rate
	   for (iPos = 1; iPos <= pNumber; iPos++) {
				tLhd_Total.at(iPos) = 0.0;
				muNew.at(iPos) = 0.0;
			for(jGene = 1; jGene <= jTot; jGene++){
				jStart = DS.geneStart.at(jGene) + jGeneStartShift;
				jEnd = DS.geneEnd.at(jGene);
					t_Gene_Elong = 0.0;
					nCodon_Elong = 0;
					gene_RPF_ln_tModel = 0.0;
				for (k = jStart; k<= jEnd - pNumber; k++) {
					kA = pAsite + k - 1; //A-site in the original data set
					indCodon = DS.iORFcodons.at(kA + iPos - pAsite);
					  tModel_k = zHAT[iPos][indCodon];
					  t_Gene_Elong = t_Gene_Elong + tModel_k;
					  nCodon_Elong = nCodon_Elong + 1;
					if (tModel_k > 0) {
						gene_RPF_ln_tModel = gene_RPF_ln_tModel +
                                                    DS.nRPF[jSet][kA] * std::log(tModel_k);
					}
				}
					 tGeneElongNew[iPos][jGene] = t_Gene_Elong;
					muNew.at(iPos) = muNew.at(iPos) + uGeneElong.at(jGene) / t_Gene_Elong;
					 tLhd_Gene[iPos][jGene] = gene_RPF_ln_tModel -
                                uGeneElong.at(jGene) * std::log(t_Gene_Elong / nCodon_Elong);
					 tLhd_Total.at(iPos) = tLhd_Total.at(iPos) + tLhd_Gene[iPos][jGene];
			}
				muNew.at(iPos) = muNew.at(iPos);
				tLhd_Total.at(iPos) = tLhd_Total.at(iPos);
		}

		// re-normalize  zHAT FingerPrints
			tDiff2 = 0;
		for (iPos = 1; iPos <= pNumber; iPos++) {
				muCheck = 0;
				coffB = muNew.at(iPos) / muInitial;
			for (indCodon = 1; indCodon <=  cNumber; indCodon++) {
				tPosCodon = zHAT[iPos][indCodon] * coffB;
			  zHAT[iPos][indCodon] = tPosCodon;
				tPosCodonDiff = (tPosCodon - zHAT_Old[iPos][indCodon]);
				tDiff2 = tDiff2 + tPosCodonDiff * tPosCodonDiff;
			  zHAT_Old[iPos][indCodon] = zHAT[iPos][indCodon];
			}

			for(jGene = 1; jGene <= jTot; jGene++){
				tGeneElong[iPos][jGene] = tGeneElongNew[iPos][jGene] * coffB;
				muCheck = muCheck + uGeneElong.at(jGene) / tGeneElong[iPos][jGene];
			}
		}
			tDiff2 = tDiff2 / rowLength;

			if (iPrint == 1){
				ssLine<< Iter<< " " << tLhd_Total.at(pAsite)<< " " << gradActiveNorm << "; "
				<< std::sqrt(tDiff2)<< " " << coffB<< std::endl;
				reportOut<<ssLine.str();
                std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
			}

			if (tDiff2 < 1E-22) {break;}  // exit Iter for
	}
	// Main Iter Ends ==============================
	//===============
	//
	// Finalize the Output
	// Calculate the Hessian at the maximum
	//   Prepare active codon array
	vector<int> jCodActive(cNumber1,0);
	int nPosCodActive1;
			pShift = 0;
        MatrixDouble mHM, mHM1, mHS;
        int iFlagInvert;
        double det_mHS;
        vector<double> diagHM1_Full(cNumber1,0.0);
	for (iPos = 1; iPos <= pNumber; iPos++) {

				nPosCodActive = 0; //get number of active codons at iPos position
			for (indCodon = 1; indCodon <=  cNumber; indCodon++) {
				if (zHAT[iPos][indCodon] > 0) {
					nPosCodActive = nPosCodActive + 1;
					jCodActive.at(nPosCodActive) = indCodon;
				}
			}

		//construct  Hessians for each iPos position
            nPosCodActive1=nPosCodActive+1;
        mHM.clear(); mHM1.clear(); mHS.clear();
        mHM.resize(nPosCodActive1, vector<double> (nPosCodActive1,0.0));
        mHM1.resize(nPosCodActive, vector<double> (nPosCodActive,0.0));
        mHS.resize(nPosCodActive, vector<double> (nPosCodActive,0.0));

		// Calculate singular Hessian
		for(jGene = 1; jGene<= jTot; jGene++) {
			wtFetha = uGeneElong.at(jGene) / tGeneElong[iPos][jGene]; // =U(jGene)/T(jGene)
			coff_Fy = wtFetha / tGeneElong[iPos][jGene];  // =U(jGene)/(T(jGene)*T(jGene))
			jStart = DS.geneStart.at(jGene) + jGeneStartShift;
			jEnd = DS.geneEnd.at(jGene);

			// Zero nPosCodInGene
				for (indCodon = 1; indCodon <=  cNumber; indCodon++) {
					nPosCodInGene[iPos][indCodon] = 0;
				}

			// Obtain nPosCodInGene for a Gene
				for (k = jStart; k<= jEnd - pNumber; k++) {
					indCod =  DS.iORFcodons.at(iPos + k - 1);
				  nPosCodInGene[iPos][indCod] = nPosCodInGene[iPos][indCod] + 1;
				}

			for (i = 1; i<= nPosCodActive; i++) {
					indCodon = jCodActive.at(i); // the current indCodon index
				for (j = 1; j<= nPosCodActive; j++) {
					iCod = jCodActive.at(j); 	// the current indCodon index
				mHM[i][j] = mHM[i][j] + coff_Fy *
					nPosCodInGene[iPos][indCodon] * nPosCodInGene[iPos][iCod];
				}
			}
		}

		// Finalize the diagonal
		for (i = 1; i<= nPosCodActive; i++) {
			  indCodon = jCodActive.at(i); //the current indCodon index
			  tPosCodon = zHAT[iPos][indCodon];
		 mHM[i][i] = mHM[i][i]  - rpfOmegaExper[iPos][indCodon] / (tPosCodon * tPosCodon);
		}

		// Prepare the diagonal of mHM-1
		for (indCodon = 1; indCodon <=  cNumber; indCodon++) {
				diagHM1_Full.at(indCodon) = 0.0;
		}

	 // Exclude k-th row/column to match it with the parameter number of the L-function (k-Par fixed)
     for(k = 1;  k<= 2; k++) {
                iA = 0;
            for (i = 1; i<= nPosCodActive; i++) {
                if (!(i == k)) {
                        iA = iA + 1;
                        jA = 0;
                    for (j = 1; j<= nPosCodActive; j++) {
                        if (!(j == k)) {
                                jA = jA + 1;
                            mHS[iA][jA] = mHM[i][j];
                        }
                    }
                }
            }

        // invert the hessian the curtailed Hessian
            VB_LUPA_Invert(1, mHS,  mHM1);

        // Record the diagonal of the inverted Hessian
                iA = 0;
            for (i = 1; i<= nPosCodActive; i++) {
                if (!(i == k)) {
                        iA = iA + 1;
                    indCodon = jCodActive.at(i); //the current active  indCodon index
                  diagHM1_Full.at(indCodon) = diagHM1_Full.at(indCodon) + mHM1[iA][iA];
                }
            }
        }

    //Finalyze ML tFP
        for (indCodon = 1; indCodon <=  cNumber; indCodon++) {
            tML_long.at(pShift + indCodon) = zHAT[iPos][indCodon];
            tML_Sigma_long.at(pShift + indCodon) = std::sqrt(std::abs(diagHM1_Full.at(indCodon) / 2));
           //wML_long(pShift + indCodon) = kRow * muCoff[iPos][indCodon] / muCoff_Cod_Sum.at(iPos);
            wML_long.at(pShift + indCodon) = kRow * muCoff[iPos][indCodon] / elongRPF_Total;
        }
		pShift = pShift + cNumber;
			muNew.at(iPos) = muNew.at(iPos);
			t_Gene_Elong = tGeneElong[iPos][2];
	}

	// Case of printing out
	// Identify left upper corner of the Print Out

        ssLine.clear(); ssLine.str("");
	if(iPrint == 1) {
			ssLine<< "Maximum Likelyhood when model depends only on the codon at iPos"<<std::endl;
			reportOut<<ssLine.str();
            std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
		for (iPos = 1; iPos <= pNumber; iPos++) {
		   ssLine<< tLhd_Total.at(iPos)<< "; ";
		}
		ssLine<<std::endl;
		reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

		ssLine<<"Absolute Max of log-Likelihood for the model defined by codons in pA and pP positions= "
			<< tLhd_Total_UB <<std::endl;
        reportOut<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

	}
}

// ====================================================================================
//
void Get_Log_Likelihood(const long& jSet, RPFdataSet& DS, MatrixDouble& zPC_Matrix,
    vector<double>&tModel, vector<double>& tGeneElong, double& tLhd_Total,
    vector<double>&tLhd_Gene, vector<double>& tLhd_Gene_UB, double& tLhd_Total_UB) {

  //Calculates likelihood function (returns tLhd_Total)
  //
  //  jSet:   input, specifies the subset of the data set to be used in calculation
  //  DS:     input, the source of information about data set, like RPFs for A-site
  //                      codons, pNumber, pAsite, etc.
  //  zFP:    input, the Matrix of zPC parameters of our model
  //  tModel: output, the model times per codon calculated in the sub from zPF Matrix
  //
  //  tGeneElong: 	output, the model time for elongation through the inner gene region
  //  tLhd_Total: 	output, the main output: likelihood of the Model for the data subset
  //  tLhd_Gene : 	output, Likelihood of the Model for each gene
  //  tLhd_Gene_UB:		output, Upper Bound of Likelihood for each gene
  //  tLhd_Total_UB:	output, Upper Bound of Likelihood for the data subset
  //
  // Implemented by Michael Pavlov
  //
	long jGene, iPos, cNumber, pNumber, pAsite;
	long  jStart, jEnd, jGeneStartShift, jTot, mRow, indCodon;
	long i, j, k, kA;

	double gene_RPF_ln_tModel, gene_RPF_ln_RPF;

        jTot = DS.geneEnd.size()-1; // 'number of genes in data set
        mRow = DS.geneEnd.at(jTot); // 'number of codons in fata set
		long mRow1=mRow+1;
        jGeneStartShift = DS.jGeneStartShift;
        pAsite = DS.pAsite;
        pNumber = DS.pNumber;
        cNumber = 64; // Number of codons

    long nCodon_Elong, nRPF_kA;
	double geneRPF_Elong, t_Gene_Elong, gene_RPF_per_Codon, gene_Time_per_Codon;
	double dRPF_kA, t_pAsite;

        tLhd_Total = 0; 	// total likelyhood function
        tLhd_Total_UB = 0;  // upper Bound of the total likelyhood function
	for(jGene = 1; jGene <= jTot; jGene++)  {
            jStart = DS.geneStart.at(jGene) + jGeneStartShift;
            jEnd = DS.geneEnd.at(jGene);
             nCodon_Elong = 0;
             geneRPF_Elong = 0;
             t_Gene_Elong = 0;
         gene_RPF_ln_tModel = 0; // Likelyhood calculation
         gene_RPF_ln_RPF = 0;    // Likelyhood Upper Bound Calculation

		for(k = jStart;  k <= (jEnd - pNumber); k++) {
				t_pAsite = 1; // time in the A-site
            for (iPos = 1; iPos <= pNumber; iPos++) {
                indCodon = DS.iORFcodons.at(iPos + k - 1);
                t_pAsite = t_pAsite * zPC_Matrix[iPos][indCodon];
            }
                kA = pAsite + k - 1; // A-site  index in the original data set
			  tModel.at(kA) = t_pAsite;

               nCodon_Elong = nCodon_Elong + 1;
               dRPF_kA = DS.nRPF[jSet][kA];
               nRPF_kA = dRPF_kA;
               geneRPF_Elong = geneRPF_Elong + dRPF_kA;
               t_Gene_Elong = t_Gene_Elong + t_pAsite;

            if (nRPF_kA > 0) {
                    gene_RPF_ln_RPF = gene_RPF_ln_RPF + dRPF_kA * std::log(dRPF_kA);
				if (t_pAsite > 0) {
                    gene_RPF_ln_tModel = gene_RPF_ln_tModel + dRPF_kA * std::log(t_pAsite);
				}
			}
		}
		   tGeneElong.at(jGene) = t_Gene_Elong;
			gene_Time_per_Codon = t_Gene_Elong / nCodon_Elong;
			tLhd_Gene.at(jGene) = gene_RPF_ln_tModel -
                        geneRPF_Elong * std::log(gene_Time_per_Codon);
			tLhd_Total = tLhd_Total + tLhd_Gene.at(jGene);

		// Upper bounds
			tLhd_Gene_UB.at(jGene) = 0;
		if (geneRPF_Elong > 0) {
			gene_RPF_per_Codon = geneRPF_Elong / nCodon_Elong;
		  tLhd_Gene_UB.at(jGene) = gene_RPF_ln_RPF -
                            geneRPF_Elong * std::log(gene_RPF_per_Codon);
		}
		  tLhd_Total_UB = tLhd_Total_UB + tLhd_Gene_UB.at(jGene);
	}
}

// ===============================================================================
void Get_Log_Likelihood_Only(const long& jSet, const RPFdataSet& DS, const MatrixDouble& zPC_Matrix,
													double& tLhd_Total) {
  //Calculates likelihood function (returns tLhd_Total)
  //
  //  jSet:   input, specifies the subset of the data set to be used in calculation
  //  DS:     input, the source of information about data set, like RPFs for A-site
  //                      codons, pNumber, pAsite, etc.
  //  zFP:    input, the Matrix of zPC parameters of our model
  //  tLhd_Total: 	output, the main output: likelihood of the Model for the data subset
  //
  // Implemented by Michael Pavlov
  //

	long jGene, iPos, cNumber, pNumber, pAsite;
	long jStart, jEnd, jGeneStartShift, jTot, mRow, indCodon;
	long i, j, k, kA;

	double gene_RPF_ln_tModel, gene_RPF_ln_RPF;

        jTot = DS.geneEnd.size()-1; // 'number of genes in data set
        jGeneStartShift = DS.jGeneStartShift;
        pAsite = DS.pAsite;
        pNumber = DS.pNumber;
        cNumber = 64; // Number of codons

    long nCodon_Elong, nRPF_kA;
	double geneRPF_Elong, t_Gene_Elong,  gene_Time_per_Codon;
	double dRPF_kA, t_pAsite;

        tLhd_Total = 0; 	// total likelyhood function
	for(jGene = 1; jGene <= jTot; jGene++)  {
            jStart = DS.geneStart.at(jGene) + jGeneStartShift;
            jEnd = DS.geneEnd.at(jGene);
             nCodon_Elong = 0;
             geneRPF_Elong = 0;
             t_Gene_Elong = 0.0;
         gene_RPF_ln_tModel = 0.0; // Likelyhood calculation
         gene_RPF_ln_RPF = 0.0;    // Likelyhood Upper Bound Calculation

      for(k = jStart;  k <= (jEnd - pNumber); k++) {
		  		t_pAsite = 1; // time in the A-site
            for (iPos = 1; iPos <= pNumber; iPos++) {
                indCodon = DS.iORFcodons.at(iPos + k - 1);
                t_pAsite = t_pAsite * zPC_Matrix[iPos][indCodon];
            }
                kA = pAsite + k - 1; // A-site  index in the original data set
               dRPF_kA = DS.nRPF[jSet][kA];
               nRPF_kA = dRPF_kA;
               geneRPF_Elong = geneRPF_Elong + dRPF_kA;
               t_Gene_Elong = t_Gene_Elong + t_pAsite;
			   nCodon_Elong = nCodon_Elong + 1;
            if (nRPF_kA > 0) {
                    gene_RPF_ln_RPF = gene_RPF_ln_RPF + dRPF_kA * std::log(dRPF_kA);
				if (t_pAsite > 0) {
                    gene_RPF_ln_tModel = gene_RPF_ln_tModel + dRPF_kA * std::log(t_pAsite);
				}
			}
		}
			gene_Time_per_Codon = t_Gene_Elong / nCodon_Elong;
			tLhd_Total = tLhd_Total + gene_RPF_ln_tModel -
                        geneRPF_Elong * std::log(gene_Time_per_Codon);
	}
}

//=============================
void RefineStepSize_New(const int& jSet, const RPFdataSet& DS, const MatrixDouble& zFP,
                                                const MatrixDouble& zShift, double& stepSize) {

	// The Sub refines the stepSize by maximizing the likelihood function along the vDir_Short
	//										direction
	//
	//  jSet:   input, specifies the subset of the data set to be used in calculation
	//  DS:     input, the source of information about data set, like RPFs for A-site
	//                      codons, pNumber, pAsite, etc.
	//  zFP:    input, the current Matrix of zPC parameters
	//  tModel: output, the model times per codon calculated in the sub from zPF Matrix
	//  tGeneElong: 	output, the model time for elongation through the inner gene region
	//  tLhd_Total: 	output, the main output: likelihood of the Model for the data subset
	//  stepSize:   main output, the step size that maximizes the likelihood along  vDir
	//  zShift_Long:   main input, the suggested direction of zFP shift
	//
	// Implemented by Michael Pavlov

	// Standard block
	int i, j, k, iStepSize;
	int iPos, pNumber, indCodon, cNumber;

	// Quadratic Interpolation Var
	double y1, y2, y3, v1, v2, v3;
	double d32, d31, d21;
	double y31_Length, y32_Length, y21_Length, yTol;

	double aCoff, bCoff, yI;
	int iDummy, iFlag_Zero, iStepLast, iA34, iFlag;

	double stepSizeMin,  stepSizeCurr;
	double tLhd_Total_Max, tLhd_Total, vI, dTol;

    int nStepMax = 19;
	int nStepMax1 = nStepMax+1;
    stepSizeCurr=stepSize;
    stepSizeMin = std::abs(stepSize) / 10;

    yTol = stepSizeMin * 0.001;
    dTol = 0.000001;

	vector<double> stepSizeSequence(nStepMax1), fValueSequence(nStepMax1);

        stepSizeSequence.at(1) = 0.0;
     Get_Log_Likelihood_Only(jSet, DS, zFP, tLhd_Total);
        fValueSequence.at(1) = tLhd_Total;
        tLhd_Total_Max = tLhd_Total;

	// do the shift to Get new z-factor and estimate the convergence
    pNumber = DS.pNumber;
	int pNumber1=pNumber+1;
    cNumber = 64; // Number of codons
    int cNumber1=cNumber+1;

	MatrixDouble zFP_New(pNumber1, vector<double>(cNumber1,0.0));


    // safeguard from the negatives in the new, post-shift parameter set zFP_New
		iFlag_Zero = 0;
	for(iA34 = 1; iA34 <= 12; iA34++) {
             iFlag = 0;
        for (iPos = 1; iPos <= pNumber; iPos++) {
                for(indCodon = 1; indCodon <= cNumber; indCodon++) {
                     zFP_New[iPos][indCodon] = zFP[iPos][indCodon] + stepSizeCurr * zShift[iPos][indCodon];
                   if (zFP_New[iPos][indCodon] < 0) { // we are in troubles: shift results in negative
                     iFlag = 1;
                     iFlag_Zero = 1;
                    break; // exit indCodon loop
                   }
                }
            if (iFlag == 1) {break;} // exit iPos loop
        }
      if(iFlag == 0) { // then we are good, exit iA34
                break;}
        else{//half the step size and re-try to step positively
                stepSizeCurr = stepSizeCurr / 2;}
	}

         iFlag = 0;
    //Check out the largest step
       Get_Log_Likelihood_Only(jSet, DS, zFP_New, tLhd_Total);
           stepSizeSequence.at(nStepMax) = stepSizeCurr; // maximal allowed step size
           fValueSequence.at(nStepMax) = tLhd_Total; // the corresponding Likelihood
       if(tLhd_Total > tLhd_Total_Max) { // got better LH then with step zero
            iFlag = iFlag + 1;
            tLhd_Total_Max = tLhd_Total;
        }

	// Reduce step size twice each time until the minimum is found

	for(iStepSize = nStepMax - 1; iStepSize >= 2; iStepSize=iStepSize-1) {
        iStepLast = iStepSize;
        stepSizeCurr = stepSizeCurr / 2; // Reduce the stepSize Twice

            for (iPos = 1; iPos <= pNumber; iPos++) {
                for(indCodon = 1; indCodon <= cNumber; indCodon++) {
					zFP_New[iPos][indCodon] = zFP[iPos][indCodon] + stepSizeCurr * zShift[iPos][indCodon];
                }
            }

        //Calculate a new  value of the L-function
		Get_Log_Likelihood_Only(jSet, DS, zFP_New, tLhd_Total);
            stepSizeSequence.at(iStepSize) = stepSizeCurr;
            fValueSequence.at(iStepSize) = tLhd_Total;

        if (tLhd_Total > tLhd_Total_Max) { // got better?
			iFlag = iFlag + 1;
			tLhd_Total_Max = tLhd_Total;
        }
         d21 = fValueSequence.at(iStepSize+1) - fValueSequence.at(iStepSize);
        if (iFlag > 0 && tLhd_Total < tLhd_Total_Max) { // function declines after the maximum?
              break; // exit iStepSize loop
		}
	}

            d21 = fValueSequence.at(nStepMax) - fValueSequence.at(nStepMax - 1);
            if (d21 >= dTol) {
                stepSize = stepSizeSequence.at(nStepMax);
                if (std::abs(d21) < dTol){
                   stepSize = stepSizeSequence.at(nStepMax-1);
                }
              return;
			}

           if (std::abs(d21) < dTol) {
               stepSize = stepSizeSequence.at(nStepMax - 1);
                return;
           }

        // Use 3-point bracket, v2 is max

        y1 = stepSizeSequence.at(iStepLast);
        y2 = stepSizeSequence.at(iStepLast + 1);
        y3 = stepSizeSequence.at(iStepLast + 2);
        v1 = fValueSequence.at(iStepLast);
        v2 = fValueSequence.at(iStepLast + 1);
        v3 = fValueSequence.at(iStepLast + 2);

	// Refine step using quadratic interpolation
	for(iDummy = 1; iDummy <= 6; iDummy++) {
		d32 = (v3 - v2) / (y3 - y2);
		d21 = (v2 - v1) / (y2 - y1);
		d31 = (v3 - v1) / (y3 - y1);
			aCoff = (d32 - d21) / (y3 - y1);
			bCoff = d31 - aCoff * (y3 + y1);
			yI = -bCoff / (2 * aCoff);

	    // Get Log-Likelihood (=vI) with yI step size
			for (iPos = 1; iPos <= pNumber; iPos++) {
                for(indCodon = 1; indCodon <= cNumber; indCodon++) {
					zFP_New[iPos][indCodon] = zFP[iPos][indCodon] + yI * zShift[iPos][indCodon];
				}
			}
		  Get_Log_Likelihood_Only(jSet, DS, zFP_New, vI);

		  //Four major cases
			iFlag = 0;

			if((yI < y2) && (vI < v2)) { // replace (1) with (I)
				y1 = yI; v1 = vI; iFlag = 1;
			}

			if((yI < y2) && (vI > v2)) { // replace (3) with (2) and (2) with (I); (1) the same
				y3 = y2; v3 = v2; y2 = yI; v2 = vI; iFlag = 1;
			}

			if((yI > y2) && (vI > v2)) { // replace (1) with (2) and (2) with (I); (3) the same
				y1 = y2; v1 = v2; y2 = yI; v2 = vI; iFlag = 1;
			}

			if((yI > y2) && (vI < v2)) {  // replace (3) with (I)
				y3 = yI; v3 = vI; iFlag = 1;
			}

			if(iFlag == 0) { // none of the above pertains: break from dummy loop
				break;
			}
			  y31_Length = y3 - y1;
			  y32_Length = y3 - y2;
			  y21_Length = y2 - y1;

            if ((y31_Length < yTol) || (y32_Length < yTol) || (y21_Length < yTol)) {
                break; //we are at the bracket boundary
            }
	 }
        stepSize = y2;
        tLhd_Total = v2;
        v2 = v2;
}

// ============== Rapid zFP refinement
void Hes_Pos_ML_Refine_zFactors(const string strReportFilePath, long jSet, long pAsite, long pNumber,
        long jGeneStartShift, RPFdataSet& DS, vector<double>& zFP_long, vector<double>& zFP_Sigma_long,
		vector<double>& wFP_long, vector<double>& zFP_Refined_long,
		vector<double>& zFP_Sigma_Refined_long,vector<double>& wFP_Refined_long) {

	// Refines z-Factors using Position-Diagonal Hessian Approximation
	// by Marquardt-Like approach
	// Implemented by Michael Pavlov

	int iPrint = 1;


    int nFus = DS.nRPF.size();
    long mRow1 = DS.iORFcodons.size();
	long mRow=mRow1-1;
    int jTot1 = DS.geneEnd.size();
	int jTot=jTot1-1;   //number of genes in data set
    int cNumber = 64;  //Number of codons
	int cNumber1=cNumber+1;
    long rowLength = pNumber * cNumber;
	long rowLength1=rowLength+1;

	// resize main output arrays
	zFP_Refined_long.resize(rowLength1,0.0);
	zFP_Sigma_Refined_long.resize(rowLength1,0.0);
	wFP_Refined_long.resize(rowLength1,0.0);

        DS.pAsite = pAsite;
        DS.pNumber = pNumber;
        DS.jGeneStartShift = jGeneStartShift;

	int pNumber1=pNumber+1;

	long i, j, k, kA, iA, jA;
	string strCodon;
	long  nActive, nRPF_kA=0;
	long jStart =0, jEnd =0;

	int jGene, iPos , nPos, nPos1;
	double tModel_kA=0.0;


	int  indCodon, iCod, indCod, iC, nCodon_Elong;
	double wtFetha, coff_Fy, geneRPF_Elong ;

	vector<long> nCodGeneElong(jTot1,0);
	vector<double> uGeneElong(jTot1,0.0);
	vector<double> tGeneElong(jTot1,0.0);
	vector<double> tGeneElongNew(jTot1,0.0);

    DS.rpfOmegaExper.clear();
	DS.rpfOmegaExper.resize(pNumber1, vector<double>(cNumber1,0.0));

	//
	double t_Gene_Elong;
	double tLhd_delta,tLhd_Total_UB, tLhd_Total;

	vector<double> tLhd_Gene_UB(jTot1);
	vector<double> tLhd_Gene(jTot1);

	double gradActiveNorm2, gradActiveNorm;
	double t_pAsite, stepSize;

	double zPosCodDiff, zDiff2, zFP_PC;
	double alfa, beta;

	//Directions, Gradients  and Progections
	 double vDir_Norm2, gradActive_Norm2, vDir_Norm, gradActive_Norm;
	 double vDir_LogLkh_Derivative, vDir_LogLkh_Derivative_Min;
	 double cos_vDir_Grad, cos_vDirOld_GradNew, vDirGradProjection;

		vector<double> tModel(mRow1,0.0);
		vector<double> rpfModel(mRow1,0.0);
	//
	// Initialize main arrays

    MatrixLong nPosCodInGene(pNumber1, vector<long>(cNumber1,0));
    MatrixDouble zHAT(pNumber1, vector<double>(cNumber1, 0.0));
    MatrixDouble zHAT_Old(pNumber1, vector<double>(cNumber1, 0.0));

    MatrixDouble zFP(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble zFP_Sigma(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble wFP(pNumber1, vector<double>(cNumber1,0.0));

	MatrixDouble zFP_Weight(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble zShift(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble zFP_Old(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble zFP_R(pNumber1, vector<double>(cNumber1,0.0));


    vector<long> indCod1(pNumber1);


    // Hessian Arrays
	MatrixDouble diag_mH_Full(pNumber1, vector<double>(cNumber1,0.0));
	D3_VectorDouble mH_Full(pNumber1, MatrixDouble(cNumber1,vector<double>(cNumber1,0.0)));
	D3_VectorDouble mH_Full_Modified(pNumber1, MatrixDouble(cNumber1,vector<double>(cNumber1,0.0)));
	D3_VectorDouble fyHessian(pNumber1, MatrixDouble(cNumber1,vector<double>(cNumber1,0.0)));
	MatrixLong ijH_A(pNumber1, vector<long>(cNumber1,0));


	int iFlagInvert;
	double  det_Hes;
	vector<double> diag_mH_PS;
	// Load current z-Factors (or  ML Fingerprints) and Transform to matrixes
        k = 0;
    for (iPos = 1; iPos <= pNumber; iPos++) {
        for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
                k = k + 1;
            zFP[iPos][indCodon] = zFP_long.at(k);
            zFP_Sigma[iPos][indCodon] = zFP_Sigma_long.at(k);
            wFP[iPos][indCodon] = wFP_long.at(k);
        }
    }

	// Convert z-factors to v-factors by normalyzing each p-position by z(p,1)
        i = 0; double z1;
    for (iPos = 1; iPos <= pNumber; iPos++) {
             z1 = zFP[iPos][1];
        for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
            zFP[iPos][indCodon] = zFP[iPos][indCodon] / z1;
            zFP_Old[iPos][indCodon] = zFP[iPos][indCodon];
            ijH_A[iPos][indCodon] = 0;
          if (zFP[iPos][indCodon] > 0) { // the position/codon  is active
            i = i + 1;
            ijH_A[iPos][indCodon] = i; // Pos/Codon => active Vector Position Index
          }
        }
    }

    // Active movement space
        nActive = i; int nActive1=nActive+1;
    vector<double> vDir(nActive1,0.0);
	vector<double> zShift_Long(rowLength1);
	vector<double> gradActive(nActive1,0.0);
	vector<long> iGradActive(nActive1,0.0);
	vector<long> ijH_Active(nActive1,0.0);
        double tLhd_Total_Old;

	// Prepare arrays for subspaces
      int nActiveShort = nActive - pNumber; // The dimention of active subspace
	  int nActiveShort1=nActiveShort+1;
	vector<double> gradActive_Short(nActiveShort1,0.0);
	vector<double> vDir_Short(nActiveShort1,0.0);
	vector<double> diag_mH_Short(nActiveShort1,0.0);
	vector<double> grad_PS, vDir_PS;
	MatrixDouble mH_PS;

	// Non-singular Hessian
	MatrixDouble mH_Short(nActiveShort1, vector<double>(nActiveShort1,0.0));

	// Calculate model elongation times from normalized zPC-factors
    for (jGene = 1; jGene <= jTot; jGene++) {
            jStart = DS.geneStart.at(jGene) + DS.jGeneStartShift;
            jEnd = DS.geneEnd.at(jGene);
        for (k = jStart; k <= jEnd - pNumber; k++) {
                t_pAsite = 1; // quazi-time
            for (iPos = 1; iPos <=  pNumber; iPos++) {
                indCodon = DS.iORFcodons.at(iPos + k - 1);
                t_pAsite = t_pAsite * zFP[iPos][indCodon];
            }
                kA = pAsite + k - 1; // the A-site in the original data set
            tModel.at(kA) = t_pAsite;
        }
    }

	double muInitial , muNew , muStandard, coffB, coffZ;

	// Initiate principle  arrays
	MatrixDouble tFetha(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble gradFull(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble rpfOmegaExper(pNumber1, vector<double>(cNumber1,0.0));
	MatrixDouble rpfOmegaModel(pNumber1, vector<double>(cNumber1,0.0));
	MatrixLong nCodPos(pNumber1, vector<long>(cNumber1,0));

	// Prepare the main rpfOmega(p,c) array of RPF fingerprint and additional arrays
        muInitial = 0;
        muStandard = 0;
	for(jGene = 1; jGene <= jTot; jGene++) {
            jStart = DS.geneStart.at(jGene) + DS.jGeneStartShift;
            jEnd = DS.geneEnd.at(jGene);
				nCodon_Elong = 0;
				geneRPF_Elong = 0;
				t_Gene_Elong = 0;
		for(k = jStart; k <= jEnd - pNumber; k++) {
                kA = pAsite + k - 1; // relative index from the A-site
                   nCodon_Elong = nCodon_Elong + 1;
                   nRPF_kA = DS.nRPF[jSet][kA];
                   geneRPF_Elong = geneRPF_Elong + nRPF_kA;
                 t_Gene_Elong = t_Gene_Elong + tModel.at(kA);
			for (iPos = 1; iPos <=  pNumber; iPos++) {
                   indCodon = DS.iORFcodons.at(iPos + k - 1);
                rpfOmegaExper[iPos][indCodon] = rpfOmegaExper[iPos][indCodon] + nRPF_kA;
                nCodPos[iPos][indCodon] = nCodPos[iPos][indCodon] + 1;
			}
		}
        uGeneElong.at(jGene) = geneRPF_Elong;  // RPFs in the inner gene region
        nCodGeneElong.at(jGene) = nCodon_Elong; //Codon Number in the inner gene region
            muInitial = muInitial + geneRPF_Elong / t_Gene_Elong;
            muStandard = muStandard + geneRPF_Elong / nCodon_Elong;
	}

	//Get Log-Likelihood function
	Get_Log_Likelihood(jSet, DS, zFP, tModel,
                        tGeneElong, tLhd_Total, tLhd_Gene, tLhd_Gene_UB, tLhd_Total_UB);

        tLhd_Total_Old = tLhd_Total;
        tLhd_delta = 0.01 * tLhd_Total_Old;

        std::stringstream ssLine;
        std::ofstream ref_Log(strReportFilePath);
	if (iPrint == 1) {
        // OUTPUT for DATASET CODONS
     ssLine.clear(); ssLine.str("");
	 string strLine;

		ssLine<< "DataSet: "<< DS.dataSetName<<std::endl;
		ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "#_of_Genes= "<< jTot<<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "subSet_#= " << jSet<<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "subSet_Name= " << DS.dataSubSetName.at(jSet)<<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "RPF_Total= " << DS.dataSubSet_RPF_Total.at(jSet)<<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "Doubling_Time= " << DS.doublingTime.at(jSet)<<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "pAsite= "<< DS.pAsite << " ;pNumber= " << DS.pNumber<<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "Upper bound of Likelihood= " << tLhd_Total_UB <<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "Initial Likelihood= " << tLhd_Total <<std::endl;
        ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine<< "Iterations to Log-Likelihood Maximum Using Position-Hessian Approximation"
		<<std::endl;
		ref_Log<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

		ssLine << "Iter;   Lhd;    Grag_Norm;   z_Shift;   (vDir,vGrad);  (eDir,vGrad);  cos_vD_vG;" <<
		"stepSize;  (eD_Old,vG_New;  cos_vD_Old_vG_New" <<std::endl;
		ref_Log<<ssLine.str(); ssLine.clear(); ssLine.str("");

		// to console
		ssLine << "Iter; Lhd; Grag_Norm; z_Shift;   stepSize;  cos(vD,vG); vD_der; " <<std::endl;
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
	}

  //Main iterations start HERE =========================================================
        alfa = 3;  // Marquard parameter
        beta = 0;  // Marquard parameter

		int iTerNumber=19;
  for (int Iter = 1; Iter <= iTerNumber; Iter++) {

	// Calculate "tFetha(iPos, indCodon)", "rpfOmegaModel(iPos, indCodon)"
    //    and "fyHessian(iPos, indCodon, iCod)"

    // Zero them first
    for (iPos = 1; iPos <= pNumber; iPos++) {
        for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
              tFetha[iPos][indCodon] = 0.0;
              rpfOmegaModel[iPos][indCodon] = 0.0;
              gradFull[iPos][indCodon] = 0.0;
            for (iCod = 1; iCod <= cNumber; iCod++)  {
                fyHessian[iPos][indCodon][iCod] = 0.0;
            }
		}
    }

    // Fill "tFetha(iPos, indCodon)", "rpfOmegaModel(iPos, indCodon)"
    //      and "fyHessian(iPos, indCodon, iCod)"
        muNew = 0;
	for(jGene = 1; jGene <= jTot; jGene++) {
            jStart = DS.geneStart.at(jGene) + DS.jGeneStartShift;
            jEnd = DS.geneEnd.at(jGene);
		wtFetha = uGeneElong.at(jGene) / tGeneElong.at(jGene); // =U(j)/T(j)
		coff_Fy = wtFetha / tGeneElong.at(jGene); // =U(j)/(T(j)*T(j))
		muNew = muNew + wtFetha;

		// Get "tFetha" and model counts "rpfModel" for gene j
       for (k = jStart; k <= jEnd - pNumber; k++) {
                kA = pAsite + k - 1;
                tModel_kA = tModel.at(kA);
			for (iPos = 1; iPos <= pNumber; iPos++) {
				indCodon = DS.iORFcodons.at(iPos + k - 1);
                tFetha[iPos][indCodon] = tFetha[iPos][indCodon] + tModel_kA;
			}
              rpfModel.at(kA) = wtFetha * tModel_kA; // Expected model RPFs
		}

		// Get model RPF count fingerprint rpfOmegaModel and fyHessian
		for (iPos = 1; iPos <= pNumber; iPos++) {
			for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
                    rpfOmegaModel[iPos][indCodon] = rpfOmegaModel[iPos][indCodon] +
                                wtFetha * tFetha[iPos][indCodon];
                for (iCod = 1; iCod <= cNumber; iCod++)  {
					fyHessian[iPos][indCodon][iCod] = fyHessian[iPos][indCodon][iCod] +
                        coff_Fy * tFetha[iPos][indCodon] * tFetha[iPos][iCod];
				}
            }
		}

		// Zero tFetha for the next gene
        for (iPos = 1; iPos <= pNumber; iPos++) {
			for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
				tFetha[iPos][indCodon] = 0;
            }
        }
	}

	//Get Full and Active Gradient Vectors
        k = 0;
		gradActiveNorm2 = 0;
	for (iPos = 1; iPos <= pNumber; iPos++) {
		for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
					k = k + 1;
					zShift_Long.at(k) = 0.0;
					zShift[iPos][indCodon] = 0.0;
				gradFull[iPos][indCodon] =
					rpfOmegaExper[iPos][indCodon] - rpfOmegaModel[iPos][indCodon];
			if (zFP[iPos][indCodon] > 0) {
				gradFull[iPos][indCodon] = gradFull[iPos][indCodon] / zFP[iPos][indCodon];
			}
				i = ijH_A[iPos][indCodon];
			if (i > 0) { // it is an active vector position
					iGradActive.at(i) = k;
					ijH_Active.at(i) = i;
					gradActive.at(i) = 0;
				if (zFP[iPos][indCodon] > 0) {
					gradActive.at(i) = gradFull[iPos][indCodon];
					gradActiveNorm2 = gradActiveNorm2 + gradActive.at(i) * gradActive.at(i);
				}
			}
		}
	}

    gradActiveNorm = std::sqrt(gradActiveNorm2);
		i = i;
	//gradActive.resize(nActive1);
	//ijH_Active.resize(nActive1);

	//	Check the Projection of the New Gradient on the Old seach directin vDir:
	// 		Note that It Should be close to Zero
        vDirGradProjection = 0;
        vDir_Norm2 = 0;
        gradActive_Norm2 = 0;
    for (iA = 1; iA <= nActive; iA++) {
        vDirGradProjection = vDirGradProjection + vDir.at(iA) * gradActive.at(iA);
      vDir_Norm2 = vDir_Norm2 + vDir.at(iA) * vDir.at(iA);
      gradActive_Norm2 = gradActive_Norm2 + gradActive.at(iA) * gradActive.at(iA);
    }
       vDir_Norm = std::sqrt(vDir_Norm2);
       gradActive_Norm = std::sqrt(gradActive_Norm2);
        vDir_LogLkh_Derivative_Min = 0;
    if (vDir_Norm > 0) {
		vDir_LogLkh_Derivative_Min = vDirGradProjection / vDir_Norm;
	}
	if (gradActive_Norm > 0) {
       cos_vDirOld_GradNew = vDir_LogLkh_Derivative_Min / gradActive_Norm;
	}
    vDir_Norm2 = vDir_Norm2;

	// Construct a Full Block-Hessian Approximation
	for (iPos = 1; iPos <= pNumber; iPos++) {
		for (iC = 1;  iC <= cNumber; iC++) {
				zFP_PC = zFP[iPos][iC];
			if (zFP_PC > 0) {
				for(iCod = iC + 1; iCod <= cNumber; iCod++) {
						mH_Full[iPos][iC][iCod] = 0.0;
					if(zFP[iPos][iCod] > 0) { mH_Full[iPos][iC][iCod] =
						fyHessian[iPos][iC][iCod] / (zFP_PC * zFP[iPos][iCod]);
					  mH_Full_Modified[iPos][iC][iCod] = mH_Full[iPos][iC][iCod];
					  mH_Full[iPos][iCod][iC] = mH_Full[iPos][iC][iCod];
					  mH_Full_Modified[iPos][iCod][iC] = mH_Full[iPos][iC][iCod];
					}
				}
				mH_Full[iPos][iC][iC] =
					(fyHessian[iPos][iC][iC] - rpfOmegaExper[iPos][iC]) / (zFP_PC * zFP_PC);
				diag_mH_Full[iPos][iC] = mH_Full[iPos][iC][iC];
				mH_Full_Modified[iPos][iC][iC] = mH_Full[iPos][iC][iC];
			}
		}
	}

	// Increase the Hessian diagonal to tilt the shift direction towards the gradient
	for (iPos = 1; iPos <= pNumber; iPos++) {
		for (iC = 1;  iC <= cNumber; iC++) {
			mH_Full_Modified[iPos][iC][iC] = mH_Full[iPos][iC][iC] +
				alfa * (diag_mH_Full[iPos][iC]-1) - beta * gradActiveNorm;
		}
	}

		// Make most frequent codon (indexed as  1) inactive at  all positions
		//   of the local context
		for (iPos = 1; iPos <= pNumber; iPos++) {
				i = ijH_A[iPos][1];
			ijH_Active.at(i) = 0;
		}

	// Accelerated solution using the block-diagonal structure of Hessian
		k = 0;
		det_Hes=1;
	for (iPos = 1; iPos <= pNumber; iPos++) {
		// Get Hessian and Gradient for each p-position
		// Get Gradient as one-dimention vector from its pos/codon matrix
		grad_PS.resize(cNumber1,0.0);
			iA = 0;
		for(iCod =2; iCod <= cNumber; iCod++) {
			i = ijH_A[iPos][iCod];
		   if(i > 0) {
				iA = iA + 1;
				grad_PS.at(iA) = gradFull[iPos][iCod];
		   }
		}
			nPos = iA;
			nPos1=nPos+1;
		grad_PS.resize(nPos1);

		// Construct Actibe Hessian for position iPos according to grad_PS
		mH_PS.clear(); vDir_PS.clear();
		mH_PS.resize(nPos1, vector<double>(nPos1,0.0));
		vDir_PS.resize(nPos1,0.0);

			iA = 0;
		for (iC = 2; iC <= cNumber; iC++) {
                i = ijH_A[iPos][iC];
            if (i > 0) {
                    iA = iA + 1;
                    jA = 0;
                for (iCod = 2; iCod <= cNumber; iCod++) {
                        j = ijH_A[iPos][iCod];
                    if (j > 0) {
                            jA = jA + 1;
                        mH_PS[iA][jA] = mH_Full_Modified[iPos][iC][iCod];
                    }
                }
            }
		}

		diag_mH_PS.clear();
		diag_mH_PS.resize(nPos1,0.0);

		// Solve  (mH_PS+(alfa-1)*diag_mH_PS)*vDir=gradActive
            VB_Gauss_Solve(1,  mH_PS, vDir_PS, grad_PS);

		// restore gradActive short and vDir_Short
		// Change  vDir_Short into the oposite, assent direction
		//      in the full active space
		for(i = 1; i <= nPos; i++) {
				k = k + 1;
			vDir_Short.at(k) = -vDir_PS.at(i);
			gradActive_Short.at(k) = grad_PS.at(i);
		}
	}
	// Get vDir_Short projection on gradient gradActive_Short
        vDirGradProjection = 0.0;
        vDir_Norm2 = 0.0;
        gradActive_Norm2 = 0.0;
    for(iA = 1; iA <= nActiveShort; iA++) {
        vDirGradProjection = vDirGradProjection + vDir_Short.at(iA) * gradActive_Short.at(iA);
      vDir_Norm2 = vDir_Norm2 + vDir_Short.at(iA) * vDir_Short.at(iA);
      gradActive_Norm2 = gradActive_Norm2 + gradActive_Short.at(iA) * gradActive_Short.at(iA);
    }
        vDir_Norm = std::sqrt(vDir_Norm2);
        gradActive_Norm = std::sqrt(gradActive_Norm2);
       vDir_LogLkh_Derivative = vDirGradProjection / vDir_Norm;
       cos_vDir_Grad = vDir_LogLkh_Derivative / gradActive_Norm;

    if(cos_vDir_Grad >= 0) {
            stepSize = 0.5 * tLhd_delta / vDir_LogLkh_Derivative;
          if (std::abs(stepSize) > 0.8) {stepSize = 0.8;}
				stepSize = 2;
        if(std::abs(vDirGradProjection) < 10000) { // Start reducing M-parameters
            alfa = alfa * 0.1;
            beta = beta * 0.5;
        }
    }
    else { // vDir direction is bad, use gradActive direction instead
			vDirGradProjection = 0;
			vDir_Norm2 = 0;
		for(iA = 1; iA <= nActiveShort; iA++) {
			vDir_Short.at(iA) = gradActive_Short.at(iA);
			vDirGradProjection = vDirGradProjection +
								vDir_Short.at(iA) * gradActive_Short.at(iA);
			vDir_Norm2 = vDir_Norm2 + vDir_Short.at(iA) * vDir_Short.at(iA);
		}
			vDir_Norm =std::sqrt(vDir_Norm2);
			vDir_LogLkh_Derivative = vDirGradProjection / vDir_Norm;

			//Try to increase Likelihood in the gradient direction
			stepSize = tLhd_Total / vDirGradProjection;
			stepSize = stepSize;
	}

	// Expand  vDir and zShift_Long to the full zPC-parameter space
        iA = 0;
    for (i = 1; i <= nActive; i++) {
            vDir.at(i) = 0;
        if (ijH_Active.at(i) > 0) {
            iA = iA + 1;
            vDir.at(i) = vDir_Short.at(iA);
        }
            k = iGradActive.at(i);
        zShift_Long.at(k) = vDir.at(i);
    }

	// record zFP_Old and transform vector zShift_long into matrix zShift
        k = 0;
    for (iPos = 1; iPos <= pNumber; iPos++) {
		for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
            zFP_Old[iPos][indCodon] = zFP[iPos][indCodon];
                k = k + 1;
            zShift[iPos][indCodon] = zShift_Long.at(k);
        }
    }

	// Refine stepSize by minimizing log-likelihood along vDir
    RefineStepSize_New(jSet, DS, zFP, zShift, stepSize);

	// Do the actual zFP-parameter shift
    for (iPos = 1; iPos <= pNumber; iPos++) {
		for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
           zFP[iPos][indCodon] = zFP_Old[iPos][indCodon] + stepSize * zShift[iPos][indCodon];
        }
	}

	// Get new Log-Likelihood function, new model times and gene times
	Get_Log_Likelihood(jSet, DS, zFP, tModel,
                        tGeneElong, tLhd_Total, tLhd_Gene, tLhd_Gene_UB, tLhd_Total_UB);
        tLhd_Total_Old = tLhd_Total;

	// Calculate the shift NORM
        zDiff2 = 0;
    for (iPos = 1; iPos <= pNumber; iPos++) {
		for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
         zPosCodDiff = (zFP[iPos][indCodon] - zFP_Old[iPos][indCodon]);
            zDiff2 = zDiff2 + zPosCodDiff * zPosCodDiff;
        }
    }
        zDiff2 = zDiff2 / rowLength;

	//Case of printing out
	iPrint=1;
	if (iPrint == 1) {
		ssLine<< Iter<< "; " << tLhd_Total<< "; " <<  gradActive_Norm << "; "
			<< std::sqrt(zDiff2)<< "; " << vDirGradProjection <<"; " << vDir_LogLkh_Derivative << "; "
			<< cos_vDir_Grad << "; " << stepSize <<"; " << vDir_LogLkh_Derivative_Min << "; "
			<< cos_vDirOld_GradNew << "; "<< std::endl;
        ref_Log<<ssLine.str(); ssLine.clear(); ssLine.str("");

        // to console
		ssLine<< Iter<< "; " << tLhd_Total<< "; " << gradActive_Norm << "; "
			<< std::sqrt(zDiff2)<< "; " << stepSize <<"; "
			<< cos_vDir_Grad << "; "<< cos_vDirOld_GradNew << "; "<< std::endl;
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
	}
	if (stepSize == 0) {break;}
  }
//Main iteration ends  HERE ==================================================

	// Finalyze z-factors
    long kTest, muWeight;
    double zPosAver, zAverTot, wCodPos;
    vector<double> zPosAverage(pNumber1);

    MatrixDouble muCoff(pNumber1, vector<double>(cNumber1, 0.0));
    vector<double> muCoff_Cod_Sum(pNumber1);
    MatrixDouble diagHM1_Full(pNumber1, vector<double>(cNumber1,0.0));
    MatrixDouble mHM1;

	// Zero muCoff=Effective weights of Codons
	for (iPos = 1; iPos <= pNumber; iPos++) {
		for (indCod = 1;  indCod <= cNumber; indCod++) {
				muCoff[iPos][indCod] = 0;
		}
	}
		muWeight = 0;
		muNew = 0;
		muStandard = 0;
	for(jGene = 1; jGene <= jTot; jGene++) {
            jStart = DS.geneStart.at(jGene) + jGeneStartShift;
            jEnd = DS.geneEnd.at(jGene);
		wtFetha = uGeneElong.at(jGene) / tGeneElong.at(jGene); // =U(j)/T(j)
		coff_Fy = wtFetha / tGeneElong.at(jGene); // =U(j)/(T(j)*T(j))
		muNew = muNew + wtFetha;
		muWeight = muWeight + nCodGeneElong.at(jGene) * wtFetha;
		muStandard = muStandard + uGeneElong.at(jGene) / nCodGeneElong.at(jGene);

		// Zero nPosCodInGene
		for (iPos = 1; iPos <= pNumber; iPos++) {
			for (indCod = 1;  indCod <= cNumber; indCod++) {
				nPosCodInGene[iPos][indCod] = 0;
			}
		}

		// Obtain nPosCodInGene
		for(k = jStart; k <= (jEnd - pNumber); k++) {
			for (iPos = 1; iPos <= pNumber; iPos++) {
				indCod = DS.iORFcodons.at(iPos + k - 1);
			  nPosCodInGene[iPos][indCod] = nPosCodInGene[iPos][indCod] + 1;
			}
		}

		// Use nPosCodInGene to get muCoff=Effective weights of Codons
	   for (iPos = 1; iPos <= pNumber; iPos++) {
			for (indCod = 1;  indCod <= cNumber; indCod++) {
				muCoff[iPos][indCod]= muCoff[iPos][indCod] +
							nPosCodInGene[iPos][indCod] * wtFetha;
			}
		}
	}

		// Finalize  muCoff=Effective weights of nucleotides
		for (iPos = 1; iPos <= pNumber; iPos++) {
			for (indCod = 1;  indCod <= cNumber; indCod++) {
				muCoff[iPos][indCod] =  muCoff[iPos][indCod] / muWeight;
			}
		}

	// Get position Hessian and Gradient
		k = 0;
		zAverTot = 1;

	// get the dimentions of the position Hessian
	for (iPos = 1; iPos <= pNumber; iPos++) {
			iA = 0;
		for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
			i = ijH_A[iPos][indCodon];
		  if (i > 0) {
				iA = iA + 1;
		  }
		}
			nPos = iA - 1;

		// get Possition Hessian with one position excluded
		mH_PS.clear(); mHM1.clear();
		mH_PS.resize(iA, vector<double> (iA,0.0));
        mHM1.resize(iA, vector<double> (iA,0.0));

		// prepare the diagonal
		for (indCodon = 1;  indCodon <= cNumber; indCodon++) {
				diagHM1_Full[iPos][indCodon] = 0;
		}

		for(kTest = 1; kTest <=2; kTest++) {
				iA = 0;
			for(indCodon = 1;  indCodon <= cNumber; indCodon++) {
			   if(!(indCodon == kTest)) {
						i = ijH_A[iPos][indCodon];
					if(i > 0) {
						   iA = iA + 1;
						   jA = 0;
						for(iCod = 1;  iCod <= cNumber; iCod++) {
							if (!(iCod == kTest)) {
									j = ijH_A[iPos][iCod];
								if (j > 0) {
									jA = jA + 1;
								   mH_PS[iA][jA] = mH_Full[iPos][indCodon][iCod];
								}
							}
						}
					}
				}
			}

			//invert the curtailed position Hessian
			VB_LUPA_Invert(1, mH_PS,  mHM1);

				iA = 0;
			for(indCodon = 1;  indCodon <= cNumber; indCodon++) {
				if(!(indCodon == kTest)) {
						i = ijH_A[iPos][indCodon];
					if(i > 0) {
						   iA = iA + 1;
						diagHM1_Full[iPos][indCodon] =
							diagHM1_Full[iPos][indCodon] + mHM1[iA][iA];
						iA = iA;
					}
				}
			}
		}

		for(kTest = 1; kTest <=2; kTest++) {
			diagHM1_Full[iPos][kTest] = 2 * diagHM1_Full[iPos][kTest];
		}

		// Get position averaged z-factors
		  zPosAver = 0;
		  wCodPos = 0;
		for(indCod = 1;  indCod <= cNumber; indCod++) {
			wCodPos = wCodPos + muCoff[iPos][indCod];
			zPosAver = zPosAver + zFP[iPos][indCod] * muCoff[iPos][indCod];
		}
		  zPosAverage.at(iPos) = zPosAver / wCodPos;
		  zAverTot = zAverTot * zPosAverage.at(iPos);
	}

	// Normalization Coffs
		coffB = muNew / muStandard;
		coffB = zAverTot * coffB;
		coffZ = exp(std::log(coffB) / pNumber);

	// Finalize the refined z-factors
		k = 0;
	for (iPos = 1; iPos <= pNumber; iPos++) {
			zPosAver = zPosAverage.at(iPos);
			//coffZ = 1
		for(indCod = 1;  indCod <= cNumber; indCod++) {
				k = k + 1;
			zFP_Refined_long.at(k) = coffZ * zFP[iPos][indCod] / zPosAver;
			zFP_Sigma_Refined_long.at(k) = coffZ * std::sqrt(std::abs(diagHM1_Full[iPos][indCod])) / zPosAver;
			wFP_Refined_long.at(k) = muCoff[iPos][indCod];
		}
	}
	ref_Log.close();
}
// ===================================================================================================
void Print_vFP_Full_As_MatrixB(const string strPrintOutputPath, string strText, const string strNorm,
    const int iPrint_rpfOmega, int iPrint_Col_Corr, int pC_First, int pC_Last, RPFdataSet& DS,
    const long jSet, vector<double>& vFP_Full, vector<double>& vFP_Sigma_Full, vector<double>& wFP_Full) {
	//
	// Prints  vector type vFP (vector FingerPrint)  as  zFP(iPos, jCod) Matrix
	// Rows correspond to position (iPos) and indexed by codon index (jCod) fron 1 to 64
	//
	// strText: the header to be printed as is
	// strNorm: takes three values: "NATIVE", "NORMALIZED", "ZF(P,1)=1" ,
	//            the later means that first element of the row is 1
	// iPrint_rpfOmega: if = 1 then prints the dataset experimental FingerPrint
	// iPrint_Col_Corr: if = 1 then report column correlations in zFP(iPos, jCod) Matrix
	//  pC_First: 	print the rows up to this as 1
	//  pC_Last :	print the rows afer this as 1
	//  DS: 		data set. Used to print additional information about the dataset
	//  jSet:   	specifies the subset of the dataset we are working with
	//  vFP_Full:   		zFP(iPos, jCod) zFP coefficients compresed to vector
	//  vFP_Signma_Full:   	Errors zFP coefficients compresed to vector
	//  wFP_Full:   		Frequency of codon jCod at local position iPos as
	//							seen by the translating ribosome
	// Implemented by Michael Pavlov
	//

	long i, j, k;
	long pNumber, pAsite, nCodon, nLength;
	long iPos, jCod, jTot;

	string  strAA, strCodon, strWeighted;

	int iFlagWeighted=0;
    nLength = vFP_Full.size()-1;
    nCodon = 64;
	int nCodon1 = nCodon+1;
    pNumber = nLength / nCodon;
    pAsite = DS.pAsite;
	int pNumber1=pNumber+1;
    jTot = DS.jTot;

	if(pC_First == 0) {pC_First = 1;}
	if(pC_Last == 0) {pC_Last = pNumber;}

    strWeighted = "Non-Weighted ";
	if (iFlagWeighted == 1) {strWeighted = "Weighted ";}
   strWeighted = "; Weighted Row Average and Dispersion";

	MatrixDouble mFP(pNumber1, vector<double>(nCodon1,0.0));
	MatrixDouble mFP_Sigma(pNumber1, vector<double>(nCodon1,0.0));
	MatrixDouble mFPW(pNumber1, vector<double>(nCodon1,0.0));

	MatrixDouble mFPexp(pNumber1, vector<double>(nCodon1,0.0));




	vector<double> mFPW_CodonSum(nCodon1,0.0);

	// convert long vectors to  matrices
        k = 0;
    for (iPos = 1; iPos <= pNumber; iPos++) {
        for (jCod = 1; jCod <= nCodon; jCod++) {
                k = k + 1;
			mFP[iPos][jCod] = vFP_Full.at(k);
            mFP_Sigma[iPos][jCod] = vFP_Sigma_Full.at(k);
            mFPW[iPos][jCod] = wFP_Full.at(k);
        }
    }

	// modify mFP matrices to be 1 at positions outside [pFirst, pLast]
	for (iPos = 1; iPos <= pC_First - 1; iPos++) {
        for (jCod = 1; jCod <= nCodon; jCod++) {
            if (mFP[iPos][jCod]> 0.0) {mFP[iPos][jCod] = 1;}
            mFP_Sigma[iPos][jCod] = 0.0;
		}
    }

    for (iPos =pC_Last + 1; iPos <= pNumber; iPos++) {
        for (jCod = 1; jCod <= nCodon; jCod++) {
            if (mFP[iPos][jCod]> 0.0) {mFP[iPos][jCod] = 1;}
            mFP_Sigma[iPos][jCod] = 0.0;
		}
    }

	// find zFP column averages, dispersions and correlations
	vector<double> mFP_ColAver, mFP_ColSigma;
	MatrixDouble mPearsonColCorr;
    Get_PCorr_Cols_mA_New(pC_First, pC_Last, mFP, mFPW, mFP_ColAver, mFP_ColSigma, mPearsonColCorr);

	// find zFP column with minimal variation
	int jCodFlat;
	double sigmaMin=0.0;
        jCodFlat = 1;
        sigmaMin = mFP_ColSigma.at(1);
    for (jCod = 1; jCod <= 40; jCod++) {
        if ((mFP_ColSigma.at(jCod) < sigmaMin) && (mFP_ColSigma.at(jCod) > 0)) {
            sigmaMin = mFP_ColSigma.at(jCod);
            jCodFlat = jCod;
        }
    }

	// calculate weighted row averages, sigmas and correlations
        vector<double> mFP_RowAver, mFP_RowSigma;
        MatrixDouble mPearsonRowCorr;
        long iCod_First=1;
	Get_PCorr_Rows_mA_New(iCod_First, nCodon, mFP, mFPW, mFP_RowAver, mFP_RowSigma,
                       mPearsonRowCorr);

	//Print  matrix and the codon information
    // Print codon Info
	std::ofstream zFP_Record(strPrintOutputPath);
	std::stringstream ssLine;
	ssLine.clear(); ssLine.str("");  //Clear the String Stream Buffer

    ssLine << "DataSet: ;" << DS.dataSetName << " ; #_Genes= ;" << jTot << "; subSet# ;" << jSet
		<< "; subSetName: ;" << DS.dataSubSetName.at(jSet)
		<< "; dbl_Time= ;" << DS.doublingTime.at(jSet)
		<< "; RPF_Tot= ;" << DS.dataSubSet_RPF_Total.at(jSet) << std::endl;
	zFP_Record<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine << strNorm << " :" << strText << strWeighted<< std::endl; // the title
	zFP_Record<<ssLine.str();
	std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine << "Codon_Index: " << "; " << "; ";
	for (j = 1; j <= nCodon; j++) {ssLine << j << "; ";}
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

		ssLine << "Number: "<< "; " << "; ";
	for (j = 1; j <= nCodon; j++) {ssLine << wFP_Full.at(j) << "; ";}
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

        ssLine << "AA_name: " << " ;" << " ;";
	for (j = 1; j <= nCodon; j++) {
		ssLine <<  DS.geneCodeTableOrdered.at(j) << " ;";
		}
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str(" ");

		ssLine << "Cod_name: " << "; Aver " << "; Sigma ";
	for (j = 1; j <= nCodon; j++) {
		ssLine <<  DS.geneCodeTableOrdered.at(j).substr(0,5) << " ;";
		}
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

	// print mFP and its averages
    for(i = 1; i <= pNumber; i++) {
			k=i;
			if(i== pAsite){k=-i;}
		if(strNorm == "NORMALIZED"){
			ssLine<< k << " ;"
			<< (mFP_RowAver.at(i) / mFP_RowAver.at(i)) <<" ;"
			<< (mFP_RowSigma.at(i) / mFP_RowAver.at(i)) << " ;";
			}
		if(strNorm == "NATIVE"){
			ssLine<< k << " ;"
			<< mFP_RowAver.at(i)  <<" ;" << mFP_RowSigma.at(i)<< " ;";
			}
		if(strNorm == "ZF(P,1)=1"){
			ssLine<< k << " ;"
			<< (mFP_RowAver.at(i)/mFP[i][1]) <<" ;"
			<< (mFP_RowSigma.at(i)/mFP[i][1]) << " ;";
			}

		for(j = 1; j <= nCodon; j++) {
			if(strNorm == "NORMALIZED"){
				ssLine<< (mFP[i][j] / mFP_RowAver.at(i)) << " ;";}
			if(strNorm == "NATIVE"){ssLine<< mFP[i][j]  << " ;";}
			if(strNorm == "ZF(P,1)=1"){ssLine<< (mFP[i][j] /mFP[i][1]) << " ;";}
        }
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
	}

		ssLine << "colAver:  "<< " ;" << " ;";
	for (j = 1; j <= nCodon; j++) {
            ssLine << mFP_ColAver.at(j) << " ;";
    }
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

		ssLine << "colSigma:  "<< " ;" << " ;";
	for (j = 1; j <= nCodon; j++) {ssLine << mFP_ColSigma.at(j) << " ;";}
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    // print Matrix of Codon Statistics
        ssLine<< "Effective Number of Codons Used= Codon Weights" << std::endl;
    for(i = 1; i <= pNumber; i++) {
			ssLine<< i <<  " ;" << " ;";
        for(j = 1; j <= nCodon; j++) {
            ssLine<< mFPW[i][j]<<" ;";
        }
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
    }

    // print Sigma Matrix of tFP Statistics
        ssLine<< "Sigma Errors in tFP/zFP factors: " << std::endl;
    for(i = 1; i <= pNumber; i++) {
			ssLine<< i <<  " ;" << " ;";
        for(j = 1; j <= nCodon; j++) {
			if(strNorm == "NORMALIZED"){
				ssLine<< (mFP_Sigma[i][j] / mFP_RowAver.at(i)) << " ;";}
			if(strNorm == "NATIVE"){ssLine<< mFP_Sigma[i][j]  << " ;";}
			if(strNorm == "ZF(P,1)=1"){
				ssLine<< (mFP_Sigma[i][j] /mFP[i][1]) << " ;";}
        }
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
    }

    // print rpfOmega fingerprint counts
    if (iPrint_rpfOmega == 1) {
        ssLine << "rpfOmega fingerprint counts" << std::endl;
            zFP_Record<<ssLine.str();
            std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
		for(i = 1; i <= pNumber; i++) {
				ssLine<< i <<  " ;" << " ;";
			for(j = 1; j <= nCodon; j++) {
				ssLine << DS.rpfOmegaExper[i][j] << " ;";
			}
			ssLine<< std::endl;
		  zFP_Record<<ssLine.str();
		  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
		}
	}
	// print correlations of mFP matrix rows
	ssLine<<" Pearson row correlation matrix " << std::endl;
		zFP_Record<<ssLine.str();
        std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
	for(i = 1; i <= pNumber; i++) {
			ssLine<< i <<  " ;" << " ;";
		for(j = 1; j <= pNumber; j++) {
            ssLine << mPearsonRowCorr[i][j]<< " ;";
        }
		ssLine<< std::endl;
	  zFP_Record<<ssLine.str();
	  std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
    }
}

// ==============================================================================
void Get_PCorr_Rows_mA_New(long& jCol_Start, long& jCol_End, const MatrixDouble& mA,
								const MatrixDouble& wA, vector<double>& mA_RowAver,
									vector<double>& mA_RowSigma, MatrixDouble& mPC) {
	//
	// It Gets Pearson correlation between rows of matrix A
	//
	// Implemented by Michael Pavlov
	//
	long i, j, k,  mRow, nCol, mRow1, nCol1;
		mRow1 = mA.size();
		nCol1 = mA[1].size();
		mRow=mRow1-1;
		nCol=nCol1-1;
	if (jCol_Start == 0) {jCol_Start = 1;}
	if (jCol_End == 0) {jCol_End = nCol;}

		mA_RowAver.clear(); mA_RowSigma.clear(); mPC.clear();
		mA_RowAver.resize(mRow1,0.0); mA_RowSigma.resize(mRow1,0.0);
		mPC.resize(mRow1, vector<double>(mRow1,0.0));
		double rowAver, rowAver2, wRow;

	// get row averages and their sigmas for matrix mA
		for(i = 1; i <= mRow; i++) {
				rowAver = 0.0;
				rowAver2 = 0.0;
				wRow = 0.0;
			for(j = jCol_Start; j <= jCol_End; j++) {
				wRow = wRow + wA[i][j]; 			//weighted sum for row iRow
				rowAver = rowAver + mA[i][j] * wA[i][j];
				rowAver2 = rowAver2 + mA[i][j] * mA[i][j] * wA[i][j];
			}
				mA_RowAver.at(i) = rowAver / wRow; //row average
			mA_RowSigma.at(i) = rowAver2 / wRow;
			mA_RowSigma.at(i) =
				std::sqrt(std::abs(mA_RowSigma.at(i) - mA_RowAver.at(i) * mA_RowAver.at(i)));
		}
	//
	// calculate row  correlations between  matrices mA and mB
	   double pCorr, wGT;
	for(i = 1; i <= mRow; i++) {
		for(k = i; k<= mRow; k++) {
				pCorr = 0;
				wRow = 0;
			for(j = jCol_Start; j <= jCol_End; j++) {
				wGT = std::sqrt(wA[i][j] * wA[k][j]); //compozed weight
				wRow = wRow + wGT;
				pCorr = pCorr +
					wGT * (mA[i][j] - mA_RowAver.at(i)) * (mA[k][j]  - mA_RowAver.at(k));
			}
					mPC[i][k] = pCorr / wRow;
				if(std::abs(mA_RowSigma.at(i)) > 0 && std::abs(mA_RowSigma.at(k)) > 0) {
					mPC[i][k] = mPC[i][k] / (mA_RowSigma.at(i) * mA_RowSigma.at(k));
				}
					mPC[k][i] = mPC[i][k];
		}
	}
}

// =============================================================================
void Get_PCorr_Cols_mA_New(long iRow_Start, long iRow_End, const MatrixDouble& mA,
		const MatrixDouble& wA, vector<double>& mA_ColAver,
								vector<double>& mA_ColSigma, MatrixDouble& mPC) {
	//
	// It Gets Pearson correlation between rows of matrix A
	//
	// Implemented by Michael Pavlov
	//
	long i, j, k, iRow, jCol, mRow, nCol, mRow1, nCol1;
    mRow1 = mA.size();
    nCol1 = mA[1].size();
    mRow=mRow1-1;
    nCol=nCol1-1;
		mA_ColAver.clear(); mA_ColSigma.clear(); mPC.clear();
		mA_ColAver.resize(nCol1,0.0); mA_ColSigma.resize(nCol1,0.0);
		mPC.resize(nCol1, vector<double>(nCol1,0.0));
	double colAver, colAver2, wCol;

	//get row averages and their sigmas for matrix mA
    for(j = 1; j <= nCol; j++) {
            colAver = 0;
            colAver2 = 0;
            wCol = 0;
        for(i = iRow_Start; i <= iRow_End; i++) {
            wCol = wCol + wA[i][j]; //Sum of weights
            colAver = colAver + mA[i][j] * wA[i][j];
            colAver2 = colAver2 + mA[i][j] * mA[i][j] * wA[i][j];
        }
            mA_ColAver.at(j) = 0;
            mA_ColSigma.at(j) = 0;
        if(wCol > 0) {
            mA_ColAver.at(j) = colAver / wCol; // codon averaged FP for position iPos
            mA_ColSigma.at(j) = colAver2 / wCol;  //codon averaged FP2 for position iPos
            mA_ColSigma.at(j) =
			   std::sqrt(std::abs(mA_ColSigma.at(j) - mA_ColAver.at(j) * mA_ColAver.at(j)));
        }
    }

	//calculate  matrices mA column correlations
	double pCorr, wGT;
	for(j = 1; j <= nCol; j++) {
		for(k = j; k <= nCol; k++) {
				pCorr = 0;
				wCol = 0;
			for (i = iRow_Start; i <= iRow_End; i++) {
				wGT = std::sqrt(wA[i][j] * wA[i][k]); // compozed weight
				wCol = wCol + wGT;
				pCorr = pCorr +
				  wGT * (mA[i][j] - mA_ColAver.at(j)) * (mA[i][k] - mA_ColAver.at(k));
			}
				mPC[j][k] = 0;
			if((wCol > 0) && (std::abs(pCorr) > 0)) {
				mPC[j][k] = pCorr / wCol;
				mPC[j][k] = mPC[j][k] / (mA_ColSigma.at(j) * mA_ColSigma.at(k));
			}
				mPC[k][j] = mPC[j][k];
		}
	}
}

// =============================================================================================
void Report_R2_zFP_Statistics(string strOutPutFile, long jSet, string strMode, string strMark,
			string strWeight, string strInfo, double tol_R2, int& pC_First, int& pC_Last,
            RPFdataSet& DS, vector<double>& zFP_long, vector<double>& zFP_Sigma_long) {
	//
	// Prints-out Exper-Model Correlaion statistics for Dataset Genes
	//
	// Implemented by Michael Pavlov
	//

	long i, j, k, kA, indCodon, iPos, kRow, kCodon;
	long  nExp, iFlagMode;
	string strGeneName;
	long cNumber, pAsite, pNumber;

	long mRow, jStart, jEnd, jTot, nCodon_Elong, nCodon_Elong1, nRPF;
	vector<double> xExperGene, yModelGene, wExperGene;
	double R2, corrME, corrME_Pearson, yModelAver, xExperAver, sigmaYM2, sigmaXE2;
	double tModelCodon, tModelGene;

	// Likelyhood variables
	vector<double>  geneRPFlnT;
	double tLhd_Total, gene_RPF_ln_tModel, tLhd_Total_UB, gene_RPF_ln_RPF;

	long geneRPF_Elong, jGene, nLength;
	double tModelGenePerCodon, tExperGenePerCodon, tExperCodon, tExperGene;
	double densGeneExper, densGeneModel;

	double t_pAsite, t_Sigma, t_Sigma2;
	double wGene, wTotal, averPRcorrW, averPRcorr;

    long jTot1 = DS.geneEnd.size();
	jTot=jTot1-1;				//number of genes in data set
    long mRow1 = DS.iORFcodons.size();
	mRow=mRow1-1;
    nExp = mRow;
    cNumber = 64; 				//Number of codons

    pNumber = DS.pNumber; 		// number of positions
    pAsite = DS.pAsite; 		//A-site position
	long pNumber1=pNumber+1;
	long cNumber1=cNumber+1;
	//
	//  Position range (from pC_First to pC_Last) to be taken into account
	//  for Model Time Calculations
	 if (pC_First == 0) {pC_First = 1;}
	 if (pC_Last == 0) {pC_Last = pNumber;}

	//  convert zFP long vectors into zFP matrix
	MatrixDouble zFP(pNumber1, vector<double>(cNumber1, 0.0));
	MatrixDouble zFP_Sigma(pNumber1, vector<double>(cNumber1, 0.0));
	MatrixDouble zFP_Sigma2(pNumber1, vector<double>(cNumber1, 0.0));

        k = 0;
    for(i = 1; i <=  pNumber; i++) {
        for(j = 1; j <= cNumber; j++) {
                k = k + 1;
            zFP[i][j] = zFP_long.at(k);
            zFP_Sigma[i][j] = zFP_Sigma_long.at(k);
          if(zFP_Sigma[i][j] > 0){
            zFP_Sigma[i][j] = zFP_Sigma[i][j] / zFP[i][j];
            zFP_Sigma2[i][j] = zFP_Sigma[i][j] * zFP_Sigma[i][j];
          }
        }
    }
	vector<double> xExper(mRow1, 0.0), yModel(mRow1, 0.0), wExper(mRow1, 0.0);


	vector<DS_Corr> CG(jTot1);
	vector<double> nGeneRPF(jTot1,0.0), tLhd_Gene(jTot1,0.0), tLhd_Gene_UB(jTot1,0.0);

	//  calculate zFactor model times
	vector<double> tModelGene_Elong(jTot1,0.0);
	vector<long> nGeneCod_Elong(jTot1,0.0);
            kRow = 0;
            tLhd_Total = 0.0;
            tLhd_Total_UB = 0.0;
    for(jGene = 1; jGene <= jTot; jGene++) {
        jStart = DS.geneStart.at(jGene) + DS.jGeneStartShift;
        jEnd = DS.geneEnd.at(jGene);

                tModelGene = 0.0;
                tExperGene = 0.0;
                geneRPF_Elong = 0.0;
                gene_RPF_ln_tModel = 0.0;
                gene_RPF_ln_RPF = 0.0;
            nCodon_Elong = jEnd - jStart + 1;
			nCodon_Elong1 =nCodon_Elong +1;
		xExperGene.clear(); yModelGene.clear();wExperGene.clear();
		xExperGene.resize(nCodon_Elong1,0.0);
		yModelGene.resize(nCodon_Elong1,0.0);
		wExperGene.resize(nCodon_Elong1,0.0);

            kCodon = 0;
        for(k = jStart; k <= jEnd - pNumber; k++) {
                kA = pAsite + k - 1; 		//A-site in the original data set
                nRPF = DS.nRPF[jSet][kA];

                t_pAsite = 1; // quasi-time  with a particular codon in the A-site
                t_Sigma2 = 0;
            for(iPos = pC_First; iPos <= pC_Last; iPos++) {
                indCodon = DS.iORFcodons.at(kA + iPos - pAsite);
                t_pAsite = t_pAsite * zFP[iPos][indCodon];
                t_Sigma2 = t_Sigma2 + zFP_Sigma2[iPos][indCodon];
            }
                t_Sigma = t_pAsite * std::sqrt(t_Sigma2);

              kCodon = kCodon + 1; // current codon in gene arrays
                yModelGene.at(kCodon) = t_pAsite;
                xExperGene.at(kCodon) = nRPF;

            if(t_pAsite > 0){gene_RPF_ln_tModel = gene_RPF_ln_tModel + nRPF * std::log(t_pAsite);}
            if(nRPF > 1) {gene_RPF_ln_RPF = gene_RPF_ln_RPF + nRPF * std::log(nRPF);}

              geneRPF_Elong = geneRPF_Elong + nRPF;
              tModelGene = tModelGene + t_pAsite;
		}
            tLhd_Gene.at(jGene) = 0.0;
            tLhd_Gene_UB.at(jGene) = 0.0;
        if(kCodon > 0) {
            if(tModelGene > 0) {
				tLhd_Gene.at(jGene) = gene_RPF_ln_tModel -
								geneRPF_Elong * std::log(tModelGene / kCodon);
			}
            if(geneRPF_Elong > 0) {
				tLhd_Gene_UB.at(jGene) = gene_RPF_ln_RPF -
							geneRPF_Elong * std::log(geneRPF_Elong / kCodon);
			}
		}
            tLhd_Total = tLhd_Total + tLhd_Gene.at(jGene);
            tLhd_Total_UB = tLhd_Total_UB + tLhd_Gene_UB.at(jGene);

        nGeneRPF.at(jGene) = geneRPF_Elong;
        tModelGene_Elong.at(jGene) = tModelGene;
        nGeneCod_Elong.at(jGene) = kCodon;

		// Transform to pausing scores
            tModelGenePerCodon = 1;
            tExperGenePerCodon = 1;
        if(kCodon > 0){
            tModelGenePerCodon = tModelGene / kCodon;
            tExperGenePerCodon = geneRPF_Elong / kCodon;
        }
            kCodon = 0;
		for(k = jStart; k <= jEnd - pNumber; k++) {
                kRow = kRow + 1;
                kCodon = kCodon + 1; // current codon in gene
              yModelGene.at(kCodon) = yModelGene.at(kCodon) / tModelGenePerCodon;
              xExperGene.at(kCodon) = xExperGene.at(kCodon) / tExperGenePerCodon;

                 wExperGene.at(kCodon) = 1;
              if(strWeight =="WEIGHTED") {wExperGene.at(kCodon) = tExperGenePerCodon;}

                yModel.at(kRow) = yModelGene.at(kCodon);
                xExper.at(kRow) = xExperGene.at(kCodon);
                wExper.at(kRow) = wExperGene.at(kCodon);
        }

        // Get Exper-Model pausing score correlations
        Get_X_Y_Correlation(1, kCodon, xExperGene, yModelGene, wExperGene,
                            xExperAver, sigmaXE2, yModelAver, sigmaYM2, corrME);
            R2 = 0;
        if(sigmaYM2 > 0 && sigmaXE2 > 0) {R2 = corrME * corrME / (sigmaYM2 * sigmaXE2);}
            corrME_Pearson = std::sqrt(R2);

        // record correlations for a gene
			CG.at(jGene).Name =DS.geneName.at(jGene);
            CG.at(jGene).setA_Aver = xExperAver;
			CG.at(jGene).setA_Sigma = std::sqrt(sigmaXE2);
			CG.at(jGene).setB_Aver = yModelAver;
			CG.at(jGene).setB_Sigma = std::sqrt(sigmaYM2);
            CG.at(jGene).setAB_Corr_Raw = corrME;
			CG.at(jGene).setAB_Corr_Pearson = corrME_Pearson;
            CG.at(jGene).setAB_R2 = R2;

        if(R2 < tol_R2) {
            if(strMark == "MARK"){ // Mark gene as badly fitting
                CG.at(jGene).Name  = "_" + DS.geneName.at(jGene);
            }
		}

        //get energy statistics
            i = 0;
        for (k = 1; k <= kCodon; k++) {
			if(yModelGene.at(k) > 0 && xExperGene.at(k) > 0) {
					i = i + 1;
				yModelGene.at(i) = std::log(yModelGene.at(k));
				xExperGene.at(i) = std::log(xExperGene.at(k));
				wExperGene.at(i) = wExperGene.at(k);
			}
        }

        Get_X_Y_Correlation(1, i, xExperGene, yModelGene, wExperGene,
                            xExperAver, sigmaXE2, yModelAver, sigmaYM2, corrME);
                R2 = 0;
        if(sigmaYM2 > 0 && sigmaXE2 > 0) {R2 = corrME * corrME / (sigmaYM2 * sigmaXE2);}
                corrME_Pearson = std::sqrt(R2);

    }

    tLhd_Total_UB = tLhd_Total_UB;
    tLhd_Total = tLhd_Total;

	//  Print statistics
	std::ofstream statistic_Output(strOutPutFile);
	 std::stringstream ssLine;
		ssLine.clear(); ssLine.str("");  //Clear Stream Buffer

    ssLine<< "DataSet:  "<< DS.dataSetName << "; #_Genes=" << jTot
		<< "; subSet# =" << jSet << "; subSet Name:" <<std::endl;
    //<< DS.dataSubSetName.at(jSet)
	statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine << strInfo << std::endl;
	statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine <<  "# ;  Gene  ; AAs ; Ui ; Ui/Cod ; R2_Time ; crPrsn ; ExperAver ; sgmExper"
    << " ; ModelAver ; sgmModel ;  ML value ;  Max ML  ;  ML/Max ML ;  Ui/Ti ; Ti/Cod"
	<<std::endl;
	statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

		wTotal = 0.0;
		averPRcorrW = 0.0;
		averPRcorr = 0.0;
	for(j = 1; j <= jTot; j++) {
			jStart = DS.geneStart.at(j) + DS.jGeneStartShift;
			jEnd = DS.geneEnd.at(j);
			nCodon_Elong = jEnd - jStart + 1;
		 wGene = nGeneRPF.at(j) / nGeneCod_Elong.at(j);
		 tModelGenePerCodon = tModelGene_Elong.at(j) / nGeneCod_Elong.at(j);
		 densGeneModel = nGeneRPF.at(j) / tModelGene_Elong.at(j);

		ssLine << j << " ; " << CG.at(j).Name << " ; "
		 << nGeneCod_Elong.at(j) << " ; " << nGeneRPF.at(j) << " ; " << wGene << " ; "
		 << CG.at(j).setAB_R2 << " ; " << CG.at(j).setAB_Corr_Pearson << " ; "
		 << CG.at(j).setA_Aver << " ; " << CG.at(j).setA_Sigma << " ; "
		 << CG.at(j).setB_Aver << " ; " << CG.at(j).setB_Sigma << " ; "
		 << tLhd_Gene.at(j) << " ; " << tLhd_Gene_UB.at(j) << " ; "
		 << (tLhd_Gene.at(j) / tLhd_Gene_UB.at(j)) << " ; "
			<< densGeneModel <<" ; " << tModelGenePerCodon << " ; " << std::endl;
		statistic_Output<<ssLine.str();
		std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");


		averPRcorrW = averPRcorrW + wGene *  CG.at(j).setAB_Corr_Pearson;
			   averPRcorr = averPRcorr + CG.at(j).setAB_Corr_Pearson;
			   wTotal = wTotal + wGene;
	}
	averPRcorrW = averPRcorrW / wTotal;
    averPRcorr = averPRcorr / jTot;

    ssLine <<  "# ;  Gene  ; AAs ; Ui ; Ui/Cod ; R2_Time ; crPrsn ; ExperAver ; sgmExper"
    << " ; ModelAver ; sgmModel ;  ML value ;  Max ML  ;  ML/Max ML ;  Ui/Ti ; Ti/Cod"
	<<std::endl;
	statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine <<  ";   ;  ;  ;  ; AverPRcorr= ;" << averPRcorr << " ; " << std::endl;
    statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine <<  ";   ;  ;  ;  ; AverPRcorrW= ;" << averPRcorrW << " ; " << std::endl;
    statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");


    // get All-Gene statistics
	Get_X_Y_Correlation(1, kRow, xExper, yModel, wExper,
                            xExperAver, sigmaXE2, yModelAver, sigmaYM2, corrME);
        R2 = corrME * corrME / (sigmaYM2 * sigmaXE2);
        corrME_Pearson = std::sqrt(R2);

    // print All-Gene statistics
    ssLine <<  ";   ;  ;  ;  ; R2_Time ; crPrsn ; ExperAver ; sgmExper"
    << " ; ModelAver ; sgmModel ;  ML value ;  Max ML  ;  ML/Max ML "
	<<std::endl;
	statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");

    ssLine << " ;   ;     ;   ;     ; " << R2 << " ; "
    << corrME_Pearson << " ; " << xExperAver << " ; " << std::sqrt(sigmaXE2) << " ; "
    << yModelAver << " ; " << std::sqrt(sigmaYM2) << " ; " << tLhd_Total << " ; "
    << tLhd_Total_UB << " ; " << (tLhd_Total / tLhd_Total_UB) << std::endl;
	statistic_Output<<ssLine.str();
    std::cout<<ssLine.str(); ssLine.clear(); ssLine.str("");
}

// =======================================================================
void Get_X_Y_Correlation(long kStart, long kEnd, vector<double>& X, vector<double>& Y,
	vector<double>& W, double& averX, double& sigmaX2, double& averY, double& sigmaY2,
																		double& corrXY) {
	//
	//  Calculates correlation  between two data sets X and Y
	//
	long i, j, k, mL, mL1;
	double vX, vY, vW, totW, averX2, averY2;

	mL1 = X.size();
	mL=mL1-1;
	// Get Averages and Sigmas
        averX = 0.0;
        averY = 0.0;
        averX2 = 0.0;
        averY2 = 0.0;
        totW = 0.0;
    for(k = kStart; k <= kEnd; k++) {
            vW = W.at(k);
            vX = X.at(k);
            vY = Y.at(k);
        averX = averX + vW * vX;
        averX2 = averX2 + vW * vX * vX;
        averY = averY + vW * vY;
        averY2 = averY2 + vW * vY * vY;
            totW = totW + vW;
    }
        sigmaX2 = 0.0;
        sigmaY2 = 0.0;
        corrXY = 0.0;
    if(totW > 0) {
            averX = averX / totW;
            averX2 = averX2 / totW;
         sigmaX2 = (averX2 - averX * averX);
            averY = averY / totW;
            averY2 = averY2 / totW;
         sigmaY2 = (averY2 - averY * averY);

		// Get correlations
            corrXY = 0;
            totW = 0;
		for(k = kStart; k <= kEnd; k++) {
				vW = W.at(k);
			corrXY = corrXY + vW * (X.at(k) - averX) * (Y.at(k) - averY);
				totW = totW + vW;
        }
            corrXY = corrXY / totW;
    }
}

//=============================================================================
void VB_Gauss_Solve(int iP, const MatrixDouble& mA,
								vector<double>& vX, const vector<double>& vB){

	// Runs Guess elimination with row permutation on extended matrix mA:vB
	// Note that the algorithm with permutations is described in Sprang
	//
	// iP:  input,  iP=0 then skip pivoting; iP>0 then do pivoting
	// mA:  input,  mA is the matrix to be LUP factorized
	//  vX:     output, solution vector
	//  vB:     input, the right side vector
	//
	// Implemented by Michael Pavlov
	//
	long i, j, k, jR, j1, iFlag;
	long mRow, nCol, mRow1, nCol1, nCol2;
	long jP =0;
	double vj, zj;
	double rRow, rRowMax, rRowMaxAbs, rRowAbs, mLij, mUjj, zTol;
	long iRowMax ;
	zTol= 0.0000000000000001;
    mRow1 = mA.size();
	nCol1 = mA[1].size();
        mRow=mRow1-1;
        nCol=nCol1-1;
		nCol2=nCol+2;
     MatrixDouble mU(mRow1,vector<double>(nCol2,0.0));
     vector<long> vP(mRow1,0);
	 vector<double> v(nCol1,0.0);
		vector<double> z(nCol1,0.0);
		vector<double> W(nCol1,0.0);
		vX.clear();
		vX.resize(nCol1,0.0);
	if(!(mRow==nCol)){
        std::cout<< "mA-matrix in not square: mRow="<<mRow
        <<"; nCol=" <<nCol<<std::endl;
	}
	//std::cout<< "mRow=" << mRow<<" ; nCol="<<nCol<<std::endl;

   // prepare mU extended matrix
	for (i = 1; i <= mRow; i++) {
		for (j = 1; j <= nCol; j++) {
			mU[i][j] = mA[i][j];
		}
			vP.at(i) = i;
			mU[i][nCol1] = vB.at(i);
	}

	// main factorization cycle
	for (j = 1; j <= nCol; j++) {
		//at step j find the largest pivot in column j below j
			iRowMax = j;
			rRowMax = mU[j][j];
			rRowMaxAbs = std::abs(rRowMax);
		for (i = j; i <= mRow; i++) {
				rRowAbs = std::abs(mU[i][j]);
			if (rRowMaxAbs < rRowAbs) {// a larger pivot is found
				rRowMaxAbs = rRowAbs;
				iRowMax = i;
			}
		}
		if (iP == 0) { iRowMax = j;} // Run without pivoting
		if (iRowMax > j) { // Swap  rows *iRowMax* and *j* in mL+mR
				// record the row swap in A in the vector vP
				jR = vP.at(j);
				vP.at(j) = vP.at(iRowMax);
				vP.at(iRowMax) = jR; // record the swap
			// actually swap the rows j and iRowMax (both in L and in the remaining of A)
			for (k = 1; k <= nCol1; k++) {
				rRow = mU[j][k];
				mU[j][k] = mU[iRowMax][k];
				mU[iRowMax][k] = rRow;
			}
		}
		// calculate  new column j of a low triangular mL matrix
		// the k element of column j mLkj contains the elimination coefficient
		// for subtracting row j from from row k creating a zero subcolumn in
		// modified A
			mUjj = mU[j][j]; // get the current pivot
			j1 = j + 1;
		if(std::abs(mUjj) >= zTol) { // run the elimination
				for (i = j1; i <= mRow; i++) {
						mLij = mU[i][j] / mUjj; // the new component of the column mL
						mU[i][j] = mLij; // save elimination coefficients in a low part of mU
				  // subtract row j multiplied by Lij from row i of modified A
					for (k = j; k <= nCol1; k++) {
						mU[i][k] = mU[i][k] - mLij * mU[j][k];
					}
				}
		}
		else {//put real zeroes to stress that mUjj=0 and skip the elimination step
					iFlag = 2;  // note the matrix mA singularity
				for (i = j; i<=mRow; i++) {
					mU[i][j] = 0;
				}
			}
		v.at(j) = mU[j][nCol1];
	}

    // Solve mU*z=v by backsubstitutions
    for (j = mRow; j>= 1; j=j-1) {
            mUjj = mU[j][j];
            zj = 8;
        if (std::abs(mUjj) > zTol) {
			zj = v[j] / mUjj;
			}
            z[j] = zj;
            vX[j] = zj;
        for (i = 1; i <= j; i++){
            v[i] = v[i] - zj * mU[i][j];
        }
	}
}

//=============================================================================
void VB_LUPA_Invert(int iP, const MatrixDouble& mA, MatrixDouble& mAM1){

	// Factorizes the matrix as P*A=L*U where L is a low and U an upper triangular
	// The row rearrangements in the A matrix  are kept in permutation vector vP)
	//
	// It then Inverts matrix A that has been LUPA factorized
	// Note that the algorithm with permutations is described in Sprang
	// Implemented by Michael Pavlov
	//
	// iP:  input,  iP=0 then skip pivoting; iP>0 then do pivoting
	// mA:  input,  mA is the matrix to be LUP factorized
	// mAM1:  output, inverted matrix mA
	//
	long i, j, k, jR, j1, iFlag;
	long mRow, nCol, mRow1, nCol1;
	long jP =0;
	double vj, zj;
	double rRow, rRowMax, rRowMaxAbs, rRowAbs, mLij, mUjj, zTol;
	long iRowMax ;
	zTol= 0.0000000000000001;
    mRow1 = mA.size();
	nCol1 = mA[1].size();
        mRow=mRow1-1;
        nCol=nCol1-1;
     MatrixDouble mL(mRow1,vector<double>(nCol1,0.0));
     MatrixDouble mU(mRow1,vector<double>(nCol1,0.0));
     vector<long> vP(mRow1,0);
	if(!(mRow==nCol)){
        std::cout<< "LUPA_Solvwe problem in not square: mRow="<<mRow
        <<"; nCol=" <<nCol<<std::endl;
	}
	//std::cout<< "mRow=" << mRow<<" ; nCol="<<nCol<<std::endl;

   // Fill mU and vP
	for (i = 1; i <= mRow; i++) {
		for (j = 1; j <= nCol; j++) {
			mU[i][j] = mA[i][j];
		}
			vP.at(i) = i;
	}

	// main factorization cycle

	for (j = 1; j <= nCol; j++) {
		//at step j find the largest pivot in column j below j
			iRowMax = j;
			rRowMax = mU[j][j];
			rRowMaxAbs = std::abs(rRowMax);
		for (i = 1; i <= mRow; i++) {
				rRowAbs = std::abs(mU[i][j]);
			if (rRowMaxAbs < rRowAbs) {// a larger pivot is found
				rRowMaxAbs = rRowAbs;
				iRowMax = i;
			}
		}
		if (iP == 0) { iRowMax = j;} // Run without pivoting
		if (iRowMax > j) { // Swap  rows *iRowMax* and *j* in mL+mR
				// record the row swap in A in the vector vP
				jR = vP.at(j);
				vP.at(j) = vP.at(iRowMax);
				vP.at(iRowMax) = jR; // record the swap
			// actually swap the rows j and iRowMax (both in L and in the remaining of A)
			for (k = 1; k <= nCol; k++) {
				rRow = mU[j][k];
				mU[j][k] = mU[iRowMax][k];
				mU[iRowMax][k] = rRow;
			}
		}
		// calculate  new column j of a low triangular mL matrix
		// the k element of column j mLkj contains the elimination coefficient
		// for subtracting row j from from row k creating a zero subcolumn in
		// modified A
			mUjj = mU[j][j]; // get the current pivot
			j1 = j + 1;
		if(std::abs(mUjj) >= zTol) { // run the elimination
				for (i = j1; i <= nCol; i++) {
						mLij = mU[i][j] / mUjj; // the new component of the column mL
						mU[i][j] = mLij; // save elimination coefficients in a low part of mU
				  // subtract row j multiplied by Lij from row i of modified A
					for (k = j1; k <= mRow; k++) {
						mU[i][k] = mU[i][k] - mLij * mU[j][k];
					}
				}
		}
		else {//put real zeroes to stress that mUjj=0 and skip the elimination step
					iFlag = 2;  // note the matrix mA singularity
				for (i = j; i<=nCol; i++) {
					mU[i][j] = 0;
				}
			}
	}
	// separate mU into mL and mR=new mU matrices
	if (nCol <= mRow) {
		for (j = 1; j <= nCol; j++) {
			for (i = j + 1; i<= mRow; i++) {
				mL[i][j] = mU[i][j];
				mU[i][j] = 0;
			}
				mL[j][j] = 1;
		}
	}
	  else {
			mL[1][1] = 0;
		for (i = 2; i <= mRow; i++) {
			for (j = 1; j<= i - 1; j++) {
				mL[i][j] = mU[i][j];
				mU[i][j] = 0;
			}
			for (j = i; j <=  nCol; j++) {
				mL[i][j] = 0;
			}
				mL[i][i] = 1;
		}
	  }

//
	// Inverts matrix A that has been LUPA factorised as P*A=L*U
	// it does it by solving n-equations with ek-right sides
	// i.e it solves A*AM1(k)=e(k) where AM1(k)
	// is a k-column of A inverse

	vector<double> v(nCol1,0.0);
	vector<double> z(nCol1,0.0);
	vector<double> W(nCol1,0.0);
	mAM1.clear();
	mAM1.resize(mRow1, vector<double>(nCol1,0.0));
	long kP=0;
	//prepare intermediate vectors
	//
	for(k = 1; k <= nCol; k++) {
		for(j = 1; j <= nCol; j++){
			v.at(j) = 0;
			z.at(j) = 0;
			W.at(j) = 0;
			jP = vP.at(j); //permutate e(k) components
		  if (jP == k) {kP = k;}
		  mAM1[j][k] = 0;
		}
			W.at(kP) = 1.0;

		// Solve L*v=W for v by forward substitutions
		for(j = kP; j <= nCol; j++){
				vj = 0;
			if (std::abs(W[j]) > 0) {
				vj = W[j] / mL[j][j];
			}
				v[j] = vj;
			if (std::abs(vj) > 0) {
				for (i = j; i <= nCol; i++) {
					W[i] = W[i] - vj * mL[i][j];
				}
			}
		}

		//then solve mU*z=v by back-substitutions
		for (j = nCol; j >= 1 ;j=j-1) {
				mUjj = mU[j][j];
				zj = 8;
			if (std::abs(mUjj) > zTol) {
				zj = v[j] / mUjj;
			}
			else{iFlag=2;}
				z[j] = zj;
				mAM1[j][k] = zj;
			for(i = 1; i<= j; i++){
				v[i] = v[i] - zj * mU[i][j];
			}
		}
	}
}


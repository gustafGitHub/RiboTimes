# R code

The code is written in R and C++ (run in R through Rcpp).

### Required software

This code requires R (for instance RStudio) with the Rcpp package installed. The Rcpp package is installed by typing in R
```install.packages("Rcpp")```

If there is a problem installing Rcpp, the original code ```main.cpp``` witten in C++ is available in the folder ./C.

### Preparation of data

The code uses input data in the form of text file with a table of counts per codon. In order to get this type of data this code by Alexander Bartholomaeus may be useful.

https://github.com/AlexanderBartholomaeus/MiMB_ribosome_profiling

### Clone git repository

Go to the directory where you want to install RiboTimes.

Type

git clone https://github.com/gustafGitHub/RiboTimes.git

### Run code

1. Load the R source code and compile the C++ code by typing

```
source("riboTimesMle.R")
```

2. Set the input and output files

```
InputFile = "/Path to RiboTimes folder/RiboTimes/data/Demo_input_file.txt"  
PathToOutput = "/Path to RiboTimes folder/RiboTimes/data/"
```
Those paths have default values as
```
InputFile = "./data/Demo_input_file.txt"
PathToOutput = "./output/"
```
Assuming that the home folder in R is set as the home folder for the repository
```
setwd("/Path to RiboTimes/RiboTimes")
```


### Input and output files

These files are described in the document User Guide.docx found in the ./documentation folder in 
the repository.


3. Run the code

Most of the code is written in C++, so the code written in R is a wrapper to run this code.
The code is documented in the document User Guide.docx in the /documentation folder included in 
the repository. To run the code in R, type
```
outputList <- MleAlgorithm(InputFile, PathToOutput)
```

### Output variables

The output from the code, outputList, is of the type "List" and contains the following fields

| Field Name | Description |
| ---------- | ----------- |
| dataSetName | Name of the dataset |
| dataSubSet_RPF_Total | Total number of RPFs in dataset |
| doublingTime | Doubling time of Cell culture in h |
| pAsite | Position of the A-site in the local context |
| pNumber | The length of the local context in codons |
| cNumber | Number of codons = 64 |
| geneCodeTable | Standard Genetic Code Table |
| geneCodeTableOrdered | Re-Ordered Standard Genetic Code Table |
| rCodonSeq | Codon sequence of Data Set |
| geneName | Gene names |
| geneStart | Gene Start in the DataSet in codons |
| geneEnd | Gene End in the DataSet in codons |
| jTot | Number of Genes in the Data set |
| geneCodons | Total Number of codons in a gene |
| geneRPFtotal | Total Number of RPFs in a gene |
| geneRPFdensity | RPF density for ORF |
| gene_Elong_AA | Number of Codons in the inner Gene part |
| gene_Ci_Exper | Number of RPFs in the inner Gene part |
| gene_di_Exper | RPFs/codon=Ui/nCodons in the inner Gene part |
| gene_fi_Model | gene_fi_Model=Ci/Ti = corrected relative gene expression level |
| gene_Gi_Model | sum of gij_Model for a gene |
| gene_Gi_Model_Time | sum of gij_Model_Time for a gene= corrected gene elongation time |
| gene_Ti_Model_Time_Abs | estimate of absolute time of gene elongation |
| global_Time_Factor | global time factor used to get absolute elong cycle time in ms |
| rpfOmegaExper | Experimental RPF fingerprint for this DataSubSet |
| iORFcodons | Indexed (from 1 to 64) codons in the dataset |
| indCodonOrder | Codon indexing (from 1 to 64) in re-ordered Genetic Table |
| nRPF | RPFs ascribed to the codon in the dataset (A-site) |
| sij_Exper | experimental RPF scores in the inner gene regions |
| sij_Model | model RPF scores in the inner gene regions |
| sij_Model_Sigma | sigmas of model RPF scores |
| sij_Model_Time | model pausing scores in the inner gene regions |
| sij_Model_Time_Sigma | sigmas of model pausing scores |
| pC_First | first position of the Narrow context |
| pC_Last | last position of the Narrow context |
| gij_Model | model gij-values in the inner gene regions |
| gij_Model_Sigma | sigmas of gij-values |
| gij_Model_Time | corrected model gij-values in the inner gene regions |
| gij_Model_Time_Sigma | sigmas of corrected model gij-values |
| tij_Model_Time_Abs | estimates of absolute cycle times at codon j in gene i |
| tij_Model_Time_Abs_Sigma | sigmas of estimates of absolute cycle times |
| zFP | z(p,c) factors estimated by the maximum likelihood algorithm |
| zFP_Sigma | Sigmas of zFP |

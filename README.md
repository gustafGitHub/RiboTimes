# R code

The code is written in R and C++ (run in R through Rcpp).

### Required software

This code requires R (for instance RStudio) and a C++ compiler.

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
InputFile = "/Path to RiboTimes folder/RiboTimes/data/ZI_30000_FA0.txt"  
PathToOutput = "/Path to RiboTimes folder/RiboTimes/data/"
```

3. Run the code

```
outputList <- MleAlgorithm(InputFile, PathToOutput)
```

### Output variables

The output from the code, outputList, is of the type "List" and contains the following fields

| Field Name | Description |
| ---------- | ----------- |
| dataSetName | Name of the dataset |
| dataSubSet_RPF_Total | Total number of RPFs in dataset |


# Code in Visual Basic in Excel

The original code in Visual Basic in an Excel macro is found in RiboTimes/excel

Instructions how to run the code is found in

Instructions_for_VBA_Activation_On_Your_Window_Computer.pdf

and

Instructions_to_Run_RPF_to_zPC_Transformer.pdf

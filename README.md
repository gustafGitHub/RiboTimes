# Required software

This code requires R and a C++ compiler.

# Preparation of data

The code uses input data in the form of text file with a table of counts per codon. In order to get this type of data this code by Alexander Bartholomaeus may be useful.

https://github.com/AlexanderBartholomaeus/MiMB_ribosome_profiling

# Clone git repository

Go to the directory where you want to install RiboTimes.

Type

git clone https://github.com/gustafGitHub/RiboTimes.git

# Run code

In an R command window, set your selected path containing the RiboTimes repository as working directory

setwd("~/RiboTimes")

To compile Rcpp code and initialise the R functions, type

source("riboTimesMle.R")

To load the data, type

C <- CountsObj()

To initiate the Mle algorithm and calculate the initial guess for the Mle algorithm, type

C <- MleInitiate(C)

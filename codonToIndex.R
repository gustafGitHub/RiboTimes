# codonToIndex
codonToIndex <- function(gene){

codonIndex = integer(0)   
codons = c ("AAA",
            "AAC",
            "AAG",
            "AAT",
            "ACA",
            "ACC",
            "ACG",
            "ACT",
            "AGA",
            "AGC",
            "AGG",
            "AGT",
            "ATA",
            "ATC",
            "ATG",
            "ATT",
            "CAA",
            "CAC",
            "CAG",
            "CAT",
            "CCA",
            "CCC",
            "CCG",
            "CCT",
            "CGA",
            "CGC",
            "CGG",
            "CGT",
            "CTA",
            "CTC",
            "CTG",
            "CTT",
            "GAA",
            "GAC",
            "GAG",
            "GAT",
            "GCA",
            "GCC",
            "GCG",
            "GCT",
            "GGA",
            "GGC",
            "GGG",
            "GGT",
            "GTA",
            "GTC",
            "GTG",
            "GTT",
            "TAA",
            "TAC",
            "TAG",
            "TAT",
            "TCA",
            "TCC",
            "TCG",
            "TCT",
            "TGA",
            "TGC",
            "TGG",
            "TGT",
            "TTA",
            "TTC",
            "TTG",
            "TTT")
  lengthGene = length(gene)
  for (i in 1:lengthGene){
    codonIndex[i]=as.integer(charmatch(gene[i],codons,0))
  }
  return(codonIndex)
}

codonToIndex2 <- function(gene){
  
  codonIndex = integer(0)   
  codons = c ("AAA",
              "AAC",
              "AAG",
              "AAU",
              "ACA",
              "ACC",
              "ACG",
              "ACU",
              "AGA",
              "AGC",
              "AGG",
              "AGU",
              "AUA",
              "AUC",
              "AUG",
              "AUU",
              "CAA",
              "CAC",
              "CAG",
              "CAU",
              "CCA",
              "CCC",
              "CCG",
              "CCU",
              "CGA",
              "CGC",
              "CGG",
              "CGU",
              "CUA",
              "CUC",
              "CUG",
              "CUU",
              "GAA",
              "GAC",
              "GAG",
              "GAU",
              "GCA",
              "GCC",
              "GCG",
              "GCU",
              "GGA",
              "GGC",
              "GGG",
              "GGU",
              "GUA",
              "GUC",
              "GUG",
              "GUU",
              "UAA",
              "UAC",
              "UAG",
              "UAU",
              "UCA",
              "UCC",
              "UCG",
              "UCU",
              "UGA",
              "UGC",
              "UGG",
              "UGU",
              "UUA",
              "UUC",
              "UUG",
              "UUU")
  lengthGene = length(gene)
  for (i in 1:lengthGene){
    codonIndex[i]=as.integer(charmatch(gene[i],codons,0))
  }
  return(codonIndex)
}

codonToIndexVar <- function(gene, codons){
  
  codonIndex = integer(0)
  
  lengthGene = length(gene)
  for (i in 1:lengthGene){
    codonIndex[i]=as.integer(charmatch(gene[i],codons,0))
  }
  return(codonIndex)
  
}

codonsIndexCount <- function(codonIndex){
  codonCount = integer(0)
  nCodons = 64 
  for (i in 1:nCodons) {
    codonCount[i] = sum(codonIndex == i)
  }
  return(codonCount)
}

nucleotideToIndex <- function(gene){
  
  nucleotideIndex = integer(0)   
  nucleotides = c ("A",
              "C",
              "G",
              "T")
  lengthGene = length(gene)
  for (i in 1:lengthGene){
    nucleotideIndex[i]=as.integer(charmatch(gene[i],nucleotides,0))
  }
  return(nucleotideIndex)
}

getNucleotides <- function(){
  nucleotides = c ("A",
              "C",
              "G",
              "T")
  
  return(nucleotides)
  
}

getCodons <- function(){
  codons = c ("AAA",
              "AAC",
              "AAG",
              "AAT",
              "ACA",
              "ACC",
              "ACG",
              "ACT",
              "AGA",
              "AGC",
              "AGG",
              "AGT",
              "ATA",
              "ATC",
              "ATG",
              "ATT",
              "CAA",
              "CAC",
              "CAG",
              "CAT",
              "CCA",
              "CCC",
              "CCG",
              "CCT",
              "CGA",
              "CGC",
              "CGG",
              "CGT",
              "CTA",
              "CTC",
              "CTG",
              "CTT",
              "GAA",
              "GAC",
              "GAG",
              "GAT",
              "GCA",
              "GCC",
              "GCG",
              "GCT",
              "GGA",
              "GGC",
              "GGG",
              "GGT",
              "GTA",
              "GTC",
              "GTG",
              "GTT",
              "TAA",
              "TAC",
              "TAG",
              "TAT",
              "TCA",
              "TCC",
              "TCG",
              "TCT",
              "TGA",
              "TGC",
              "TGG",
              "TGT",
              "TTA",
              "TTC",
              "TTG",
              "TTT")
  
  return(codons)
  
}

getCodonsRNA <- function(){
  codons = c ("AAA",
              "AAC",
              "AAG",
              "AAU",
              "ACA",
              "ACC",
              "ACG",
              "ACU",
              "AGA",
              "AGC",
              "AGG",
              "AGU",
              "AUA",
              "AUC",
              "AUG",
              "AUU",
              "CAA",
              "CAC",
              "CAG",
              "CAU",
              "CCA",
              "CCC",
              "CCG",
              "CCU",
              "CGA",
              "CGC",
              "CGG",
              "CGU",
              "CUA",
              "CUC",
              "CUG",
              "CUU",
              "GAA",
              "GAC",
              "GAG",
              "GAU",
              "GCA",
              "GCC",
              "GCG",
              "GCU",
              "GGA",
              "GGC",
              "GGG",
              "GGU",
              "GUA",
              "GUC",
              "GUG",
              "GUU",
              "UAA",
              "UAC",
              "UAG",
              "UAU",
              "UCA",
              "UCC",
              "UCG",
              "UCU",
              "UGA",
              "UGC",
              "UGG",
              "UGU",
              "UUA",
              "UUC",
              "UUG",
              "UUU")
  
  return(codons)
  
}

tRNA_abundance_old <- function(){
  
  codons = c (#"AAA" 
    0.0297,
              #"AAC" 
    0.0185,
              #"AAG" 
    0.0297,
              #"AAT" 
    0.0185,
              #"ACA" 
    0.0142,
              #"ACC" 
    0.0016,
              #"ACG" 
    0.0084,
              #"ACT" 
    0.0016,
              #"AGA" 
    0.0134,
              #"AGC" 
    0.0218,
              #"AGG" 
    0.0065,
              #"AGT" 
    0.0218,
              #"ATA" 
    0,
              #"ATC" 
    0.0539,
              #"ATG" 
    0.0188 + 0.0111 + 0.0109,
              #"ATT" 
    0.0539,
              #"CAA" 
    0.0118,
              #"CAC" 
    0.0099,
              #"CAG" 
    0.0136,
              #"CAT" 
    0.0099,
              #"CCA" 
    0.0090,
              #"CCC" 
    0.0111,
              #"CCG" 
    0.0138,
              #"CCT" 
    0.0111,
              #"CGA" 
    0.0737,
              #"CGC" 
    0.0737,
              #"CGG" 
    0.0099,
              #"CGT" 
    0.0737,
              #"CTA" 
    0.0103,
              #"CTC" 
    0.0146,
              #"CTG" 
    0.0694,
              #"CTT" 
    0.0146,
              #"GAA" 
    0.0732,
              #"GAC" 
    0.0372,
              #"GAG" 
    0.0732,
              #"GAT" 
    0.0372,
              #"GCA" 
    0.0504,
              #"GCC" 
    0.0095,
              #"GCG" 
    0.0504,
              #"GCT" 
    0.0504,
              #"GGA" 
    0.0331,
              #"GGC" 
    0.0676,
              #"GGG" 
    0.0331,
              #"GGT" 
    0.0676,
              #"GTA" 
    0.0596,
              #"GTC" 
    0.0097,
              #"GTG" 
    0.0596,
              #"GTT" 
    0.0596,
              #"TAA" 
    0,
              #"TAC" 
    0.0119,
              #"TAG" 
    0,
              #"TAT" 
    0.0119,
              #"TCA" 
    0.0201,
              #"TCC" 
    0.0218,
              #"TCG" 
    0.0201,
              #"TCT" 
    0.0201,
              #"TGA" 
    0.0034,
              #"TGC" 
    0.0246,
              #"TGG" 
    0.0146,
              #"TGT" 
    0.0246,
              #"TTA" 
    0.0160,
              #"TTC" 
    0.0160,
              #"TTG" 
    0.0297,
              #"TTT" 
    0.0160
    )
  
}

tRNA_abundance2 <- function(){
  
  codons = c (
    #"AAA" 
    2.26,
    #"AAC" 
    1.01,
    #"AAG" 
    0.71,
    #"AAT" 
    0.84,
    #"ACA" 
    0.33,
    #"ACC" 
    1.34,
    #"ACG" 
    1.51,
    #"ACT" 
    0.94,
    #"AGA" 
    1.34,
    #"AGC" 
    1.41,
    #"AGG" 
    0.65,
    #"AGT" 
    0.77,
    #"ATA" 
    0.49,
    #"ATC" 
    1.54,
    #"ATG" 
    4.08,
    #"ATT" 
    3.34,
    #"CAA" 
    1.18,
    #"CAC" 
    0.6,
    #"CAG" 
    1.36,
    #"CAT" 
    0.39,
    #"CCA" 
    0.2,
    #"CCC" 
    0.53,
    #"CCG" 
    1.78,
    #"CCT" 
    0.99,
    #"CGA" 
    0.61,
    #"CGC" 
    2.73,
    #"CGG" 
    0.99,
    #"CGT" 
    4.03,
    #"CTA" 
    0.05,
    #"CTC" 
    0.83,
    #"CTG" 
    7.96,
    #"CTT" 
    0.63,
    #"GAA" 
    6.13,
    #"GAC" 
    1.39,
    #"GAG" 
    1.19,
    #"GAT" 
    2.33,
    #"GCA" 
    1.48,
    #"GCC" 
    0.95,
    #"GCG" 
    2.45,
    #"GCT" 
    1.12,
    #"GGA" 
    1.39,
    #"GGC" 
    3.68,
    #"GGG" 
    1.92,
    #"GGT" 
    3.08,
    #"GTA" 
    1.18,
    #"GTC" 
    0.89,
    #"GTG" 
    2.81,
    #"GTT" 
    1.97,
    #"TAA" 
    0,
    #"TAC" 
    1.88,
    #"TAG" 
    0,
    #"TAT" 
    1.26,
    #"TCA" 
    0.45,
    #"TCC" 
    0.67,
    #"TCG" 
    1.24,
    #"TCT" 
    1.36,
    #"TGA" 
    0.0034,
    #"TGC" 
    1.69,
    #"TGG" 
    1.46,
    #"TGT" 
    0.77,
    #"TTA" 
    0.54,
    #"TTC" 
    0.79,
    #"TTG" 
    4.03,
    #"TTT" 
    0.81
  )
  
  sumCodons = sum(codons)
  
  codons_n = codons/sumCodons
  
  return(codons)
  
}


tRNA_abundance3 <- function(){
  
  codons = c (
    #"AAA" 
    10.43,
    #"AAC" 
    7.29,
    #"AAG" 
    10.43,
    #"AAT" 
    7.29,
    #"ACA" 
    6.89,
    #"ACC" 
    0.67 + 5.54,
    #"ACG" 
    6.89 + 3.12,
    #"ACT" 
    6.89 + 0.67 + 5.54,
    #"AGA" 
    1.34,
    #"AGC" 
    1.41,
    #"AGG" 
    0.65,
    #"AGT" 
    0.77,
    #"ATA" 
    24.74,
    #"ATC" 
    24.74,
    #"ATG" 
    4.43,
    #"ATT" 
    24.74,
    #"CAA" 
    4.38,
    #"CAC" 
    4.38,
    #"CAG" 
    4.38 + 6.27,
    #"CAT" 
    4.38,
    #"CCA" 
    2.56,
    #"CCC" 
    2.56 + 3.75,
    #"CCG" 
    2.56 + 2.67,
    #"CCT" 
    2.56 + 3.75,
    #"CGA" 
    25.57,
    #"CGC" 
    25.57,
    #"CGG" 
    2.30,
    #"CGT" 
    25.57,
    #"CTA" 
    3.17,
    #"CTC" 
    5.93,
    #"CTG" 
    3.17 + 22.20,
    #"CTT" 
    3.17 + 5.93,
    #"GAA" 
    29.35,
    #"GAC" 
    15.46,
    #"GAG" 
    29.35,
    #"GAT" 
    15.46,
    #"GCA" 
    20.97,
    #"GCC" 
    20.97 + 3.57,
    #"GCG" 
    20.97,
    #"GCT" 
    20.97 + 3.57,
    #"GGA" 
    11.08,
    #"GGC" 
    24.96,
    #"GGG" 
    11.08,
    #"GGT" 
    24.96,
    #"GTA" 
    20.39,
    #"GTC" 
    2.79 + 4.42,
    #"GTG" 
    20.39,
    #"GTT" 
    20.39+2.79+4.42,
    #"TAA" 
    0,
    #"TAC" 
    4.19 + 5.04,
    #"TAG" 
    0,
    #"TAT" 
    4.19 + 5.04,
    #"TCA" 
    7.36,
    #"TCC" 
    4.03,
    #"TCG" 
    7.36 + 1.45,
    #"TCT" 
    7.36 + 4.03,
    #"TGA" 
    0.0034,
    #"TGC" 
    7.07,
    #"TGG" 
    5.02,
    #"TGT" 
    7.07,
    #"TTA" 
    3.78,
    #"TTC" 
    5.11,
    #"TTG" 
    3.78 + 9.30,
    #"TTT" 
    5.11
  )
  
  sumCodons = sum(codons)
  
  codons_n = codons/sumCodons
  
  return(codons_n)
  
}

tRNAconcentrationPerAA <- function(){
  #Table 3. The molar ratio of tRNA/ribosome at different
  #growth rates (doublings/hour)
  #tRNA 0.4 0.7 1.07 1.6 2.5
  Ala1B = c(0.65, 0.54, 0.46, 0.44, 0.40)
  M = matrix(Ala1B, ncol = 5, nrow = 1)
  rownames(M)<-"Ala1B"
  Ala2 = c(0.12, 0.10, 0.08, 0.08, 0.07)
  M<-rbind(M, Ala2)
  Arg2 = c(0.95, 0.67, 0.50, 0.60, 0.49)
  M<-rbind(M, Arg2)
  Arg3 = c(0.13, 0.12, 0.05, 0.06, 0.04)
  M<-rbind(M, Arg3)
  Arg4 = c(0.17, 0.11, 0.09, 0.08, 0.07)
  M<-rbind(M, Arg4)
  Arg5 = c(0.08, 0.07, 0.05, 0.06, 0.04)
  M<-rbind(M, Arg5)
  Asn = c(0.24, 0.18, 0.14, 0.15, 0.14)
  M<-rbind(M, Asn)
  Asp1 = c(0.48, 0.37, 0.27, 0.30, 0.29)
  M<-rbind(M, Asp1)
  Cys = c(0.32, 0.22, 0.17, 0.18, 0.14)
  M<-rbind(M, Cys)
  Gln1 = c(0.15, 0.12, 0.12, 0.08, 0.08)
  M<-rbind(M, Gln1)
  Gln2 = c(0.18, 0.14, 0.11, 0.13, 0.12)
  M<-rbind(M, Gln2)
  Glu2 = c(0.94, 0.71, 0.54, 0.61, 0.56)
  M<-rbind(M, Glu2)
  Gly1 = c(0.43, 0.33, 0.25, 0.28, 0.21)
  M<-rbind(M, Gly1)
  Gly2 = c(0.43, 0.33, 0.25, 0.28, 0.21)
  M<-rbind(M, Gly2)
  Gly3 = c(0.87, 0.70, 0.55, 0.50, 0.48)
  M<-rbind(M, Gly3)
  His = c(0.13, 0.10, 0.09, 0.08, 0.08)
  M<-rbind(M, His)
  Ile1 = c(0.69, 0.54, 0.43, 0.47, 0.47)
  M<-rbind(M, Ile1)
  Ile2 = c(0.69, 0.54, 0.43, 0.47, 0.47)
  M<-rbind(M, Ile2)
  Leu1 = c(0.89, 0.68, 0.55, 0.54, 0.42)
  M<-rbind(M, Leu1)
  Leu2 = c(0.19, 0.16, 0.13, 0.12, 0.11)
  M<-rbind(M, Leu2)
  Leu3 = c(0.13, 0.11, 0.09, 0.08, 0.06)
  M<-rbind(M, Leu3)
  Leu4 = c(0.38, 0.29, 0.23, 0.24, 0.18)
  M<-rbind(M, Leu4)
  Leu5 = c(0.23, 0.16, 0.13, 0.09, 0.07)
  M<-rbind(M, Leu5)
  Lys = c(0.38, 0.31, 0.24, 0.22, 0.20)
  M<-rbind(M, Lys)
  Met_f1 = c(0.24, 0.22, 0.19, 0.16, 0.19)
  M<-rbind(M, Met_f1)
  Met_f2 = c(0.14, 0.10, 0.08, 0.09, 0.07)
  M<-rbind(M, Met_f2)
  Met_m = c(0.14, 0.12, 0.09, 0.10, 0.09)
  M<-rbind(M, Met_m)
  Phe = c(0.21, 0.17, 0.14, 0.12, 0.10)
  M<-rbind(M, Phe)
  Pro1 = c(0.18, 0.11, 0.11, 0.07, 0.05)
  M<-rbind(M, Pro1)
  Pro2 = c(0.14, 0.12, 0.07, 0.10, 0.07)
  M<-rbind(M, Pro2)
  Pro3 = c(0.12, 0.09, 0.07, 0.06, 0.05)
  M<-rbind(M, Pro3)
  SelCys = c(0.04, 0.04, 0.03, 0.03, 0.02)
  M<-rbind(M, SelCys)
  Ser1 = c(0.26, 0.25, 0.18, 0.17, 0.14)
  M<-rbind(M, Ser1)
  Ser2 = c(0.07, 0.05, 0.04, 0.03, 0.03)
  M<-rbind(M, Ser2)
  Ser3 = c(0.28, 0.20, 0.15, 0.14, 0.11)
  M<-rbind(M, Ser3)
  Ser5 = c(0.15, 0.12, 0.09, 0.09, 0.08)
  M<-rbind(M, Ser5)
  Thr1 = c(0.02, 0.02, 0.02, 0.01, 0.01)
  M<-rbind(M, Thr1)
  Thr2 = c(0.11, 0.09, 0.07, 0.07, 0.06)
  M<-rbind(M, Thr2)
  Thr3 = c(0.22, 0.17, 0.13, 0.12, 0.11)
  M<-rbind(M, Thr3)
  Thr4 = c(0.18, 0.14, 0.11, 0.12, 0.13)
  M<-rbind(M, Thr4)
  Trp = c(0.19, 0.13, 0.11, 0.10, 0.10)
  M<-rbind(M, Trp)
  Tyr1 = c(0.15, 0.11, 0.09, 0.12, 0.08)
  M<-rbind(M, Tyr1)
  Tyr2 = c(0.25, 0.18, 0.12, 0.13, 0.10)
  M<-rbind(M, Tyr2)
  Val1 = c(0.77, 0.55, 0.36, 0.48, 0.39)
  M<-rbind(M, Val1)
  Val2A = c(0.13, 0.09, 0.08, 0.07, 0.05)
  M<-rbind(M, Val2A)
  Val2B = c(0.13, 0.11, 0.09, 0.09, 0.08)
  M<-rbind(M, Val2B)
  #4.5 S RNA 0.08 0.07 0.06 0.05 0.05
  return(M)
  
}

getAminoAcids<-function(){
  
  amino<-c(
    
    "Ala1B",
    "Ala2",
    "Arg2",
    "Arg3",
    "Arg4",
    "Arg5",
    "Asn",
    "Asp1",
    "Cys",
    "Gln1",
    "Gln2",
    "Glu2",
    "Gly1",
    "Gly2",
    "Gly3",
    "His",
    "Ile1",
    "Ile2",
    "Leu1",
    "Leu2",
    "Leu3",
    "Leu4",
    "Leu5",
    "Lys",
    "Met_f1",
    "Met_f2",
    "Met_m",
    "Phe",
    "Pro1",
    "Pro2",
    "Pro3",
    "SelCys",
    "Ser1",
    "Ser2",
    "Ser3",
    "Ser5",
    "Thr1",
    "Thr2",
    "Thr3",
    "Thr4",
    "Trp",
    "Tyr1",
    "Tyr2",
    "Val1",
    "Val2A",
    "Val2B")
    
    return(amino)
    
  
}


getAminoAcids2<-function(){
  
  amino<-c(
    
    "Phe",
    "Leu",
    "Ile",
    "Met",
    "Val",
    "Ser",
    "Pro",
    "Thr",
    "Ala",
    "Tyr",
    "His",
    "Gln",
    "Asn",
    "Lys",
    "Asp",
    "Glu",
    "Cys",
    "Trp",
    "Arg",
    "Gly"
    
  )
 
  return(amino) 
}

AminoAcid2codon2<-function(){
  
  codonList = list()
  codons <- getCodonsRNA()
  
  codonList[[1]] = charmatch(c("UUU", "UUC"), codons) # Phe
  codonList[[2]] = charmatch(c("UUA", "UUG", "CUU", "CUC", "CUA", "CUG"), codons) # Leu
  codonList[[3]] = charmatch(c("AUU", "AUC", "AUA"), codons) # Ile
  codonList[[4]] = charmatch("AUG", codons) # Met
  codonList[[5]] = charmatch(c("GUU", "GUC", "GUA", "GUG"), codons) # Val
  codonList[[6]] = charmatch(c("UCU", "UCC", "UCA", "UCG", "AGU", "AGC"), codons) # Ser
  codonList[[7]] = charmatch(c("CCU", "CCC", "CCA", "CCG"), codons) # Pro
  codonList[[8]] = charmatch(c("ACU", "ACC", "ACA", "ACG"), codons) # Thr
  codonList[[9]] = charmatch(c("GCU", "GCC", "GCA", "GCG"), codons) # Pro
  codonList[[10]] = charmatch(c("UAU", "UAC"), codons) # Tyr
  codonList[[11]] = charmatch(c("CAU", "CAC"), codons) # His
  codonList[[12]] = charmatch(c("CAA", "CAG"), codons) # Gln
  codonList[[13]] = charmatch(c("AAU", "AAC"), codons) # Asn
  codonList[[14]] = charmatch(c("AAA", "AAG"), codons) # Lys
  codonList[[15]] = charmatch(c("GAU", "GAC"), codons) # Asp
  codonList[[16]] = charmatch(c("GAA", "GAG"), codons) # Glu
  codonList[[17]] = charmatch(c("UGU", "UGC"), codons) # Cys
  codonList[[18]] = charmatch("UGG", codons) # Trp
  codonList[[19]] = charmatch(c("CGU", "CGC", "CGA", "CGG", "AGA", "AGG"), codons) # Arg
  codonList[[20]] = charmatch(c("GGU", "GGC", "GGA", "GGG"), codons) # Gly
  
  return(codonList)

}

aminoAcid2Index<-function(aminoAcid){
  
  aminoAcidAll<-getAminoAcids()
  aaIndex=as.integer(charmatch(aminoAcid,aminoAcidAll,0))
  return(aaIndex)
  
}

AminoAcid2codon<-function(){
  
  codonList = list()
  codons <- getCodonsRNA()
  
  codonList[[1]] = charmatch(c("GCU", "GCA", "GCG"), codons) #"Ala1B"
  codonList[[2]] = charmatch("GCC", codons) #"Ala2"
  codonList[[3]] = charmatch(c("CGU", "CGC", "CGA"), codons) #"Arg2"
  codonList[[4]] = charmatch("CGG", codons) #"Arg3"
  codonList[[5]] = charmatch("AGA", codons) #"Arg4"
  codonList[[6]] = charmatch("AGG", codons) #"Arg5"
  codonList[[7]] = charmatch(c("AAC", "AAU"), codons) #"Asn"
  codonList[[8]] = charmatch(c("GAC", "GAU"), codons) #"Asp1"
  codonList[[9]] = charmatch(c("UGC", "UGU"), codons) #"Cys"
  codonList[[10]] = charmatch("CAA", codons) #"Gln1"
  codonList[[11]] = charmatch("CAG", codons) #"Gln2"
  codonList[[12]] = charmatch(c("GAA", "GAG"), codons) #"Glu2"
  codonList[[13]] = charmatch("GGG", codons) #"Gly1"
  codonList[[14]] = charmatch(c("GGA", "GGG"), codons) #"Gly2"
  codonList[[15]] = charmatch(c("GGC", "GGU"), codons) #"Gly3"
  codonList[[16]] = charmatch(c("CAC", "CAU"), codons) #"His"
  codonList[[17]] = charmatch(c("AUC", "AUU"), codons) #"Ile1"
  codonList[[18]] = charmatch("AUA", codons) #"Ile2"
  codonList[[19]] = charmatch("CUG", codons) #"Leu1"
  codonList[[20]] = charmatch(c("CUC", "CUU"), codons) #"Leu2"
  codonList[[21]] = charmatch(c("CUA", "CUG"), codons) #"Leu3"
  codonList[[22]] = charmatch("UUG", codons) #"Leu4"
  codonList[[23]] = charmatch(c("UUA", "UUG"), codons) #"Leu5"
  codonList[[24]] = charmatch(c("AAA", "AAG"), codons) #"Lys"
  codonList[[25]] = charmatch("AUG", codons) #"Met_f1"
  codonList[[26]] = charmatch("AUG", codons) #"Met_f2"
  codonList[[27]] = charmatch("AUG", codons) #"Met_m"
  codonList[[28]] = charmatch(c("UUC", "UUU"), codons) #"Phe"
  codonList[[29]] = charmatch("CCG", codons) #"Pro1"
  codonList[[30]] = charmatch(c("CCC", "CCU"), codons) #"Pro2"
  codonList[[31]] = charmatch(c("CCA", "CCU", "CCG"), codons) #"Pro3"
  codonList[[32]] = charmatch("UGA", codons) #"SelCys"
  codonList[[33]] = charmatch(c("UCA", "UCU", "UCG"), codons) #"Ser1"
  codonList[[34]] = charmatch("UCG", codons) #"Ser2"
  codonList[[35]] = charmatch(c("AGC", "AGU"), codons) #"Ser3"
  codonList[[36]] = charmatch(c("UCC", "UCU"), codons) #"Ser5"
  codonList[[37]] = charmatch(c("ACC", "ACU"), codons) #"Thr1"
  codonList[[38]] = charmatch("ACG", codons) #"Thr2"
  codonList[[39]] = charmatch(c("ACC", "ACU"), codons) #"Thr3"
  codonList[[40]] = charmatch(c("ACA", "ACU", "ACG"), codons) #"Thr4"
  codonList[[41]] = charmatch("UGG", codons) #"Trp"
  codonList[[42]] = charmatch(c("UAC", "UAU"), codons) #"Tyr1"
  codonList[[43]] = charmatch(c("UAC", "UAU"), codons) #"Tyr2"
  codonList[[44]] = charmatch(c("GUA", "GUG", "GUU"), codons) #"Val1"
  codonList[[45]] = charmatch(c("GUC", "GUU"), codons) #"Val2A"
  codonList[[46]] = charmatch(c("GUC", "GUU"), codons) #"Val2B"
  
  return(codonList)
  
}


codon2AminoAcid<-function(codonIndex){
  
  aa<-getAminoAcids()
  cc<-AminoAcid2codon()
  ccu<-unlist(lapply(cc, function(x) x == codonIndex))
  aai <- AminoAcid2codonIndex()
  
  aaName = aa[aai[ccu]]
  
  return(aaName)
  
}

codon2AminoAcidNumber<-function(codonIndex){
  
  aa<-getAminoAcids()
  cc<-AminoAcid2codon()
  ccu<-unlist(lapply(cc, function(x) x == codonIndex))
  aai <- AminoAcid2codonIndex()
  
  aaNumber = aai[ccu]
  
  return(aaNumber)
  
}

AminoAcid2codonIndex<-function(){
  
  cc<-AminoAcid2codon()
  ucc = unlist(lapply(cc, length))
  aaIndex = NULL
  for(i in 1:length(ucc)){
    aaIndex = c(aaIndex, rep(i, ucc[i]))
    
  }
  
  return(aaIndex)
  
}

AminoAcid2codonIndex2<-function(){
  
  cc<-AminoAcid2codon2()
  ucc = unlist(lapply(cc, length))
  aaIndex = NULL
  for(i in 1:length(ucc)){
    aaIndex = c(aaIndex, rep(i, ucc[i]))
    
  }
  
  return(aaIndex)
  
}

getCharge<-function(){
  
  charge<-c(
    
    0, #"Ala1B",
    0, #"Ala2",
    1, #"Arg2",
    1, #"Arg3",
    1, #"Arg4",
    1, #"Arg5",
    0, #"Asn",
    -1, #"Asp1",
    0, #"Cys",
    0, #"Gln1",
    0, #"Gln2",
    -1, #"Glu2",
    0, #"Gly1",
    0, #"Gly2",
    0, #"Gly3",
    0.1, #"His",
    0, #"Ile1",
    0, #"Ile2",
    0, #"Leu1",
    0, #"Leu2",
    0, #"Leu3",
    0, #"Leu4",
    0, #"Leu5",
    1, #"Lys",
    0, #"Met_f1",
    0, #"Met_f2",
    0, #"Met_m",
    0, #"Phe",
    0, #"Pro1",
    0, #"Pro2",
    0, #"Pro3",
    0, #"SelCys",
    0, #"Ser1",
    0, #"Ser2",
    0, #"Ser3",
    0, #"Ser5",
    0, #"Thr1",
    0, #"Thr2",
    0, #"Thr3",
    0, #"Thr4",
    0, #"Trp",
    0, #"Tyr1",
    0, #"Tyr2",
    0, #"Val1",
    0, #"Val2A",
    0 #"Val2B"
    
  )
  
  return(charge)
  
}

getHydropathy<-function(){
  
  hydro<-c(
    
    rep(1.8, 4), #"Ala"
    rep(-4.5, 6), #"Arg"
    rep(-3.5, 2), #"Asn"
    rep(-3.5, 2), #"Asp"
    rep(2.5, 2), #"Cys"
    rep(-3.5, 2), #"Gln"
    rep(-3.5, 2), #"Glu"
    rep(-0.4, 5), #"Gly"
    rep(-3.2, 2), #"His"
    rep(4.5, 3), #"Ile"
    rep(3.8, 8), #"Leu"
    rep(-3.9, 2), #"Lys"
    rep(1.9, 3), #"Met"
    rep(2.8, 2), #"Phe"
    rep(-1.6, 6), #"Pro"
    0, #"SelCys"
    rep(-0.8, 8), #"Ser"
    rep(-0.7, 8), #"Thr"
    -0.9, #"Trp",
    rep(-1.3, 4), #"Tyr",
    rep(4.2, 7) #"Val"
    
  )
  
  return(hydro)
  
}

codon2Charge <- function(codonIndex){
  
  aa<-getAminoAcids()
  cc<-AminoAcid2codon()
  ccu<-unlist(lapply(cc, function(x) x == codonIndex))
  aai <- AminoAcid2codonIndex()
  chi = getCharge()
  
  charge = chi[aai[ccu]]
  
  return(charge)
  
}

codon2Hydro <- function(codonIndex){
  
  cc<-AminoAcid2codon()
  ccu<-unlist(lapply(cc, function(x) x == codonIndex))
  #aai <- AminoAcid2codonIndex()
  h <- getHydropathy()
  
  charge = h[ccu]
  
  return(charge)
  
}


tRNA_abundance4 <- function(ConcentrationIndex){
  
  tRNAconc <- tRNAconcentrationPerAA()
  tRNAa = numeric()
  for(i in 1:64){
    
    tRNAa[i] = sum(tRNAconc[codon2AminoAcidNumber(i),ConcentrationIndex])
      
  }
  
  return(tRNAa)
  
  
}
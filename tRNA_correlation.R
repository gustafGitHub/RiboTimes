# tRNA correlation
source("codonToIndex.R")

CodonIndexes <- read.table("./data/CodonIndexes.txt")
codonNameOrder = substr(CodonIndexes[[2]], 1,3)
codons <- getCodonsRNA()
codonIndexOrder <- codonToIndexVar(codons, codonNameOrder)
load("savedMle")

z_a_site = Mle@z_initial[8,]

tRNA_anbundance = tRNA_abundance4(3)

z_a_site_sorted = Mle@z_initial[7,codonIndexOrder]

is_stop = get_is_stop()

corTest = cor.test(tRNA_anbundance[!is_stop], z_a_site_sorted[!is_stop])
corTest2 = corr_function(tRNA_anbundance[!is_stop], z_a_site_sorted[!is_stop])

pdf(file = "plots/correlation.pdf", width = 6, height = 5)

plot(tRNA_anbundance[!is_stop], z_a_site_sorted[!is_stop], xlab = "tRNA abundance", ylab = "z factor p site", main = "Correlation between tRNA abundance and z score")

dev.off()

corr_function <- function(X, Y){
  
  r_out = list()
  
  Xm = mean(X)
  Ym = mean(Y)
  
  n = length(X)
  n2 = length(Y)
  
  if(n != n2){stop('X and Y are not the same size')}
  
  s1 = 0
  s2 = 0
  s3 = 0
  
  for(i in 1:n){
    
    s1 = s1 + (X[i] - Xm)*(Y[i] - Ym)
    s2 = s2 + (X[i] - Xm)^2
    s3 = s3 + (Y[i] - Ym)^2
    
  }
  
  r = s1/(sqrt(s2)*sqrt(s3))
  
  Fr = atanh(r)
  z = Fr*sqrt(n - 3)
  p = 2*pnorm(-z)
  
  r_out$r = r
  r_out$p = p
  
  return(r_out)
  
}
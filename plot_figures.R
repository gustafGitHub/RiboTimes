# Figure 2

rpsQ_index = which(outputList$geneName == "rpsQ")
atpE_index = which(outputList$geneName == "atpE")

geneStart_rpsQ = outputList$geneStart[rpsQ_index]
geneEnd_rpsQ = outputList$geneEnd[rpsQ_index]
geneStart_atpE = outputList$geneStart[atpE_index]
geneEnd_atpE = outputList$geneEnd[atpE_index]

counts_rpsQ = outputList$nRPF[geneStart_rpsQ: geneEnd_rpsQ]
counts_atpE= outputList$nRPF[geneStart_atpE: geneEnd_atpE]

pdf(file = "./plots/Figure2_RpsQ_Counts.pdf", width = 30, height = 11.5)

plot(1:length(counts_rpsQ), counts_rpsQ, type = 'h')

dev.off()

s_exp_rpsQ = outputList$sij_Exper[geneStart_rpsQ: geneEnd_rpsQ]
s_mod_rpsQ = outputList$sij_Model[geneStart_rpsQ: geneEnd_rpsQ]

pdf(file = "./plots/Figure3A_RpsQ_S_exp_vs_mod.pdf", width = 30, height = 11.5)

plot(8:(length(counts_rpsQ) - 8), s_exp_rpsQ[8:(length(s_mod_rpsQ) - 8)], type = 'b', col = 'red', ylim = c(0, 12))
points(8:(length(counts_rpsQ) - 8), s_mod_rpsQ[8:(length(s_mod_rpsQ) - 8)], type = 'b', col = 'blue')

dev.off()

s_exp_atpE = outputList$sij_Exper[geneStart_atpE: geneEnd_atpE]
s_mod_atpE = outputList$sij_Model[geneStart_atpE: geneEnd_atpE]

pdf(file = "./plots/Figure3B_atpE_S_exp_vs_mod.pdf", width = 30, height = 11.5)

plot(8:(length(counts_atpE) - 8), s_exp_atpE[8:(length(s_mod_atpE) - 8)], type = 'b', col = 'red', ylim = c(0, 14))
points(8:(length(counts_atpE) - 8), s_mod_atpE[8:(length(s_mod_atpE) - 8)], type = 'b', col = 'blue')

dev.off()
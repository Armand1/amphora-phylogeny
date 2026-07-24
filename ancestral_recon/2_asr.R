suppressMessages({library(treeio); library(ape)})
setwd("/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree")
phy <- as.phylo(read.beast("runs/srvf38pc/output/simple_run_srvf38PC_mcc.tree"))
S <- read.csv("prep/srvf38_tip_scores.csv", row.names=1)
S <- as.matrix(S[phy$tip.label,])                       # order to tips
anc <- sapply(1:ncol(S), function(j) ace(S[,j], phy, type="continuous", method="pic")$ace)
colnames(anc) <- colnames(S)
rownames(anc) <- (Ntip(phy)+1):(Ntip(phy)+phy$Nnode)    # internal node numbers
write.csv(anc, "prep/srvf38_node_scores.csv")
cat("ASR done:", nrow(anc), "internal nodes x", ncol(anc), "PCs\n")

library(phylosim)
set.seed(833)

setwd("AU/meta_strains_project/evo_sim")

tree<-read.tree("tree_B_root.nwk")
plot(tree)

library("Biostrings")

fastaFile <- readDNAStringSet("my.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)


#root.seq <- NucleotideSequence(length=5000000)
root.seq <- NucleotideSequence(length=100000)
p <- JC69()
attachProcess(root.seq, p)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
sim$rootSeq <- root.seq
Simulate(sim)

args=(commandArgs(TRUE))
start = as.real(args[1])
step = as.real(args[2])
iter = as.integer(args[3])
length =as.integer( args[4])
library(phylosim)
file=file(paste('./data/', args[5], sep=""), "w")
for(i in 0:(iter-1)){
	branchLength = (start + i * step)/2
	cat(paste("(t1:",branchLength,",t2:",branchLength,");", sep=""), file = "./data/2taxa.nwk")
	sim = Simulate(PhyloSim(root.seq = sampleStates(NucleotideSequence(len = length,proc = list(list(JC69())))), phylo = read.tree("3taxa.nwk")))
	writeLines(paste(">t1:", branchLength, "|", length, sep=""), file)
	writeLines(paste(sim$alignment["t1",],collapse=''), file)
	writeLines(paste(">t2:", branchLength, "|", length, sep=""), file)
	writeLines(paste(sim$alignment["t2",],collapse=''), file)
	
}


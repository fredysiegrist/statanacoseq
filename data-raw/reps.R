reps <- rep(NA, times=65*(length(GtRNAdb2species)+1))
dim(reps) <- c(65,length(GtRNAdb2species)+1)
reps[,1] <- readtRNAout()$reps
colnames(reps) <- rep("", times=length(GtRNAdb2species)+1)
colnames(reps)[1] <- "Etef"
rownames(reps) <- names(aa_ac)


for (species in 1:length(GtRNAdb2species)) {
  reps[,species+1] <- readfasta(species)$reps
  colnames(reps)[species] <- GtRNAdb2species[species]
}
heatmap(reps, scale='column', col=topo.colors(100))

for (species in 1:length(GtRNAdb2species)) {
  if (!all(rownames(readfasta(species))==names(aa_ac))) print(species)
}

library(statanacoseq)
context("All species have all codons defined")

for (species in 1:length(GtRNAdb2species)) {
  test_that(paste(species, "has all codons present"), {
    expect_equal(rownames(readfasta(species)), names(aa_ac))
  })
  gc()
}

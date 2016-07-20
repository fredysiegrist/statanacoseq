library(statanacoseq)
context("Frequency of optimal codons are calculated (with added opal stop codon) as initially set")

x <-  c(0.4238506, 0.3888255, 0.3604938, 0.4260870, 0.3328059, 0.4127660, 0.3576799, 0.4117647, 0.3787276, 0.3541430, 0.3304094, 0.4055944, 0.3810445, 0.3392857, 0.3554422, 0.3421053, 0.3851852, 0.3325876, 0.4074074, 0.3966942, 0.3662420, 0.3571429, 0.3950512, 0.3770492, 0.3447368, 0.3473684, 0.3447368, 0.3447368, 0.3900227, 0.3825243, 0.4026316, 0.3888889, 0.3407407, 0.4523810, 0.3696172, 0.3280255, 0.3787234)
for (n in 1:length(mylist(whatout=1))) {
  test_that("Sample CDS has correct Fop calculated also with opal stop codon",
            expect_equal(ComputeFop(toupper(paste(c2s(mylist(whatout=1)[[n]]), "tga", sep="")),
                ref='tRNAgene', codonusagec=Tef), x[n],  tolerance = .0000001)
                         )}

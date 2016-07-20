library(statanacoseq)
context("Check for the correct calculation of effective numbers of codons")

  test_that("Check if Nc calculation form Fuglsang 2004 equals the one from Wright 1990",
            expect_equal(vhica::CUB(sequence=mylist(whatout=1)),  unlist(lapply(mylist(whatout=1), function(x) {ComputeNEC(c2s(x))})))
  )




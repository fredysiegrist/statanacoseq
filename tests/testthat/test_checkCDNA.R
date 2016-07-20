library(statanacoseq)
context("Sample CDS throws warning and missed stop codon can be rescued")

for (testcdna in mylist(whatout=1)) {
  test_that(paste(attr(testcdna, "name"), "throws missing stop codon warning"), {
    expect_warning(checkCDS(c2s(testcdna)), "CDS does not end with stop codon")
  })
  test_that(paste(attr(testcdna, "name"), "don't throw warning with stop codon"), {
    expect_true(checkCDS(paste(c2s(testcdna), "tag", sep="")))
  })
  gc()
}

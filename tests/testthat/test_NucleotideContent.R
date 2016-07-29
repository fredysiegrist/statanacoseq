library(statanacoseq)
context("Nucleotide composition of gene, compare nucleotide content with probabilities")
library(seqinr)

letters <- c("A", "C", "G", "T")
probabilities <- c(0.25348595, 0.08518904, 0.22941381, 0.43191120)
# Make a random sequence of length n letters, using the multinomial model with probabilities "probabilities"
seq <- sample(letters, 1000, rep=TRUE, prob=probabilities) # Sample with replacement
expec <- as.numeric(table(seq)/rep(length(seq), times=length(table(seq))))
names(expec) <- names(table(seq))
  test_that("Nucleotide content is calculating correctly if all positions given",
            expect_equal(NucleotideContent(c2s(seq)), expec,
                         tolerance = .00001)
                         )

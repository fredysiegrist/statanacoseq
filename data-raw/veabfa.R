setwd("E:/R/viruses")
v <- dir(pattern=".stats.fa$")

setwd("E:/R/eukariots")
e <- dir(pattern=".stats.fa$")

setwd("E:/R/archaeas")
a <- dir(pattern=".stats.fa$")

setwd("E:/R/bacterias")
b <- dir(pattern=".stats.fa$")

veabfa <- c(v, e, a, b)
setwd("E:/R")
save(veabfa, file="veabfa.RData")
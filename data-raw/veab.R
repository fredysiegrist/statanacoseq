setwd("E:/R/viruses")
v <- dir(pattern=".stats.html$")

setwd("E:/R/eukariots")
e <- dir(pattern=".stats.html$")

setwd("E:/R/archaeas")
a <- dir(pattern=".stats.html$")

setwd("E:/R/bacterias")
b <- dir(pattern=".stats.html$")

veab <- c(v, e, a, b)
setwd("E:/R")
save(veab, file="veab.RData")
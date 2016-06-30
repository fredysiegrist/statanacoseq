setwd("E:/R/gtrnadb.ucsc.edu/GtRNAdb2/genomes")
directories <- dir(getwd())        # c("archaea" ,  "bacteria" , "eukaryota" ,"viruses" ) )

eukariots <- dir(paste(getwd(), directories[3], sep="/"))
eukariots <-  gsub("#", "%23", eukariots)
archeas <- dir(paste(getwd(), directories[1], sep="/"))
virus <- dir(paste(getwd(), directories[4], sep="/"))
bacis <- dir(paste(getwd(), directories[2], sep="/"))
bacis <-  gsub("#", "%23", bacis)
bacis <-  gsub("Pseu_poae_RE_1_1_14", "Pseu_poae_RE*1_1_14", bacis)



GtRNAdb2species <- c(virus, eukariots, archeas, bacis)
setwd("E:/R")
save(GtRNAdb2species, file="GtRNAdb2species.RData")
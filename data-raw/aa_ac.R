aa_ac <- vir
for (x in names(aa_ac)) aa_ac[x] <- 0
names(aa_ac) <- paste(substr(names(aa_ac), 1, 3), substr(names(aa_ac), nchar(names(aa_ac))-2, nchar(names(aa_ac))), sep="_")
names(aa_ac)[64] <- "SeC_TCA"
aa_ac <- c(aa_ac, "Und_???"=0)
save(aa_ac, file="aa_ac.RData")
a <- rep(0, times=37)
b <- a
for (i in 1:37) {
a[i] <- mean(ComputeGC12syn(toupper(c2s(mylist(whatout=1)[[i]]))))
b[i] <- ComputeGC3syn(toupper(c2s(mylist(whatout=1)[[i]])))
}
plot(b, a, xlab="GC content at 3rd pos", ylab="GC content at 1st and 2nd pos", pch=16)
abline(coefficients(lm(a~b)), col="blue4", lwd=2)

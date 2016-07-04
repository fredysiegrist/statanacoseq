for (n in 1:length(mylist(whatout=1))) {print(paste(attr(mylist(whatout=1)[[n]], 'name'), round(ComputeFop(toupper(c2s(mylist(whatout=1)[[n]])), ref='tRNAgene', codonusagec=Tef),3)), sep=': ')}

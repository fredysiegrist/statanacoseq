for (n in 1:length(mylist(whatout=1))) {print(paste(attr(mylist(whatout=1)[[n]], 'name'), round(ComputeNEC(toupper(c2s(mylist(whatout=1)[[n]]))),3)), sep=': ')}

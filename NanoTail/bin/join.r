g<-read.table("test2ko.assigned", sep="\t", header=FALSE)
a<-read.table("test2ko.findr.len", sep="\t", header=FALSE)
b<-read.table("test2ko.nanopol.len", sep="\t", header=FALSE)

c<-merge(a,b, all.x=T, by.x="V1", by.y="V1", sort=FALSE)
d<-merge(c,g, all.x=T, by.x="V1", by.y="V1", sort=FALSE)i

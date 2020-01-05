findr	<- "test2ko.findr.len"
nanopol	<- "test2ko.nanopol.len"
assigned <- "test2ko.assigned"
prefix = "test2ko"

# R --slave --args findr nanopol assigned prefix < join.r

args<-commandArgs(TRUE)
findr	<- args[1]
nanopol	<- args[2]
assigned <- args[3]
prefix <- args[4]

a<-read.table(findr, sep="\t", header=FALSE)
b<-read.table(nanopol, sep="\t", header=FALSE)
g<-read.table(assigned, sep="\t", header=FALSE)

c<-merge(a,b, all.x=T, by.x="V1", by.y="V1", sort=FALSE)
d<-merge(c,g, all.x=T, by.x="V1", by.y="V1", sort=FALSE)

colnames(d)<-c("Read name", "Tailfindr", "Nanopolish", "Gene Name")

e<-na.omit(d)
e[,2]<-log(e[,2]+1,base=2)
e[,3]<-log(e[,3]+1,base=2)

write.table(d, file = paste0(prefix, "_joined.txt"), sep="\t", row.names = FALSE)


x<-e[,2]
y<-e[,3]

data.lm<-lm(y ~ x)
rsquared<-round(summary(data.lm)$r.squared, 2)

subtitle=bquote(R^2 == .(rsquared))
pdf(paste0(prefix, "_corr.pdf"))
plot(x, y, xlab="TailfindR", ylab="Nanopolish", ylim=c(0,10), xlim=c(0,10), pch=19, frame=FALSE, main=prefix, sub=subtitle)
abline(data.lm, col = "blue")
dev.off()

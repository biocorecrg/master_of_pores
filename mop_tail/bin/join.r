# R --slave --args findr nanopol assigned prefix < join.r

args<-commandArgs(TRUE)
findr	<- args[1]
nanopol	<- args[2]
assigned <- args[3]
prefix <- args[4]

#Arguments
zinput1<-gzfile(findr, "rt")
zinput2<-gzfile(nanopol, "rt")
zinput3<-gzfile(assigned, "rt")

a<-read.delim(zinput1, sep="\t", header=FALSE)
b<-read.delim(zinput2, sep="\t", header=FALSE)
g<-read.delim(zinput3, sep="\t", header=FALSE)

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

rsquare=bquote(R^2 == .(rsquared))

subtitle=expression('PolyA tail length estimation log'[2]*'(n+1)')

png(paste0(prefix, "_corr.png"), width = 1600, height = 1200)
op <- par(mar=c(8, 8, 8, 8) + 0.1)
plot(x, y, xlab="TailfindR", ylab="Nanopolish", ylim=c(0,10), xlim=c(0,10), cex=.2, cex.lab=2, cex.sub=2,, cex.main=3, cex.axis=2, pch=19, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5), frame=FALSE, main=prefix)
mtext(side=3, line=0, at=6.5, adj=1, cex=2, subtitle)
mtext(side=3, line=-65, at=10, adj=1, cex=2, rsquare)

par(op)
abline(data.lm, col = "blue")
dev.off()

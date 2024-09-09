# R --slave --args epi tombo prefix < join.r

args<-commandArgs(TRUE)
epi	<- args[1]
tombo	<- args[2]
prefix <- args[3]

a<-read.table(epi, sep="\t", header=FALSE)
b<-read.table(tombo, sep="\t", header=FALSE)

a$V2<-"epinano"
b$V2<-"tombo"

c<-merge(a,b, all=T, by.x="V1", by.y="V1", sort=FALSE)
colnames(c)<-c("positions", "epinano", "tombo")

c[c=="epinano"]<-1
c[c=="tombo"]<-1
c[is.na(c)]<-0

write.table(c, file = "RNA_modifications.txt", sep="\t", row.names = FALSE)

library(VennDiagram)
library(RColorBrewer)
myCol <- c("#B3E2CD","#FDCDAC")

venn.diagram(
    x = list(a$V1, b$V1),
    category.names = c("Epinano", "Tombo"),
    filename = 'venn_diagram.png',
    imagetype = "png",
    output=TRUE,
    height = 1024, 
    width = 1024 , 
    resolution = 300,
    compression = "lzw",
    main.pos = c(0.5,0.7),
    
	# Circles
	lwd = 2,
	lty = 'blank',
	fill = myCol,
	margin = 0.6,
    main = prefix,
    main.fontfamily = "sans",
     
	# Numbers
	cex = .6,
	fontface = "bold",
	fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-135, 135),
    cat.fontfamily = "sans"
)



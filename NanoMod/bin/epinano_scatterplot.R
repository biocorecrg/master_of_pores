#Scatter plots 
#Rscript epinano_scatterplot.R input1 label1 input2 label2 feature
#Libraries needed
library(plyr)
library(ggplot2)
library(ggrepel) 
library(MASS)
library(reshape2)

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
zinput1<-gzfile(args[1], "rt")
zinput2<-gzfile(args[3], "rt")

input1 <- read.delim(zinput1,sep=",")  #1st variable
label1 <- as.character(args[2])  #1st label
input2 <-read.delim(zinput2,sep=",") #2nd variable
label2 <- as.character(args[4]) #2nd label
feature<- as.character(args[5]) #Feature



#Cleanup 
cleanup <- function(input, label) {
	#Filter low coverage reads
	input <- subset(input, cov>30)
	#Filter read starts
	input <- subset(input, pos>20)
	#Add summed errors column
	input$sum <- input$mis + input$del + input$ins 
	#Add a column with position 
	input$position<- paste(input$X.Ref,input$pos)
	#Change column names 
	input <- input[, c("X.Ref","pos","position", "base", feature)]
	colnames(input)<- c("Chr","Position","chr_pos","base",feature )
	data_melted<- melt(data = input, id.vars = c("Chr", "Position", "chr_pos", "base"))
    colnames(data_melted)[which(names(data_melted) == "value")] <- paste(label, "value", sep="_")
	return(data_melted)
}


#Cleanup and process the data
data1 <- cleanup(input1, label1)
data2 <- cleanup(input2, label2)

merged <- join(data1,data2, by="chr_pos")
merged$Chr <- NULL
merged$Position <- NULL
merged$base <- NULL
merged$variable <- NULL
			

plot<- function(data)
for (chr in unique(data$Chr)) {
	subs <- subset(data,  Chr==chr)
        if(nrow(subs)>0){
		res<- rlm(subs[,c(paste(label1, "value", sep="_"))] ~ subs[,c(paste(label2, "value", sep="_"))]) #linear model  
		res_vec <- res$residuals#this contains residuals 
		threshold <-  5 * sd(res_vec) #The threshold
		subs$score<- abs(subs[,c(paste(label1, "value", sep="_"))] - subs[,c(paste(label2, "value", sep="_"))])
		pdf(file=paste(chr,feature, label1, label2, "scatter.pdf", sep="_"),height=5,width=5,onefile=FALSE)
			print(ggplot(subs, aes_string(x=paste(label1, "value", sep="_"), y=paste(label2, "value", sep="_"))) +
				geom_point(size=2, color="grey")+
				geom_abline(slope=1, intercept=0,linetype="dashed")+
     				geom_point(data=subset(subs, score>threshold), size=2, color="red")+
     				geom_text_repel(data=subset(subs, score>threshold), aes(label=chr_pos), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
     				ggtitle(feature)+
				xlab(label1)+
				ylab(label2) +
				theme_bw()+
				theme(axis.text.x = element_text(face="bold", color="black",size=11),
					 axis.text.y = element_text(face="bold", color="black", size=11),
				plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
				axis.title.x = element_text(color="black", size=15, face="bold"),
				axis.title.y = element_text(color="black", size=15, face="bold"),
				panel.background = element_blank(),
				axis.line = element_line(colour = "black", size=0.5),
				legend.title = element_text(color = "black", size = 15,face="bold"),
				legend.text = element_text(color = "black", size=15),
				panel.grid.major = element_blank(), panel.grid.minor = element_blank())
				)
			dev.off()
		}
	}

plot(merged)

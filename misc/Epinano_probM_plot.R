#!/usr/bin/env Rscript

suppressMessages (library(reshape2))
suppressMessages (library(ggplot2))
suppressMessages (library(ggrepel))
suppressMessages (library(tidyverse))

library(stringr)



args <- commandArgs (trailingOnly=TRUE)

df <- read.csv (args[1], header=T)[,c("Window","Ref","Strand","ProbM","prediction")]

split <- data.frame (str_split_fixed (df$Window, '-', 2))

split$middle <- as.integer (as.character(split$X1)) +  2

df$Position <- split$middle

df<- df[,c("Ref", "Position","Strand", "ProbM","prediction")]

barplot <- function (df, out_pdf) {
   mod = df[df$prediction=="mod",]
   p <- ggplot(df, aes_string(x="Position", y="ProbM")) +
	     	ggtitle ("Probability of modification") +
   		geom_bar(stat = "identity", width=0.1, fill="#2a7886") +
	        xlab("Positions")+
        	ylab("ProbM") +
        	theme_bw()+
	        theme(axis.text.x = element_text(face="bold", color="black",size=11),
        	axis.text.y = element_text(face="bold", color="black", size=11),
	        plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
        	axis.title.x = element_text(color="black", size=15, face="bold"),
	        axis.title.y = element_text(color="black", size=15, face="bold"),
	        panel.background = element_blank(),
        	legend.position = "none",
	        axis.line = element_line(colour = "black", size=0.5)
	)
    p <- p + geom_text_repel(data=mod, aes(x=Position, y=ProbM, label=Position), color="red", box.padding=3, point.padding=3, segment.size  = 0.1, segment.color = "green")+
    pdf(file=out_pdf,height=5,width=20,onefile=FALSE)
    print (p)
    dev.off()
}

for (chr in unique(df$Ref)) {
	sub <- df [df$Ref == chr,]
	pdf_out = paste (chr, ".ProbM.pdf" ,sep="")
	barplot (sub, pdf_out)
}

#!/usr/bin/env Rscripts 

args = commandArgs (trailingOnly=TRUE)

suppressMessages (library(plyr))
suppressMessages (library(dplyr))
suppressMessages (library(tidyr))
suppressMessages (library(ggplot2)) 
suppressMessages (library(ggrepel))
suppressMessages (library(MASS))
suppressMessages (library(reshape2))
suppressMessages (library (tidyverse))

library (ggpubr)

input <- read.csv(args[1], header=T) #, comment.char="#")

input$IDX<- paste(input$Ref, input$Window,  input$Strand,input$X.Kmer, sep=";")

#head (input)

input <- input %>% separate (Window, c("pos1","pos2")) 

input$Pos <- as.numeric (input$pos1) + 2

#input <- input[, c("IDX", "Ref", "Pos", "sample", "ProbM")]
input <- input[, c("IDX", "Ref", "Pos", "prediction", "ProbM")]

#data_melted<- melt(data = input, id.vars = c("IDX", "sample","Ref", "Pos"))
data_melted<- melt(data = input, id.vars = c("IDX", "prediction","Ref", "Pos"))

for (chr in unique(data_melted$Ref)) {
  #pdf(file=paste(chr,  "probability.dot-density.pdf", sep="_"), height=5, width=20, onefile=FALSE)
  subs = subset(data_melted, Ref=chr)
  subs$ProbM = subs$value
  sp <- ggscatter(subs, x = "Pos", y = "ProbM", color = "prediction")
  border ()
  xplot <- ggdensity(subs, "value", fill="prediction", palette="jc")
  yplot <- ggdensity(subs, "value", fill="prediction", palette="jc") + rotate()
  yplot <- yplot + clean_theme()
  xplot <- xplot + clean_theme()
  ## Arranging the plot
  ggarrange (xplot, NULL, sp, yplot, ncol=2, nrow=2, align='hv', widths=c(2,1), heights=c(1,2), common.legend = T)
  ggsave(file=paste(chr,"probability.dot-density.pdf",sep="."))
  dev.off()
}

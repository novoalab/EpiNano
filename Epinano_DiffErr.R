#!/usr/bin/env Rscript

args = commandArgs (trailingOnly=TRUE)
	if (length(args) == 0 ) {
		stop ("\n\ntry --help|-h to display help msg!!\n\n")
}

suppressMessages (library(optparse))
option_list <- list (
	make_option (c("-c","--coverage"), type="integer", default=30, 
		help="minimum coverage/depth; default: 30"),
	make_option (c("-t", "--threshold"), type="double", default=3, 
		help="minimum z-score (i.e., number of standard deviation from mean) to determine modified sites; default: 3"),
	make_option (c("-d", "--deviance"), type="double", default=0.1, 
		help="minimum deviance of selected feature between two samples; default: 0.1"),
	make_option (c("-f","--feature"), type="character", 
		help="the feature (column name(s) in from input file) to use to predict modifications"),
	make_option (c("-k", "--ko_sample"), type="character", 
		help="knockout/unmodified sample"),
	make_option (c("-w","--wt_sample"), type="character", 
 		help="wildtype/modified sample"), 
	make_option (c("-o","--out_prefix"), type="character", 
		help="output prefix"),
	make_option (c("-p", '--plot'), type="logical", default=0, action="store_true",
		help = "whether or not generate plots;  default: no plots will be generated because Epinano_Plot.R can do the job")
)

parser <- parse_args (OptionParser (option_list=option_list, usage="
	DiffErr.R v0.1 compares a given feature between two samples. It predict potential modified sites mainly through two methods:
		1. compute deviance of selected featuers between samples and then calculate z-scores. Outliers or potential modified sites will then 
		be determined based on user-defined threshold. Note that this is not suit for our published curlcakes construct data becos they are full of modifications. 
		2. fit a linear regression model between two samples.
			1) detect residuals outliers of the given linear model. 
			2) compute z-scores of residuals for each observation in turn and determine outliers using user-defined z-score threshold.
		Examples:	
			1 compare sum_err between two samples
			Rscript Epinano_DiffErr.R -k ko.csv -w wt.csv -t 3 -o Test -c 30 -f sum_err  -d 0.1 
			2 same as above, but generate plots, one for each reference. 
			Rscript Epinano_DiffErr.R -k ko.csv -w wt.csv -t 3 -o Test -c 30 -f sum_err  -d 0.1 -p 
	   "))

suppressMessages (library (outliers))
suppressMessages (library(reshape2))
suppressMessages (library(ggplot2))
suppressMessages (library(car))
suppressMessages (library(ggrepel))
suppressMessages (library(tidyverse))

if (!is.na(parser$feature)) { 
	feature <- parser$feature
} else { 
	stop ('please provide the feature you would like to use to detect modification' ) 
}

if (!is.na(parser$out_prefix)) { 
	prefix <- parser$out_prefix
} else { 
	prefix <- 'DirrErrOut'
}

if (!is.na(parser$coverage)) { 
	coverage <- parser$coverage
} else { 
	coverage <- 30
}

if (!is.na(parser$deviance)) { 
	deviance <- parser$deviance
} else { 
	deviance <- 0.1
}

if (!is.na(parser$threshold)) { 
	threshold <- parser$threshold
} else { 
	threshold <- 3
}

if (!is.na(parser$plot)) { 
	plot <- parser$plot
} else { 
	plot <- 0
}

if (!is.na(parser$ko_sample)) { 
	ko <- parser$ko_sample
} else { 
	stop ('please provide the 1st Error/Variants features table you would like to use to detect modification' ) 
}

if (!is.na(parser$wt_sample)) { 
	wt <- parser$wt_sample
} else { 
	stop ('please provide the 2nd Error/Variants features table you would like to use to detect modification' ) 
}

out1 = paste (prefix,".","delta-",feature,".prediction.csv", sep="")
out2 = paste (prefix,"linear-regression","prediction.csv", sep=".")

ko <- read.csv (ko, header = T)
wt <- read.csv (wt, header = T)

Chrs <- unique(unique(ko$X.Ref), unique(wt$X.Ref))

cleanup <- function(input, label, coverage, feature) {
	input <- input[input$cov>coverage,]
	#Filter read starts
	input <- input[input$pos>20, ]
	#Add a column with position 
	input$position <- paste(input$X.Ref,input$pos, input$base, input$strand)
	input$sum_err <- rowSums(input[,c("mis", "ins", "del")])
	#Change column names 
	input <- input[, c("X.Ref","pos","position", "base", "strand", feature)]
	colnames(input) <- c("Chr","Position","chr_pos","base","strand", feature)
	data_melted <- melt(data = input, id.vars = c("Chr", "Position", "chr_pos", "base", "strand"))
	colnames(data_melted)[which(names(data_melted) == "value")] <- paste(label, "value", sep="_")
  	to_drop <- c("Chr","Position", "base", "strand", "variable")
  	data_melted <- data_melted[,!(colnames(data_melted) %in% to_drop)]
	return(data_melted)
}

dat1  <- cleanup(ko, 'ko', coverage, feature)
dat2  <- cleanup(wt, 'wt', coverage, feature)
#write.table (combine,file="ko.csv",sep=",",quote=FALSE, row.names=FALSE)
#write.table (combine,file="wt.csv",sep=",",quote=FALSE, row.names=FALSE)
combine <- merge(dat1, dat2, by="chr_pos")
#write.table (combine,file="combine.csv",sep=",",quote=FALSE, row.names=FALSE)


#primary_filt <- function (combine, feature, feature_deviance) {
	#delta = paste("delta_", feature, sep="")
	#combine$delta <- abs (combine$wt_value - combine$ko_value)
	#names(combine)[ncol(combine)] <- delta
	#print (head(combine	))
	#stop()
	#combine <- combine[combine[, ncol(combine)] > feature_deviance, ]
#}

univariate_outlier <- function (combine, Threshold, deviance, feature) { 
	delta = paste("delta_", feature, sep="")
	combine$delta <- abs (combine$wt_value - combine$ko_value)
	names(combine)[ncol(combine)] <- delta
	combine$z_scores <- scores (combine[,ncol(combine)], type="z")  # aka, analyze delta feature 
	combine$z_score_prediction <- ifelse(combine$z_scores > Threshold & combine[delta] > deviance, "mod", "unm")
	colnames (combine)[which (names(combine) == "ko_value")] <- paste ("ko","feature",sep="_")
	colnames (combine)[which (names(combine) == "wt_value")] <- paste ("wt","feature",sep="_")
	return (combine)
}

multi_variate_outlier <- function (combine, deviance, feature) {
	lmFit <- lm (wt_value ~ ko_value, data= combine)
	test<-outlierTest(lmFit, cutoff=0.05, n.max=ncol(combine))
	outlier_names <- names(test$rstudent)
	combine$lm_Bonferroni_outlier_test <- ifelse (rownames(combine) %in% names(test$rstudent) & combine$wt_value-combine$ko_value > deviance, "mod","unm")
	combine$lm_residuals <- lmFit$residuals
	combine$lm_residuals_z_score <- scores (combine$lm_residuals, type='z')
	combine$lm_residuals_z_scores_prediction <- ifelse (combine$lm_residuals_z_score > threshold & combine$wt_value - combine$ko_value>deviance, "mod","unm")
	colnames (combine)[which (names(combine) == "ko_value")] <- paste ("ko",feature,sep="_")
	colnames (combine)[which (names(combine) == "wt_value")] <- paste ("wt",feature,sep="_")
	return (combine)
}

scatter_plot <- function (df, feature, out_pdf) {
	pdf(file=out_pdf,height=5,width=20,onefile=FALSE)
	ko_feature = paste ("ko_",feature,sep='')
	wt_feature = paste ("wt_",feature,sep='')
	#df$tmpX  <- df[, colnames(df) %in% ko_feature]
	#df$tmpY  <- df[, colnames(df) %in% wt_feature]
	mod = df[df$lm_residuals_z_scores_prediction=="mod",]
	title = paste ('ko_',feature,' ~ ','wt_',feature,sep="")
    print(ggplot(df, aes_string(x=ko_feature, y=wt_feature)) +
                geom_point(size=2, color="grey")+
                geom_abline(slope=1, intercept=0, linetype="dashed")+
                geom_point(data=mod, size=2, color="red")+ 
                geom_text_repel (data=mod, aes(label=Position), color='red',
                	box.padding = 0.35,
                	point.padding = 0.5,
                	segment.color = 'green') + 
                ggtitle(title) +
                xlab(ko_feature) +
                ylab(wt_feature) +
                theme_bw() +
                xlim (0,1) + ylim (0,1) +
                theme(axis.text.x = element_text(face="bold", color="black",size=11),
                                axis.text.y = element_text(face="bold", color="black", size=11),
                                plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
                                axis.title.x = element_text(color="black", size=15, face="bold"),
                                axis.title.y = element_text(color="black", size=15, face="bold"),
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black", size=0.5),
                                legend.title = element_text(color = "black", size = 15,face="bold"),
                                legend.text = element_text(color = "black", size=15),
                        		panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                      ) + coord_fixed() 
         		) 
	dev.off()
    	#ggsave (out_pdf)
}

bar_plot <- function (df, feature, out_pdf) {
	pdf(file=out_pdf,height=5,width=20,onefile=FALSE)
	df$tmp_feature <- df[,4] # in order to pass name to ggplot
	mod = df[df$z_score_prediction=="mod",]
	#write.table(df, file = paste(out_pdf, '.csv',sep=""),sep=",", quote=FALSE, row.names=FALSE )	
	print(ggplot(df, aes_string(x="Position", y="tmp_feature")) +
      geom_bar(stat = "identity", width=0.1, fill="#2a7886") +
      geom_text_repel(data=mod, aes_string("Position", "tmp_feature", label="Position"), size=3, color="red", segment.size  = 1, segment.color = "black")+
      ggtitle(paste(chr, feature, sep="_"))+
      xlab("Positions")+
      ylab(feature) +
      theme_bw()+
      theme(axis.text.x = element_text(face="bold", color="black",size=11),
        axis.text.y = element_text(face="bold", color="black", size=11),
        plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black", size=0.5)))

	dev.off()
#	ggsave (out_pdf)
}

n <- 0 
for (chr in Chrs) {
	sub <- combine[grepl(chr, combine$chr_pos, fixed=TRUE), ]
	#sub <- primary_filt (sub, feature, deviance)

	nrows = nrow (sub)
	if (nrows>0) {  
		#sub_out = paste (chr,"sub.out.csv",sep="")
		#print (chr)
		#write.table(sub, file=sub_out, sep=",", quote=FALSE, row.names=FALSE)
		delta_feature <-univariate_outlier(sub, threshold, deviance, feature)
		lm_feature <- multi_variate_outlier(sub, deviance, feature)

		if (n==0) {
			write.table (delta_feature, file=out1, sep=",", quote=FALSE, row.names=FALSE)
			write.table (lm_feature, file=out2, sep=",", quote=FALSE, row.names=FALSE)
		} else {
			write.table (delta_feature, file=out1, sep=",", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
			write.table (lm_feature, file=out2, sep=",", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)	
		}
		pos<-c()
		for (x in strsplit(delta_feature$chr_pos,' ', 4)) {pos <- c(pos, x[2])}
		delta_feature$Position <- as.numeric (pos)
	#delta_feature <- delta_feature [order(delta_feature$Position),]

		pos<-c()
		for (x in strsplit(lm_feature$chr_pos,' ', 4)) {pos <- c(pos, x[2])}
		lm_feature$Position <- as.numeric (pos)
	#lm_feature <- lm_feature [order(lm_feature$Position),]
	
		n <- n + 1
		if (plot) {
			barplot <- paste (chr,".",prefix,".","delta-",feature,".bar.pdf", sep="")
			bar_plot (delta_feature, paste ("delta_",feature, sep=""), barplot)
			xyplot <- paste (chr,prefix,"linear-regression","scatter.pdf", sep=".")
			scatter_plot (lm_feature, feature, xyplot)
		}
	}
}

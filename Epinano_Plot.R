#!/usr/bin/env Rscript

args <- commandArgs (trailingOnly=TRUE)
if (length(args) == 0 ) {
        stop ("\n\nUSAGE:\n------------------------------------------------------------------------------------------------------
Rscript Epinano_Plot.R <modification prediction results from Epinano_DiffErr.R and Epinano_Predict.py> \n\n")
}

type <<- "pdf"

if (length (args) == 2 & args[2] == "png") {
    type <<- 'png'
}

suppressMessages (library(reshape2))
suppressMessages (library(ggplot2))
suppressMessages (library(ggrepel))
suppressMessages (library(tidyverse))
suppressMessages (library(stringr))

diff_err_scatter_plot <- function (df, feature, out_fig) {
    if (type == "png") { 
	   png (filename = out_fig, height=5,width=20)
    } else {
        pdf(file = out_fig, height=5,width=20,onefile=FALSE)
    }
	ko_feature <- colnames(df)[2]
    wt_feature <- colnames(df)[3]
	mod = df[df$lm_residuals_z_scores_prediction=="mod",]
	title = paste ('ko_',feature,' ~ ','wt_',feature,sep="")
    print(ggplot(df, aes_string(x=ko_feature, y=wt_feature)) +
                geom_point(size=2, color="grey")+
                geom_abline(slope=1, intercept=0, linetype="dashed")+
                geom_point(data=mod, size=2, color="red")+ 
                geom_text_repel (data=mod, aes(label=Position), color='red', 
                	box.padding = 1,
                	point.padding = 1,
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
    	#ggsave (out_fig)
}


diff_err_bar_plot <- function (df, feature, out_fig) {
    if (type == "png") { 
       png (filename = out_fig, height=5,width=20)
    } else {
        pdf(file = out_fig, height=5,width=20,onefile=FALSE)
    }
    df$tmp_feature <- df[,4] # in order to pass name to ggplot
    df$Position <- as.integer (as.character (separate(data=df, col=chr_pos, into =c("chr", "pos"),sep=" ")$pos))
    mod = df[df$z_score_prediction=="mod",]
    #write.table(df, file = paste(out_pdf, '.csv',sep=""),sep=",", quote=FALSE, row.names=FALSE )	
	print(ggplot(df, aes_string(x="Position", y="tmp_feature")) +
      geom_bar(stat = "identity", width=0.1, fill="#2a7886") +
      geom_text_repel(data=mod, aes_string("Position", "tmp_feature", label="Position"), color="red", segment.size  = 1, segment.color = "black")+
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
        axis.line = element_line(colour = "black", size=0.5)) 
      ) 
	dev.off()
	#ggsave (out_fig)
}

svm_probm_plot <- function (inputfile, out_prefix) {
    df <- read.csv (inputfile, header=T)[,c("Window","Ref","Strand","ProbM","prediction")]
    split <- data.frame (str_split_fixed (df$Window, '-', 2))
    split$middle <- as.integer (as.character(split$X1)) +  2
    df$Position <- split$middle
    df<- df[,c("Ref", "Position","Strand", "ProbM","prediction")]
    
    barplot <- function (df, out_fig) {
        if (type %in% "png") { 
            png (file = out_fig, height=5,width=20)
        } else {
            pdf(file = out_fig, height=5,width=20,onefile=FALSE)
        }
        mod = df[df$prediction=="mod",]
        print ( ggplot (df, aes_string(x="Position", y="ProbM")) +
                    ggtitle ("Probability of modification") +
                    geom_bar (stat = "identity", width=0.1, fill="#2a7886") +
                    xlab ("Positions")+
                    ylab ("ProbM") +
                    theme_bw ()+
                    theme (axis.text.x = element_text(face="bold", color="black",size=11),
                           axis.text.y = element_text(face="bold", color="black", size=11),
                            plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
                            axis.title.x = element_text(color="black", size=15, face="bold"),
                            axis.title.y = element_text(color="black", size=15, face="bold"),
                            panel.background = element_blank(),
                            legend.position = "none",
                            axis.line = element_line(colour = "black", size=0.5)) + 
                    geom_text_repel (data=mod, aes(x=Position, y=ProbM, label=Position), color="red", box.padding=3, point.padding=3, segment.size  = 0.6, segment.color = "green") 
                    
            )
        dev.off()
    }

    for (chr in unique(df$Ref)) {
            sub <- df [df$Ref == chr,]
            out_fig = paste (chr, out_prefix, ".ProbM.", type, sep="")
            barplot (sub, out_fig)  
    }
}

## main 
df <- read.csv (args[1], header = TRUE) 

if ("z_score_prediction" %in% colnames (df)) {
    feature <- gsub ('delta_','',colnames (df)[4])
    Chrs <- unique (separate(data=df, col=chr_pos, into =c("chr", "pos"),sep=" ")$chr)
    for (chr in Chrs) {
        sub <- df[grepl(chr, df$chr_pos, fixed=TRUE), ]
        out_fig = paste (chr,".delta_",feature,".",type, sep="")
        diff_err_bar_plot (sub, feature, out_fig)
    }
} else if ("lm_residuals_z_scores_prediction" %in% colnames(df)) {
    feature <- gsub ('ko_','',colnames(df)[2])
    Chrs <- unique (separate(data=df, col=chr_pos, into =c("chr", "pos"),sep=" ")$chr)
    df$Position <- as.integer (as.character (separate(data=df, col=chr_pos, into =c("chr", "pos"),sep=" ")$pos))
    for (chr in Chrs) {
        sub <- df[grepl(chr, df$chr_pos, fixed=TRUE), ]
        out_fig = paste (chr,".",feature,".lm_regression.",type, sep="")
        diff_err_scatter_plot (sub, feature, out_fig)
    }
} else if ('dist' %in% colnames(df) & 'ProbM' %in% colnames (df)) {
    prefix <- str_split (args[1], "MODEL", 2)[[1]][2]
    svm_probm_plot (args[1], prefix)
} 

#add current intensity values analysis 
#("SVM_Predict_delta_features.mis3.del3.q3.MODEL.rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump.csv","MODEL",2)[[1]][2]

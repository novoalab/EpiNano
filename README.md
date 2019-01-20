# EpiNano
Detection of RNA modifications from Oxford Nanopore direct RNA sequencing reads

### About EpiNano
EpiNano is a tool to identify RNA modifications present in direct RNA sequencing reads. The current algorithm has been trained and tested on detecting m6A RNA modifications. 

EpiNano will extract a set of 'features' from direct RNA sequencing reads, which will be in turn used to predict whether the 'error' is caused by the presence of an RNA modification or not. 
Features extracted include: 
- k-mer current intensity
- read quality
- per-base quality
- per-base mismatch frequency
- per-base deletion frequency
- per-base insertion frequency

The software has been trained and tested upon a set of 'unmodified' and 'modified' sequences containing m6A at known sites or A. Its use to detect other RNA modifications has not yet been tested. 

## What's included
- Scripts to extract features from FAST5 files
- Scripts to process mapped BAM files into kmer pileups (similar to samtools mpileup format but for 5mer sequences)
- Support Vector Machine training (SVM) & testing to predict m6A RNA modifications

## Getting Started and pre-requisites
In addition to python2.7 and python3, minimap2, jvarkit-sam2tsv.jar, and samtools1.4+ should be installed.
Python modules including h5py, seaborn, numpy, pandas, and scikitlearn should also be installed.

## Running the software
* To extract features from basecalled FASTQ files: 
```
#1 trim the first and last few bad quality bases from raw fastq with NanoFilt

cat MOD.fastq|./NanoFilt -q 0 --headcrop 5 --tailcrop 3 --readtype 1D --logfile mod.nanofilt.log > mod.h5t3.fastq
cat UNM.fastq|./NanoFilt -q 0 --headcrop 5 --tailcrop 3 --readtype 1D --logfile unm.nanofilt.log > unm.h5t3.fastq

#2 'U' to 'T' conversion

awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' mod.h5t3.fastq > mod.U2T.fastq
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' unm.h5t3.fastq > unm.U2T.fastq

#3 mapping to reference using minimap2

minimap2 -ax map-ont ref.fasta mod.h5t3.fastq | samtools view -bhS - | samtools sort -@ 6 - mod  && samtools index mod.bam
minimap2 -ax map-ont ref.fasta unm.h5t3.fastq | samtools view -bhS - | samtools sort -@ 6 - unm  && samtools index unm.bam

#4 calling variants for each single read-to-reference alignment

java -jar sam2tsv.jar -r  ref.fasta mod.bam > mod.bam.tsv
java -jar sam2tsv.jar -r  ref.fasta unm.bam > unm.bam.tsv

#5 convert results from step 4 and generate per_read variants information; the input file can be splitted based on read into smaller files to speed this step up.

per_read_var.stats.py mod.bam.tsv > mod.per_read.var.csv
per_read_var.stats.py unm.bam.tsv > unm.per_read.var.csv

#6 sumarize results from step 4 and generate variants information according the reference sequences (i.e., per_site variants); the input file can be splitted based on ref into smaller ones to speed this step up. 

per_site_var.py mod.bam.tsv > mod.per_site.var.csv
per_site_var.py unm.bam.tsv > unm.per_site.var.csv

#7 slide per_site variants with window size of 5, so that fast5 event table information can be combined

slide_per_site_var.py mod.ref.per_site.var.csv > mod.per_site.var.sliding.win.csv
slide_per_site_var.py unm.ref.per_site.var.csv > unm.per_site.var.sliding.win.csv

```


* To extract features from FAST5 files: 
``` 
#1 extract event table from fast5 files; 

fast5ToEventTbl.py input.fast5 > output.event.tbl

#2 extract features needed (esp. current intensity) for downstream analyses

event_tbl_feature_extraction.py output.event.tbl > output.event.tbl.features

#3 combine extracted features with per_read and per_site variants information

fastq_len.py h5t3.fastq > h5t3.fastq.len
adjust_read_base_positions.py  h5t3.fastq.len output.event.tbl.features number_of_chopped_leading_bases number_of_chopped_end_bases > output.event.tbl.features.readposition.adj.csv
assign_current_to_per_read_kmer.py output.event.tbl.features.readposition.adj.csv  > per_read.var.current.csv
per_read_kmer_intensity_to_per_site_kmer_intensity.py per_read.var.current.csv per_site.varsliding.win.csv > per_site.var.current.csv
```
* To build SVM and get predictions:
```
This step includes SVM training, prediction and performance assessment using single and multiple parameters.
The analyses are coded and performed in python3 programming environment. 
Uses: Pandas (0.23.4), sklearn (0.19.2) and numpy (1.15.1). 
The details can be found in "SVM_analyses" python jupyter-notebook. 

```

* To visulize results:
```
Python matplotlib.pyplot, seaborn 
R ggplot
```
## Citing this work:
If you find this work useful, please cite:

Huanle Liu, Oguzhan Begik, Morghan Lucas, Christopher E. Mason, Schraga Schwartz, John S. Mattick, Martin A. Smith and Eva Maria Novoa. Accurate detection of m6A RNA modifications in native RNA sequences. bioRxiv 2018

### License 
See LICENSE.md for details

### Contact
Please email us at huanle.liu@crg.eu or eva.novoa@crg.eu if you have any doubts/concerns/suggestions.
Thanks!

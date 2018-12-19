# EpiNano
Detection of RNA modifications from Oxford Nanopore direct RNA sequencing reads

### About EpiNano
EpiNano is a tool to identify RNA modifications present in direct RNA sequencing reads. 

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

## Running the software
* To extract features from basecalled FASTQ files: 
```
1. 'U' to 'T' conversion
2. nanofilt to trim leading and tailing bad quality bases
3. mapping to reference using minimap2
4. calling variants for each single read-to-reference alignment
5. slide results from step 4 with a window size of 5 and generate per_read variants information 
6. assign current intensity information from fast5 event table to per_read variants.
7. sumarize results from step 4 and generate variants information according the reference sequences (i.e., per_site variants)
8. slide per_site variants with window size of 5, so that fast5 event table information can be combined
```



* To extract features from FAST5 files: 
``` 
1. extract event table from fast5 files
fast5ToEventTbl.py input.fast5 > output.event.tbl
2. extract features needed (esp. current intensity) for downstream analyses
extract_feature_from_event_tbl.py output.event.tbl > output.event.tbl.features
3. combine extracted features with per_read and per_site variants information

```
* To build SVM and get predictions:
```
This step includes SVM training, perdiction and prediction performance assessment using single and multiple parameters.
The analyses are coded and performed in python3 programming environment. 
Pandas (0.23.4), sklearn (0.19.2) and numpy (1.15.1). 
The details can be found in "SVM_analyses" python jupyter-notebook. 
```

## Citing this work:
If you find this work useful, please cite:

Huanle Liu, Oguzhan Begik, Morghan Lucas, Christopher E. Mason, Schraga Schwartz, John S. Mattick, Martin A. Smith and Eva Maria Novoa. Accurate detection of m6A RNA modifications in native RNA sequences. bioRxiv 2018

### License 
See LICENSE.md for details

### Contact
Please email us at huanle.liu@crg.eu or eva.novoa@crg.eu if you have any doubts/concerns/suggestions.
Thanks!

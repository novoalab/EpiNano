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
- Support Vector Machine training & testing to predict m6A RNA modifications

## Getting Started and pre-requisites
These scripts are written in XXX. 

HUANLE PLEASE ADD DETAILS HERE OF PRE-INSTALLED SOFTWARE THAT is REQUIRED

## Running the software

* To extract features from FAST5 files: 
``` 
HUANLE PLEASE ADD DETAILS HERE 

```
* To extract features from basecalled FASTQ files: 
```
HUANLE PLEASE ADD DETAILS HERE

```

* To build SVM and get predictions:
```
HUANLE PLEASE ADD DETAILS HERE 

```

## Citing this work:
If you find this work useful, please cite:

Huanle Liu, Oguzhan Begik, Morghan Lucas, Christopher E. Mason, Schraga Schwartz, John S. Mattick, Martin A. Smith and Eva Maria Novoa. Accurate detection of m6A RNA modifications in native RNA sequences. bioRxiv 2018

### License 
See LICENSE.md for details

### Contact
Please email us at huanle.liu@crg.eu or eva.novoa@crg.eu if you have any doubts/concerns/suggestions.
Thanks!

# EpiNano
Detection of RNA modifications from Oxford Nanopore direct RNA sequencing reads

### About EpiNano
EpiNano is a tool to identify RNA modifications present in direct RNA sequencing reads. 
It is based on the detection of distinct 'features' that occur due to the presence of the modification, including changes in read quality, per-base quality, increased mismatch errors, increased deletion frequencies

## Testing
The software has been tested upon a set of 'unmodified' and 'modified' sequences containing m6A at known sites or A.
Fasta sequences used in the training/testing are included in the release.

### What's included
- Scripts to extract features from FAST5 files (e.g. current intensity)
- Scripts to process basecalled FASTQ files into kmer pileups (similar to samtools mpileup format but for 5mer sequences). 
- Support Vector Machine training/testing to predict m6A RNA modifications

### Getting Started
These scripts are written in XXX. 

HUANLE PLEASE ADD DETAILS HERE OF SOFTWARE REQUIRED

### Prerequisites

HUANLE PLEASE ADD DETAILS HERE OF PRE-INSTALLED SOFTWARE THAT is REQUIRED

### Running the demo

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

### Citing this work:
If you find this work useful, please cite:

Huanle Liu, Oguzhan Begik, Morghan Lucas, Christopher E. Mason, Schraga Schwartz, John S. Mattick, Martin A. Smith and Eva Maria Novoa. Accurate detection of m6A RNA modifications in native RNA sequences. bioRxiv 2018


### License 
See LICENSE.md for details

### Contact
Please email us at huanle.liu@crg.eu or eva.novoa@crg.eu if you have any doubts/concerns/suggestions.
Thanks!

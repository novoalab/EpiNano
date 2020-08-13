<img src="image/Epinano_logo_v1.jpg" width='300' height='120'>

Detection of RNA modifications from Oxford Nanopore direct RNA sequencing reads  

## Upgrades

**EpiNano 1.2** - latest version, includes pretrained m6A models base-called with *Guppy* v 3.1.5

* NOTE: WE ARE CURRENTLY UPGRADING THE REPO FROM EPINANO 1.1 for 1.2 -- AT THE MOMENT YOU WILL DOWNLOAD 1.1

* This version of EpiNano allows to make predictions using two different strategies: *EpiNano-Error* and *EpiNano-SVM*. 

* This version now includes modules for visualization of your RNA modification predictions


<img src="image/EpiNano_modes.png" width='800' height='400'>

*EpiNano-Error* can only be run in pairwise mode (e.g. WT and KO or KD). It combines the different types of base-calling errors that appear in a given dataset (mismatches, deletions, insertions) as well as alterations in per-base-calling qualities. RNA modification predictions are based on the differences in error patterns observed in two matched samples. This strategy can be used with FASTQ data base-called with any given base-calling algorithm version

*EpiNano-SVM* can be run using either pre-trained models for a given RNA modification, or by building your own models. The first versions of Epinano (1.1. and 1.0) allowed the user to run the code in standalone mode (i.e. only 1 condition). However, we should note that using a matched control (e.g. KO or KD) is still highly recommended, due to the noisy nature of direct RNA sequencing reads, which are 'error'-rich. Therefore, EpiNano-SVM (version 1.2) now requires a matched control to be provided in the command line. 


**EpiNano 1.1** - a new and slim version, written in python3 has been  released, which is available [here](https://github.com/enovoa/EpiNano/releases).

* This version is the one currently implemented in [MasterOfPores](https://www.frontiersin.org/articles/10.3389/fgene.2020.00211), a workflow to analyze direct RNA sequencing data. 

* The major differences with EpiNano 1.0 are (i) it is much faster (ii) Uses python3 instead of python2 (iii) Does not extract current intensity in the feature table, as this feature was not used to train the final models.

*  Includes pre-trained m6A models base-called with *Albacore* version 2.1.7. 

* Works both with Guppy and Albacore, but the SVM predictions will be only accurate if your data has been base-called using Albacore 2.1.7.

* If using Guppy base-called data,  EpiNano can still be used as a toolkit to extract per 5-mer base-calling 'errors' (mismatch, insertion, deletion, quality), which are a proxy of RNA modifications present in a given dataset. We recommend running EpiNano in paired mode, i.e. computing the features in two datasets (WT-KO) to then accurately predict the RNA modified sites (i.e. those showing largest differences in their base-calling 'error' features). 


**EpiNano 1.0** - original code used in [Liu, Begik et al., Nature Comm 2019](https://www.nature.com/articles/s41467-019-11713-9), which is available [here](https://github.com/enovoa/EpiNano/releases).

* Includes pre-trained m6A models base-called with *Albacore* version 2.1.7. 

* It extracts both base-calling 'errors' (mismatch, insertion, seletion, per-base quality) as well as current intensity values

* Current intensity information is extracted from the base-called Albacore FAST5 files.

* Does not work with Guppy base-called datasets.


## About EpiNano

EpiNano is a tool to identify RNA modifications present in direct RNA sequencing reads.

EpiNano will extract a set of 'features' from direct RNA sequencing reads, which will be in turn used to predict whether the 'error' is caused by the presence of an RNA modification o r not. Features extracted include:
- k-mer current intensity
- read quality
- per-base quality
- per-base mismatch frequency
- per-base deletion frequency
- per-base insertion frequency

EpiNano version 1.2 can predict RNA-modified sites in two different ways:

1. **EpiNano-Error** 
* Base-calling algorithm independent. 
* Applicable to any given RNA modification that causes an effect in the base-calling features. 

2. **EpiNano-SVM** 
* Base-calling algorithm dependent (data must be base-called with Guppy 3.1.5) 
* Can use both base-calling error features as well as current intensity features
* It can be used to train your own models as well as be applied to datasets for which a pre-trained model is available (m6A)
* The available m6A SVM models has been trained and tested upon a set of 'unmodified' and 'modified' sequences containing m6A at known sites or A. 


## Considerations when using EpiNano

* EpiNano relies on the use of base-calling 'errors' to detect RNA modifications; however, direct RNA sequencing base-calling produces a significant amount of 'errors' in unmodified sequences. Therefore, to obtain higher confidence m6A-modified sites, we recommend to sequence both modified and unmodified datasets (e.g. treated with demethylase, or comparing a wild-type vs knockout/knockdown). Coupling a "control" (KD/KO) is not required in earlier Epinano versions, but is highly recommended.
* You can use EpiNano as a feature extractor to predict RNA modifications based on alterations in base-called features (as used [here](https://www.biorxiv.org/content/10.1101/2020.07.06.189969v2)), as well as use the pre-trained SVMs to detect m6A RNA modifications (as used [here](https://www.nature.com/articles/s41467-019-11713-9)). In the latest version of EpiNano (1.2) we provide scripts to predict RNA modifications using both modes (EpiNano-Error and EpiNano-SVM). 
* EpiNano does not have per-read resolution. We are currently working on an improved version of EpiNano to obtain predictions at per-read level.
* The performance of the algorithm is dependent on the stoichiometry of the site (i.e. sites with very low stoichiometry will be often missed by the algorithm)
* Pre-trained models to predict m6A sites are included in each release. Please note that if you use pre-trained m6A models, your data should be base-called with the SAME base-calling algorithm and version (i.e. Guppy 3.1.5 if you use EpiNano 1.2, and Albacore 2.1.7 if you use EpiNano 1.0 or 1.1).
*  If you are using a different base-calling algorithm version, you can still use EpiNano to extract features (i.e. 'errors'), but the SVM predictions (ProbM) will not be accurate.



## Pre-requisites
The following softwares and modules were used by EpiNano

| Software  | Version |
| ------------- | ------------- |
| NanoFilt  | 2.2.0  |
| minimap2  | 2.14-r886  |
| samtools  | 0.1.19  |
| sam2tsv  | a779a30d6af485d9cd669aa3752465132cf21eec   |
| python  | 3.6.7   |
| java openjdk   | 1.8.0    |
| h5py  | 2.8.0    |
| numpy  | 1.15.4     |
| pandas  | 0.23.4        |
| sklearn  | 0.20.2     |
| nanopolish | 0.12.4 |


## Getting the code

To download the latest version EpiNano version 1.2, you just need to clone the repo:

```
git clone git@github.com:enovoa/EpiNano.git 
```

You can choose to download EpiNano 1.1 [HERE](https://github.com/enovoa/EpiNano/releases)

You can choose to download EpiNano 1.0 [HERE](https://github.com/enovoa/EpiNano/releases)

### What's included in the releases
- Scripts to extract features from FAST5 files
- Scripts to process mapped BAM files into kmer pileups (similar to samtools mpileup format but for 5mer sequences)
- Support Vector Machine training (SVM) & testing to predict m6A RNA modifications

## Running the code: EpiNano 1.2

### STEP 1. Extract base-calling error features
*Epinano_Variants.py*, outputs two feature tables: 
(i) sample.per.site.var.csv, which contains base-calling ‘error’ information for each reference position, and 
(ii) sample.per_site.5mer.csv, which contains the same base-called features organized in 5-mer windows.

```
$ python $EPINANO_HOME/Epinano_Variants.py -h

Usage: Epinano_Variants.py [-h] [-r READS] [-R REFERENCE] [-b BAM] [-f FILE]
                           [-t THREADS] [-s SAM2TSV] [-T TYPE] [-p]

optional arguments:
  -h, --help            show this help message and exit
  -r READS, --reads READS
                        fastq(a) reads input
  -R REFERENCE, --reference REFERENCE
                        samtools faidx indexed reference file
  -b BAM, --bam BAM     bam file; if given; no need to offer reads file;
                        mapping will be skipped
  -f FILE, --file FILE  tsv file generated by sam2tsv.jar; if given, reads
                        mapping and sam2tsv conversion will be skipped
  -t THREADS, --threads THREADS
                        number of threads
  -s SAM2TSV, --sam2tsv SAM2TSV
                        /path/to/sam2tsv.jar; needed unless a sam2tsv.jar
                        produced file is already given
  -T TYPE, --type TYPE  reference types, which is either g(enome) or
                        t(ranscriptome);
  -p, --per_read_variants
                        compute per reads variants statistics
```


The '--type' flag indicates the type of reference that was used to obtain the bam file. If the reads were mapped to a genome reference with splicing mapping, '--type g' should be specified, and EpiNano will distinguish and process separately the reads mapped to forward and reverse strands. Otherwise, by default, the script assumes the bam file was generated by mapping the reads to reference transcriptome and the reads were only mapped to the forward strand. 

Example: 
```
python $EPINANO_HOME/Epinano_Variants.py -t 6 -R reference.fasta -b sample.reads.bam -s /path/to/sam2tsv/sam2tsv.jar --type t 
```

Note 1: it is possible to organize the variants in any kmer length using *$EPINANO_HOME/misc/Slide_Variants.py*.  
Note 2: the users should split the computations for each  reference sequences if the reference genome is large.

### STEP 2. Extraction of current intensity values (optional - only needed if using EpiNano-SVM)

**Epinano_Current.py** uses [Nanopolish]((https://nanopolish.readthedocs.io/en/latest/quickstart_eventalign.htmlhttps://) to extract current signal level information and then collapses it on a single position basis. 


Note 1: Please add the /path/to/nanopolish to environmental **$PATH** variable, otherwise the script will fail.

```
sh EPINANO_HOME/Epinano_Current.sh -h

Epinano_Current.sh [-h] [-b bam -r reads -f genome/transriptome reference -d fast5dir -t threads -m bam_file_mapping_type]

        it runs nanopolish eventalign; aggreagets current intensity values associated with single positions

        nanopolish, samtools have to be installed and added to environmental paths!!

        -h [this help text]
        -b [bam file]
        -r [fastq reads file]
        -f [reference file used to produce bam file]
        -t [number of threads]
        -m [t: reads mapping performed using reference transcriptome; g: reads mapping performed with reference genome]


```

Example:

```
sh $EPINANO_HOME/Epinano_Current.sh -b sample.reads.bam -r sample.reads.fastq -f reference.fasta  -t 6  -m t -d fast5_folder/ 
```


### STEP 3. Predict RNA modifications 
EpiNano offers two alternative methods to predict RNA modifications: 

i) **EpiNano-Error** uses of the variants/error frequencies computed above 

ii) **EpiNano-SVM** uses an SVM algorithm to train models and predict modifications. 

#### a) Predicting RNA modifications using EpiNano-Error:

 *Epinano_DiffErr.R* uses  ordinary least square (OLS) to fit a linear regression model between paired unmodified and modified samples, and then detects outliers, i.e., observations with large residuals, which tend to be underlined by base modifications. 

Note 1: different types of RNA base modification show distinct biases toward the spefic types of errors. Thus, offered *Epinano_sumErr.py* to combine mismatches, indels and even quality scores. Just like the independent types of errors, the combined error is internally performed when running *Epinano_ErrDiff.R*. 

```
$ Rscript $EPINANO_HOME/misc/Epinano_DiffErr.R

Usage:
   Epinano_Diff_Err.R [-k <sample1> -w <ssample2> -x <Label1> -y <label2> -o <output prefix> -s <feature> -l 0 -p <type of plot> -c <covergae/depth> -t <threshold>]

Examples:
        #1 scatter plot without marking candidates
        Rscript Epinano_Diff_Err.R -k ko1.plus_strand.per.site.var.csv -w wt1.plus_strand.per.site.var.csv -x ko -y wt -s mis -o scatter -l 0 -p c -c 30 -t 5

        #2 bar-plot with candidates marked
        Rscript Epinano_Diff_Err.R -k ko1.plus_strand.per.site.var.csv -w wt1.plus_strand.per.site.var.csv -x ko -y wt -s mis -o scatter -l 0 -p c -c 30 -t 5

        #3 scatter-plot using sum_err
        Rscript Epinano_Diff_Err.R -k ko1.sum_err.csv -w wt1.sum_err.csv -x ko -y wt -s sum_err -o scatter_sum_err -l 0 -p c -c 30 -t 5

Options:
   -k sample1
   -w sample2
   -x label1
   -y label2
   -s feature (e.g, mis, ins, del, sum_err) to use for plot
   -o output prefix
   -l 0: do NOT label outliers; 1: label outliers [default:0]
   -p c: scatter-plot; b; barplot
   -c minimum coverage/depth; default: 30
   -t threshold: if a number *N* is used, N * sd (i.e., standard deviation [Ordinary Least Squares regression (OLS)] of ) will be used as threshold to determine candidates; default: 5
```

#### b) Predicting RNA modifications using EpiNano-SVM:

EpiNano v 1.2 includes pre-trained models (found in *$EPINANO_HOME/models/*), which have been trained using synthetic molecules (curlcakes) with and without introduced m6A modiifcations. However, the user can train their own models using *EpiNano_Predict.py*, employing  the features generated with *EpiNano_Variants.py* and/or *EpiNano_Current.py* as shown in the previous steps.

```
python $EPINANO_HOME/Epinano_Predict.py -h

usage: Epinano_Predict.py [-h] [-k KERNEL] [-o OUT_PREFIX] [-a] [-M MODEL]
                          [-t TRAIN] [-mc MODIFICATION_STATUS_COLUMN] -p
                          PREDICT -cl COLUMNS

optional arguments:
  -h, --help            show this help message and exit
  -k KERNEL, --kernel KERNEL
                        kernel used for training SVM, choose any one from
                        'linear', 'poly', 'rbf', 'sigmoid'; if no choice made,
                        all 4 kernels will be used
  -o OUT_PREFIX, --out_prefix OUT_PREFIX
                        output file prefix
  -a, --accuracy_estimation
                        '-a' performs accuracy estimation with known modified
                        status from --predict file; only feasible when there
                        are prior knolwdge of modiifcaton status indiciated by
                        --modification_status_column
  -M MODEL, --model MODEL
                        pre-trained model that can ben used for prediction
  -t TRAIN, --train TRAIN
                        file name of feature table used for training; can be
                        gzipped
  -mc MODIFICATION_STATUS_COLUMN, --modification_status_column MODIFICATION_STATUS_COLUMN
                        column number from (input file1, i.e, traing file)
                        that contains modification status information

required arguments:
  -p PREDICT, --predict PREDICT
                        file name of feature table used for making predictions
                        or testing accuracy; can be gzipped. when this file is
                        the same the one used for training, half of the data
                        will be chosen for training.
  -cl COLUMNS, --columns COLUMNS
                        comma seperated column number(s) that contain features
                        used for training and prediciton

```

Example: 

```
python $EPINANO_HOME/Epinano_Predict.py 
          --train ko_wt_combined.per_site_raw_feature.rrach.5mer.csv  
          --predict ko_wt_combined.per_site_raw_feature.rrach.5mer.csv  
          --accuracy_estimation --out_prefix train_and_test 
          --columns 8,13,23 --modification_status_column 26
```

While the user can choose to train the algorithm with one sample (--train) and test it on an independent sample (--predict),  it is also possible to use the same input file both for training and testing the model, as depicted in the example above. In this scenario, Epinano_Predict will train the models with 50% of the input data, and make predictions with the remaining 50% of the data.

In the above command, ‘--columns’ denotes the column numbers of features that are used for training models (in this case, corresponding to ‘q3’, ’mis3’ and ‘del3’), while ‘--modification_status_column’ indicates the prior knowledge of the modification statuses, i.e., the labels ‘mod’ and ‘unm’. Switching on --accuracy_estimation will report the accuracy of the trained model(s).  Unless ‘--kernel’ is used, Epinano_Predict will train models with multiple kernels. 

With the trained models, the user can make predictions of modifications.
```
python $EPINANO_HOME/Epinano_Predict.py 
--model q3.mis3.del3.MODEL.linear.model.dump 
--predict some_sample.per_site.5mer.csv 
--columns 8,13,23
--out_prefix some_sample.modification 
```
In the command above, we employ a previously trained model ‘q3.mis3.del3.MODEL.linear.model.dump’ that will predict m6A modifications in RRACH k-mers on a dataset that is specified with ‘--predict’. Please note that if using this pre-trained model, the csv datasets should be filtered to include only RRACH k-mers, which were the k-mers used for training this model.

## Running the code: EpiNano 1.1

* Build feature table (on which predictions will be made)

    For step-by-step instructions to build a feature table, please take a look at the [Wiki](https://github.com/enovoa/EpiNano/wiki)

* To train SVM and perform predictions:
```
This step includes SVM training, prediction and performance assessment using single and multiple features.
$ python3 SVM.py -h

Commad:  scripts/SVM.py -h
usage: SVM.py [-h] [-k KERNEL] [-o OUT_PREFIX] [-a] [-M MODEL] [-t TRAIN]
              [-mc MODIFICATION_STATUS_COLUMN] -p PREDICT -cl COLUMNS

optional arguments:
  -h, --help            show this help message and exit
  -k KERNEL, --kernel KERNEL
                        kernel used for training SVM, choose any one from
                        'linear', 'poly', 'rbf', 'sigmoid'; if no choice made,
                        all 4 kernels will be used
  -o OUT_PREFIX, --out_prefix OUT_PREFIX
                        ouput file prefix
  -a, --accuracy_estimation
                        -a perform accuracy estimation with known modified
                        status from --predict file
  -M MODEL, --model MODEL
                        pre-trained model that can ben used for prediction; if
                        this is not available SVM model will be trained and
                        dumped; there can be multiple models, which should be
                        in the same order as kernels applied
  -t TRAIN, --train TRAIN
                        file name of feature table used for training
  -mc MODIFICATION_STATUS_COLUMN, --modification_status_column MODIFICATION_STATUS_COLUMN
                        column number from (input file1, i.e, traing file)
                        that contains modification status information

required arguments:
  -p PREDICT, --predict PREDICT
                        file name of feature table used for making predictions
                        or testing accuracy. when this file is the same the
                        one used for training, half of the data will be chosen
                        for training.
  -cl COLUMNS, --columns COLUMNS
                        comma seperated column number(s) that contain features
                        used for training and prediciton

```                        

### Example

With the example svm input files from example/svm_input folder:

* training with the example feature tables
```
	# the command below will train models using all quality scores are all positions from sample1 and make prediction on sample2
	python3 SVM.py -a -t sample1.csv -p sample2.csv -cl 1-5 -mc 11 -o test

	# while this command will do the same thing except choosing a 'linear' kernel for SVM training
	python3 SVM.py -a -k linear -cl 1-5 -t sample1.csv -p sample2.csv -mc 11

	# this 3rd command uses base quality and mismatch frequencies of the centred bases for SVM training
	python3 SVM.py -t sample1.csv -p sample2.csv -cl 3,7 -mc 11
```
* predict modifications   
```
	#use previously trained model and epinano-scripts-generated site-wise feature table to make predictions
	python3 SVM.py -a -M M6A.mis3.del3.q3.poly.dump -p test.csv -cl 7,12,22 -mc 28 -o pretrained.prediction
```


## Further Documentation

Please check the [Wiki](https://github.com/enovoa/EpiNano/wiki) for additional information on usage on EpiNano 1.1

## Citing this work:
If you find this work useful, please cite:

Huanle Liu*, Oguzhan Begik*, MorghanC  Lucas, Jose Miguel Ramirez, Christopher E. Mason, David Wiener, Schraga Schwartz, John S. Mattick, Martin A. Smith and Eva Maria Novoa. Accurate detection of m6A RNA modifications in native RNA sequences . Nature Communications 2019, 10:4079.

Link to paper: https://www.nature.com/articles/s41467-019-11713-9

### License
See LICENSE.md for details

### Contact
Please read the [Wiki](https://github.com/enovoa/EpiNano/wiki) before opening an issue. Also, please go over other [issues](https://github.com/enovoa/EpiNano/issues) that may have been previously resolved (check out "closed" issues).
If you still have doubts/concerns/suggestions, please open a new Issue.
Thanks!

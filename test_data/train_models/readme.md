### Raw data to train models

Raw data used to train the SVM models (i.e. m6A-modified and unmodified 'curlcakes') can be found at https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP174366

Reference sequences can be found in cc.fasta

### Processed data to train models 

In this folder you will also find error and electric signal features extracted with Epinano_Variants and Epinano_Current_Intensity scripts.

### How to train my own SVM model?

In this folder you fill find the folder **train_test.sh** with commands to train SVM models and test their accuracies.

### How to train and test my own models?

You can also find in ../make_predictions/run.sh how to prepare feature table for training and making predictions. 

### Additional observations

- *middleAs* means [GCT][GCT]A[GCT][GCT] kmers
- All A-bases are modified from the mod samples (i.e. curlcakes are modified at 100% stoichiometry).

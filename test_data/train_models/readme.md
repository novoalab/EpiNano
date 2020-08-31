## Additional documentation and examples:  EpiNano 1.2 

### How to train my own SVM model?

In this folder you will find the script [train_test.sh](https://github.com/enovoa/EpiNano/tree/master/test_data/train_models/train_test.sh), with commands illustrating how to train SVM models and test their accuracies.

### How to prepare the input data to train and test my own models?

Please take a look at script [run.sh](https://github.com/enovoa/EpiNano/tree/master/test_data/make_predictions/run.sh), to see examples on different commands showcasing how to prepare feature tables for training models, as well as making predictions. 

### Raw data used to train EpiNano 1.2 SVM models 

Raw data used to train the SVM models (i.e. m6A-modified and unmodified 'curlcakes') can be found at https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP174366

Reference sequences can be found in cc.fasta

### Processed data used to train EpiNano 1.2 SVM models 

In this folder, you will also find error and electric signal features extracted with *Epinano_Variants* and *Epinano_Current* scripts.

### Additional observations regarding  EpiNano 1.2 SVM models 

- *middleAs* means [GCT][GCT]A[GCT][GCT] kmers
- All A-bases are modified from the mod samples (i.e. curlcakes are modified at 100% stoichiometry).

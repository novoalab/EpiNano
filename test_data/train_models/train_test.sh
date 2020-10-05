#!/usr/bin/env bash

echo "in the follwoing examples we have prior knowlwdgw of modification status"
echo "please refer to  ../make_predicitons/run.sh for how to prepare feature tables for training or making predicitons !"
echo "combine and label mod and unm samples "

sh ../../misc/Epinano_LabelSamples.sh -m rep1.mod.middleAs.5mer.csv -u rep1.unm.middleAs.5mer.csv -o rep1.raw.5mer.csv
sh ../../misc/Epinano_LabelSamples.sh -m rep2.mod.middleAs.5mer.csv -u rep2.unm.middleAs.5mer.csv -o rep2.raw.5mer.csv

echo "create delta features; This is only feasible when you have paired modified and unmodified samples!"
python ../../misc/Epinano_make_delta.py rep1.mod.middleAs.5mer.csv rep1.unm.middleAs.5mer.csv 30 5 > rep1.delta.mod.5mer.csv
python ../../misc/Epinano_make_delta.py rep2.mod.middleAs.5mer.csv rep2.unm.middleAs.5mer.csv 30 5 > rep2.delta.mod.5mer.csv
python ../../misc/Epinano_make_delta.py rep1.unm.middleAs.5mer.csv rep2.unm.middleAs.5mer.csv 30 5 > delta.unm.5mer.csv

sh ../../misc/Epinano_LabelSamples.sh -m rep1.delta.mod.5mer.csv -u delta.unm.5mer.csv -o rep1.delta.5mer.csv
sh ../../misc/Epinano_LabelSamples.sh -m rep2.delta.mod.5mer.csv -u delta.unm.5mer.csv -o rep2.delta.5mer.csv


echo "make sumErr feature"
python3.6 ../../misc/Epinano_sumErr.py --file rep1.mod.middleAs.5mer.csv --kmer 5 --out rep1.mod.middleAs.sumErr.csv
python3.6 ../../misc/Epinano_sumErr.py --file rep2.mod.middleAs.5mer.csv --kmer 5 --out rep2.mod.middleAs.sumErr.csv
python3.6 ../../misc/Epinano_sumErr.py --file rep1.unm.middleAs.5mer.csv --kmer 5 --out rep1.unm.middleAs.sumErr.csv
python3.6 ../../misc/Epinano_sumErr.py --file rep2.unm.middleAs.5mer.csv --kmer 5 --out rep2.unm.middleAs.sumErr.csv

sh ../../misc/Epinano_LabelSamples.sh -m rep1.mod.middleAs.sumErr.csv -u rep1.unm.middleAs.sumErr.csv -o rep1.middleAs.sumErr.csv
sh ../../misc/Epinano_LabelSamples.sh -m rep2.mod.middleAs.sumErr.csv -u rep2.unm.middleAs.sumErr.csv -o rep2.middleAs.sumErr.csv

echo "make delta-sumErr feature"
python ../../misc/Epinano_delta_sumErr.py  rep1.mod.middleAs.sumErr.csv rep1.unm.middleAs.sumErr.csv  rep1.delta_sumErr.5mer.csv
python ../../misc/Epinano_delta_sumErr.py rep2.mod.middleAs.sumErr.csv rep2.unm.middleAs.sumErr.csv  rep2.delta_sumErr.5mer.csv
python ../../misc/Epinano_delta_sumErr.py rep1.unm.middleAs.sumErr.csv rep2.unm.middleAs.sumErr.csv  control.delta_sumErr.5mer.csv
sh  ../../misc/Epinano_LabelSamples.sh -m rep2.delta_sumErr.5mer.csv -u control.delta_sumErr.5mer.csv -o delta_sumErr.5mer.for_train.csv 
awk -F',' '{if (/^#/)print $0",sample"; else print $0",mod"}'  rep1.delta_sumErr.5mer.csv >  delta_sumErr.5mer.for_test.csv
echo "train models and test it with delta-sumErr feature"
echo "train with single site delta_sum_err"
python3.6 ../../Epinano_Predict.py  -mc 10 -cl 7 -t delta_sumErr.5mer.for_train.csv  -p delta_sumErr.5mer.for_test.csv  -o single_delta_sumErr -a 
echo "train with 5 consecutive sites delta_sum_err"
python3.6 ../../Epinano_Predict.py  -mc 10 -cl 5-9 -t delta_sumErr.5mer.for_train.csv  -p delta_sumErr.5mer.for_test.csv  -o 5mer_delta_sumErr -a 


echo "train SVM models and assess prediction accuracies using raw error features"
python3.6 ../../Epinano_Predict.py -o tr2pr1 -mc 26 -cl 8,13,23 -t rep2.raw.5mer.csv -p rep1.raw.5mer.csv -a
echo "tr1pr2.q3.mis3.del3.MODEL.linear.model.dump is the trained model using mis3, q3 and del3 from rep2(mod2 + unm2) "
echo  "tr1pr2.q3.mis3.del3.MODEL.kernel.linear.csv contains predicitons on rep1(mod1 + unm1) with the trained model"
echo "now you can use this model to make predictions."

echo "train with delta features"
python3.6 ../../Epinano_Predict.py -o tr2pr1 -mc 26 -cl 8,13,23 -t rep2.raw.5mer.csv -p rep1.raw.5mer.csv  -a 

echo "train with sumErr feature"
python3.6 ../../Epinano_Predict.py -o tr2pr1 -mc 11 -cl 8 -p rep1.middleAs.sumErr.csv -t rep2.middleAs.sumErr.csv -a  

echo "It is possible to use any kind of feature(s) to train models and make predictions"



echo "the *accray file contains accuracies of different models, which should help you choose the *best* one!"

#the 6 models are 
rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump: model was trained with difference between Q3,Mis3, and Del3 of modified and un-modified samples.

rrach.q3.mis3.del3.linear.dump: modle was trained with Q3, Mis3 and Del3

tr2pr1.sum_err.MODEL.rbf.model.dump: model trained with a combination of mis, ins, and del 

tr2pr1.sum_err1.sum_err2.sum_err3.sum_err4.sum_err5.MODEL.linear.model.dump: model trained with sum_err in 5-mer model
 
tr2pr1.delta_sum_err3.MODEL.rbf.model.dump: model trained with single-site delta_sum_err feature

tr2pr1.delta_sum_err1.delta_sum_err2.delta_sum_err3.delta_sum_err4.delta_sum_err5.MODEL.linear.model.dump: model trained with delta_sum_err at consecutive 5-sites

# ~~~~~~~~~~~~~~
sum_err feature can be derived using *Epinano_sumErr.py*, with --kmer 0 if input is per.site.var.csv or --kmer 5 if input is per.site.var.5mer.csv
delta_sum_err can be computed with *Epinano_delta_sumErr.py*

## feel free to train models on your preferred motifs and features !!! 

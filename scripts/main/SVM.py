#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn import datasets, svm
from sklearn.model_selection import train_test_split
from collections import defaultdict
import sys
import os
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-k","--kernel",  help="kernel used for training SVM, choose any one from 'linear', 'poly', 'rbf', 'sigmoid'; if no choice made, all 4 kernels will be used")
requiredGrp = ap.add_argument_group('required arguments')
requiredGrp.add_argument("-f1","--train", required=True, help="file name of feature table used for training")
requiredGrp.add_argument("-f2","--test", required=True, help="file name of feature table used for testing")
requiredGrp.add_argument('-c','--columns',required=True,help = "comma seperated column number(s) that contain features used for training and prediciton")

args = vars (ap.parse_args())
'''
print (args['kernel'])
print (args['train'])
print (args['test'])
print (args['columns'])
usgae: 
python3 SVM_genral.py rep1_mod_unm_per_site_slided_current.csv  8,11 0
input1 MOD_UNM (file contain labeled mod and unm information)
input2 column index of parametes used for train and prediction
input3  test_dataset, which has to be mod and unm combined (input is 0 or file name)
input4   kernel signified by user
'''
def evaluate_on_test_data (test, predicitons):
    correct_classifications = 0
    for i in range (len(y_test)):
        if predictions[i] == test[i]:
            correct_classifications += 1
    accuracy = correct_classifications * 100 / len (test)
    return accuracy

def plot_ROC (y_test, probas,fig_out,kn, parameter):
    fpr, tpr, thresholds = roc_curve(y_test,probas[:,1])
    roc_auc = auc(fpr,tpr)
#    plt.clf()
    plt.plot (fpr, tpr, label = 'ROC Curve for ' + 'kernel '+kn +' (Area under ROC = %0.2f)'% roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC of Modification Prediction based On ' + parameter)
    plt.legend(loc="lower right")
    plt.savefig(fig_out)

m_u_var =args['train'] 
m_u_var_df = pd.read_csv(m_u_var)

columns_to_consider = map (int, args['columns'].split(',')) ### 1-based column index e.g.  8,11 (q, mis)
cols = [i -1 for i in columns_to_consider]

if args['test'] != args['train']:
    test2 = pd.read_csv (args['test'])
names = list (m_u_var_df.columns[cols])
out_prefix = '.'.join (names) + '.SVM'
df_tmp = m_u_var_df.dropna (subset=names)
df_tmp = df_tmp.reset_index(drop=True) # do not insert index as extra column
test_tmp = test2.dropna(subset=names)
test_tmp = test_tmp.reset_index(drop=True)

X = df_tmp.iloc[:,cols]
indices = df_tmp.index
Y = df_tmp[['sample']]

old_header = ",".join (df_tmp.columns.values)
if args['test'] == args['train']:
    X_train, X_test, y_train, y_test, indices_train, indices_test = train_test_split(X,Y.values.ravel(), indices, test_size=0.25, random_state= 100)
else:
    X_train, X_test, y_train, y_test, indices_train, indices_test = train_test_split(X,Y.values.ravel(), indices, test_size=0, random_state= 100)

if args['test'] != args['train']:
    X_test = test_tmp.iloc[:,cols]
    indices_test =  X_test.index
    y_test = test_tmp[['sample']]
    y_test = y_test.values

kernels = ('linear', 'poly', 'rbf', 'sigmoid')  if args['kernel'] is None else (args['kernel'],)

accuracies = {}
probabilities = defaultdict (list)

for _, kn in enumerate (kernels):
    model = svm.SVC(kernel=kn, probability=True)
    model_fit = model.fit (X_train, y_train)
    predictions = model_fit.predict (X_test) #if sys.argv[5] == '0' else model_fit.predict()
    outh = open( out_prefix+'.kernel.' + kn + '.csv','w')
    print (old_header + ',prediction'+',dist,ProbM,ProbU',file=outh)
    for t in range (len(indices_test)):
        idx = indices_test[t]  
        original_line = ",".join(map(str, df_tmp.iloc[idx].values)) if sys.argv[3] == '0' else ",".join (map (str, test_tmp.iloc[idx].values))
        dist = map (str, model.decision_function([X_test.iloc[t]]))
        probM, probU = map (str, model.predict_proba ([X_test.iloc[t]])[0])
        print (original_line + ',' + predictions[t] +',' + ",".join (dist) + ',' +  probM + ',' + probU, file = outh)
    accuracies[kn] = evaluate_on_test_data(y_test, predictions)
acc_sort = sorted (accuracies.items(), key = lambda kv:kv[1])
best_kn  = acc_sort[-1][0]
best_acc = acc_sort[-1][1]
best_prediciton = out_prefix + '.best-kernel.' + best_kn + '.accuracy' 
print ("Best accuracy {} %  obtained with kernel = {}".format(best_acc,best_kn), file = open (best_prediciton,'w'))
del accuracies[best_kn]
for k,v in accuracies.items ():
		print (" {} % accuracy obtained with kernel = {}".format(v,k), file = open (best_prediciton,'a'))

#####################################################################################
############# remove prdiciton results based on kernels that are less accurate ######
#####################################################################################
for _, kn in enumerate(kernels):
    if kn != best_kn:
        os.remove(out_prefix +'.kernel.'+kn + '.csv')
        #os.rename (out_prefix +'.kernel.'+kn + '.csv', out_prefix + '.best-kernel.' + kn + '.csv')
####################################

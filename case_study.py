# -*- coding: utf-8 -*-
"""
Created on 2019-7-16

@author: Jun Zhang
"""
import os
import argparse
import numpy as np
from Bio import SeqIO
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score, precision_score
from sklearn.metrics import matthews_corrcoef, auc, roc_curve, classification_report, precision_recall_curve
from sklearn.metrics import coverage_error
from sklearn.metrics import label_ranking_loss
from sklearn.metrics import hamming_loss
from sklearn.metrics import f1_score
import numpy.random as r
import math
def lauc_score(true_labels, pred_scores, fpr_c):
    fpr, tpr, tresholds = roc_curve(true_labels, pred_scores, sample_weight=None)
    # fpr_c = 0.05
    lowc = 0
    for i in range(len(fpr)):
        if fpr[i] >= fpr_c:
            lowc = i
            break
    lauc = auc(fpr[:lowc], tpr[:lowc])#*2/(fpr_c*fpr_c)
    return lauc

def ratio_recall_curve(labels1, labels2, scores):
    print len(labels1), len(labels2), len(scores)
    tpr = []
    cro = []
    thr = []
    for c in np.arange(1.0, -0.001, -0.001):
        pred_labels = []
        for e in scores:
            if e >= c:
                pred_labels.append(1.0)
            else:
                pred_labels.append(0.0)
        tpr.append(recall_score(labels1, pred_labels))
        cro.append(recall_score(labels2, pred_labels))
        thr.append(c)
    print len(tpr), len(cro)
    return tpr, cro, thr

label_dna = []
label_rna = []
label_disorder = []
prots = list(SeqIO.parse('data/pdb_example_8.txt', 'fasta'))
for prot in prots:
    with open('examples/'+prot.id+'_label.txt', 'r') as f:
        for line in f:
            label_dna.append(float(line.strip().split()[0]))
            label_rna.append(float(line.strip().split()[1]))
            label_disorder.append(float(line.strip().split()[2]))
rna_pred1 = []
dna_pred1 = []
rna_pred2 = []
dna_pred2 = []
rna_pred3 = []
dna_pred3 = []
rna_scores1 = []
dna_scores1 = []
rna_scores2 = []
dna_scores2 = []
rna_scores3 = []
dna_scores3 = []
disorder_scores = []
for prot in prots:
    if not os.path.isfile('result/nucbind/' + prot.id + '_DNA.txt'):
        for i in range(len(prot.seq)):
            dna_scores3.append(0.0)
            dna_pred3.append(0.0)
            rna_scores3.append(0.0)
            rna_pred3.append(0.0)
    else:
        with open('result/nucbind/' + prot.id + '_DNA_svm.txt') as f:
            for line in f.readlines()[1:]:
                dna_scores3.append(float(line.strip().split()[3]))
                dna_pred3.append(float(line.strip().split()[2]))
        with open('result/nucbind/' + prot.id + '_RNA_svm.txt') as f:
            for line in f.readlines()[1:]:
                rna_scores3.append(float(line.strip().split()[3]))
                rna_pred3.append(float(line.strip().split()[2]))

with open('result/pdb_example_8_results_DRNApred.txt', 'r') as f:
    for line in f.readlines():
        if '>' in line or 'Amino' in line:
            continue
        dna_pred1.append(float(line.strip().split()[2]))
        rna_pred1.append(float(line.strip().split()[4]))
        dna_scores1.append(float(line.strip().split()[1]))
        rna_scores1.append(float(line.strip().split()[3]))

with open('result/pdb_example_8_results-nowind-148.txt', 'r') as f:
    for line in f.readlines():
        if '>' in line or 'Amino' in line:
            continue
        dna_pred2.append(float(line.strip().split()[2]))
        rna_pred2.append(float(line.strip().split()[4]))
        dna_scores2.append(float(line.strip().split()[1]))
        rna_scores2.append(float(line.strip().split()[3]))
        disorder_scores.append(float(line.strip().split()[5]))

disorder_pred = []
for s in disorder_scores:
    if s >= 0.32:
        disorder_pred.append(1.0)
    else:
        disorder_pred.append(0.0)
r1 = 0
d1 = 0
r2 = 0
d2 = 0
dis = 0
for i in range(len(label_dna)):
    if label_disorder[i] > 0:
        dis += 1
        if dna_pred1[i] > 0:
            d1 += 1
        if rna_pred1[i] > 0:
            r1 += 1
        if dna_pred2[i] > 0:
            d2 += 1
        if rna_pred2[i] > 0:
            r2 += 1
print dis
print d1,d2,r1,r2

mcc1 = roc_auc_score(label_dna, dna_scores1)
mcc3 = roc_auc_score(label_dna, dna_scores3)
mcc2 = roc_auc_score(label_dna, dna_scores2)
print 'DNA AUC:', mcc1,mcc3, mcc2

mcc1 = roc_auc_score(label_rna, rna_scores1)
mcc3 = roc_auc_score(label_rna, rna_scores3)
mcc2 = roc_auc_score(label_rna, rna_scores2)
print 'RNA AUC:', mcc1,mcc3, mcc2

mcc1 = lauc_score(label_dna, dna_scores1, 0.054)
mcc3 = lauc_score(label_dna, dna_scores3, 0.054)
mcc2 = lauc_score(label_dna, dna_scores2, 0.054)
print 'DNA LAUC:', mcc1,mcc3, mcc2

mcc1 = lauc_score(label_rna, rna_scores1, 0.045)
mcc3 = lauc_score(label_rna, rna_scores3, 0.045)
mcc2 = lauc_score(label_rna, rna_scores2, 0.045)
print 'RNA LAUC:', mcc1, mcc3,mcc2

mcc2 = roc_auc_score(label_disorder, disorder_scores)
print 'Disorder AUC:', mcc2

mcc1 = matthews_corrcoef(label_dna, dna_pred1)
mcc3 = matthews_corrcoef(label_dna, dna_pred3)
mcc2 = matthews_corrcoef(label_dna, dna_pred2)
print 'DNA MCC:', mcc1,mcc3, mcc2

mcc1 = matthews_corrcoef(label_rna, rna_pred1)
mcc3 = matthews_corrcoef(label_rna, rna_pred3)
mcc2 = matthews_corrcoef(label_rna, rna_pred2)
print 'RNA MCC:', mcc1, mcc3, mcc2

mcc1 = recall_score(label_dna, dna_pred1)
mcc3 = recall_score(label_dna, dna_pred3)
mcc2 = recall_score(label_dna, dna_pred2)
print 'DNA Recall:', mcc1, mcc3,mcc2

mcc1 = recall_score(label_rna, rna_pred1)
mcc3 = recall_score(label_rna, rna_pred3)
mcc2 = recall_score(label_rna, rna_pred2)
print 'RNA Recall:', mcc1,mcc3, mcc2

mcc1 = recall_score(label_rna, dna_pred1)
mcc3 = recall_score(label_rna, dna_pred3)
mcc2 = recall_score(label_rna, dna_pred2)
print 'DNA Ratio:', mcc1, mcc3, mcc2

mcc1 = recall_score(label_dna, rna_pred1)
mcc3 = recall_score(label_dna, rna_pred3)
mcc2 = recall_score(label_dna, rna_pred2)
print 'RNA Ratio:', mcc1,mcc3, mcc2

mcc1 = precision_score(label_dna, dna_pred1)
mcc3 = precision_score(label_dna, dna_pred3)
mcc2 = precision_score(label_dna, dna_pred2)
print 'DNA Precision:', mcc1, mcc3, mcc2

mcc1 = precision_score(label_rna, rna_pred1)
mcc3 = precision_score(label_rna, rna_pred3)
mcc2 = precision_score(label_rna, rna_pred2)
print 'RNA Precision:', mcc1, mcc3, mcc2

mcc = matthews_corrcoef(label_disorder, disorder_pred)
rec = recall_score(label_disorder, disorder_pred)
pre = precision_score(label_disorder, disorder_pred)
print 'Disorder:', mcc, rec, pre
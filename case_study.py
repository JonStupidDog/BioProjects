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

label_dna = []
label_rna = []
label_disorder = []
with open('examples/6joo_label.txt', 'r') as f:
    for line in f:
        label_dna.append(float(line.strip().split()[0]))
        label_rna.append(float(line.strip().split()[1]))
        label_disorder.append(float(line.strip().split()[2]))
rna_pred1 = []
dna_pred1 = []
rna_pred2 = []
dna_pred2 = []
disorder_pred = []
with open('result/6JOO_results_DRNApred.txt','r') as f:
    for line in f.readlines()[2:]:
        dna_pred1.append(float(line.strip().split()[2]))
        rna_pred1.append(float(line.strip().split()[4]))

with open('result/6JOOA_results.txt','r') as f:
    for line in f.readlines()[2:929]:
        dna_pred2.append(float(line.strip().split()[2]))
        rna_pred2.append(float(line.strip().split()[4]))

mcc1 = matthews_corrcoef(label_dna, dna_pred1)
mcc2 = matthews_corrcoef(label_dna, dna_pred2)
print mcc1, mcc2

mcc1 = matthews_corrcoef(label_rna, rna_pred1)
mcc2 = matthews_corrcoef(label_rna, rna_pred2)
print mcc1, mcc2
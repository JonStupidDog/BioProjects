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
import datetime
from time import time
import matplotlib.pyplot as plt
import numpy.random as r

def random_data():
    size = 100000
    y_true = np.array([1 if i <= 0.2974 else 0 for i in r.random(size)], dtype=np.float32)
    y_pred = r.random(size)
    fpr, tpr, th = roc_curve(y_true, y_pred)
    precision, recall, th = precision_recall_curve(y_true, y_pred)
    return precision, recall, fpr, tpr

def lauc_score(true_labels, pred_scores):
    fpr, tpr, tresholds = roc_curve(true_labels, pred_scores, sample_weight=None)
    fpr_c = 0.054
    lowc = 0
    for i in range(len(fpr)):
        if fpr[i] >= fpr_c:
            lowc = i
            break
    lauc = auc(fpr[:lowc], tpr[:lowc], reorder=True)*2/(fpr_c*fpr_c)
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
    print len(tpr),len(cro)
    return tpr, cro, thr


test_file_path = 'drnapred_data/TEST.fasta_seq.txt'
# test_dir = './data/2017_'+tp+'_test_data_seq'
# test_file_path = './data/test_all_seq.txt'
# test_file_path = './data/New_data_seq.txt'
# test_file_path = './data/YK16_test_seq3.txt'
test_prots = list(SeqIO.parse(test_file_path, 'fasta'))

true_labels1 = []
true_labels2 = []
true_labels3 = []
pred_scores11 = []
pred_scores12 = []
pred_scores21 = []
pred_scores22 = []
pred_scores23 = []
spot_scores = []
with open('result/netsurfp_results.txt', 'r') as f:
    netsurfp_scores = [float(e.strip()) for e in f]
names = {}
with open('result/spotd_results/list.desc', 'r') as f:
    for line in f:
        tmp = line.strip().split()
        names[tmp[1]] = tmp[0]
for prot in test_prots:
    with open('./result/test_label/'+prot.id+'_d.txt') as f:
        labels1 = [float(label.strip()) for label in f]
    with open('./result/test_label/'+prot.id+'_r.txt') as f:
        labels2 = [float(label.strip()) for label in f]
    with open('./result/YK17/' + prot.id + '.txt') as f:
        scores1 = []
        scores2 = []
        for line in f:
            scores1.append(float(line.strip().split()[0]))
            scores2.append(float(line.strip().split()[1]))

    with open('./result/spotd_results/' + names[prot.id] + '.spotd') as f:
        lines = f.readlines()
        for line in lines[1:]:
            spot_scores.append(float(line.strip().split()[2]))

    with open('./result/YK17_test_with_disorder/' + prot.id + '.txt') as f:
        lines = f.readlines()
        for i in range(len(prot.seq)):
            pred_scores23.append(float(lines[i].strip().split()[8]))
            if labels1[i] < 0 or labels2[i] < 0:
                pass
                # pred_scores21.append(float(lines[i].strip().split()[4]))
                # pred_scores22.append(float(lines[i].strip().split()[6]))
            else:
                # pass
                pred_scores21.append(float(lines[i].strip().split()[4]))
                pred_scores22.append(float(lines[i].strip().split()[6]))

    for i in range(len(prot.seq)):
        if labels1[i] < 0 or labels2[i] < 0:
            true_labels3.append(1.0)
            # pass
            # if labels1[i] < 0:
            #     true_labels1.append(0.0)
            # else:
            #     true_labels1.append(labels1[i])
            # if labels2[i] < 0:
            #     true_labels2.append(0.0)
            # else:
            #     true_labels2.append(labels2[i])
            # pred_scores11.append(scores1[i])
            # pred_scores12.append(scores2[i])
        else:
            true_labels3.append(0.0)
            # pass
            true_labels1.append(labels1[i])
            true_labels2.append(labels2[i])
            pred_scores11.append(scores1[i])
            pred_scores12.append(scores2[i])

# aucv = roc_auc_score(true_labels, pred_scores)
# acc = accuracy_score(true_labels, pred_labels)
# rec = recall_score(true_labels, pred_labels)
# cro = recall_score(true_labels, pred_labels)
# pre = precision_score(true_labels, pred_labels)
# mcc = matthews_corrcoef(true_labels, pred_labels)

# lowc = 0
# for i in range(len(fpr)):
#     if fpr[i] >= fpr_c:
#         lowc = i
#         break
# lauc = auc(fpr[:lowc], tpr[:lowc], reorder=True)*2/(fpr_c*fpr_c)
#
# print 'ACC:%.4f\tAUC:%.4f\tREC:%.4f\tPRE:%.4f\tMCC:%.4f\tAULC:%.6f\tRc:%.4f' % (acc, aucv, rec, pre, mcc, lauc, cro)

plt.figure(1)
plt.title('P-R Curve')  # give plot a title
plt.xlabel('True Positive Rate')  # make axis labels False Positive Rate
plt.ylabel('Cross-predicted Rate')  # True Positive Rate

auc1 = roc_auc_score(true_labels1, pred_scores11)
auc2 = roc_auc_score(true_labels1, pred_scores21)
print 'DNA AUC:', auc1, auc2
auc1 = roc_auc_score(true_labels2, pred_scores12)
auc2 = roc_auc_score(true_labels2, pred_scores22)
print 'RNA AUC:', auc1, auc2
lauc1 = lauc_score(true_labels1, pred_scores11)
lauc2 = lauc_score(true_labels1, pred_scores21)
print 'DNA LAUC:', lauc1, lauc2
lauc1 = lauc_score(true_labels2, pred_scores12)
lauc2 = lauc_score(true_labels2, pred_scores22)
print 'RNA LAUC:', lauc1, lauc2

auc1 = roc_auc_score(true_labels3, spot_scores)
auc2 = roc_auc_score(true_labels3, netsurfp_scores)
auc3 = roc_auc_score(true_labels3, pred_scores23)
print 'Disorder AUC:', auc1, auc2, auc3

# fpr1, tpr1, tresholds1 = roc_curve(true_labels2, pred_scores12)
# fpr2, tpr2, tresholds2 = roc_curve(true_labels2, pred_scores22)

# fpr1, tpr1, tresholds1 = roc_curve(true_labels3, spot_scores)
# fpr2, tpr2, tresholds2 = roc_curve(true_labels3, netsurfp_scores)
# fpr3, tpr3, tresholds3 = roc_curve(true_labels3, pred_scores23)
#
# fpr1, tpr1, tresholds1 = precision_recall_curve(true_labels3, spot_scores)
# fpr2, tpr2, tresholds2 = precision_recall_curve(true_labels3, netsurfp_scores)
# fpr3, tpr3, tresholds3 = precision_recall_curve(true_labels3, pred_scores23)

# tpr1, cro1, tresholds1 = ratio_recall_curve(true_labels2, true_labels1, pred_scores12)
# tpr2, cro2, tresholds2 = ratio_recall_curve(true_labels2, true_labels1, pred_scores22)

# pre1, rec1, tresholds1 = precision_recall_curve(true_labels1, pred_scores11)
# pre2, rec2, tresholds2 = precision_recall_curve(true_labels1, pred_scores21)

# precision_rand, recall_rand, fpr, tpr = random_data()




# plt.plot(fpr1, tpr1, color='#0066CC', label='SPOT-disorder')
# plt.plot(fpr2, tpr2, color='blue', label='NetSurfP-2.0')
# plt.plot(fpr3, tpr3, color='red', label='DRBRPred')

# plt.plot(tpr1, cro1, color='#0066CC', label='DRNAPred')
# plt.plot(tpr2, cro2, color='red', label='DRBRPred')

# plt.plot(pre1, rec1, color='#0066CC', label='DRNAPred')
# plt.plot(pre2, rec2, color='red', label='DRBRPred')

# plt.plot(fpr, tpr, color='gray', label='Random', linestyle="--")
# plt.plot(precision4, recall4, color='red', label='iDRBP_MMC')
# plt.legend()
# plt.grid()
# plt.show()

# -*- coding: utf-8 -*-
"""
Created on 2019-7-16

@author: Jun Zhang
"""
import os
import argparse
import numpy as np
from Bio import SeqIO
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score, precision_score, f1_score
from sklearn.metrics import matthews_corrcoef, auc, roc_curve, classification_report, precision_recall_curve
from sklearn.metrics import coverage_error
from sklearn.metrics import label_ranking_loss
from sklearn.metrics import hamming_loss
from sklearn.metrics import f1_score
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
    print len(tpr), len(cro)
    return tpr, cro, thr

def predict_label(true_labels, pred_scores):
    threshold1 = 0.0
    threshold2 = 0.0
    for c in np.arange(0.5, 0.0, -0.005):
        pred_labels = []
        for i in range(len(pred_scores)):
            if pred_scores[i][0] >= c:
                pred_labels.append(1.0)
            else:
                pred_labels.append(0.0)

        if sum(pred_labels) >= sum(true_labels[:, 0]):
            threshold1 = c
            break
    for c in np.arange(0.5, 0.0, -0.005):
        pred_labels = []
        for i in range(len(pred_scores)):
            if pred_scores[i][1] >= c:
                pred_labels.append(1.0)
            else:
                pred_labels.append(0.0)

        if sum(pred_labels) >= sum(true_labels[:, 1]):
            threshold2 = c
            break

    pred_labels = []
    for i in range(len(pred_scores)):
        if pred_scores[i][0] >= threshold1 and pred_scores[i][1] >= threshold2:
            pred_labels.append([1, 1])
            # print i
        elif pred_scores[i][0] >= threshold1:
            pred_labels.append([1, 0])
        elif pred_scores[i][1] >= threshold2:
            pred_labels.append([0, 1])
        else:
            pred_labels.append([0, 0])
    print 'Threshold: %.4f, %.4f' % (threshold1, threshold2)
    return pred_labels


test_file_path = 'drnapred_data/TEST.fasta_seq.txt'
# test_dir = './data/2017_'+tp+'_test_data_seq'
# test_file_path = './data/test_all_seq.txt'
# test_file_path = './data/New_data_seq.txt'
# test_file_path = './data/YK16_test_seq3.txt'
test_prots = list(SeqIO.parse(test_file_path, 'fasta'))

true_labels = []
true_labels3 = []
pred_scores1 = []
pred_scores2 = []
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
            pred_scores2.append([float(lines[i].strip().split()[4]), float(lines[i].strip().split()[6])])

    for i in range(len(prot.seq)):
        pred_scores1.append([scores1[i], scores2[i]])
        if labels1[i] < 0 or labels2[i] < 0:
            true_labels3.append(1.0)
            if labels1[i] < 0 and labels2[i] < 0:
                true_labels.append([0.0, 0.0])
            elif labels1[i] < 0:
                true_labels.append([0.0, labels2[i]])
            elif labels2[i] < 0:
                true_labels.append([labels1[i], 0.0])
            else:
                true_labels.append([labels1[i], labels2[i]])

        else:
            true_labels3.append(0.0)
            true_labels.append([labels1[i], labels2[i]])



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
true_labels = np.asarray(true_labels)
pred_scores1 = np.asarray(pred_scores1)
pred_scores2 = np.asarray(pred_scores2)
print true_labels.shape
print pred_scores2.shape
print pred_scores1.shape

auc2 = roc_auc_score(true_labels, pred_scores2, average='micro')
auc1 = roc_auc_score(true_labels, pred_scores1, average='micro')
print 'Micro AUC:', auc1, auc2
auc1 = roc_auc_score(true_labels, pred_scores1, average='macro')
auc2 = roc_auc_score(true_labels, pred_scores2, average='macro')
print 'Macro AUC:', auc1, auc2
pred_labels1 = predict_label(true_labels, pred_scores1)
pred_labels2 = predict_label(true_labels, pred_scores2)
# loss1 = hamming_loss(true_labels, pred_labels1)
# loss2 = hamming_loss(true_labels, pred_labels2)
# print 'Hamming loss:', loss1, loss2
f1_1 = f1_score(true_labels, pred_labels1, average='macro')
f1_2 = f1_score(true_labels, pred_labels2, average='macro')
print 'Macro F1-score:', f1_1, f1_2
f1_1 = f1_score(true_labels, pred_labels1, average='micro')
f1_2 = f1_score(true_labels, pred_labels2, average='micro')
print 'Micro F1-score:', f1_1, f1_2


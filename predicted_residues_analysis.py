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
    return np.asarray(pred_labels)

def distance_labels(labels, distance):
    temp_labels = []
    true_labels = []
    for e in labels:
        if e >= 0:
            true_labels.append(e)
        else:
            true_labels.append(0)
    for i in range(len(true_labels)):
        if i < distance:
            if sum(true_labels[:i+1+distance]) > 0:
                temp_labels.append(1)
            else:
                temp_labels.append(0)
        elif i >= len(true_labels)-distance-1:
            if sum(true_labels[i-distance:]) > 0:
                temp_labels.append(1)
            else:
                temp_labels.append(0)
        else:
            if sum(true_labels[i-distance:i+distance]) > 0:
                temp_labels.append(1)
            else:
                temp_labels.append(0)
    return temp_labels

def get_perfomance(labels, scores, distance, threshold):
    dTP = 0.0
    TP = 0.0
    TN = 0.0
    FP = 0.0
    FN = 0.0
    for true_labels, pred_scores in zip(labels, scores):
        if len(true_labels) != len(pred_scores):
            print 'ERROR!'
        else:
            for i in range(len(true_labels)):
                if true_labels[i] > 0 and pred_scores[i] >= threshold:
                    TP += 1.0
                if true_labels[i] > 0 and pred_scores[i] < threshold:
                    FN += 1.0
                if true_labels[i] == 0 and pred_scores[i] >= threshold:
                    FP += 1.0
                if true_labels[i] == 0 and pred_scores[i] < threshold:
                    TN += 1.0

                if i < distance:
                    if sum(true_labels[:i+1+distance]) > 0:
                        if pred_scores[i] >= threshold:
                            dTP += 1.
                elif i >= len(true_labels)-distance-1:
                    if sum(true_labels[i-distance:]) > 0:
                        if pred_scores[i] >= threshold:
                            dTP += 1.
                else:
                    if sum(true_labels[i-distance:i+distance+1]) > 0:
                        if pred_scores[i] >= threshold:
                            dTP += 1.

    ca = dTP - TP
    TP = dTP
    FP = FP - ca
    print TP, TN, FP, FN
    SN = TP / (TP + FN)  # Sensitivity = TP/P  and P = TP + FN
    SP = TN / (FP + TN)  # Specificity = TN/N  and N = TN + FP
    MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    ACC = (TP+TN) / len(labels)
    return SN, SP, MCC, ACC

def get_labels_scores(test_prots):
    true_labels1 = []
    true_labels2 = []
    pred_scores11 = []
    pred_scores12 = []
    pred_scores21 = []
    pred_scores22 = []
    for prot in test_prots:
        with open('./result/test_label/' + prot.id + '_d.txt') as f:
            labels1 = []
            for label in f:
                if float(label.strip()) < 0:
                    labels1.append(0.0)
                else:
                    labels1.append(float(label.strip()))
            true_labels1.append(labels1)
        with open('./result/test_label/' + prot.id + '_r.txt') as f:
            labels2 = []
            for label in f:
                if float(label.strip()) < 0:
                    labels2.append(0.0)
                else:
                    labels2.append(float(label.strip()))
            true_labels2.append(labels2)

        with open('./result/YK17/' + prot.id + '.txt') as f:
            scores1 = []
            scores2 = []
            for line in f:
                scores1.append(float(line.strip().split()[0]))
                scores2.append(float(line.strip().split()[1]))
            pred_scores11.append(scores1)
            pred_scores12.append(scores2)
        with open('./result/YK17_test_with_disorder/' + prot.id + '.txt') as f:
            lines = f.readlines()
            scores1 = []
            scores2 = []
            for i in range(len(prot.seq)):
                scores1.append(float(lines[i].strip().split()[4]))
                scores2.append(float(lines[i].strip().split()[6]))
            pred_scores21.append(scores1)
            pred_scores22.append(scores2)
    return true_labels1, true_labels2, pred_scores11, pred_scores12, pred_scores21, pred_scores22

test_file_path = 'drnapred_data/TEST.fasta_seq.txt'
# test_dir = './data/2017_'+tp+'_test_data_seq'
# test_file_path = './data/test_all_seq.txt'
# test_file_path = './data/New_data_seq.txt'
# test_file_path = './data/YK16_test_seq3.txt'
test_prots = list(SeqIO.parse(test_file_path, 'fasta'))

true_labels1, true_labels2, pred_scores11, pred_scores12, pred_scores21, pred_scores22 = get_labels_scores(test_prots)

# Threshold: 0.5000, 0.1500
# Threshold: 0.4950, 0.1100
print 'DNA'
SN, SP, MCC, ACC = get_perfomance(true_labels1, pred_scores11, 4, 0.50)
print 'SN, MCC:', SN, MCC
SN, SP, MCC, ACC = get_perfomance(true_labels1, pred_scores21, 4, 0.495)
print 'SN, MCC:', SN, MCC
print 'RNA'
SN, SP, MCC, ACC = get_perfomance(true_labels2, pred_scores12, 4, 0.15)
print 'SN, MCC:', SN, MCC
SN, SP, MCC, ACC = get_perfomance(true_labels2, pred_scores22, 4, 0.11)
print 'SN, MCC:', SN, MCC

# true_labels = np.asarray(true_labels)
# pred_scores1 = np.asarray(pred_scores1)
# pred_scores2 = np.asarray(pred_scores2)
#
# pred_labels1 = predict_label(true_labels, pred_scores1)
# pred_labels2 = predict_label(true_labels, pred_scores2)
#
# true_labels = []
# for prot in test_prots:
#     with open('./result/test_label/' + prot.id + '_d.txt') as f:
#         labels1 = [float(label.strip()) for label in f]
#         labels1 = distance_labels(labels1, 1)
#     with open('./result/test_label/' + prot.id + '_r.txt') as f:
#         labels2 = [float(label.strip()) for label in f]
#         labels2 = distance_labels(labels2, 1)
#     for i in range(len(prot.seq)):
#         true_labels.append([labels1[i], labels2[i]])
# true_labels = np.asarray(true_labels)
# # print true_labels[:, 0]
# # print true_labels[:, 1]
# #
# # exit()
#
# f1_1 = matthews_corrcoef(true_labels[:, 0], pred_labels1[:, 0])
# f1_2 = matthews_corrcoef(true_labels[:, 0], pred_labels2[:, 0])
# print 'DNA MCC:', f1_1, f1_2
# f1_1 = matthews_corrcoef(true_labels[:, 1], pred_labels1[:, 1])
# f1_2 = matthews_corrcoef(true_labels[:, 1], pred_labels2[:, 1])
# print 'RNA MCC:', f1_1, f1_2
#
# f1_1 = recall_score(true_labels[:, 0], pred_labels1[:, 0])
# f1_2 = recall_score(true_labels[:, 0], pred_labels2[:, 0])
# print 'DNA Recall:', f1_1, f1_2
# f1_1 = recall_score(true_labels[:, 1], pred_labels1[:, 1])
# f1_2 = recall_score(true_labels[:, 1], pred_labels2[:, 1])
# print 'RNA Recall:', f1_1, f1_2
#

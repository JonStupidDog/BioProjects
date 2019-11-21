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

def lauc_score(true_labels, pred_scores, fpr_c):
    fpr, tpr, tresholds = roc_curve(true_labels, pred_scores, sample_weight=None)
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

def evaluate(true_label, pred_label):
    tp = 0
    fp = 0
    tn = 0
    fn = 0

    for tl, pl in zip(true_label, pred_label):
        if tl > 0:
            if pl > 0:
                tp += 1
            else:
                fn += 1
        else:
            if pl > 0:
                fp += 1
            else:
                tn += 1
    return float(tp),float(fp),float(tn),float(fn)


test_file_path = 'data/Disorder_binding_Test_36_seq.txt'
# test_dir = './data/2017_'+tp+'_test_data_seq'
# test_file_path = './data/test_all_seq.txt'
# test_file_path = './data/New_data_seq.txt'
# test_file_path = './data/YK16_test_seq3.txt'
test_prots = list(SeqIO.parse(test_file_path, 'fasta'))

true_labels1 = []
true_labels2 = []

pred_scores11 = []
pred_scores12 = []
pred_scores21 = []
pred_scores22 = []
pred_scores31 = []
pred_scores32 = []
pred_scores41 = []
pred_scores42 = []

pred_labels11 = []
pred_labels12 = []
pred_labels21 = []
pred_labels22 = []
pred_labels31 = []
pred_labels32 = []
pred_labels41 = []
pred_labels42 = []


with open('result/Disorder_binding_Test_36_results.txt') as f:
    lines = f.readlines()
    for line in lines:
        if '>' in line or 'Amino' in line:
            pass
        else:
            pred_scores21.append(float(line.strip().split()[1]))
            pred_scores22.append(float(line.strip().split()[3]))
            pred_labels21.append(float(line.strip().split()[2]))
            pred_labels22.append(float(line.strip().split()[4]))

with open('result/DRNApred_disorder_test_36_results.txt') as f:
    lines = f.readlines()
    for line in lines:
        if '>' in line or 'Amino' in line:
            pass
        else:
            pred_scores11.append(float(line.strip().split()[1]))
            pred_scores12.append(float(line.strip().split()[3]))
            pred_labels11.append(float(line.strip().split()[2]))
            pred_labels12.append(float(line.strip().split()[4]))

with open('result/DisoRDPbind_36_results.txt') as f:
    lines = f.readlines()
    for i in range(len(lines)):
        if '>' in lines[i]:
            RNA_pred = lines[i + 2].strip().split(':')[1]
            RNA_score = lines[i + 3].strip().split(':')[1].strip(',').split(',')
            DNA_pred = lines[i + 4].strip().split(':')[1]
            DNA_score = lines[i + 5].strip().split(':')[1].strip(',').split(',')
            for j in range(len(RNA_score)):
                pred_scores41.append(float(DNA_score[j]))
                pred_scores42.append(float(RNA_score[j]))
                pred_labels41.append(float(DNA_pred[j]))
                pred_labels42.append(float(RNA_pred[j]))



for prot in test_prots:
    with open('./data/DisDRBP_test_36_label/'+prot.id+'_d.txt') as f:
        labels1 = [float(label.strip()) for label in f]
    with open('./data/DisDRBP_test_36_label/'+prot.id+'_r.txt') as f:
        labels2 = [float(label.strip()) for label in f]

    with open('result/NucBind_disorder_test_36/' + prot.id + '_DNA_svm.txt') as f:
        scores31 = []
        labels31 = []
        for line in f.readlines()[1:]:
            scores31.append(float(line.strip().split()[3]))
            labels31.append(float(line.strip().split()[2]))
    with open('result/NucBind_disorder_test_36/' + prot.id + '_RNA_svm.txt') as f:
        scores32 = []
        labels32 = []
        for line in f.readlines()[1:]:
            scores32.append(float(line.strip().split()[3]))
            labels32.append(float(line.strip().split()[2]))

    true_labels1.extend(labels1)
    true_labels2.extend(labels2)
    pred_scores31.extend(scores31)
    pred_scores32.extend(scores32)
    pred_labels31.extend(labels31)
    pred_labels32.extend(labels32)


true1_labels1 = []
true1_labels2 = []

pred1_scores11 = []
pred1_scores12 = []
pred1_scores21 = []
pred1_scores22 = []
pred1_scores31 = []
pred1_scores32 = []
pred1_scores41 = []
pred1_scores42 = []

pred1_labels11 = []
pred1_labels12 = []
pred1_labels21 = []
pred1_labels22 = []
pred1_labels31 = []
pred1_labels32 = []
pred1_labels41 = []
pred1_labels42 = []

n1=0
n2=0
for i in range(len(true_labels1)):
    if true_labels1[i] < 0 or true_labels2[i] < 0:
        if true_labels1[i] < 0 and true_labels2[i] < 0:
            pass
        else:
            if true_labels1[i] < 0:
                n2+=1
            else:
                n1+=1
    else:
        true1_labels1.append(true_labels1[i])
        true1_labels2.append(true_labels2[i])

        pred1_scores11.append(pred_scores11[i])
        pred1_scores12.append(pred_scores12[i])
        pred1_scores21.append(pred_scores21[i])
        pred1_scores22.append(pred_scores22[i])
        pred1_scores31.append(pred_scores31[i])
        pred1_scores32.append(pred_scores32[i])
        pred1_scores41.append(pred_scores41[i])
        pred1_scores42.append(pred_scores42[i])

        pred1_labels11.append(pred_labels11[i])
        pred1_labels12.append(pred_labels12[i])
        pred1_labels21.append(pred_labels21[i])
        pred1_labels22.append(pred_labels22[i])
        pred1_labels31.append(pred_labels31[i])
        pred1_labels32.append(pred_labels32[i])
        pred1_labels41.append(pred_labels41[i])
        pred1_labels42.append(pred_labels42[i])
print n1, n2
print sum(true1_labels1), sum(true1_labels2)
tp1,fp1,tn1,fn1 = evaluate(true_labels1, pred1_labels11)
tp2,fp2,tn2,fn2 = evaluate(true_labels1, pred1_labels21)
tp3,fp3,tn3,fn3 = evaluate(true_labels1, pred1_labels31)
tp4,fp4,tn4,fn4 = evaluate(true_labels1, pred1_labels41)
print 'tp1,fp1,tn1,fn1:',tp1,fp1,tn1,fn1,fp1/(fp1+tn1)
print 'tp2,fp2,tn2,fn2:',tp2,fp2,tn2,fn2,fp2/(fp2+tn2)
print 'tp3,fp3,tn3,fn3:',tp3,fp3,tn3,fn3,fp3/(fp3+tn3)
print 'tp4,fp4,tn4,fn4:',tp4,fp4,tn4,fn4,fp4/(fp3+tn4)
mcc1 = matthews_corrcoef(true1_labels1, pred1_labels11)
mcc2 = matthews_corrcoef(true1_labels1, pred1_labels21)
mcc3 = matthews_corrcoef(true1_labels1, pred1_labels31)
mcc4 = matthews_corrcoef(true1_labels1, pred1_labels41)
print 'DNA MCC:', mcc1, mcc2, mcc3, mcc4
tp1,fp1,tn1,fn1 = evaluate(true1_labels2, pred1_labels12)
tp2,fp2,tn2,fn2 = evaluate(true1_labels2, pred1_labels22)
tp3,fp3,tn3,fn3 = evaluate(true1_labels2, pred1_labels32)
tp4,fp4,tn4,fn4 = evaluate(true1_labels2, pred1_labels42)
print 'tp1,fp1,tn1,fn1:',tp1,fp1,tn1,fn1,fp1/(fp1+tn1)
print 'tp2,fp2,tn2,fn2:',tp2,fp2,tn2,fn2,fp2/(fp2+tn2)
print 'tp3,fp3,tn3,fn3:',tp3,fp3,tn3,fn3,fp3/(fp3+tn3)
print 'tp4,fp4,tn4,fn4:',tp4,fp4,tn4,fn4,fp4/(fp4+tn4)
mcc1 = matthews_corrcoef(true1_labels2, pred1_labels12)
mcc2 = matthews_corrcoef(true1_labels2, pred1_labels22)
mcc3 = matthews_corrcoef(true1_labels2, pred1_labels32)
mcc4 = matthews_corrcoef(true1_labels2, pred1_labels42)
print 'RNA MCC:', mcc1, mcc2, mcc3, mcc4

# aucv = roc_auc_score(true_labels, pred_scores)
# acc = accuracy_score(true_labels, pred_labels)
# rec = recall_score(true_labels, pred_labels)
# cro = recall_score(true_labels, pred_labels)
# pre = precision_score(true_labels, pred_labels)
# mcc = matthews_corrcoef(true_labels, pred_labels)
# print mcc
# lowc = 0
# for i in range(len(fpr)):
#     if fpr[i] >= fpr_c:
#         lowc = i
#         break
# lauc = auc(fpr[:lowc], tpr[:lowc], reorder=True)*2/(fpr_c*fpr_c)
#
# print 'ACC:%.4f\tAUC:%.4f\tREC:%.4f\tPRE:%.4f\tMCC:%.4f\tAULC:%.6f\tRc:%.4f' % (acc, aucv, rec, pre, mcc, lauc, cro)

plt.figure(1)
plt.title('ROC Curve')  # give plot a title
plt.xlabel('False Positive Rate')  # make axis labels False Positive Rate
plt.ylabel('True Positive Rate')  # True Positive Rate

auc1 = roc_auc_score(true1_labels1, pred1_scores11)
auc2 = roc_auc_score(true1_labels1, pred1_scores21)
auc3 = roc_auc_score(true1_labels1, pred1_scores31)
auc4 = roc_auc_score(true1_labels1, pred1_scores41)
print 'DNA AUC:', auc1, auc2, auc3, auc4
auc1 = roc_auc_score(true1_labels2, pred1_scores12)
auc2 = roc_auc_score(true1_labels2, pred1_scores22)
auc3 = roc_auc_score(true1_labels2, pred1_scores32)
auc4 = roc_auc_score(true1_labels2, pred1_scores42)
print 'RNA AUC:', auc1, auc2, auc3, auc4

# lauc1 = lauc_score(true_labels1, pred_scores11, 0.054)
# lauc2 = lauc_score(true_labels1, pred_scores21, 0.054)
# lauc3 = lauc_score(true_labels1, pred_scores31, 0.054)
# lauc4 = lauc_score(true_labels1, pred_scores41, 0.054)
# print 'DNA LAUC:', lauc1, lauc2, lauc3, lauc4
# lauc1 = lauc_score(true_labels2, pred_scores12, 0.045)
# lauc2 = lauc_score(true_labels2, pred_scores22, 0.045)
# lauc3 = lauc_score(true_labels2, pred_scores32, 0.045)
# lauc4 = lauc_score(true_labels2, pred_scores42, 0.045)
# print 'RNA LAUC:', lauc1, lauc2, lauc3, lauc4

# mcc1 = recall_score(true_labels2, dna_pred1)
# mcc2 = recall_score(true_labels2, dna_pred2)
# print 'DNA Ratio:', mcc1, mcc2
#
# mcc1 = recall_score(true_labels1, rna_pred1)
# mcc2 = recall_score(true_labels1, rna_pred2)
# print 'RNA Ratio:', mcc1, mcc2

# auc1 = roc_auc_score(true_labels3, netsurfp_scores)
# auc4 = roc_auc_score(true_labels3, espritz_scores)
# auc2 = roc_auc_score(true_labels3, spot_scores)
# auc3 = roc_auc_score(true_labels3, pred_scores23)
# print 'Disorder AUC:', auc1, auc2, auc4, auc3


# tpr1, cro1, t1 = ratio_recall_curve(true_labels1, true_labels2, pred_scores11)
# tpr2, cro2, t2 = ratio_recall_curve(true_labels1, true_labels2, pred_scores21)
# tpr3, cro3, t3 = ratio_recall_curve(true_labels1, true_labels2, pred_scores31)
# tpr4, cro4, t4 = ratio_recall_curve(true_labels1, true_labels2, pred_scores41)
# aurc1 = auc(tpr1, cro1)
# aurc2 = auc(tpr2, cro2)
# aurc3 = auc(tpr3, cro3)
# aurc4 = auc(tpr4, cro4)
# print 'DNA AURC:', aurc1, aurc2, aurc3, aurc4

# fpr_c = 0.5
# lowc = 0
# for i in range(len(tpr1)):
#     if tpr1[i] >= fpr_c:
#         lowc = i
#         break
# lauc1 = auc(tpr1[:lowc], cro1[:lowc])
# lowc = 0
# for i in range(len(tpr2)):
#     if tpr2[i] >= fpr_c:
#         lowc = i
#         break
# lauc2 = auc(tpr2[:lowc], cro2[:lowc])
# lowc = 0
# for i in range(len(tpr3)):
#     if tpr3[i] >= fpr_c:
#         lowc = i
#         break
# lauc3 = auc(tpr3[:lowc], cro3[:lowc])
# lowc = 0
# for i in range(len(tpr4)):
#     if tpr4[i] >= fpr_c:
#         lowc = i
#         break
# lauc4 = auc(tpr4[:lowc], cro4[:lowc])
# print 'DNA AURLC:', lauc1, lauc2, lauc3, lauc4

# tpr1, cro1, t1 = ratio_recall_curve(true_labels2, true_labels1, pred_scores12)
# tpr2, cro2, t2 = ratio_recall_curve(true_labels2, true_labels1, pred_scores22)
# tpr3, cro3, t3 = ratio_recall_curve(true_labels2, true_labels1, pred_scores32)
# tpr4, cro4, t4 = ratio_recall_curve(true_labels2, true_labels1, pred_scores42)
# aurc1 = auc(tpr1, cro1)
# aurc2 = auc(tpr2, cro2)
# aurc3 = auc(tpr3, cro3)
# aurc4 = auc(tpr4, cro4)
# print 'RNA AURC:', aurc1, aurc2, aurc3, aurc4

# fpr_c = 0.5
# lowc = 0
# for i in range(len(tpr1)):
#     if tpr1[i] >= fpr_c:
#         lowc = i
#         break
# lauc1 = auc(tpr1[:lowc], cro1[:lowc])
# lowc = 0
# for i in range(len(tpr2)):
#     if tpr2[i] >= fpr_c:
#         lowc = i
#         break
# lauc2 = auc(tpr2[:lowc], cro2[:lowc])
# lowc = 0
# for i in range(len(tpr3)):
#     if tpr3[i] >= fpr_c:
#         lowc = i
#         break
# lauc3 = auc(tpr3[:lowc], cro3[:lowc])
# lowc = 0
# for i in range(len(tpr4)):
#     if tpr4[i] >= fpr_c:
#         lowc = i
#         break
# lauc4 = auc(tpr4[:lowc], cro4[:lowc])
# print 'RNA AURLC:', lauc1, lauc2, lauc3, lauc4

fpr1, tpr1, tresholds1 = roc_curve(true1_labels2, pred1_scores12)
fpr2, tpr2, tresholds2 = roc_curve(true1_labels2, pred1_scores22)
fpr3, tpr3, tresholds3 = roc_curve(true1_labels2, pred1_scores32)
fpr4, tpr4, tresholds4 = roc_curve(true1_labels2, pred1_scores42)

# fpr1, tpr1, tresholds1 = roc_curve(true1_labels1, pred1_scores11)
# fpr2, tpr2, tresholds2 = roc_curve(true1_labels1, pred1_scores21)
# fpr3, tpr3, tresholds3 = roc_curve(true1_labels1, pred1_scores31)
# fpr4, tpr4, tresholds4 = roc_curve(true1_labels1, pred1_scores41)
#
# fpr1, tpr1, tresholds1 = roc_curve(true_labels3, netsurfp_scores)
# fpr2, tpr2, tresholds2 = roc_curve(true_labels3, spot_scores)
# fpr4, tpr4, tresholds4 = roc_curve(true_labels3, espritz_scores)
# fpr3, tpr3, tresholds3 = roc_curve(true_labels3, pred_scores23)
# #
# fpr1, tpr1, tresholds1 = precision_recall_curve(true_labels3, spot_scores)
# fpr2, tpr2, tresholds2 = precision_recall_curve(true_labels3, netsurfp_scores)
# fpr3, tpr3, tresholds3 = precision_recall_curve(true_labels3, pred_scores23)
#
# pre1, rec1, tresholds1 = precision_recall_curve(true_labels2, pred_scores12)
# pre2, rec2, tresholds2 = precision_recall_curve(true_labels2, pred_scores22)
# pre3, rec3, tresholds3 = precision_recall_curve(true_labels2, pred_scores32)

precision_rand, recall_rand, fpr, tpr = random_data()

plt.plot(fpr4, tpr4, color='#8968CD', label='DisoRDPbind')
plt.plot(fpr1, tpr1, color='#424242', label='DRNAPred')
plt.plot(fpr3, tpr3, color='#0066CC', label='NucBind')
plt.plot(fpr2, tpr2, color='red', label='iDRBR_MSL')

#
# plt.plot(fpr2, tpr2, color='#8968CD', label='SPOT-disorder')
# plt.plot(fpr1, tpr1, color='#8B4726', label='NetSurfP-2.0')
# plt.plot(fpr4, tpr4, color='#CDAD00', label='Espritz')
# plt.plot(fpr3, tpr3, color='red', label='DRBRPred')

# plt.plot(tpr1, cro1, color='#424242', label='DRNAPred')
# plt.plot(tpr3, cro3, color='#0066CC', label='DRNAPred')
# plt.plot(tpr2, cro2, color='red', label='DRBRPred')

# plt.plot(pre1, rec1, color='#424242', label='DRNAPred')
# plt.plot(pre3, rec3, color='#0066CC', label='DRNAPred')
# plt.plot(pre2, rec2, color='red', label='DRBRPred')

plt.plot(fpr, tpr, color='#D3D3D3', label='Random', linestyle="--")
# plt.plot(precision4, recall4, color='red', label='iDRBP_MMC')
plt.legend()
plt.grid()
plt.show()

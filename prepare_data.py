from Bio import SeqIO
import numpy as np
import os
def get_seq():
    proteins = list(SeqIO.parse('data/2017_DNA_test_data.txt', 'fasta'))
    print proteins[0].seq

    with open('data/2017_RNA_test_data.txt','r') as f:
        lines=f.readlines()
        with open('data/2017_RNA_test_data_seq.txt','w') as wf:
            for i in range(len(lines)):
                if int(i+1) % 3 == 0:
                    continue
                else:
                    wf.write(lines[i])
def count_protein_len():
    proteins = list(SeqIO.parse('data/2017_DNA_train_data_seq.txt', 'fasta'))
    proteins.extend(list(SeqIO.parse('data/2017_RNA_train_data_seq.txt', 'fasta')))
    tmax=0
    tmin=9999
    L=[]
    for prot in proteins:
        L.append(len(prot.seq))
        print len(prot.seq)
        if len(prot.seq)>tmax:
            tmax=len(prot.seq)
        if len(prot.seq)<tmin:
            tmin=len(prot.seq)
    print tmin,tmax
    L.sort()
    print L

def get_label():
    with open('data/2017_DNA_train_data.txt', 'r') as f:
        if not os.path.isdir('data/2017_DNA_train_data_label/labels_masking'):
            os.makedirs('data/2017_DNA_train_data_label/labels_masking')
        lines = f.readlines()
        for i in range(len(lines)):
                if '>' in lines[i]:
                    with open('data/2017_DNA_train_data_label/'+lines[i].strip().split('>')[-1]+'.txt', 'w') as wf:
                        for e in lines[i+2].strip():
                            if e == '2':
                                wf.write('-1\n')
                            else:
                                wf.write(e + '\n')
                else:
                    pass
def get_bioseq_analysis_style():
    # proteins = list(SeqIO.parse('data/2017_DNA_test_data.txt', 'fasta'))
    with open('data/2017_RNA_test_data.txt', 'r') as f:
        lines = f.readlines()
        with open('data/2017_RNA_test_data_label_bioseq_analysis.txt', 'w') as wf:
            with open('data/2017_RNA_test_data_seq_bioseq_analysis.txt', 'w') as df:
                for i in range(len(lines)):
                    if int(i + 1) % 3 == 0:
                        for j in range(len(lines[i].strip())):
                            if lines[i].strip()[j] == '2':
                                pass
                            else:
                                wf.write(lines[i].strip()[j]+' ')
                                df.write(lines[i-1].strip()[j])
                        wf.write('\n')
                        df.write('\n')
                    elif '>' in lines[i]:
                        wf.write(lines[i].strip() + '\n')
                        df.write(lines[i].strip() + '\n')
                    else:
                        pass

if __name__ == '__main__':
    get_bioseq_analysis_style()
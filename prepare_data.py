from Bio import SeqIO
import numpy as np
import os
def get_seq():
    # proteins = list(SeqIO.parse('data/2017_DNA_test_data.txt', 'fasta'))
    # print proteins[0].seq
    with open('data/YK16_3.5A_DNA_training.fasta','r') as f:
        lines = f.readlines()[4:]
        with open('data/YK16_3.5A_DNA_training_seq.txt','w') as wf:
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
        with open('data/2017_RNA_test_data_label.txt', 'w') as wf:
            with open('data/2017_RNA_test_data_seq.txt', 'w') as df:
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

def get_seq_and_label():
    file_name = 'YK16_5A_DNA_training'
    with open('data/'+file_name+'.txt', 'r') as f:
        lines = f.readlines()[4:]
        with open('data/'+file_name+'_seq.txt', 'w') as wf:
            for i in range(len(lines)):
                if int(i+1) % 3 == 0:
                    continue
                else:
                    if '>' in lines[i]:
                        wf.write(lines[i].split()[0]+'\n')
                    else:
                        wf.write(lines[i].strip()+'\n')
        if not os.path.isdir('data/'+file_name+'_label'):
            os.makedirs('data/'+file_name+'_label')
        for i in range(len(lines)):
                if '>' in lines[i]:
                    with open('data/'+file_name+'_label/'+lines[i].split()[0].split('>')[-1]+'.txt', 'w') as wf:
                        for e in lines[i+2].strip():
                            if e == '2':
                                wf.write('-1\n')
                            else:
                                wf.write(e + '\n')
                else:
                    pass

def get_multi_label():
    file_name = 'New_D31'
    rfile_name = 'New_R15'
    with open('data/'+file_name+'.txt', 'r') as f:
        lines = f.readlines()[4:]
        if not os.path.isdir('data/'+file_name+'_label'):
            os.makedirs('data/'+file_name+'_label')
        for i in range(len(lines)):
                if '>' in lines[i]:
                    with open('data/'+file_name+'_label/'+lines[i].split()[0].split('>')[-1]+'_d.txt', 'w') as wf:
                        for e in lines[i+2].strip():
                            if e == '2':
                                wf.write('-1\n')
                            else:
                                wf.write(e + '\n')
                    with open('data/'+file_name+'_label/'+lines[i].split()[0].split('>')[-1]+'_r.txt', 'w') as wf:
                        if os.path.isfile('data/'+rfile_name+'_label/'+lines[i].split()[0].split('>')[-1]+'.txt'):
                            with open('data/'+rfile_name+'_label/'+lines[i].split()[0].split('>')[-1]+'.txt') as rf:
                                for e in rf:
                                    wf.write(e.strip()+'\n')
                        else:
                            for e in lines[i+2].strip():
                                if e == '2':
                                    wf.write('-1\n')
                                else:
                                    wf.write('0\n')
                else:
                    pass

def get_seq_and_label_combine():
    file_name = 'TRAINING.fasta'
    with open('drnapred_data/'+file_name+'.txt', 'r') as f:
        lines = f.readlines()
        with open('drnapred_data/'+file_name+'_seq.txt', 'w') as wf:
            for i in range(len(lines)):
                if '>' in lines[i]:
                    wf.write(lines[i].strip()+'\n')
                    wf.write(lines[i+1].strip()+'\n')
                else:
                    pass
        if not os.path.isdir('drnapred_data/'+file_name+'_label'):
            os.makedirs('drnapred_data/'+file_name+'_label')
        for i in range(len(lines)):
                if '>' in lines[i]:
                    with open('drnapred_data/'+file_name+'_label/'+lines[i].strip().split('>')[-1]+'_d.txt', 'w') as wf:
                        for e in lines[i+2].strip():
                            if e == '2':
                                wf.write('-1\n')
                            else:
                                wf.write(e + '\n')
                    with open('drnapred_data/'+file_name+'_label/'+lines[i].strip().split('>')[-1]+'_r.txt', 'w') as wf:
                        for e in lines[i+3].strip():
                            if e == '2':
                                wf.write('-1\n')
                            else:
                                wf.write(e + '\n')
                else:
                    pass
if __name__ == '__main__':
    # get_multi_label()
    prots5 = list(SeqIO.parse('data/train_all_seq.txt', 'fasta'))
    idd = []
    for prot in prots5:
        idd.append(prot.id)
    tmp = []
    prots5 = list(SeqIO.parse('data/test_all_seq.txt', 'fasta'))
    for prot in prots5:
        if prot.id not in idd:
            tmp.append(prot)
            idd.append(prot.id)

    # prots5 = list(SeqIO.parse('data/YK16_test_seq.txt', 'fasta'))
    # for prot in prots5:
    #     if prot.id not in idd:
    #         tmp.append(prot)
    #         idd.append(prot.id)
    #
    # prots5 = list(SeqIO.parse('data/New_data_seq.txt', 'fasta'))
    # for prot in prots5:
    #     if prot.id not in idd:
    #         tmp.append(prot)
    #         idd.append(prot.id)
    #
    print len(tmp)
    SeqIO.write(tmp, 'data/test_all_seq.txt', 'fasta')




    # with open('drnapred_data/lists_DNA.txt') as f:
    #     idd = [e.strip() for e in f]
    # with open('drnapred_data/lists_RNA.txt') as f:
    #     idr = [e.strip() for e in f]
    # print len(set(idd) & set(idr))
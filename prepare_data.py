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

def split_DRNApred_results(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if '>' in lines[i]:
                with open('result/YK16_3.5/' + lines[i].strip().strip('>') + '.txt', 'w') as wf:
                    n = 2
                    while n > 0:
                        wf.write(lines[i+n].strip().split()[1]+'\t'+lines[i+n].strip().split()[3]+'\n')
                        n += 1
                        if i+n >= len(lines):
                            break
                        elif '>' in lines[i+n]:
                            n = 0
                        else:
                            pass


if __name__ == '__main__':
    prots = list(SeqIO.parse('data/examples_seq_A_output.txt', 'fasta'))
    ids = []
    tmp_prots = []
    for prot in prots:
       if len(prot.seq) > 50 and 'U' not in prot.seq and 'X' not in prot.seq:
           tmp_prots.append(prot)

    # prots5 = list(SeqIO.parse('data/New_data_seq.txt', 'fasta'))
    # id1 = []
    # for prot in prots5:
    #     id1.append(prot.id)
    # id2 = []
    # prots5 = list(SeqIO.parse('data/YK16_test_seq3.txt', 'fasta'))
    # for prot in prots5:
    #     id2.append(prot.id)
    #
    # prots5 = list(SeqIO.parse('drnapred_data/TRAINING.fasta_seq.txt', 'fasta'))
    # id3 = []
    # for prot in prots5:
    #     id3.append(prot.id)
    #
    # prots5 = list(SeqIO.parse('drnapred_data/TEST.fasta_seq.txt', 'fasta'))
    # id4 = []
    # for prot in prots5:
    #     id4.append(prot.id)
    #
    # print len(tmp)
    SeqIO.write(tmp_prots, 'data/examples_seq.txt', 'fasta')
    # split_DRNApred_results('result/results.txt')
    #
    #
    #
    # with open('drnapred_data/lists_DNA.txt') as f:
    #     idd = [e.strip() for e in f]
    # with open('drnapred_data/lists_RNA.txt') as f:
    #     idr = [e.strip() for e in f]
    # print len(set(id3) & set(id1))
    # print len(set(id4) & set(id1))
    # print len(set(id2) & set(id1))
    # A = [0, 319, 320, 529, 182, 353, 538, 545, 183, 316, 319, 348, 350, 353, 448, 449, 395, 525, 529, 28,
    #      0, 319, 320, 322, 525, 529, 353, 401, 416, 474, 510, 528, 538, 545, 183, 316, 348, 449, 453, 36,
    #      0, 319, 320, 529, 182, 353, 538, 545, 183, 316, 319, 348, 350, 353, 448, 449, 395, 525, 529, 28]
    # B = [12,  13,  14, 153, 183, 254, 316, 319, 348, 350, 353, 357, 448, 449, 450, 451, 452, 453, 495, 496,
    #      0,  37,  70,  88,  90,  96, 133, 148, 183, 216, 319, 322, 358, 373, 495, 519, 520, 524, 525, 529,
    #      187, 193, 196, 201, 219, 223, 231, 251, 256, 263, 289, 293, 352, 388, 426, 477, 481, 491, 520, 537]
    # # C = set(A)&set(B)
    # # print len(C)
    # # print C
    #
    # tmp1 = []
    # for e in A:
    #     if e not in tmp1:
    #         tmp1.append(e)
    # print len(tmp1)
    # print tmp1
    #
    # tmp2 = []
    # for e in B:
    #     if e not in tmp2:
    #         tmp2.append(e)
    # print len(tmp2)
    # print tmp2
    #
    # C = set(tmp2)|set(tmp1)
    # print len(C)
    # print C

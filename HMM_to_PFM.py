import numpy as np
import os

hmm_dir='/home/zhangjun/PycharmProjects/iDRBP_Residues/data/2017_DNA_train_data_seq/hh_profile'
pfm_dir='/home/zhangjun/PycharmProjects/iDRBP_Residues/data/2017_DNA_train_data_seq/PFM'
if not os.path.isdir(pfm_dir):
    os.makedirs(pfm_dir)
file_list=os.listdir(hmm_dir)
for fn in file_list:
    print fn
    hmm_path=hmm_dir+'/'+str(fn)
    pfm_path=pfm_dir+'/'+str(fn).split('.')[0]+'.pfm'
    with open(hmm_path, 'r') as f:
        lines=f.readlines()
        i=0
        flg=False
        with open(pfm_path, 'w') as wf:
            while i < len(lines):
                if flg:
                    tmp=lines[i].strip().split()
                    if len(tmp)>20:
                        wf.write(tmp[0]+'\t')
                        for e in tmp[2:22]:
                            if e=='*':
                                wf.write('0.0000\t')
                            else:
                                wf.write('%.4f\t'%2.0**(-float(e)/1000))
                    elif len(tmp)>1:
                        for e in tmp:
                            if e == '*':
                                wf.write('0.0000\t')
                            else:
                                wf.write('%.4f\t' % 2.0 ** (-float(e) / 1000))
                        wf.write('\n')
                    else:
                        pass
                else:
                    pass
                if '#' in lines[i]:
                    flg=True
                    i+=5
                else:
                    i+=1





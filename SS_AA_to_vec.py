import os
ss_file_path='/home/zhangjun/PycharmProjects/iDRBP_Residues/data/2017_DNA_test_data_seq/2017_DNA_test_ss_aa.out.ss'
aa_file_path='/home/zhangjun/PycharmProjects/iDRBP_Residues/data/2017_DNA_test_data_seq/2017_DNA_test_ss_aa.out.acc'
pfm_dir='/home/zhangjun/PycharmProjects/iDRBP_Residues/data/2017_DNA_test_data_seq/PFM'
vec_dir='/home/zhangjun/PycharmProjects/iDRBP_Residues/data/2017_DNA_test_data_seq/Vectors'
if not os.path.isdir(vec_dir):
    os.makedirs(vec_dir)
with open(ss_file_path, 'r') as f1:
    ss_lines=f1.readlines()
with open(aa_file_path, 'r') as f2:
    aa_lines=f2.readlines()

for i in range(len(ss_lines)):
    ss_name=ss_lines[i].strip().split('>')[1]
    ss=ss_lines[i+1].strip().split('>')[0]
    aa_name = aa_lines[i].strip().split('>')[1]
    aa = aa_lines[i + 1].strip().split('>')[0]
    if ss_name==aa_name and len(aa)==len(ss):
        pfm_path=pfm_dir+'/'+aa_name+'.pfm'
        with open(pfm_path,'r') as f3:
            pfm_lines=f3.readlines()
        vec_path = vec_dir + '/' + aa_name + '.pfm'
        if len(pfm_lines) == len(aa):
            with open(vec_path, 'w') as f4:
                for j in range(len(aa)):
                    tmp_line = pfm_lines[j].strip().split()
                    if aa[j] == 'e':
                        tmp_line.append('1.0000')
                    else:
                        tmp_line.append('0.0000')
                    if ss[j] == 'C':
                        tmp_line.append('1.0000')
                        tmp_line.append('0.0000')
                        tmp_line.append('0.0000')
                    elif ss[j] == 'E':
                        tmp_line.append('0.0000')
                        tmp_line.append('1.0000')
                        tmp_line.append('0.0000')
                    else:
                        tmp_line.append('0.0000')
                        tmp_line.append('0.0000')
                        tmp_line.append('1.0000')
                    f4.write('\t'.join(tmp_line)+'\n')
        else:
            print 'ERROR!!'
            print 'ss_len/aa_len,pfm_len:', len(aa), len(pfm_lines)
    else:
        print 'ERROR!!'
        print 'ss_name,aa_name:', ss_name, aa_name
        print 'ss_len,aa_len:', len(ss), len(aa)

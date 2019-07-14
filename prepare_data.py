from Bio import SeqIO
import numpy as np
# proteins = list(SeqIO.parse('data/2017_DNA_test_data.txt', 'fasta'))
# print proteins[0].seq

# with open('data/2017_RNA_test_data.txt','r') as f:
#     lines=f.readlines()
#     with open('data/2017_RNA_test_data_seq.txt','w') as wf:
#         for i in range(len(lines)):
#             if int(i+1) % 3 == 0:
#                 continue
#             else:
#                 wf.write(lines[i])

# proteins = list(SeqIO.parse('data/2017_DNA_train_data.txt', 'fasta'))
# tmp=0
# for prot in proteins:
#     print len(prot.seq)
#     if len(prot.seq)>tmp:
#         tmp=len(prot.seq)
# print tmp
A=[2,1,3,4]
print [e**2 for e in A]
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import SeqIO
from Bio.PDB import *
import os
import random

AA_dic = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','PHE':'F','TRP':'W','TYR':'Y','ASP':'D','ASN':'N',
          'GLU':'E','LYS':'K','GLN':'Q','MET':'M','SER':'S','THR':'T','CYS':'C','PRO':'P','HIS':'H','ARG':'R'}
AA = ['GLY','ALA','VAL','LEU','ILE','PHE','TRP','TYR','ASP','ASN',
      'GLU','LYS','GLN','MET','SER','THR','CYS','PRO','HIS','ARG']
DA = ['DA', 'DC', 'DT', 'DG']
RA = ['A', 'C', 'T', 'U', 'G']
def download_pdb_ent(ids):
    for pid in ids:
        pdb1 = PDBList()
        pdb1.retrieve_pdb_file(pid, pdir='./examples/', file_format='pdb')
        if os.path.isfile('examples/pdb' + pid + '.ent'):
            parser = PDBParser()
            structure = parser.get_structure(pid, 'examples/pdb' + pid + '.ent')
            print pid, structure.header['head']

def get_pdb_ids():
    pids = []
    prots = SeqIO.parse('drnapred_data/TEST.fasta_seq.txt', 'fasta')
    for prot in prots:
        pids.append(prot.id[:4])
    download_pdb_ent(pids)
    return pids

def calculate_distance(pid):
    parser = PDBParser()
    structure = parser.get_structure(pid, 'examples/pdb' + pid + '.ent')
    atom_list = []
    ids = []
    residues = []
    for model in structure:
        n = 0
        for chain in model:
            if n == 0:
                for residue in chain:
                    if str(residue.get_resname()).strip() in AA:
                        # print residue.get_id()[1]
                        # print residue.get_resname()
                        ids.append(residue.get_id()[1])
                        residues.append(AA_dic[residue.get_resname()])
                        for atom in residue:
                            atom_list.append(atom)
            else:
                for residue in chain:
                    # print str(residue.get_resname()).strip()
                    if (str(residue.get_resname()).strip() in DA) or (str(residue.get_resname()).strip() in RA):
                        for atom in residue:
                            atom_list.append(atom)
                    else:
                        pass
            n += 1

    searcher = NeighborSearch(atom_list=atom_list, bucket_size=10)
    contac = searcher.search_all(radius=3.5, level='R')
    return ids, residues, contac

def check_residue(pid):
    ids = []
    parser = PDBParser()
    structure = parser.get_structure(pid, 'examples/pdb' + pid + '.ent')
    for model in structure:
        for chain in model:
            n = 0
            for residue in chain:
                # print residue
                if n == 0:
                    print residue.get_resname(), residue.get_id()[1]
                if residue.get_resname() in AA:
                    # print residue.get_id()
                    ids.append(residue.get_id()[1])
                    # print AA_dic[residue.get_resname()], residue.get_id()[1]
                    n += 1
                else:
                    pass
                    # print residue.get_resname(), residue.get_id()[1]
            print chain, n
    return ids

def chain_a_labels(prot):
    print prot.id
    ids, residues, contac = calculate_distance(prot.id)
    print len(ids), max(ids)
    for i in range(len(prot.seq)):
        f = 0
        for j in range(i, i+2):
            if list(prot.seq)[j] == residues[j-i]:
                f += 1
        x = random.randint(i, len(residues)-1)
        c = i-int(ids[0])
        print i, x, ids[x], c, len(prot.seq), len(residues)
        if f == 2 and list(prot.seq)[int(ids[x])+c] == residues[x]:
            print i
            break
    contac_labels = {}
    for e in contac:
        if (str(e[0].get_resname()).strip() in AA) and (str(e[1].get_resname()).strip() in RA):
            # print e[0].get_id()[1], str(e[1].get_resname()).strip()
            if int(e[0].get_id()[1]) + c not in contac_labels.keys():
                contac_labels[int(e[0].get_id()[1]) + c] = [0, 1]
            else:
                contac_labels[int(e[0].get_id()[1]) + c][1] = 1
        if (str(e[0].get_resname()).strip() in RA) and (str(e[1].get_resname()).strip() in AA):
            # print e[1].get_id()[1], str(e[0].get_resname()).strip()
            if int(e[1].get_id()[1]) + c not in contac_labels.keys():
                contac_labels[int(e[1].get_id()[1]) + c] = [0, 1]
            else:
                contac_labels[int(e[1].get_id()[1]) + c][1] = 1
        if (str(e[0].get_resname()).strip() in AA) and (str(e[1].get_resname()).strip() in DA):
            # print e[0].get_id()[1], str(e[1].get_resname()).strip()
            if int(e[0].get_id()[1]) + c not in contac_labels.keys():
                contac_labels[int(e[0].get_id()[1]) + c] = [1, 0]
            else:
                contac_labels[int(e[0].get_id()[1]) + c][0] = 1
        if (str(e[0].get_resname()).strip() in DA) and (str(e[1].get_resname()).strip() in AA):
            # print e[1].get_id()[1], str(e[0].get_resname()).strip()
            if int(e[1].get_id()[1]) + c not in contac_labels.keys():
                contac_labels[int(e[1].get_id()[1]) + c] = [1, 0]
            else:
                contac_labels[int(e[1].get_id()[1]) + c][0] = 1
    with open('examples/' + prot.id + '_label.txt', 'w') as f:
        for i in range(len(prot.seq)):
            if i-c not in ids:
                if i not in contac_labels.keys():
                    f.write('0\t0\t1\n')
                else:
                    f.write('%d\t%d\t1\n' % tuple(contac_labels[i]))
            else:
                if i not in contac_labels.keys():
                    f.write('0\t0\t0\n')
                else:
                    f.write('%d\t%d\t0\n' % tuple(contac_labels[i]))

if __name__ == '__main__':
    prots = list(SeqIO.parse('data/5ZQ0.txt', 'fasta'))
    for prot in prots:
        # check_residue(prot.id)
        chain_a_labels(prot)
    # ids = check_residue('5THE')
    # print ids
    # ids = []
    # with open('examples/pdb_ids.txt', 'r') as f:
    #     lines = f.readline().split(', ')
    #     for e in lines:
    #         if e.strip() not in ids:
    #             ids.append(e.strip())
    # print len(ids)
    # download_pdb_ent(['5ZQ0'])#['6J9B','6JBX', '6ON0', '5ZJQ', '6AJK', '6DCC', '6DU5', '5ZQ8']
    # ids, residues, contac = calculate_distance('6DU5')
    # print ids
    # print residues
    # print contac





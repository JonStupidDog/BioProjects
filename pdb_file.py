from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PPBuilder
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import SeqIO
from Bio.PDB import *
import os

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
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_list.append(atom)
    searcher = NeighborSearch(atom_list=atom_list, bucket_size=10)
    contac = searcher.search_all(radius=3.5, level='R')
    return contac

def check_residue(pid):
    parser = PDBParser()
    structure = parser.get_structure(pid, 'examples/pdb' + pid + '.ent')
    for model in structure:
        for chain in model:
            n = 0
            for residue in chain:
                if residue.get_resname() in AA:
                    print AA_dic[residue.get_resname()], residue.get_id()[1]
                    n += 1
            print chain, n



if __name__ == '__main__':
    # check_residue('6OMF')
    contac = calculate_distance('6OMF')
    for e in contac:
        # print e[0].get_resname(), e[1].get_resname()
        if (str(e[0].get_resname()).strip() in AA) and (str(e[1].get_resname()).strip() in RA):
            print e[0].get_id()[1], str(e[1].get_resname()).strip()
        if (str(e[0].get_resname()).strip() in RA) and (str(e[1].get_resname()).strip() in AA):
            print e[1].get_id()[1], str(e[0].get_resname()).strip()
        if (str(e[0].get_resname()).strip() in AA) and (str(e[1].get_resname()).strip() in DA):
            print e[0].get_id()[1], str(e[1].get_resname()).strip()
        if (str(e[0].get_resname()).strip() in DA) and (str(e[1].get_resname()).strip() in AA):
            print e[1].get_id()[1], str(e[0].get_resname()).strip()

    # download_pdb_ent(['6OMF'])





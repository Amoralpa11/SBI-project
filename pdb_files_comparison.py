from Bio.PDB import *
from Bio import pairwise2
import numpy
import gzip
import os
import re
import sys

fasta_p = re.compile(".pdb")

# Creating a list with the file/s passed:
# if options.infile:
#     if os.path.isfile(options.infile):
#         pdb_files.append(options.infile)
#     else:
#         files = filter(fasta_p.search, os.listdir(options.infile))
#         for i in files:
#             pdb_files.append(''.join([options.infile, i]))
# else:
#     pdb_files = filter(fasta_only.search, os.listdir(os.getcwd()))

pdb_files = ["PAIR_HG.pdb", "PAIR_HHGG.pdb", "PAIR_IH.pdb", "PAIR_JC.pdb", "PAIR_JG.pdb", "PAIR_JI.pdb", "PAIR_KH.pdb", "PAIR_LE.pdb",
             "PAIR_LG.pdb", "PAIR_LK.pdb"]
pairwise_interact = {}


def str_comparison_list(str1, str2):
    ls = 0
    ppb = PPBuilder()
    ls_count = [[0, 0], [0, 0]]
    for pp1 in ppb.build_peptides(str1):
        seq1 = pp1.get_sequence()

        if ls == 0:
            ls += 1
        else:
            ls = 0
        print(ls)
        ls2 = 0
        for pp2 in ppb.build_peptides(str2):
            seq2 = pp2.get_sequence()
            alignment = pairwise2.align.globalxx(seq1, seq2)
            score = alignment[0][2]
            ident_perc = score/len(seq1)  # to look at, choose longest one?

            if ident_perc > 0.95:
                ls_count[ls][ls2] = 1
            else:
                ls_count[ls][ls2] = 0

            ls2 += 1
    return ls_count


def str_comparison_superimpose(str1, str2):
    res = 0
    sup = Superimposer()
    print(list(str1.get_atoms()))
    print(list(str2.get_atoms()))
    print("superimposition")
    sup.set_atoms(list(str1.get_atoms()), list(str2.get_atoms()))
    print(str1)
    print(str2)
    print(numpy.abs(sup.rms))
    if numpy.abs(sup.rms) > 2.5:
        res += 1
    return res


def dict_filler(pdb_list, pdb_interact_dict):
    for pdb_file in pdb_list:
        pdb_id, ext = pdb_file.split('.')
        structure = PDBParser().get_structure(pdb_id, pdb_file)
        count_ls = []
        counter = 0
        if pdb_interact_dict:
            count_ls = 0

            for pdb_struct in pdb_interact_dict.values():
                res_ls = str_comparison_list(pdb_struct, structure)
                print(res_ls)
                if 1 in res_ls[0] and 1 in res_ls[1]:
                    counter += str_comparison_superimpose(pdb_struct, structure)

                else:
                    counter += 1

        else:
            pdb_interact_dict[pdb_file] = structure

        if counter == len(pdb_interact_dict):
            pdb_interact_dict[pdb_file] = structure

dict_filler(pdb_files, pairwise_interact)

print(pairwise_interact)

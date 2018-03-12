from Bio.PDB import *
import numpy
import gzip
import os
import re
import sys

fasta_p = re.compile(".pdb")

# Creating a list with the file/s passed:
if options.infile:
    if os.path.isfile(options.infile):
        pdb_files.append(options.infile)
    else:
        files = filter(fasta_p.search, os.listdir(options.infile))
        for i in files:
            pdb_files.append(''.join([options.infile, i]))
else:
    pdb_files = filter(fasta_only.search, os.listdir(os.getcwd()))

pairwise_interact = {}


def str_comparison(str1, str2):
    ls = 0
    ppb = PPBuilder()
    for pp1 in ppb.build_peptides(str1):
        seq1 = pp1.get_sequence()

        if ls == 0:
            ls += 1
        else:
            ls = 0
            ls_count = [[0, 0], [0, 0]]

        ls2 = 0
        for pp2 in ppb.build_peptides(str2):
            seq2 = pp2.get_sequence()
            alignment = pairwise2.align.globalxx(seq1, seq2)
            score = alignment[0][2]
            ident_perc = score / len(seq1) #### to look at, choose longest one?

            if ident_perc > 0.95:
                ls_count[ls][ls2] = 1
            else:
                ls_count[ls][ls2] = 0

            ls2 += 1
    return ls_count


def dict_filler(pdb_list, pdb_interact_dict):
    for i in pdb_list:
        name = i.split('.')
        structure = parser.get_structure(name, i)

        if pdb_interact_dict:
            count_ls = 0
            for pdb_struct in pdb_interact_dict:
                res_ls = str_comparison(pdb_struct, structure)

                if 1 in res_ls[0] and 1 in res_ls[1]:
                    sup = Superimposer()
                    if sup.set_atoms(pdb_struct.get_atoms(), structure.get_atoms()) < 5:
                        break

                    else:
                        count_ls += 1 ## finish the counter to save the str to the dict

        else:
            pdb_interact_dict[i] = structure

    if count_ls == 0:
        pdb_interact_dict[i] = structure

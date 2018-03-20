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




def str_comparison_list(str1, str2):
    '''This function returns a list indicating which chains
    from str1 are the same in str2 based on % of identity.
    It requires 2 arguments, which are the 2 srtuctures we
    want to compare.'''
    ls = 0
    ppb = PPBuilder()
    ls_count = [[0, 0], [0, 0]]
    for pp1 in ppb.build_peptides(str1):
        seq1 = pp1.get_sequence()

        if ls == 0:
            ls += 1
        else:
            ls = 0
        ls2 = 0
        for pp2 in ppb.build_peptides(str2):
            seq2 = pp2.get_sequence()
            alignment = pairwise2.align.globalxx(seq1, seq2)
            score = alignment[0][2]
            length = max(len(seq1), len(seq2))
            ident_perc = score/length  # to look at, choose longest one?

            if ident_perc > 0.95:
                ls_count[ls][ls2] = 1
            else:
                ls_count[ls][ls2] = 0

            ls2 += 1
    return ls_count


def str_comparison_superimpose(str1, str2):
    '''
    This function compares if 2 structures are structurally similar
    and returns a 0 if they are the same or 1 if they are different.
    It requires 2 arguments, the 2 structures we want to compare.
    '''

    chains_list = [x.get_id() for x in str1.get_chains()]
    res = 0
    sup = Superimposer()
    # print(list(str1.get_atoms()))
    # print(list(str2.get_atoms()))
    # print("superimposition")
    mean_distances = []
    for chain_id1 in chains_list:
        for chain_id2 in chains_list:
            if chain_id1.upper() != chain_id2.upper():
                continue
            for round in range(10):
                sup.set_atoms(list(str1[0][chain_id1].get_atoms()), list(str2[0][chain_id2].get_atoms()))
                sup.apply(str2)

            distance_array = []
            other_chain1 = [x for x in chains_list if x != chain_id1][0]
            other_chain2 = [x for x in chains_list if x != chain_id2][0]
            CA_other1 = [x['CA'] for x in str1[0][other_chain1].get_residues() if
                        'CA' in [y.get_id() for y in x.get_atoms()]]
            CA_other2 = [x['CA'] for x in str2[0][other_chain2].get_residues() if
                        'CA' in [y.get_id() for y in x.get_atoms()]]
            for pair in zip(CA_other1,CA_other2):
                
                distance_array.append(pair[0]-pair[1])

            mean_distances.append(sum(distance_array)/len(distance_array))

            if min(mean_distances) < 9:
                print('\t%s' % min(mean_distances))
                return 0
    print('\t%s' % min(mean_distances))
    return 1



def dict_filler(pdb_list, pdb_interact_dict):
    '''
    This functions fills a dictionary which may be empty or not with
    pdb structures which have not yet been saved in the dictionary.
    :param pdb_list: list of pdb files we want to add to the dictionary.
    :param pdb_interact_dict: dictionary which we want to fill with or add structures.
    '''
    for pdb_file in pdb_list:
        pdb_id, ext = pdb_file.split('.')
        structure = PDBParser().get_structure(pdb_id, pdb_file)
        count_ls = []
        counter = 0
        if pdb_interact_dict:
            count_ls = 0

            for pdb_struct in pdb_interact_dict.values():
                res_ls = str_comparison_list(pdb_struct, structure)
                if 1 in res_ls[0] and 1 in res_ls[1]:
                    counter += str_comparison_superimpose(pdb_struct, structure)

                else:
                    counter += 1

        else:
            pdb_interact_dict[pdb_file] = structure

        if counter == len(pdb_interact_dict):
            pdb_interact_dict[pdb_file] = structure

# if __name__ == '__main__':
#
#     pdb_files = ["PAIR_HG.pdb", "PAIR_HHGG.pdb", "PAIR_IH.pdb", "PAIR_JC.pdb", "PAIR_JG.pdb", "PAIR_JI.pdb",
#                  "PAIR_KH.pdb", "PAIR_LE.pdb",
#                  "PAIR_LG.pdb", "PAIR_LK.pdb"]
#     pairwise_interact = {}
#
#     dict_filler(pdb_files, pairwise_interact)
#
#     print(pairwise_interact)

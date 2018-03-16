from Bio.PDB import *
from Bio.PDB import Structure
from Bio.PDB import Model
from Bio.PDB import Chain
from Bio.PDB import Residue
from Bio.PDB import Atom
from Bio import pairwise2
import numpy as np
import string
import gzip
import os
import re
import sys
import argparse


def get_structure_name(filename):

    p = re.compile('(.*).pdb')
    m = p.match(filename)
    return m.group(1)

def get_min_distance(chain1, chain2):
    distance_list = []
    for res1 in chain1.get_residues():
        for res2 in chain2.get_residues():
            distance_list.append(res1['CA']-res2['CA'])

    min_distance = min(distance_list)
    return min_distance

class ChainSelect(Select):

    def __init__(self,chain_to_include1, chain_to_include2):
        self.chain_to_include1 = chain_to_include1
        self.chain_to_include2 = chain_to_include2

    def accept_chain(self,chain):
        if chain.get_id() == self.chain_to_include1 or chain.get_id() == self.chain_to_include2:
            return True
        else:
            return False



def compare_chains(chain1, chain2):

    seq1 = get_sequence_from_chain(chain1)
    seq2 = get_sequence_from_chain(chain2)


    if seq1 and seq2:

        alignment = pairwise2.align.globalxx(seq1, seq2)
        score = alignment[0][2]
        ident_perc = score / len(seq1)

        if ident_perc > 0.95:
            # print(alignment)
            return True
        else:
            return False
    else:
        return False


def get_numeric_array(seq_aln):
    numeric_array = []
    n=0
    for character in seq_aln:
        if character == "-":
            numeric_array.append("-")
        else:
            numeric_array.append(n)
        n+=1
    return numeric_array

def get_sequence_from_chain(chain):
    res_list = [x.get_resname() for x in chain.get_residues() if 'CA' in [y.get_id() for y in x.get_atoms()]]
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    res_short_list = []

    for res in res_list:
        try:
            res_short_list.append(d[res])
        except KeyError:
            return False

    return "".join(res_short_list)

def trim_to_superimpose(chain1, chain2):


    seq1 = get_sequence_from_chain(chain1)
    seq2 = get_sequence_from_chain(chain2)


    alignment = pairwise2.align.globalxx(seq1, seq2)

    score = alignment[0][2]
    ident_perc = score / len(seq1)

    # print("%s-%s" % (chain1.get_id(),chain2.get_id()))
    #
    # print(alignment[0][0])
    # print(alignment[0][1])

    if ident_perc > 0.95:
        seq1_array = list(alignment[0][0])
        seq2_array = list(alignment[0][1])
        seq1_numeric = get_numeric_array(alignment[0][0])
        seq2_numeric = get_numeric_array(alignment[0][1])
        to_delete_from_1 = []
        to_delete_from_2 = []
        pairs1 = zip(seq1_array,seq2_numeric)
        for pair in pairs1:
            if pair[0] == '-':
                to_delete_from_2.append(list(chain2.get_residues())[pair[1]].get_id())

        pairs2 = zip(seq2_array, seq1_numeric)
        for pair in pairs2:
            if pair[0] == '-':
                to_delete_from_1.append(list(chain1.get_residues())[pair[1]].get_id())

        # print(list(chain1.get_residues())[0])
        # print(list(chain2.get_residues())[0])

        for residue_to_delete in to_delete_from_1:
            chain1.__delitem__(residue_to_delete)

        for residue_to_delete in to_delete_from_2:
            chain2.__delitem__(residue_to_delete)

        # print(list(chain1.get_residues())[0])
        # print(list(chain2.get_residues())[0])


def compare_interactions(interaction1, interaction2):
    structure1 = Structure.Structure('1')
    structure2 = Structure.Structure('2')

    structure1.add(Model.Model(0))
    structure2.add(Model.Model(0))

    homodimer = False

    for chain in interaction1:
        if len(list(structure1[0].get_chains())) == 1 and compare_chains(chain,list(structure1[0].get_chains())[0] ):
            homodimer = True

        structure1[0].add(Chain.Chain(chain.get_id()))
        res_counter = 0
        for residue in chain:
            if 'CA' in [x.get_id() for x in residue.get_atoms()]:
                atom = residue['CA']
                structure1[0][chain.get_id()].add(Residue.Residue(('',res_counter,''),residue.get_resname(),residue.get_segid()))

                structure1[0][chain.get_id()][('',res_counter,'')].add(atom.copy())
                res_counter += 1

    for chain in interaction2:

        structure2[0].add(Chain.Chain(chain.get_id()))
        res_counter = 0
        for residue in chain:
            if 'CA' in [x.get_id() for x in residue.get_atoms()]:
                atom = residue['CA']
                structure2[0][chain.get_id()].add(Residue.Residue(('',res_counter,''),residue.get_resname(),residue.get_segid()))

                structure2[0][chain.get_id()][('',res_counter,'')].add(atom.copy())
                res_counter += 1

    if homodimer:
        for int in [structure1[0],structure2[0]]:
            trim_to_superimpose(list(int.get_chains())[0],list(int.get_chains())[1])

    for chain1 in structure1[0]:
        for chain2 in structure2[0]:
            if chain1.get_id() != chain2.get_id():
                continue
            trim_to_superimpose(chain1,chain2)

                # print(list(chain1.get_residues())[0])
                # print(list(chain2.get_residues())[0])




    # print(list(structure1.get_chains()))
    # print(list(structure2.get_chains()))
    result = str_comparison_superimpose(structure1,structure2)

    return result


def get_seq_dict(structure):
    ppb = PPBuilder()
    seq_dict = {}
    for chain in structure.get_chains():

        pp = ppb.build_peptides(chain)

        try:
            seq = pp[0].get_sequence()
            seq_dict[chain.get_id()] = seq
        except IndexError:
            continue
    return seq_dict


fasta_p = re.compile(".pdb")
alphabet = list(string.ascii_uppercase) + list(string.ascii_lowercase)





def str_comparison_list(str1, str2):
    '''This function returns a list indicating which chains
    from str1 are the same in str2 based on % of identity.
    It requires 2 arguments, which are the 2 srtuctures we
    want to compare.
    '''
    ls = 0
    ls_count = [[0, 0], [0, 0]]
    for chain1 in str1.get_chains():
        seq1 = get_sequence_from_chain(chain1)

        if ls == 0:
            ls += 1
        else:
            ls = 0
        ls2 = 0
        print(seq1)
        print(list(str2.get_chains()))
        for chain2 in str2.get_chains():
            seq2 = get_sequence_from_chain(chain2)
            print(seq2)
            alignment = pairwise2.align.globalxx(seq1, seq2)
            score = alignment[0][2]
            print(alignment)
            length = max(len(seq1), len(seq2))
            ident_perc = score / length  # to look at, choose longest one?

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

    chains_list1 = [x.get_id() for x in str1.get_chains()]
    chains_list2 = [x.get_id() for x in str2.get_chains()]
    sup = Superimposer()
    mean_distances = []

    for chain_id1 in chains_list1:
        for chain_id2 in chains_list2:
            # print(compare_chains(str1[0][chain_id1], str2[0][chain_id2]))
            print(chain_id1)
            print(chain_id2)
            if not compare_chains(str1[0][chain_id1], str2[0][chain_id2]):
                continue
            for round in range(10):
                sup.set_atoms(list(str1[0][chain_id1].get_atoms()), list(str2[0][chain_id2].get_atoms()))
                sup.apply(str2)

            distance_array = []
            other_chain1 = [x for x in chains_list1 if x != chain_id1][0]
            other_chain2 = [x for x in chains_list2 if x != chain_id2][0]
            CA_other1 = [x['CA'] for x in str1[0][other_chain1].get_residues() if
                        'CA' in [y.get_id() for y in x.get_atoms()]]
            CA_other2 = [x['CA'] for x in str2[0][other_chain2].get_residues() if
                        'CA' in [y.get_id() for y in x.get_atoms()]]
            for pair in zip(CA_other1, CA_other2):
                pair1 = pair[0].get_coord()
                pair2 = pair[1].get_coord()
                resta = pair1 - pair2
                squares=list(map(lambda x: pow(x, 2), resta))
                dist = np.sqrt(sum(squares))
                # dis2 = pair[0] - pair[1]
                # distance_array.append(pair[0] - pair[1])
                distance_array.append(dist)

            mean_distances.append(sum(distance_array)/len(distance_array))

            if min(mean_distances) < 9:
                print('\t%s' % min(mean_distances))
                return 0
    return 1

def dict_filler(pdb_list, pdb_interact_dict):
    '''
    This functions fills a dictionary which may be empty or not with
    pdb structures which have not yet been saved in the dictionary.
    :param pdb_list: list of pdb files we want to add to the dictionary.
    :param pdb_interact_dict: dictionary which we want to fill with or add structures.
    '''
    for pdb_file in pdb_list:
        m = 0
        n = 0
        pdb_id, ext = pdb_file.split('.')
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(pdb_id, pdb_file)
        counter = 0

        if pdb_interact_dict:

            for pdb_struct in pdb_interact_dict:
                tmp_chain_tup = []
                res_ls = str_comparison_list(pdb_interact_dict[pdb_struct], structure)
                if 1 in res_ls[0] and 1 in res_ls[1]:

                    str_comp = str_comparison_superimpose(pdb_interact_dict[pdb_struct], structure) # 1=diferentes 0=iguales
                    counter += str_comp

                    if str_comp == 1:
                        num = []
                        repeat_comp = []
                        for key in [id for id in pdb_interact_dict if list(id[:2]) == [pdb_struct[0], pdb_struct[1]]]:
                            repeat_comp.append(str_comparison_superimpose(pdb_interact_dict[key], structure))
                            num.append(key[2])
                        if sum(repeat_comp) == len(pdb_interact_dict):
                            tmp_chain_tup = [pdb_struct[0], pdb_struct[1], max(num) + 1]
                            counter = len(pdb_interact_dict)
                            break
                        else:
                            counter = 0
                        break
                    else:
                        break

                elif 1 in res_ls[0]:
                    if pdb_struct[0] not in tmp_chain_tup:
                        tmp_chain_tup.append(pdb_struct[0])
                    counter += 1

                elif 1 in res_ls[1]:
                    if pdb_struct[0] not in tmp_chain_tup:
                        tmp_chain_tup.append(pdb_struct[1])
                    counter += 1

                else:
                    counter += 1

        else:
            seq = []
            for chain1 in structure.get_chains():
                seq.append(get_sequence_from_chain(chain1))

            align = pairwise2.align.globalxx(seq[0], seq[1])
            score = align[0][2] / max([len(seq[0]), len(seq[1])])

            if score < 0.95:
                m += 1
                chain_tup = (alphabet[n], alphabet[m], 0)

            else:
                chain_tup = (alphabet[n], alphabet[m], 0)
            pdb_interact_dict[chain_tup] = structure

            if chain_tup[0] in alphabet:
                alphabet.remove(chain_tup[0])

            if chain_tup[1] in alphabet:
                alphabet.remove(chain_tup[1])

        if counter == len(pdb_interact_dict):
            seq = []

            for pp in structure.get_chains():
                seq.append(get_sequence_from_chain(pp))

            align = pairwise2.align.globalxx(seq[0], seq[1])
            score = align[0][2] / max([len(seq[0]), len(seq[1])])

            if len(tmp_chain_tup) == 0:
                if score < 0.95:
                    m += 1
                    chain_tup = (alphabet[n], alphabet[m], 0)
                else:
                    chain_tup = (alphabet[n], alphabet[m], 0)

            if len(tmp_chain_tup) == 1:
                if score < 0.95:
                    chain_tup = tuple([tmp_chain_tup[0], alphabet[n], 0])
                else:
                    chain_tup = tuple([tmp_chain_tup[0], tmp_chain_tup[0], 0])

            if len(tmp_chain_tup) == 2:
                tmp_chain_tup.append(0)
                chain_tup = tuple(tmp_chain_tup)

            if len(tmp_chain_tup) == 3:
                chain_tup = tuple(tmp_chain_tup)


            pdb_interact_dict[chain_tup] = structure

            if tmp_chain_tup and tmp_chain_tup[0] in alphabet:
                alphabet.remove(tmp_chain_tup[0])
            if tmp_chain_tup and tmp_chain_tup[1] in alphabet:
                alphabet.remove(tmp_chain_tup[1])


if __name__ == '__main__':


    # pdb_files = ["PAIR_HG.pdb", "PAIR_HHGG.pdb", "PAIR_IH.pdb", "PAIR_JC.pdb", "PAIR_JG.pdb",
    #              "PAIR_JI.pdb", "PAIR_KH.pdb", "PAIR_LE.pdb", "PAIR_LG.pdb", "PAIR_LK.pdb"]
    # pdb_files = ["PAIR_HG.pdb", "PAIR_HHGG.pdb", "PAIR_KH.pdb"]

    pdb_files = []
    fasta_p = re.compile(".pdb")

    # Creating a list with the file/s passed:
    folder = ''

    if os.path.isfile(folder):
        pdb_files.append(folder)
    else:
        files = filter(fasta_p.search, os.listdir(folder))
        for i in files:
            pdb_files.append('%s/%s' % (folder, i))


    pairwise_interact = {}
    similar_chains = {}

    dict_filler(pdb_files, pairwise_interact)

    print(pairwise_interact)



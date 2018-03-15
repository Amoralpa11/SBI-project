from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Structure
from Bio.PDB import Model
from Bio.PDB import Chain
from Bio.PDB import Residue
from Bio.PDB import Atom
from Bio.PDB import NeighborSearch
from Bio.PDB import PDBIO
from Bio.PDB import Select
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pdb_files_comparison import str_comparison_superimpose
import copy

sns.set()

"""En este archivo estarán las funciones cuyo proposito sea procesar la información del input que proporcione el usuario. Esto implica, identificar proteinas similares presentes en mas de una interacción, clasificar cuantos tipos de interacciones establece una misma proteina, teniendo en cuenta con qué protein interaccionan y si la superficie de contacto es la misma. """

parser = PDBParser(PERMISSIVE=1)

structure_id = '5vox'
filename = '5vox.pdb'
structure = parser.get_structure(structure_id,filename)

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



def compare_chains(chain1, chain2,seq_dict):

    if chain1.get_id() not in seq_dict or chain2.get_id() not in seq_dict:
        return False

    seq1 = seq_dict[chain1.get_id()]
    seq2 = seq_dict[chain2.get_id()]


    alignment = pairwise2.align.globalxx(seq1, seq2)
    score = alignment[0][2]
    ident_perc = score / len(seq1)

    if ident_perc > 0.95:
        # print(alignment)
        return True
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
    res_list = [x.get_resname() for x in chain.get_residues()]
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    res_short_list = []

    for res in res_list:
        res_short_list.append(d[res])
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


def compare_interactions(interaction1, interaction2,similar_sequences):
    structure1 = Structure.Structure('1')
    structure2 = Structure.Structure('2')

    structure1.add(Model.Model(0))
    structure2.add(Model.Model(0))

    homodimer = False

    for chain in interaction1:
        chain_id = similar_sequences[chain].get_id()
        if chain_id in [x.get_id() for x in structure1.get_chains()]:
            homodimer = True
            if chain_id.upper() == chain_id:
                chain_id = chain_id.lower()
            else:
                chain_id = chain_id.upper

        structure1[0].add(Chain.Chain(chain_id))
        res_counter = 0
        for residue in chain:
            if 'CA' in [x.get_id() for x in residue.get_atoms()]:
                atom = residue['CA']
                structure1[0][chain_id].add(Residue.Residue(('',res_counter,''),residue.get_resname(),residue.get_segid()))

                structure1[0][chain_id][('',res_counter,'')].add(atom.copy())
                res_counter += 1

    for chain in interaction2:

        chain_id = similar_sequences[chain].get_id()
        if chain_id in [x.get_id() for x in structure2.get_chains()]:
            if chain_id.upper() == chain_id:
                chain_id = chain_id.lower()
            else:
                chain_id = chain_id.upper

        structure2[0].add(Chain.Chain(chain_id))
        res_counter = 0
        for residue in chain:
            if 'CA' in [x.get_id() for x in residue.get_atoms()]:
                atom = residue['CA']
                structure2[0][chain_id].add(Residue.Residue(('',res_counter,''),residue.get_resname(),residue.get_segid()))

                structure2[0][chain_id][('',res_counter,'')].add(atom.copy())
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

def get_interaction_pairs (pdb_filename):

    parser = PDBParser(PERMISSIVE=1)

    structure_id = get_structure_name(pdb_filename)
    filename = pdb_filename
    structure = parser.get_structure(structure_id,filename)

    neighbor_chains = {}

    ns = NeighborSearch(list(structure.get_atoms()))

    for chain in structure.get_chains():
        neighbor_chains[chain] = set([])
        for atom in chain.get_atoms():
           for chain2 in ns.search(atom.get_coord(),5,level='C'):
               if chain2 != chain and chain2 not in neighbor_chains.keys():
                   neighbor_chains[chain].add(chain2)
    # print(neighbor_chains)

    similar_sequences = {}
    seq_dict = get_seq_dict(structure)

    chain_list2 = list(structure.get_chains())
    for chain in structure.get_chains():
        if chain not in similar_sequences:
            similar_sequences[chain] = chain
        chain_list2.remove(chain)
        for chain2 in chain_list2:
            cmp = compare_chains(chain,chain2,seq_dict)

            if cmp:

                similar_sequences[chain2] = similar_sequences[chain]
    print(similar_sequences)


    interaction_dict = {}


    for chain1 in neighbor_chains:
        for chain2 in neighbor_chains[chain1]:
            nr_interaction = tuple(sorted([similar_sequences[chain1].get_id(),similar_sequences[chain2].get_id()]))
            if tuple(sorted([similar_sequences[chain1].get_id(), similar_sequences[chain2].get_id()])) not in interaction_dict:
                interaction_dict[nr_interaction] = []

            interaction_dict[nr_interaction].append([chain1,chain2])



    for pair in interaction_dict:
        list_to_remove = []
        print('\n')
        print(pair)
        interaction_list1 = copy.copy(interaction_dict[pair])
        interaction_list2 = copy.copy(interaction_dict[pair])
        for interaction1 in interaction_list1:
            if interaction1 in list_to_remove:
                continue
            print(interaction1)
            for interaction2 in interaction_list2:
                if interaction1 == interaction2:
                    continue
                if interaction2 in list_to_remove:
                    continue
                print('\t%s' % interaction2)
                if not compare_interactions(interaction1,interaction2,similar_sequences):
                    list_to_remove.append(interaction2)
        for interaction in list_to_remove:
            interaction_list1.remove(interaction)
        interaction_dict[pair] = interaction_list1


    if not os.path.exists(structure_id):
        os.makedirs(structure_id)
    else:
        for the_file in os.listdir(structure_id):
            file_path = os.path.join(structure_id, the_file)
            if os.path.isfile(file_path):
                os.unlink(file_path)

    io = PDBIO()
    io.set_structure(structure)

    for pair in interaction_dict:
        for interaction in interaction_dict[pair]:
                io.save('%s/%s_%s%s.pdb' % (structure_id,structure_id,interaction[0].get_id(),interaction[1].get_id()),ChainSelect(interaction[0].get_id(),interaction[1].get_id()))


    # distance_matrix  = np.array([[distance_matrix[chain1][chain2] for chain2 in sorted(distance_matrix[chain1].keys())] for chain1 in sorted(distance_matrix.keys())])
    #
    # distance_matrix = distance_matrix.flatten()
    # # sns.heatmap(distance_matrix,annot=True)
    # sns.distplot(distance_matrix)




# for model in structure.get_list():
#     for chain in model.get_list():
#         io = PDBIO()
#         io.set_structure(structure)
#         io.save('5vox_%s.pdb' % chain.get_id(),ChainSelect(chain.get_id()))

if __name__ == '__main__':



    get_interaction_pairs('5vox.pdb')
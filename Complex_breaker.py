from Bio.PDB.PDBParser import PDBParser
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

def compare_chains(chain1, chain2):
    ppb = PPBuilder()
    pp1 = ppb.build_peptides(chain1)

    pp2 = ppb.build_peptides(chain2)

    try:
        seq1 = pp1[0].get_sequence()
        seq2 = pp2[0].get_sequence()
    except IndexError:
        return False

    alignment = pairwise2.align.globalxx(seq1, seq2)
    score = alignment[0][2]
    ident_perc = score / len(seq1)

    if ident_perc > 0.95:
        return True
    else:
        return False


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
    print(neighbor_chains)

    similar_sequences = {}
    chain_list2 = list(structure.get_chains())
    for chain in structure.get_chains():
        if chain not in similar_sequences:
            similar_sequences[chain] = chain
        chain_list2.remove(chain)
        for chain2 in chain_list2:
            cmp = compare_chains(chain,chain2)

            if cmp:

                similar_sequences[chain2] = similar_sequences[chain]
    print(similar_sequences)


    interaction_dict = {}

    for chain1 in neighbor_chains:
        for chain2 in neighbor_chains[chain1]:
            nr_interaction = tuple(sorted([similar_sequences[chain1].get_id(),similar_sequences[chain2].get_id()]))
            if tuple(sorted([similar_sequences[chain1].get_id(),similar_sequences[chain2].get_id()])) not in interaction_dict:
                interaction_dict[nr_interaction] = []

            interaction_dict[nr_interaction].append([chain1,chain2])


    for pair in interaction_dict:
        interaction_list1 = interaction_dict[pair]
        interaction_list2 = interaction_dict[pair]
        for interaction1 in interaction_list1:
            interaction_list2.remove(interaction1)
            for interaction2 in interaction_list2:

                if compare_interactions(interaction1,interaction2):
                    interaction_list1.remove(interaction2)
                    interaction_list2.remove(interaction2)



    for pair in interaction_dict:
        print(pair)
        print(interaction_dict[pair])




    if not os.path.exists(structure_id):
        os.makedirs(structure_id)
    else:
        for the_file in os.listdir(structure_id):
            file_path = os.path.join(structure_id, the_file)
            if os.path.isfile(file_path):
                os.unlink(file_path)

    for chain1 in neighbor_chains.keys():
        for chain2 in neighbor_chains[chain1]:
                io = PDBIO()
                io.set_structure(structure)
                io.save('%s/%s_%s%s.pdb' % (structure_id,structure_id,similar_sequences[chain1].get_id(),similar_sequences[chain2].get_id()),ChainSelect(similar_sequences[chain1].get_id(),similar_sequences[chain2].get_id()))


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
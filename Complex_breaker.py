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
from pdb_files_comparison import *
import copy

sns.set()

"""En este archivo estarán las funciones cuyo proposito sea procesar la información del input que proporcione el usuario. Esto implica, identificar proteinas similares presentes en mas de una interacción, clasificar cuantos tipos de interacciones establece una misma proteina, teniendo en cuenta con qué protein interaccionan y si la superficie de contacto es la misma. """

parser = PDBParser(PERMISSIVE=1)

structure_id = '5vox'
filename = '5vox.pdb'
structure = parser.get_structure(structure_id, filename)


def get_interaction_pairs(pdb_filename):
    parser = PDBParser(PERMISSIVE=1)

    structure_id = get_structure_name(pdb_filename)
    filename = pdb_filename
    structure = parser.get_structure(structure_id, filename)

    neighbor_chains = {}

    ns = NeighborSearch(list(structure.get_atoms()))

    for chain in structure.get_chains():
        neighbor_chains[chain] = set([])
        for atom in chain.get_atoms():
            for chain2 in ns.search(atom.get_coord(), 5, level='C'):
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
            cmp = compare_chains(chain,chain2)

            if cmp:
                similar_sequences[chain2] = similar_sequences[chain]
    print(similar_sequences)

    interaction_dict = {}

    for chain1 in neighbor_chains:
        for chain2 in neighbor_chains[chain1]:
            nr_interaction = tuple(sorted([similar_sequences[chain1].get_id(), similar_sequences[chain2].get_id()]))
            if tuple(sorted(
                    [similar_sequences[chain1].get_id(), similar_sequences[chain2].get_id()])) not in interaction_dict:
                interaction_dict[nr_interaction] = []

            interaction_dict[nr_interaction].append([chain1, chain2])

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
                if not compare_interactions(interaction1,interaction2):
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
            io.save('%s/%s_%s%s.pdb' % (structure_id, structure_id, interaction[0].get_id(), interaction[1].get_id()),
                    ChainSelect(interaction[0].get_id(), interaction[1].get_id()))

    # distance_matrix  = np.array([[distance_matrix[chain1][chain2] for chain2 in sorted(distance_matrix[chain1].keys())] for chain1 in sorted(distance_matrix.keys())])
    #
    # distance_matrix = distance_matrix.flatten()
    # # sns.heatmap(distance_matrix,annot=True)
    # sns.distplot(distance_matrix)


def get_all_interaction_pairs(pdb_filename):
    parser = PDBParser(PERMISSIVE=1)

    structure_id = get_structure_name(pdb_filename)
    filename = pdb_filename
    structure = parser.get_structure(structure_id, filename)

    neighbor_chains = {}

    ns = NeighborSearch(list(structure.get_atoms()))

    for chain in structure.get_chains():
        neighbor_chains[chain] = []
        for atom in chain.get_atoms():
            for chain2 in ns.search(atom.get_coord(), 5, level='C'):
                if chain2 != chain and chain2 not in neighbor_chains.keys():
                    neighbor_chains[chain].append(chain2)
    # print(neighbor_chains)

    if not os.path.exists('%s_all_interactions' % structure_id):
        os.makedirs('%s_all_interactions' % structure_id)
    else:
        for the_file in os.listdir('%s_all_interactions' % structure_id):
            file_path = os.path.join('%s_all_interactions' % structure_id, the_file)
            if os.path.isfile(file_path):
                os.unlink(file_path)

    io = PDBIO()
    io.set_structure(structure)

    for chain in neighbor_chains:
        for other_chain in neighbor_chains[chain]:
            io.save(
                '%s_all_interactions/%s_%s%s.pdb' % (structure_id, structure_id, chain.get_id(), other_chain.get_id()),
                ChainSelect(chain.get_id(), other_chain.get_id()))


# for model in structure.get_list():
#     for chain in model.get_list():
#         io = PDBIO()
#         io.set_structure(structure)
#         io.save('5vox_%s.pdb' % chain.get_id(),ChainSelect(chain.get_id()))

if __name__ == '__main__': 
    get_all_interaction_pairs('5vox.pdb')

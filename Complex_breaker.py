from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Structure
from Bio.PDB import Model
from Bio.PDB import Chain
from Bio.PDB import Residue
from Bio.PDB import NeighborSearch
from Bio.PDB import PDBIO
from Bio.PDB import Select
from Bio import pairwise2
import re
import os
from pdb_files_comparison import str_comparison_superimpose
import copy

"""In this file there are the functions that will take part in the processing of the user input. That implies 
identify similar proteins in the input pdb files, classifiying in a non redundant way the types of interactions that 
a protein can stablish. """


def get_structure_name(filename):
    """
    Takes a .pdb file name, returns the string without extension
    :param filename: name of the file we want to process
    :return: the string without the extension
    """

    p = re.compile('(.*).pdb')
    m = p.match(filename)
    return m.group(1)


class ChainSelect(Select):
    """

    """

    def __init__(self, chain_to_include1, chain_to_include2):
        self.chain_to_include1 = chain_to_include1
        self.chain_to_include2 = chain_to_include2

    def accept_chain(self, chain):
        if chain == self.chain_to_include1 or chain == self.chain_to_include2:
            return True
        else:
            return False


def compare_chains(chain1, chain2, seq_dict):
    """
    Takes two chain objects and return true if their sequences are more than 95% similar and False if not.
    :param chain1: chain we want to compare.
    :param chain2: other chain we want to compare.
    :param seq_dict: dictionary relating sequences with their chain
    :return: True or False if they are the same or not.
    """
    if seq_dict:
        seq1 = seq_dict[chain1]
        seq2 = seq_dict[chain2]

    if seq1 and seq2:  # If there is an error in the sequecing of a chain the value of seq will be false

        alignment = pairwise2.align.globalxx(seq1, seq2)

        score = alignment[0][2]
        length = max(len(seq1), len(seq2))
        ident_perc = score / length

        if ident_perc > 0.95:
            # print(alignment)
            return True
        else:
            return False
    else:
        return False


def get_numeric_array(seq_aln):
    """
    Takes a sequence alignment and returns an array with numbers instead of residues.
    :param seq_aln: alignment we want to process.
    :return: list with numbers instead of residues.
    """

    numeric_array = []
    n = 0
    for character in seq_aln:
        if character == "-":
            numeric_array.append("-")
        else:
            numeric_array.append(n)
        n += 1
    return numeric_array


def get_sequence_from_chain(chain):
    """
    Takes a chain and returns a string holding its sequence in a single character code
    :param chain: chain we want to know the sequence of.
    :return: sequence of the chain.
    """

    # The next line avoids taking residuals that are not amino acids by only taking those with an alpha carbon
    res_list = [x.get_resname() for x in chain.get_residues() if 'CA' in [y.get_id() for y in x.get_atoms()] or 'P' in [y.get_id() for y in x.get_atoms()]]

    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
         'UNK': 'X', ' DG': 'g', ' DA': 'a', ' DC': 'c', ' DT': 'u',
         'ATV': 'z', '  A': 'a', '  U': 'u', '  G': 'g', '  C': 'c'}

    res_short_list = []

    for res in res_list:

        try:
            res_short_list.append(d[res])
        except KeyError:
            return False

    return "".join(res_short_list)


def trim_to_superimpose(chain1, chain2):
    """
    Takes two chains and removes the residues that do not have a match in the sequence alignment.
    :param chain1: first chain.
    :param chain2:  second chain.
    :return: returns chain 1 and 2 trimmed so that they have the same atom length.
    """
    seq1 = get_sequence_from_chain(chain1)
    seq2 = get_sequence_from_chain(chain2)

    alignment = pairwise2.align.globalxx(seq1, seq2)

    score = alignment[0][2]
    length = max(len(seq1), len(seq2))
    ident_perc = score / length

    # print("%s-%s" % (chain1.get_id(),chain2.get_id()))
    #
    # print(alignment[0][0])
    # print(alignment[0][1])

    # in principle, there will not get here two chains that do not have more than 95% of similarity
    if ident_perc > 0.95:

        seq1_array = list(alignment[0][0])  # Storing the alignment sequences as arrays
        seq2_array = list(alignment[0][1])
        seq1_numeric = get_numeric_array(alignment[0][0])  # Storing those sequences with numbers instead of residuals
        seq2_numeric = get_numeric_array(alignment[0][1])
        to_delete_from_1 = []  # To avoid modifying a list while iterating through it we will store the elements to
        to_delete_from_2 = []  # remove in these lists
        pairs1 = zip(seq1_array, seq2_numeric)  # we pair the sequnce alignment with the numeric list of the other chain
        for pair in pairs1:
            if pair[0] == '-':  # If there is a gap in the sequence, we will remove the nth residue of the other chain
                to_delete_from_2.append(list(chain2.get_residues())[pair[1]].get_id())

        pairs2 = zip(seq2_array, seq1_numeric)
        for pair in pairs2:
            if pair[0] == '-':
                to_delete_from_1.append(list(chain1.get_residues())[pair[1]].get_id())

        # print(list(chain1.get_residues())[0])
        # print(list(chain2.get_residues())[0])

        for residue_to_delete in to_delete_from_1:  # Removing the residuals from the sequences
            chain1.__delitem__(residue_to_delete)

        for residue_to_delete in to_delete_from_2:
            chain2.__delitem__(residue_to_delete)

        # print(list(chain1.get_residues())[0])
        # print(list(chain2.get_residues())[0])


def compare_interactions(interaction1, interaction2, similar_sequences):
    """
    This function takes two structures with two chains each one and a dictionary with chains as keys and keys as
    values relating them if they have more than a 95% of similarity and returns 1 if the two interactions are
    different and 0 if they are the same interaction.
    :param interaction1: one of the interactions you want to compare.
    :param interaction2: the other interaction you want to compare.
    :param similar_sequences: dictionary which relates sequences by similiarity.
    :return: returns true if they are the same and false if they ar enot.
    """

    homodimer = False  # This variable will be true if the chains in the interaction are more than a 95% similar

    chain_list1 = []
    chain_list2 = []

    for chain in interaction1:
        chain_id = similar_sequences[chain].get_id()  # To identify similar chains in the superimposition we name them
        #  as the main chain of its type
        if chain_id in [x.get_id() for x in chain_list1]:
            homodimer = True  # if the second chain is similar to the first we change homodimer to true and
        chain_list1.append(Chain.Chain(chain_id))
        res_counter = 0
        for residue in chain:
            if 'CA' in [x.get_id() for x in residue.get_atoms()]:  # for every residue in chain that have an alpha carbon
                atom = residue['CA']  # storing the alpha carbon
                chain_list1[-1].add(
                    Residue.Residue(('', res_counter, ''), residue.get_resname(), residue.get_segid()))  # adding the
                #  residue
                chain_list1[-1][('', res_counter, '')].add(atom.copy())  # adding a copy of the atom to avoid
                #  modifiying the original ones
                res_counter += 1
            if 'P' in [x.get_id() for x in
                        residue.get_atoms()]:  # for every residue in chain that have an alpha carbon
                atom = residue['P']  # storing the alpha carbon
                chain_list1[-1].add(
                    Residue.Residue(('', res_counter, ''), residue.get_resname(),
                                    residue.get_segid()))  # adding the
                #  residue
                chain_list1[-1][('', res_counter, '')].add(atom.copy())  # adding a copy of the atom to avoid
                #  modifiying the original ones
                res_counter += 1

    for chain in interaction2:  # Doing the same for the structure 2

        chain_id = similar_sequences[chain].get_id()

        chain_list2.append(Chain.Chain(chain_id))
        res_counter = 0
        for residue in chain:
            if 'CA' in [x.get_id() for x in residue.get_atoms()]:
                atom = residue['CA']
                chain_list2[-1].add(
                    Residue.Residue(('', res_counter, ''), residue.get_resname(), residue.get_segid()))

                chain_list2[-1][('', res_counter, '')].add(atom.copy())
                res_counter += 1
            if 'P' in [x.get_id() for x in residue.get_atoms()]:
                atom = residue['P']
                chain_list2[-1].add(
                    Residue.Residue(('', res_counter, ''), residue.get_resname(), residue.get_segid()))

                chain_list2[-1][('', res_counter, '')].add(atom.copy())
                res_counter += 1

    if homodimer:  # if the chain is an homodimer we remove different residues from chains in the same interaction
        for int in [chain_list1, chain_list2]:
            trim_to_superimpose(int[0], int[1])

    for chain1 in chain_list1:  # Removing different residues betwen similar chains in different interactions
        for chain2 in chain_list2:
            if chain1.get_id() != chain2.get_id():
                continue
            trim_to_superimpose(chain1, chain2)

    result = str_comparison_superimpose(chain_list1, chain_list2)

    return result


def get_seq_dict(chain_list):
    """
    Takes a list of chains and returns a dictionary with chain objects as keys
    and strings holding their sequences as values.
    :param chain_list: list with chain objects.
    :return: dictionary with chain objects as keys and strings holding their sequences as values
    """

    seq_dict = {}
    for chain in chain_list:
        seq = get_sequence_from_chain(chain)
        seq_dict[chain] = seq

    return seq_dict


def get_neighbor_chains(structure):
    """
    Takes an structure and returns a dictionary with chains as keys and a list of chains as values holding the
    chains with alpha carbons at less than 8 amstrongs from an alpha carbon of the key chain
    :param structure: structure we want to check for clashes.
    :return: dictionary with the clashes between chains.
    """

    neighbor_chains = {}
    ns = NeighborSearch(list(structure.get_atoms()))
    for chain in structure.get_chains():
        neighbor_chains[chain] = set([])

        neighbor_dict = {}
        for atom in [atom for atom in chain.get_atoms() if
                     atom.get_id() == 'CA' or atom.get_id() == 'P']:  # For every alpha carbon in chain
            for atom2 in ns.search(atom.get_coord(), 8, level='A'):
                if atom2.get_id() == 'CA' or atom2.get_id() == 'P':  # for every alpha carbon at 8 amstrongs or less from atom
                    chain2 = atom2.get_parent().get_parent()  # Gettin to wich chain it belongs
                    if chain2 != chain and chain2 not in neighbor_chains.keys():
                        # If it is not in the same chain and it is not already a
                        if chain2 not in neighbor_dict:
                            neighbor_dict[chain2] = 0
                        neighbor_dict[chain2] += 1

        print('\n%s' % chain)
        for close_chain, contacts in neighbor_dict.items():
            print('%s: %s' % (close_chain, contacts))
            if contacts > 8:
                neighbor_chains[chain].add(close_chain)
    return neighbor_chains


def get_similar_sequences(chain_list, seq_dict):
    """
    Takes a list of chain objects and a dictionary with chains as keys
    and their sequence as values and returns a dictionary with every chain
    as a key and a chain they are similar to as value.
    :param chain_list: list with all the chains passed by the user
    :param seq_dict:
    :return:
    """

    similar_sequences = {}
    chain_list2 = copy.copy(chain_list)
    for chain in chain_list:
        if chain in similar_sequences:  # if the chain is already in the dictionary we skip it
            continue
        if chain not in similar_sequences:  # If not, it is the first time we saw it its value will be itself
            similar_sequences[chain] = chain
        chain_list2.remove(chain)  # Removing the chain from the second list to avoid comparing the same chain
        remove_list = []
        for chain2 in chain_list2:
            if compare_chains(chain, chain2, seq_dict):  # if the chains have more than a 95 % of similarity
                similar_sequences[chain2] = similar_sequences[chain]  # Conecting them in the dict
                remove_list.append(chain2)  # We will remove this chain from the list 2
        #         Todo: discutir si la linea de arriba produce un mal funcionamiento

        for chain in remove_list:
            chain_list2.remove(chain)

    return similar_sequences


def get_interaction_pairs(pdb_filename):
    """
    This function Takes a pdb file path and generates a folder with pdb files holding the unique pairwise
    interactions in the first pdb
    :param pdb_filename:
    :return: ...
    """

    parser = PDBParser(PERMISSIVE=1)

    structure_id = get_structure_name(pdb_filename)
    filename = pdb_filename
    structure = parser.get_structure(structure_id, filename)

    neighbor_chains = get_neighbor_chains(structure)

    seq_dict = get_seq_dict(structure.get_chains())

    similar_sequences = get_similar_sequences(list(structure.get_chains()), seq_dict)

    # print(similar_sequences)

    interaction_dict = {}
    # Here we organize the data in similar_sequences and neighbor_chains in a dictionary with pairs of chain types (
    # an id representing all chains with more than 95% of similarity) with all the pairwise interactions within this
    # two chain types

    for chain1 in neighbor_chains:
        for chain2 in neighbor_chains[chain1]:
            nr_interaction = tuple(sorted([similar_sequences[chain1].get_id(), similar_sequences[chain2].get_id()]))
            if tuple(sorted(
                    [similar_sequences[chain1].get_id(), similar_sequences[chain2].get_id()])) not in interaction_dict:
                interaction_dict[nr_interaction] = []

            interaction_dict[nr_interaction].append([chain1, chain2])

    clean_interaction_dict(interaction_dict, similar_sequences)

    counter = 0
    print('\n')
    for pair in interaction_dict:
        print(pair)
        for int in interaction_dict[pair]:
            print("\t%s" % int)
            counter += 1
    print(counter)

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
                    ChainSelect(interaction[0], interaction[1]))

    return structure_id


def clean_heteroatoms(interaction_dict):
    """
    Function that removes the heteroatoms from the chains in the dictionary passed .
    :param interaction_dict: dictionary with the non-redundant interactions passed by the user.
    :return: Modifies interaction_dict.
    """
    chain_set = set([])
    for interaction_list in interaction_dict.values():
        for interaction in interaction_list:
            for chain in interaction:
                chain_set.add(chain)

    for chain in chain_set:
        for residue in chain.get_residues():
            if residue.get_id()[0] != ' ':
                chain.detach_child(residue.get_id())


def clean_interaction_dict(interaction_dict, similar_sequences):
    """
    Takes an dictionary with tuples of 2 strings and a list of lists of chains and a dictionary with chains as
    keys and similar chains as values and removes chain pairs interacting in a similar way from the interaction_dict.
    :param interaction_dict:
    :param similar_sequences:
    :return:
    """

    # Todo discutir si es conveniente eliminar de la lista dos aquellas interacciones que hayan encontrado interacciones similares
    list_to_inverse = set([])
    for pair in interaction_dict:

        list_to_remove = []  # to avoid modifiying a list while looping through it we store here the elements we want
        #  to remove and do it at the end

        print('\n')
        print(pair)
        interaction_list1 = copy.copy(interaction_dict[pair])
        interaction_list2 = copy.copy(interaction_dict[pair])
        for interaction1 in interaction_list1:
            if interaction1 in list_to_remove:  # if an interaction is repeated we skip it
                continue
            print(interaction1)
            interaction_list2.remove(interaction1)
            for interaction2 in interaction_list2:

                if pair == ('ab', 'ad'):
                    print('hey')

                if interaction2 in list_to_remove:
                    continue
                print('\t%s' % interaction2)

                comparison_result = compare_interactions(interaction1, interaction2, similar_sequences)
                if len(comparison_result) == 1:
                    if True in comparison_result:
                        list_to_remove.append(interaction2)
                else:
                    if False in comparison_result:
                        if True in comparison_result:
                            list_to_remove.append(interaction2)
                            list_to_inverse.add((pair, tuple(interaction1)))
                    else:
                        list_to_remove.append(interaction2)

        for interaction in list_to_remove:
            interaction_list1.remove(interaction)
        interaction_dict[pair] = interaction_list1

    for pair, interaction in list_to_inverse:
        interaction_copy = copy.copy(list(interaction))
        interaction_dict[pair].append(interaction_copy[::-1])

    clean_heteroatoms(interaction_dict)


def get_all_interaction_pairs(pdb_filename, print_files=True):
    """
    Takes a pdb file path and generates a folder with all the pairs of interacting chains without checking if
    there is redundant content. This simulates the user input
    :param pdb_filename:  pdb file with the structure we want to break into interactions
    :param print_files: parameter indicating if we want to output the interaction pairs to a directory.
    :return: a directory with pdb files of the interactions and a list with the first element being the list of all
    interactions, ... to finish this with adri
    """

    parser = PDBParser(PERMISSIVE=1)

    structure_id = get_structure_name(pdb_filename)
    filename = pdb_filename
    structure = parser.get_structure(structure_id, filename)

    neighbor_chains = get_neighbor_chains(structure)

    if print_files:
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
                    '%s_all_interactions/%s_%s%s.pdb' % (
                    structure_id, structure_id, chain.get_id(), other_chain.get_id()),
                    ChainSelect(chain, other_chain))
    else:
        interaction_list = []
        structure_counter = 0
        for chain, neighbor in neighbor_chains.items():
            for chain2 in neighbor:
                new_str = Structure.Structure('%s_%s' % (structure_id, structure_counter))
                structure_counter += 1
                new_str.add(Model.Model(0))
                new_str[0].add(chain)
                new_str[0].add(chain2)
                interaction_list.append(new_str)

        return [interaction_list, 's%s_all_interactions' % structure_id]


def get_pdb_from_directory(directory):
    """
    Takes a directory path and return a list of strings holding the names of pdb files in the directory
    :param directory: directory where the pdb files with the interactions are
    :return: returns a list with the names of the files.
    """

    files_list = ['%s/%s' % (directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

    return files_list


def get_id_list():
    """
    :return: returns a list of all the combinations of two letters
    """

    id_list = []
    alphabets = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
                 'u', 'v', 'w', 'x', 'y', 'z']

    for letter in alphabets:
        for letter2 in alphabets:
            id_list.append('%s%s' % (letter, letter2))

    return id_list[::-1]


def get_id_dict(structure_list):
    """
    Takes a list of structures and returns a dictionary with chains as keys and unique identifiers as values
    :param structure_list: list of structures we want to process.
    :return: dictionary with chains as the keys and their identifier as the value
    """
    """"""

    id_list = get_id_list()
    id_dict = {}
    for structure in structure_list:
        for chain in structure.get_chains():
            id_dict[chain] = id_list.pop()
    return id_dict


def get_interaction_pairs_from_input(directory):
    """
    Takes the path of a directory and returns a list holding the interaction dictionary
    of the pdbs in this directory, a similar chains dictionary and a dictionary that
    relates every chain with its id.
    :param directory: directory from where the pdb files we want to process are.
    :return: list holding the interaction dictionary
    of the pdbs in this directory, a similar chains dictionary and a dictionary that
    relates every chain with its id
    """

    files_list = get_pdb_from_directory(directory)
    structure_list = []

    parser = PDBParser(PERMISSIVE=1)
    for file in files_list:
        structure_id = get_structure_name(file)
        structure = parser.get_structure(structure_id, file)
        if len(list(structure.get_chains())) == 2:
            structure_list.append(structure)
        else:
            structure_list += get_all_interaction_pairs(file, False)[0]

    id_dict = get_id_dict(structure_list)

    chain_list = []
    for structure in structure_list:
        chain_list += list(structure.get_chains())

    seq_dict = get_seq_dict(chain_list)

    similar_sequences = get_similar_sequences(chain_list, seq_dict)

    interaction_dict = {}

    for structure in structure_list:
        chains = list(structure.get_chains())
        nr_interaction = tuple(sorted([id_dict[similar_sequences[chains[0]]], id_dict[similar_sequences[chains[1]]]]))
        if nr_interaction not in interaction_dict.keys():
            interaction_dict[nr_interaction] = []

        interaction_dict[nr_interaction].append(chains)

    clean_interaction_dict(interaction_dict, similar_sequences)

    print('\n')

    counter = 0
    for pair in interaction_dict:
        print(pair)
        for int in interaction_dict[pair]:
            print("\t%s" % int)
            counter += 1
    print(counter)

    # TODO: Eliminar cadenas no utilizadas de similar sequences

    return [interaction_dict,id_dict,similar_sequences, seq_dict]

if __name__ == '__main__':
    get_all_interaction_pairs('5vox.pdb')


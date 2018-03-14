from Bio.PDB import *
import sys

struct = PDBParser().get_structure('X', '5vox.pdb') # load your molecule

def clash_identifier(structure):
    """
    This function checks if there is a clash within a structure.
    It returns 1 (T) if no clash and 0 (F) if there is a clash.
    :param str: the structure we want to evaluate
    :return: T or F depending on if there is a clash or not
    """
    clash = 0
    ns = NeighborSearch(list(structure.get_atoms()))
    for chain in structure.get_chains():
        for atom in chain.get_atoms():
            for chain2 in ns.search(atom.get_coord(), 2, level='C'):
                if chain2 != chain:
                    clash += 0
                else:
                    clash += 1
    if clash > 0:
        clash = 1

    return clash


print(clash_identifier(struct))


def add_chain(fixed_struct, mobile_struct, chain):
    """
    This function adds a new chain to an existing object using superimposition.
    :param fixed_struct: structure that we are updating.
    :param mobile_struct: structure with the chain we want to add to the model.
    :return: returns a structure with the new chain.
    """
    sup = sup.set_atoms(list(str1.get_atoms()), list(str2.get_atoms()))
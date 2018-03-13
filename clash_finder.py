from Bio.PDB import *
import sys

struct = PDBParser().get_structure('X', '5vox.pdb') # load your molecule

def clash_identifier(structure):
    """
    This function checks if there is a clash within a structure
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
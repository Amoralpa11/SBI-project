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
import string
import Complex_id
import pickle
from Complex_breaker import get_sequence_from_chain

parser = PDBParser(PERMISSIVE=1)

structure = parser.get_structure('3kuy', '3kuy.pdb')

for chain in structure.get_chains():
    res_set = set()
    for residue in chain.get_residues():
        ident = residue.id
        if ident[0] == ' ':
            res_set.add(residue.get_resname())
    ls = sorted(res_set)
    if sorted(res_set) == [' DA', ' DC', ' DG', ' DT']:
        print(res_set)

dir(chain)

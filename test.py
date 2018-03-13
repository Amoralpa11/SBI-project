from Bio.PDB import *
from Complex_breaker import *
import pdb_files_comparison

parser = PDBParser(PERMISSIVE=1)

pdb_filename = '2f1d/2f1d_AB.pdb'

structure_id = get_structure_name(pdb_filename)
filename = pdb_filename
structure1 = parser.get_structure(structure_id,filename)

pdb_filename = '2f1d/2f1d_AD.pdb'

structure_id = get_structure_name(pdb_filename)
filename = pdb_filename
structure2 = parser.get_structure(structure_id,filename)

print(pdb_files_comparison.str_comparison_superimpose(structure1,structure2))
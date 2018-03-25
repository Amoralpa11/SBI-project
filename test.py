from Bio.PDB import *
from ResidueDepth_copy import *



from Complex_breaker import *
import pdb_files_comparison


parser = PDBParser(PERMISSIVE=1)

vox = parser.get_structure('5vox','5vox.pdb')

rd = ResidueDepth(vox[0]['A'],5)

print('hey')
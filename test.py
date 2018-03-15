from Bio.PDB import *
from Complex_breaker import *
import pdb_files_comparison


structure1 = Structure.Structure('1')
structure2 = Structure.Structure('2')

structure1.add(Model.Model(0))
structure2.add(Model.Model(0))


structure1[0].add(Chain.Chain('A'))
structure2[0].add(Chain.Chain('A'))

print(hash(structure1[0]['A']))
print(hash(structure2[0]['A']))
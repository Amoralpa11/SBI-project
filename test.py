from Bio.PDB import *
from Complex_breaker import *
import pdb_files_comparison

parser = PDBParser(PERMISSIVE=1)

pdb_filename = '5vox/5vox_AF.pdb'

structure_id = get_structure_name(pdb_filename)
filename = pdb_filename
structure1 = parser.get_structure(structure_id,filename)

pdb_filename = '5vox/5vox_BC.pdb'

structure_id = get_structure_name(pdb_filename)
filename = pdb_filename
structure2 = parser.get_structure(structure_id,filename)

res = 0
sup = Superimposer()
# print(list(str1.get_atoms()))
# print(list(str2.get_atoms()))
# print("superimposition")
for round in range(100):
    sup.set_atoms(list(structure1[0]['A'].get_atoms()), list(structure2[0]['C'].get_atoms()))
    # sup.set_atoms(list(structure1.get_atoms()), list(structure2.get_atoms()))
    sup.apply(structure2)

io = PDBIO()
io.set_structure(structure1)
io.save('proba1.pdb')

iio = PDBIO()
io.set_structure(structure2)
io.save('proba2.pdb')

print (sup.rms)
distance_array = []
for res in [x.get_id() for x in  structure1[0]['F'].get_residues() if 'CA' in [y.get_id() for y in x.get_atoms()]]:
    print(structure1[0]['F'][res].get_id())
    print(structure2[0]['B'][res].get_id())
    distance_array.append(structure1[0]['F'][res]['CA']-structure2[0]['B'][res]['CA'])
    print(structure1[0]['F'][res]['CA']-structure2[0]['B'][res]['CA'])

    print(min(distance_array))
    print(max(distance_array))
    print(sum(distance_array)/len(distance_array))
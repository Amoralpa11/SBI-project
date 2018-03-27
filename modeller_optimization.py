from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
from Bio.PDB import PDBIO


def write_to_pdb(structure):
    """
    This function writes a pdb file from a structure
    :param structure: structure we want to write the pdb file from.
    :return: writes a PDB file to the working directory
    """

    code = structure.id
    io = PDBIO()
    io.set_structure(structure)
    io.save(code + '.pdb')


def structure_optimization(pdb_file):
    """
    This functions optimizes a protein structure saved in a pdb file.
    It optimizes the stereochemistry of the given model including non-bonded contacts.
    :param pdb_file: pdb file with the structure we want to optimize.
    :return: returns a pdb file with optimized structure.
    """

    env = environ()
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    path, ext = pdb_file.split('.')
    dir, code = path.split('/')

    mdl = complete_pdb(env, pdb_file)
    mdl.write(file=code + '.ini')

    # Select all atoms:
    atmsel = selection(mdl)

    # Generate the restraints:
    mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
    mdl.restraints.write(file=code + '.rsr')

    mpdf = atmsel.energy()
    print("The energy of " + code + " is: " + str(mpdf[0]))

    # Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='REPORT')
    md = molecular_dynamics(output='REPORT')

    # Open a file to get basic stats on each optimization
    trcfil = open(code + '.D00000001', 'w')

    # Run CG on the all-atom selection; write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))
    # Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
    # 10 steps during the run, and write stats every 10 steps
    md.optimize(atmsel, temperature=300, max_iterations=50,
                actions=[actions.write_structure(10, code + '.D9999%04d.pdb'),
                         actions.trace(10, trcfil)])
    # Finish off with some more CG, and write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20,
                actions=[actions.trace(5, trcfil)])

    mpdf = atmsel.energy()
    print("The energy of " + code + " is: " + str(mpdf[0]))

    mdl.write(file=dir + '/' + code + '_optimized' + '.' + 'pdb')


if __name__ == "__main__":
    structure_optimization("5vox.pdb")

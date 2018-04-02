from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
from Bio.PDB import PDBIO
from plotting import energy_profile_plot


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


def modeller_funcs(pdb_file, options):
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
    # mdl.write(file=code + '.ini')

    # Select all atoms:
    atmsel = selection(mdl)

    # Calculate initial energy and energy profile for the built structure
    mpdf_ini = atmsel.energy()
    z_score_ini = mdl.assess_normalized_dope()
    mdl_ep_ini = atmsel.get_dope_profile()
    mdl_ep_ini_smoothed = mdl_ep_ini.get_smoothed()
    energy_profile_txt_path = dir + '/' + code + '_DOPE_EnergyProfile.txt'
    mdl_ep_ini_smoothed.write_to_file(energy_profile_txt_path)
    print("\nModel energy")
    print("The unoptimized model energy of " + code + " is: " + str(mpdf_ini[0]))
    print("\nZ-score")
    print("The unoptimized Z-score of " + code + " is: " + str(z_score_ini))

    if options.optimization is None:
        energy_profile_plot(options, path, energy_profile_txt_path)

    else:
        # Create optimizer objects and set defaults for all further optimizations
        cg = conjugate_gradients(output='NO_REPORT')

        # Open a file to get basic stats on each optimization
        trcfil = open(dir + '/optimization_stats_' + code + '.txt', 'w')

        # Run CG on the all-atom selection; write stats every 5 steps
        cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))

        # Finish off with some more CG, and write stats every 5 steps
        cg.optimize(atmsel, max_iterations=20,
                    actions=[actions.trace(5, trcfil)])

        mpdf = atmsel.energy()
        # Calculate the normalized Z-score for the model after optimization
        z_score = mdl.assess_normalized_dope()
        print("\nModel energy")
        print("The final energy of " + code + " is: " + str(mpdf[0]))
        print("\nZ-score")
        print("The final z-score of " + code + " is: " + str(z_score))

        # Getting the energy profile of our optimized model
        mdl_ep_fin = atmsel.get_dope_profile()

        # Smooth the energy profile by applying a smoothing window of 50
        mdl_ep_fin_smoothed = mdl_ep_fin.get_smoothed(window=50)
        energy_profile_txt_path_opt = path + '_optimized_DOPE_EnergyProfile.txt'
        mdl_ep_fin_smoothed.write_to_file(energy_profile_txt_path_opt)
        mdl.write(file=path + '_optimized.pdb')
        energy_profile_plot(options, path, energy_profile_txt_path, energy_profile_txt_path_opt)

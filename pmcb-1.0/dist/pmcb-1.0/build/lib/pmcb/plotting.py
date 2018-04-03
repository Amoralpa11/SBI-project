import matplotlib.pyplot as plt
import pandas as pd

def energy_profile_plot(options, dir, code, data_unopt, data_opt=None):
    """
    This function recieves 1 or 2 datasets and draws the energy profiles.
    :param data_unopt: data of the not optimized structure.
    :param data_opt: data from the optimized structure
    :return: graph of the energy profiles.
    """
    plt.clf()
    data = pd.read_csv(data_unopt, sep="\s+", header=None)
    line_norm, = plt.plot(data[0], data[1], linewidth=1, label='Not Optimized')

    if options.optimize:
        data_o = pd.read_csv(data_opt, sep="\s+", header=None)
        line_opt, = plt.plot(data_o[0], data_o[1], linewidth=1, label='Optimized')
        plt.legend(handles=[line_norm, line_opt])

    plt.title(code + ' energy profile')
    plt.xlabel('Residue')
<<<<<<< HEAD:plotting.py
    plt.ylabel('DOPE energy')
=======
    plt.ylabel('DOPE_energy')
>>>>>>> b9bd86b097a10685696b32a17540b58b47d12c97:pmcb-1.0/dist/pmcb-1.0/build/lib/pmcb/plotting.py
    plt.axhline(linewidth=1, linestyle=':', color='r')

    # plt.show()
    plt.savefig(dir + '/' + code + '_EnergyProfile_plot.tiff', dpi=300)

    return

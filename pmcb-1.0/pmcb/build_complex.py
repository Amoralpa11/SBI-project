#! python

from pmcb.macrocomplex_builder import *
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="This program receives fasta files and returns an ordered list by sequence length of "
                    "id+length+molecular weight")
    parser.add_argument('-i', '--input',
                        dest="infile",
                        action="store",
                        help="Input FASTA formatted file or a directory containing fasta files")

    parser.add_argument('-k', '--subunit_limit',
                        dest="subunit_n",
                        action="store",
                        type=int,
                        help="number of subunits to form the macro-complex in case it is theoretically endless.")

    parser.add_argument('-v', "--verbose",
                        dest="verbose",
                        help="increase output verbosity",
                        action="store_true")

    parser.add_argument('-opt', "--optimize",
                        dest="optimize",
                        help="optimize the macro-complex",
                        action="store_true")

    parser.add_argument('-int', "--intensive",
                        dest="intensive",
                        help="Perform an intensive search or just return the 1st structure found.",
                        action="store_true")

    parser.add_argument('-st', "--stoichiometry",
                        dest="st",
                        help="Allows the user to pass the stoichiometry once the chains have been processed.",
                        default=False,
                        action='store_true')

    parser.add_argument('-br', '--break',
                        dest="break_complex",
                        action='store',
                        help='Indicate if you want to obtain all the pairwise interactions, "all", or just one of '
                             'each type, "unique". '
                        )
    parser.add_argument('-dir', '--directory',
                        dest="dir",
                        action='store',
                        help='Indicate a name to create a directory where the outputs will be written'
                        )

    options = parser.parse_args()


    class MyException(Exception):
        pass


    if options.infile is None:
        raise MyException("Input directory is a required argument.")

    if options.dir is None and options.break_complex is None:
        raise MyException("New output directory is a required argument.")

    if options.break_complex not in [None, 'all', 'unique']:
        raise MyException("br if set must be 'all' or 'unique'.")

    if not options.break_complex:
        result = get_interaction_pairs_from_input(options.infile)
        id_dict = result[1]
        interaction_dict = result[0]
        similar_sequences = result[2]
        seq_dict = result[3]
        macrocomplex_builder(id_dict, similar_sequences, interaction_dict, seq_dict, options)


    elif options.break_complex == "all":
        get_all_interaction_pairs(options.infile)

    elif options.break_complex == "unique":
        get_interaction_pairs(options.infile)

if __name__ == "__main__":

    main()

    


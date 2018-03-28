if __name__ == "__main__":
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
                        help="number of subunits to form the macro-complex in case it is theoretically endless.")

    parser.add_argument('-v', "--verbose",
                        dest="verbose",
                        help="increase output verbosity",
                        action="store_true")

    parser.add_argument('-opt', "--optimize",
                        dest="optimize",
                        help="optimize the macro-complex",
                        default=True,
                        action="store_true")

    parser.add_argument('-int', "--intensive",
                        dest="intensive",
                        help="Perform an intensive search or just return the 1st structure found.",
                        action="store_true")

    parser.add_argument('-br', '--break',
                        dest="break_complex",
                        action='Indicate if you want to obtain all the pairwise interactions, "all", or just one of each type, "unique".',
                        default=None,
                        help="regex to be searched for in the sequence"
                        )

    parser.add_argument('-clash_par', '--clash_parameters',
                        dest="clash_parameters",
                        action="store",
                        help="File with the clash parameters specifications."
                        )

    parser.add_argument('-inter_par', '--interaction_parameters',
                        dest="interaction_parameters",
                        action="store",
                        help="File with the interaction parameters specifications."
                        )

    parser.add_argument('-stoic', '--stoichiometry',
                        dest="stoichiometry",
                        action="store",
                        help="Parameter defining the stoichiometry of the complex, the program will assume the "
                             "different interaction pdb files passed are all the interactions forming the "
                             "macro-complex. "
                        )
    options = parser.parse_args()

    # get_all_interaction_pairs('')
    print('starting')
    if !options.break_complex:
        result = get_interaction_pairs_from_input(options.infile)
        if stoichiometry:

    elif options.break_complex == "all":
        directory = get_all_interaction_pairs(options.infile)[1]
        result = get_interaction_pairs_from_input(directory)

    elif options.break_complex == "unique":
        directory = get_interaction_pairs(options.infile)
        result = get_interaction_pairs_from_input(directory)

    else:
        # TODO raise exception saying option passed to break is not valid.
        raise Exception

    result = get_interaction_pairs_from_input('5oom_all_interactions')
    print(result)
    id_dict = result[1]
    interaction_dict = result[0]
    similar_sequences = result[2]

    macrocomplex_builder(id_dict, similar_sequences, interaction_dict, options.infile)

    structures_generated = generate_struct_from
    CI(CI_dict)

    for struct in structures_generated:
        write_to_pdb(struct)
        if optimize:
            structure_optimization(str(struct.id) + '.pdb')
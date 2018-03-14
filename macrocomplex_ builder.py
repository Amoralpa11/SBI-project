for interact in str_dict.values():
    macrocomplex_builder(interact)



def macrocomplex_builder(prot_srtuct, sup_dict, ):
    """
    This function rebuilds a complex with the interactions we obtained from the pdb files.
    :param str_dict: dictionary with all the interactions we want to build the complex with.
    :return: returns a XXXXXXX with the macro-complex built and... the middle steps??
    """

    for interact in str_dict.values():
        sup = Superimposer()
        sup.set_atoms(list(prot_srtuct.get_atoms()), list(interact.get_atoms()))

        if not clash_identifier(sup) and numpy.abs(sup.rms) < 10:
            macrocomplex_builder(sup)
        else:
            sup_dict[sup] = sup

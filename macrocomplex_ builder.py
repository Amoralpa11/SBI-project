def update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict):

    for other_CI in [ ident for ident in complex_id_dict[len(complex_id.nodes())] if ident != complex_id]:
        if complex_id.compare_with(other_CI):  # how to compare CI...
            return

    for nodes in complex_id.nodes():
        for interact in [ interaction[0] for interaction in nodes.get_interact_dict().items() if interaction[1] is None]:
            if similar_seq[interact[0]] == similar_seq[interact[1]]:
                # TODO function to copy chains without copying id of the chain
                chain_str2_copy = copy.deepcopy(interact[0])
                superimpose_fun(base_struct, interact, nodes, chain_str2_copy, complex_id)
                update_structure(base_struct, complex_id)
            else:
                for i in interact:
                    if similar_seq[chains_str_dict[nodes.chain_type()]] == similar_seq[i]:
                        superimpose_fun(base_struct, interact, nodes, i, complex_id)
                        update_structure(base_struct, complex_id.copy())
                        break
        # TODO if the structure is modified when we go back a node
        if structure modified:
            pop_
#########

def superimpose_fun(str1, str2, node, chain_str2, complex_id):
    """
    This functions superimposes 2 structure using nodes
    :param str1: fixed structure
    :param str2: mobile structure
    :param node: node in which we are working
    :param chain_str2: chain of str2 we want to superimpose
    :param complex_id: complex_id information
    :return: complex_id with new node if superimposition is feasible or a clash if it is not.
    """
    chain1 = node.get_chain()
    str2_copy = copy.deepcopy(str2)

    trim_to_superimpose(str1[0][chain1], str2[0][chain_str2])
    atoms_chain1 = list(str1[0][chain1].get_atoms())
    atoms_chain2 = list(str2[0][chain_str2].get_atoms())
    sup = Superimposer()

    sup.set_atoms(atoms_chain1, atoms_chain2)
    sup.apply(str2_copy)

    other_chain2 = [x for x in str2_copy.get_chains() if x != chain_str2]
    # CA_other2 = [x['CA'] for x in str2_copy[0][other_chain2].get_residues() if
    #             'CA' in [y.get_id() for y in x.get_atoms()]]

    distance_array = []
    # TODO get neighbour chains
    if no clashes:
        complex_id.add_node(other_chain2, node, str2)
        str1[0].add(other_chain2)
        similar_sequences[other_chain2.get_id()] = similar_sequences[chain_str2]
    else:
        node.add_interaction("clash", str2)

#########

def macrocomplex_builder(str_dict, id_dict, similar_seq, interaction_dict):
    """
    This function rebuilds a complex with the interactions we obtained from the pdb files.
    :param str_dict: dictionary with all the interactions we want to build the complex with.
    :param id_dict: dictionary with all the chains with their specific key.
    :return: returns a XXXXXXX with the macrocomplex built and... the middle steps??
    """

    # returns a set with all the chains passed by the user
    chains = set([item for sublist in [x for x in str_dict] for item in sublist])

    # unique structure selection
    unique_chains = set(similar_seq.values())

    # Create a dictionary with unique chains with chain name as key and str in the value
    chains_str_dict = {v: k for k, v in id_dict.items() if k in unique_chains}

    # initialize a complex id dictionary
    complex_id_dict = {}

    for chain in chains_str_dict:
        # initialize an empty structure
        base_struct = Structure.Structure('1')
        base_struct.add(Model.Model(0))
        copy_chain = copy.deepcopy(Chain.Chain(chains_str_dict[chain]))
        # add chain to the new structure
        base_struct[0].add(copy_chain)
        similar_seq[1] = similar_seq[chain]
        complex_id = ComplexId(interaction_dict, id_dict, similar_seq, base_struct)
        complex_id_dict[len(complex_id.nodes())] = complex_id
        update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict)

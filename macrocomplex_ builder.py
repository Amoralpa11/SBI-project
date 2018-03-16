def str_comparison_superimpose(str1, str2):
    '''
    This function compares if 2 structures are structurally similar
    and returns a 0 if they are the same or 1 if they are different.
    It requires 2 arguments, the 2 structures we want to compare.
    '''

    chains_list = [x.get_id() for x in str1.get_chains()]
    res = 0
    sup = Superimposer()
    mean_distances = []
    for chain_id1 in chains_list:
        for chain_id2 in chains_list:
            if chain_id1.upper() != chain_id2.upper():
                continue
            for round in range(10):
                sup.set_atoms(list(str1[0][chain_id1].get_atoms()), list(str2[0][chain_id2].get_atoms()))
                sup.apply(str2)

            distance_array = []
            other_chain1 = [x for x in chains_list if x != chain_id1][0]
            other_chain2 = [x for x in chains_list if x != chain_id2][0]
            CA_other1 = [x['CA'] for x in str1[0][other_chain1].get_residues() if
                         'CA' in [y.get_id() for y in x.get_atoms()]]
            CA_other2 = [x['CA'] for x in str2[0][other_chain2].get_residues() if
                         'CA' in [y.get_id() for y in x.get_atoms()]]
            for pair in zip(CA_other1, CA_other2):
                distance_array.append(pair[0] - pair[1])

            mean_distances.append(sum(distance_array) / len(distance_array))

            if min(mean_distances) < 9:
                print('\t%s' % min(mean_distances))
                return 0
    print('\t%s' % min(mean_distances))
    return 1


def update_structure(chain, str1, str2):
    """
    This functions returns the superimposed structure between str1 and str2.
    :param chain: chain in common for superimposition.
    :param str1: fixed structure.
    :param str2: mobile structure.
    :return: str1 with str2 updated.
    """
    sup = Superimposer()
    for round in range(10):
        sup.set_atoms(list(str1[0][chain].get_atoms()), list(str2[0][chain].get_atoms()))
        sup.apply(str2)

    return structure_updated


def macrocomplex_builder(str_dict, id_dict, similar_seq, complex_id_dict):
    """
    This function rebuilds a complex with the interactions we obtained from the pdb files.
    :param str_dict: dictionary with all the interactions we want to build the complex with.
    :param id_dict: dictionary with all the chains with their specific key.
    :return: returns a XXXXXXX with the macrocomplex built and... the middle steps??
    """

    # returns a set with all the chains passed by the user
    chains = set([item for sublist in [x for x in str_dict] for item in sublist])

    # Create a dictionary with structures equivalent only to the chain in the key
    chains_str_dict = {v: k for k, v in id_dict.items() if id in chains}

    for chain in chains_str_dict:
        base_struct = Structure.Structure('1')
        base_struct.add(Model.Model(0))
        copy_chain = copy.deepcopy(Chain.Chain(chains_str_dict[chain]))
        base_struct[0].add(copy_chain)
        #complex_id -- build for a chain
        update_structure(base_struct, complex_id)


##########

def update_structure(base_struct, complex_id):
    if complex_id not in complex_id_dict[complex_id.nodes()]:
        for nodes in complex_id.nodes():
            for interact in nodes.interact_list():
                if nodes.interact_list[interact] == None:
                    if similar_seq[interact[0]] == similar_seq[interact[1]]:
                        chain_str2_copy = copy.deepcopy(interact[0])
                        superimpose_fun(base_struct, interact, node, chain_str2_copy , complex_id)
                        update_structure(base_struct, complex_id)
                    else:
                        for i in interact:
                            if similar_seq[chains_str_dict[node.chain_type]] == similar_seq[i]:
                                chain_str2_copy = copy.deepcopy(i)
                                superimpose_fun(base_struct, interact, node, chain_str2_copy, complex_id)
                                update_structure(base_struct, complex_id)
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
    chain2 = chain_str2
    atoms_chain1 = list(str1[0][chain1].get_atoms())
    atoms_chain2 = copy.deepcopy(list(str2[0][chain2].get_atoms()))
    sup = Superimposer()

    sup.set_atoms(atoms_chain1, atoms_chain2)
    sup.apply(str2)

    other_chain1 = [x for x in chains_list if x != chain1][0]
    other_chain2 = [x for x in chains_list if x != chain2][0]

    CA_other1 = [x['CA'] for x in str1[0][other_chain1].get_residues() if
                 'CA' in [y.get_id() for y in x.get_atoms()]]
    CA_other2 = [x['CA'] for x in str2[0][other_chain2].get_residues() if
                 'CA' in [y.get_id() for y in x.get_atoms()]]


    distance_array = []
    for pair in zip(CA_other1, CA_other2):
        distance_array.append(pair[0] - pair[1])

    mean_distances.append( sum(distance_array)/len(distance_array) )

    if mean_distances > 2:
        complex_id.add_node(node, str2, chain_str2)
    else:
        complex_id.add_node(node, str2, "clash_royale")
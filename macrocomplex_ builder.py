form Bio.PDB import *
from  Complex_breaker import *

def get_clash_chains(structure, chain):

    """
    This function recieves a hain and a structure and calculates if there is a clash between these.
    :param structure: structure where we want to add the chain
    :param chain: chain object we want to add to the structure
    :return: True or false, True if there is clash and false if there is no clash.
    """
    # center_residues = chain.get_residues()
    chain_atoms = chain.get_atoms()
    # chain_atoms = Selection.unfold_entities(center_residues, 'A')
    atom_list = structure.get_atoms()
    ns = NeighbourSearch(atom_list)
    # clashing_chains = {res for chain_atoms in chain_atoms
    #                    for res in ns.search(chain_atoms.get_coord(), 1.2, 'C')}
    clashing_chain_ls = []
    for atom1 in chain_atoms:
        for chain in ns.search(atom1.get_coord(), 1.2, 'C'):
            # if chain != atom1.get_parent().get_parent():
            clashing_chain_ld.append(chain)

    if clashing_chain_ls:
        return True
    else:
        return False

#########

def interaction_finder(structure, ref_chain_id, complex_id):
    """
    This function recieves a hain and a structure and calculates if there is a clash between these.
    :param structure: structure where we want to add the chain
    :param chain: chain object we want to add to the structure
    :return: True or false, True if there is clash and false if there is no clash.
    """
    neighbor_chains = []
    ns = NeighborSearch(list(structure.get_atoms()))
    ref_chain = structure[0][ref_chain_id]
    for atom in [atom for atom in ref_chain.get_atoms() if atom.get_id() == 'CA']:  # For every alpha carbon in chain
        for atom2 in ns.search(atom.get_coord(), 8, level='A'):
            if atom2.get_id() == 'CA':  # for every alpha carbon at 8 armstrongs or less from atom
                chain2 = atom2.get_parent().get_parent()  # Getting to which chain it belongs
                if chain2 != ref_chain and chain2 not in neighbor_chains:
                    neighbor_chains.add(chain2)  # If it is not in the same chain and it is not already a
                    # key (we already have its interactions) we add the chain as a value

    for chain in neighbor_chains:
        for interaction in complex_id.get_interaction_dict():
            if complex_id.id_dict[complex_id.similar_sequences[chain]] in interaction and complex_id.nodes[-1].get_chain_type() in interaction:
                if not compare_interactions([structure[0][ref_chain_id], chain], interaction, complex_id.similar_sequences):
                    complex_id.nodes[-1].add_interaction(complex_id[chain.get_id()], interaction)
                    complex_id[chain.get_id()].add_interaction(complex_id.nodes[-1], interaction)
                    break

#########

def copy_chain(chain, id):
    """
    This function creates a new chain which is a copy of the passed one but with a different ID
    :param chain: chain we want to copy
    :return: copy of the passed chain
    """

    new_chain = copy.deepcopy(chain)
    new_chain.id = id
    return new_chain
#########

def superimpose_fun(str1, str2, node, chain_str2, complex_id, similar_seq):
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

    if not get_clash_chains(str1, other_chain2): ## returns T if there is a clash and F if there isn't.
        complex_id.add_node(other_chain2, node, str2)
        str1[0].add(other_chain2)
        similar_seq[other_chain2.get_id()] = similar_seq[chain_str2]
        interaction_finder(str1, other_chain2.get_id(), complex_id)
        return True
    else:
        node.add_interaction("clash", str2)

#########

def update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict):
    modified_str = 0
    for other_CI in [ ident for ident in complex_id_dict[len(complex_id.nodes())] if ident != complex_id]:
        if complex_id.compare_with(other_CI):
            return

    for nodes in complex_id.nodes():
        for interact in [ interaction[0] for interaction in nodes.get_interact_dict().items() if interaction[1] is None ]:

            if similar_seq[interact[0]] == similar_seq[interact[1]]:
                # TODO len de complex id
                chain_str2_copy = copy_chain(interact[0], len(complex_id.get_nodes()))
                modified_str = superimpose_fun(base_struct, interact, nodes, chain_str2_copy, complex_id, similar_seq)

                if modified_str:
                    update_structure(base_struct, complex_id)
                    complex_id.pop_structure(base_struct)

            else:
                for i in interact:
                    if similar_seq[chains_str_dict[nodes.chain_type()]] == similar_seq[i]:
                        chain_str2_copy = copy_chain(i, len(complex_id.get_nodes()))
                        modified_str = superimpose_fun(base_struct, interact, nodes, chain_str2_copy, complex_id, similar_seq)

                        if modified_str:
                            update_structure(base_struct, complex_id)
                            complex_id.pop_structure(base_struct)
                        break

def macrocomplex_builder(str_dict, id_dict, similar_seq, interaction_dict):
    """
    This function rebuilds a complex with the interactions we obtained from the pdb files.
    :param str_dict: dictionary with all the interactions we want to build the complex with.
    :param id_dict: dictionary with all the chains with their specific key.
    :return: returns a XXXXXXX with the macrocomplex built and... the middle steps??
    """

    # # returns a set with all the chains passed by the user
    # chains = set([item for sublist in [x for x in str_dict] for item in sublist])

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
        # copy_chain = copy.deepcopy(Chain.Chain(chains_str_dict[chain]))
        chain_copied = copy_chain(chains_str_dict[chain], 1)
        # add chain to the new structure
        base_struct[0].add(chain_copied)
        similar_seq[1] = similar_seq[chain]
        complex_id = ComplexId(interaction_dict, id_dict, similar_seq, base_struct)
        complex_id_dict[len(complex_id.nodes())] = complex_id
        update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict)

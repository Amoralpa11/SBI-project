from Bio.PDB import *
from Complex_breaker import *
from Complex_id import *
from ResidueDepth_copy import *


branch_id = [1]
pdb_counter = 1

def write_to_pdb(structure):
    """
    This function writes a pdb file from a structure
    :param structure: structure we want to write the pdb file from.
    :return: writes a PDB file to the working directory
    """
    final_structure = Structure.Structure(1)
    model_counter = 0

    id_list = []

    for chain in structure.get_chains():

        if not id_list:
            id_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                       'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B',
                       'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
                       'V', 'W', 'X', 'Y', 'Z'][::-1]
            model_counter +=1
            final_structure.add(Model.Model(model_counter))

        new_chain = chain.copy()
        new_chain.id = id_list.pop()

        final_structure[model_counter].add(new_chain)

    global pdb_counter
    code = pdb_counter
    pdb_counter += 1
    io = PDBIO()
    io.set_structure(final_structure)
    io.save('result/'+str(code) + '.pdb')

def get_clash_chains(structure, chain, prev_chain):
    """
    This function recieves a hain and a structure and calculates if there is a clash between these.
    :param structure: structure where we want to add the chain
    :param chain: chain object we want to add to the structure
    :return: True or false, True if there is clash and false if there is no clash.
    """
    # center_residues = chain.get_residues()
    chain_atoms = [x for x in list(chain.get_atoms()) if x.get_id() == 'CA']
    # chain_atoms = Selection.unfold_entities(center_residues, 'A')
    atom_list = list(structure.get_atoms())
    ns = NeighborSearch(atom_list)
    # clashing_chains = {res for chain_atoms in chain_atoms
    #                    for res in ns.search(chain_atoms.get_coord(), 1.2, 'C')}

    clash_counter = 0
    # first_clash = True
    for atom1 in chain_atoms:
        atom_produces_clash = False
        for atom in ns.search(atom1.get_coord(), 1.2, 'A'):
            # print('%s(%s)-%s(%s)' % (atom1.get_id(),atom1.get_parent().get_id(), atom.get_id(),atom1.get_parent().get_id()))
            # if chain != atom1.get_parent().get_parent():
            clashing_chain = atom.get_parent().get_parent().get_id()
            if clashing_chain != prev_chain:
                # if first_clash:
                #     first_clash = False
                #     rd = ResidueDepth(chain,4)
                # ca_dep = rd[(chain.get_id(),atom1.get_parent().get_id())][1]
                atom_produces_clash = True
                break
        if atom_produces_clash:
            clash_counter += 1
            if clash_counter < 5:
                print('more than 5 clashes found')
                # print('hey')
                return True

    return False


#########

def interaction_finder(structure, ref_chain_id, complex_id, node):
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
                if chain2 != ref_chain and chain2 not in neighbor_chains and chain2.get_id() != node.get_chain():
                    neighbor_chains.append(chain2)  # If it is not in the same chain and it is not already a
                    # key (we already have its interactions) we add the chain as a value
    print('%s interactions found' % (len(neighbor_chains)))
    for chain in neighbor_chains:
        tup = sorted([complex_id.id_dict[complex_id.similar_sequences[chain]],complex_id.nodes[-1].get_chain_type()])
        tup = tuple(tup)

        if tup in complex_id.get_interaction_dict():

            for interaction_type in complex_id.get_interaction_dict()[tup]:

                comparison_result = compare_interactions([structure[0][ref_chain_id], chain], interaction_type, complex_id.similar_sequences)

                if tup[0] == tup[1]:
                    if comparison_result == [True,False]:
                        complex_id.nodes[-1].add_interaction(complex_id[chain.get_id()], interaction_type)
                        complex_id[chain.get_id()].add_interaction(complex_id.nodes[-1], interaction_type[::-1])

                        break
                    elif comparison_result == [True,True]:
                        complex_id.nodes[-1].add_interaction(complex_id[chain.get_id()], interaction_type)
                        complex_id[chain.get_id()].add_interaction(complex_id.nodes[-1], interaction_type)

                        break


                elif True in comparison_result:

                    complex_id.nodes[-1].add_interaction(complex_id[chain.get_id()], interaction_type)
                    complex_id[chain.get_id()].add_interaction(complex_id.nodes[-1], interaction_type)

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

def superimpose_fun(str1, str2, node, i, complex_id, similar_seq, homodimer):
    """
    This functions superimposes 2 structure using nodes
    :param str1: fixed structure
    :param str2: mobile structure
    :param node: node in which we are working
    :param chain_str2: chain of str2 we want to superimpose
    :param complex_id: complex_id information
    :return: complex_id with new node if superimposition is feasible or a clash if it is not.
    """
    print('\nSuperimposing %s over chain %s' %(str2, node.get_chain()))
    chain1 = node.get_chain()
    str2_copy = copy.deepcopy(str2)
    node_chain_copy = copy.deepcopy(str1[0][chain1])
    chain_str2 = copy.deepcopy(i)

    trim_to_superimpose(node_chain_copy, chain_str2)
    atoms_chain1 = [atom for atom in list(node_chain_copy.get_atoms()) if atom.get_id() == 'CA']
    atoms_chain2 = [atom for atom in list(chain_str2.get_atoms())if atom.get_id() == 'CA']
    sup = Superimposer()

    sup.set_atoms(atoms_chain1, atoms_chain2)



    if not homodimer:
        other_chain2 = [x for x in str2_copy if x.get_id() != chain_str2.get_id()][0]
        other_chain2_original = [x for x in str2 if x.get_id() != chain_str2.get_id()][0]
    else:
        other_chain2 =str2_copy[1]
        other_chain2_original = str2[1]

    sup.apply(other_chain2)

    if not get_clash_chains(str1, other_chain2,chain1): ## returns T if there is a clash and F if there isn't.

        other_chain2.id = len(complex_id.get_nodes())+1
        similar_seq[other_chain2] = similar_seq[other_chain2_original]
        complex_id.add_node(other_chain2, node, str2)
        str1[0].add(other_chain2)

        interaction_finder(str1, other_chain2.get_id(), complex_id,node)
    
        return True
    else:
        return 'clash'


#########

def update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict):

    global branch_id
    branch_id.append(0)

    for node in complex_id.get_nodes():
        print('%s: %s' % (node.get_chain(),node))
        for interaction, value in node.get_interaction_dict().items():
            print("%s: %s " % (interaction, value))

    for other_CI in [ident for ident in complex_id_dict[len(complex_id.get_nodes())]]:

        if complex_id.compare_with(other_CI,4):
            print('Encontrada complex id repetida')

            branch_id.pop()
            print('\nvolviendo a la rama %s' % ".".join([str(x) for x in branch_id[:-1]]))

            return
    complex_id_dict[len(complex_id.get_nodes())].append(complex_id)

    for nodes in complex_id.get_nodes():
        for interact in [ interaction[0] for interaction in nodes.get_interaction_dict().items() if interaction[1] is None ]:

            branch_id[-1] += 1
            print("\nStarting new Branch: %s" % ".".join([str(x) for x in branch_id]))

            if branch_id == [1,1]:
                print('stop')


            for node in complex_id.get_nodes():
                print('%s: %s' % (node.get_chain(), node))
                for interaction, value in node.get_interaction_dict().items():
                    print("%s: %s " % (interaction, value))
            print('\n')
            print(list(base_struct.get_chains()))
            print('\n')

            complex_id_copy = complex_id.copy()
            copied_current_node = complex_id_copy[nodes.get_chain()]

            if similar_seq[interact[0]] == similar_seq[interact[1]]:
                # TODO len de complex id
                chain_str2_copy = copy_chain(interact[0], len(complex_id.get_nodes()))
                modified_str = superimpose_fun(base_struct, interact, copied_current_node, chain_str2_copy, complex_id_copy, similar_seq, True)

                if modified_str == 'clash':
                    nodes.add_interaction("clash", interact)
                    modified_str = None

                if modified_str:

                    if len(complex_id_copy.get_nodes()) not in complex_id_dict:
                        complex_id_dict[len(complex_id_copy.get_nodes())] = []

                    update_structure(base_struct, complex_id_copy, complex_id_dict, similar_seq, chains_str_dict)
                    print('Haciendo un pop')
                    complex_id_copy.pop_structure(base_struct)

            else:
                for i in interact:
                    if similar_seq[chains_str_dict[nodes.get_chain_type()]] == similar_seq[i]:

                        modified_str = superimpose_fun(base_struct, interact, copied_current_node, i, complex_id_copy, similar_seq, False)
                        if modified_str == 'clash':
                            nodes.add_interaction("clash", interact)
                            modified_str = None

                        if modified_str:

                            if len(complex_id_copy.get_nodes()) not in complex_id_dict:
                                complex_id_dict[len(complex_id_copy.get_nodes())] = []

                            update_structure(base_struct, complex_id_copy, complex_id_dict, similar_seq, chains_str_dict)
                            print('Haciendo un pop')
                            complex_id_copy.pop_structure(base_struct)
                        break
    # target_list = []
    # for node in complex_id.get_nodes():
    #     for target in node.get_interaction_dict().values():
    #         target_list.append(target)
    # if None not in target_list:
    write_to_pdb(base_struct)
    exit(0)
    branch_id.pop()
    print('\nvolviendo a la rama %s' % ".".join([str(x) for x in branch_id[:-1]]))

def macrocomplex_builder(id_dict, similar_seq, interaction_dict):
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

    if not os.path.exists('result'):
        os.makedirs('result')
    else:
        for the_file in os.listdir('result'):
            file_path = os.path.join('result', the_file)
            if os.path.isfile(file_path):
                os.unlink(file_path)

    global branch_id
    for chain in chains_str_dict:
        print("\nStarting new Branch: %s" % ".".join([str(x) for x in branch_id]))
        # initialize an empty structure
        base_struct = Structure.Structure('1')
        base_struct.add(Model.Model(0))
        # copy_chain = copy.deepcopy(Chain.Chain(chains_str_dict[chain]))
        chain_copied = copy_chain(chains_str_dict[chain], 1)
        # add chain to the new structure
        base_struct[0].add(chain_copied)
        similar_seq[chain_copied] = similar_seq[chains_str_dict[chain]]
        complex_id = ComplexId(interaction_dict, id_dict, similar_seq, base_struct)
        complex_id_dict[len(complex_id.get_nodes())] = []
        update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict)
        branch_id[-1] += 1

    print('hey')


if __name__ == '__main__':

    # get_all_interaction_pairs('')

    result = get_interaction_pairs_from_input('5vox_all_interactions')

    id_dict = result[1]
    interaction_dict = result[0]
    similar_sequences = result[2]

    macrocomplex_builder(id_dict, similar_sequences, interaction_dict)


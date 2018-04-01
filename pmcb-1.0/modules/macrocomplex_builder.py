from Bio.PDB import *
from Complex_breaker import *
from Complex_id import *
from ResidueDepth_copy import *
from Complex_breaker import trim_to_superimpose
from modeller_optimization import modeller_funcs

branch_id = [1]
pdb_counter = 1


def write_to_pdb(structure, directory):
    """
    This function writes a pdb file from a structure
    :param structure: structure we want to write the pdb file from.
    :return: writes a PDB file to the working directory
    """
    global pdb_counter

    final_structure = Structure.Structure(structure.id)
    model_counter = 0

    id_list = []

    for chain in structure.get_chains():

        if not id_list:
            id_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                       'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B',
                       'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
                       'V', 'W', 'X', 'Y', 'Z'][::-1]
            model_counter += 1
            final_structure.add(Model.Model(model_counter))

        new_chain = chain.copy()
        new_chain.id = id_list.pop()

        final_structure[model_counter].add(new_chain)

    io = PDBIO()
    io.set_structure(final_structure)
    file_name = 'result_' + directory + '/' + structure.id + str(pdb_counter) + '.pdb'
    io.save(file_name)
    pdb_counter += 1
    return file_name


def get_clash_chains(structure, chain, prev_chain, options):
    """
    This function recieves a hain and a structure and calculates if there is a clash between these.
    :param structure: structure where we want to add the chain
    :param chain: chain object we want to add to the structure
    :return: True or false, True if there is clash and false if there is no clash.
    """

    # Select the atoms we want to check if there are clashes between.
    chain_atoms = [x for x in list(chain.get_atoms()) if x.get_id() == 'CA' or x.get_id() == 'P']
    atom_list = list(structure.get_atoms())
    ns = NeighborSearch(atom_list)

    # Check if there are clashes between the selected atoms.
    for atom1 in chain_atoms:
        for atom in ns.search(atom1.get_coord(), 1.2, 'A'):
            print('%s(%s)-%s(%s)' % (atom1.get_id(), atom1.get_parent().get_id(), atom.get_id(), atom1.get_parent().get_id()))
            clashing_chain = atom.get_parent().get_parent().get_id()
            if clashing_chain != prev_chain:
                if options.verbose:
                    print('clash found')
                return True

    return False


#########

def interaction_finder(structure, ref_chain_id, complex_id, node, options):
    """
    This function receives a chain and a structure and calculates if there is a clash between these.
    :param structure: structure where we want to add the chain
    :param ref_chain_id: chain object we want to add to the structure
    :param complex_id:
    :param node:
    :param options: the arguments passed by the user
    :return: True or false, True if there is clash and false if there is no clash.
    """
    neighbor_chains = []
    ns = NeighborSearch(list(structure.get_atoms()))
    ref_chain = structure[0][ref_chain_id]
    for atom in [atom for atom in ref_chain.get_atoms() if atom.get_id() == 'CA' or atom.get_id() == 'P']:  # For every alpha carbon in chain
        for atom2 in ns.search(atom.get_coord(), 8, level='A'):
            if atom2.get_id() == 'CA' or atom2.get_id() == 'P':  # for every alpha carbon at 8 angstroms or less from atom
                chain2 = atom2.get_parent().get_parent()  # Getting to which chain it belongs
                if chain2 != ref_chain and chain2 not in neighbor_chains and chain2.get_id() != node.get_chain():
                    neighbor_chains.append(chain2)  # If it is not in the same chain and it is not already a
                    # key (we already have its interactions) we add the chain as a value
    if options.verbose:
        print('%s interactions found' % (len(neighbor_chains)))

    for chain in neighbor_chains:
        tup = sorted([complex_id.id_dict[complex_id.similar_sequences[chain]], complex_id.nodes[-1].get_chain_type()])
        tup = tuple(tup)

        if tup in complex_id.get_interaction_dict():

            for interaction_type in complex_id.get_interaction_dict()[tup]:

                comparison_result = compare_interactions([structure[0][ref_chain_id], chain], interaction_type,
                                                         complex_id.similar_sequences)

                if tup[0] == tup[1]:
                    if comparison_result == [True, False]:
                        complex_id.nodes[-1].add_interaction(complex_id[chain.get_id()], interaction_type)
                        complex_id[chain.get_id()].add_interaction(complex_id.nodes[-1], interaction_type[::-1])

                        break
                    elif comparison_result == [True, True]:
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

    new_chain = chain.copy()
    new_chain.id = id
    return new_chain


#########

def superimpose_fun(str1, str2, node, i, complex_id, similar_seq, homodimer, options):
    """
    This functions superimposes 2 structure using nodes
    :param str1: fixed structure
    :param str2: mobile structure
    :param node: node in which we are working
    :param chain_str2: chain of str2 we want to superimpose
    :param complex_id: complex_id information
    :return: complex_id with new node if superimposition is feasible or a clash if it is not.
    """
    if options.verbose:
        print('\nSuperimposing %s over chain %s' % (str2, node.get_chain()))
    chain1 = node.get_chain()
    str2_copy = [x.copy() for x in str2]
    node_chain_copy = str1[0][chain1].copy()
    chain_str2 = i.copy()

    # Trim the chains so that they have the same number of atoms to superimpose
    trim_to_superimpose(node_chain_copy, chain_str2)
    atoms_chain1 = [atom for atom in list(node_chain_copy.get_atoms()) if atom.get_id() == 'CA' or atom.get_id() == 'P']
    atoms_chain2 = [atom for atom in list(chain_str2.get_atoms()) if atom.get_id() == 'CA' or atom.get_id() == 'P']
    sup = Superimposer()

    # Superimpose the chains
    sup.set_atoms(atoms_chain1, atoms_chain2)

    # Select the chains that haven't been superimposed
    if not homodimer:
        other_chain2 = [x for x in str2_copy if x.get_id() != chain_str2.get_id()][0]
        other_chain2_original = [x for x in str2 if x.get_id() != chain_str2.get_id()][0]
    else:
        other_chain2 = str2_copy[1]
        other_chain2_original = str2[1]

    # Apply the rotation matrix to the chain we want to add to the complex
    sup.apply(other_chain2)

    # Assess if there is a clash between the chain we want to add and the others
    if not get_clash_chains(str1, other_chain2, chain1, options):

        # add the chain to the macrocomplex
        other_chain2.id = len(complex_id.get_nodes()) + 1
        similar_seq[other_chain2] = similar_seq[other_chain2_original]
        complex_id.add_node(other_chain2, node, str2)
        str1[0].add(other_chain2)

        if options.intensive:
            interaction_finder(str1, other_chain2.get_id(), complex_id, node, options)

        return True
    else:
        return 'clash'


#########

def update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict, stoichiometry_dict,
                     directory, options):
    """
    This is a recursive function which tries to add new chains until it is impossible to add any other due to clashes or
    to the specifications set by the user.
    :param base_struct: structure on to which we want to add new chains.
    :param complex_id: complex id of the base structure.
    :param complex_id_dict: dictionary of all the complex id set so far.
    :param similar_seq: dictionary relating sequences that are similar.
    :param chains_str_dict: dictionary with unique chains with chain name as key and str in the value
    :param stoichiometry_dict: XXXXXXX
    :param directory: name of the directory where the interactions files are.
    :param options: parameters passed by the user using the command line.
    :return: base structure updated with the new chain if it can add it.
    """
    global branch_id

    # Assess if the subunit limit has been reached and act accordingly
    if options.subunit_n or options.subunit_n == 0:
        if options.subunit_n > 0:
            options.subunit_n -= 1
        else:
            file_name = write_to_pdb(base_struct, directory)
            modeller_funcs(file_name, options)
            if options.subunit_n == 0:
                exit(0)
    branch_id.append(0)

    # Iterate for each chain the macro-complex for
    if options.verbose:
        for node in complex_id.get_nodes():
            print('%s: %s' % (node.get_chain(), node))
            for interaction, value in node.get_interaction_dict().items():
                    print("%s: %s " % (interaction, value))

    # Compare if we have already obtained this complex id
    for other_CI in [ident for ident in complex_id_dict[len(complex_id.get_nodes())]]:

        if complex_id.compare_with(other_CI, 4):
            if options.verbose:
                print('Repeated Complex id found')
            # Go back a branch since it has found that the ID has already been done before
            branch_id.pop()
            if options.verbose:
                print('\nReturning to branch %s' % ".".join([str(x) for x in branch_id[:-1]]))

            return

    complex_id_dict[len(complex_id.get_nodes())].append(complex_id)

    # iterate for the chains forming the complex at that point and the interactions that chain can do
    for nodes in complex_id.get_nodes():
        for interact in [interaction[0] for interaction in nodes.get_interaction_dict().items() if
                         interaction[1] is None]:

            branch_id[-1] += 1
            if options.verbose:
                print("\nStarting new Branch: %s" % ".".join([str(x) for x in branch_id]))

            for node in complex_id.get_nodes():
                if options.verbose:
                    print('%s: %s' % (node.get_chain(), node))
                for interaction, value in node.get_interaction_dict().items():
                    if options.verbose:
                        print("%s: %s " % (interaction, value))
            print('\n')
            print(list(base_struct.get_chains()))
            print('\n')

            complex_id_copy = complex_id.copy()
            copied_current_node = complex_id_copy[nodes.get_chain()]

            # if the sequences interacting are the same
            if similar_seq[interact[0]] == similar_seq[interact[1]]:
                if not options.st or (similar_seq[interact[0]] not in stoichiometry_dict or stoichiometry_dict[
                    similar_seq[interact[0]]]):
                    chain_str2_copy = copy_chain(interact[0], len(complex_id.get_nodes()))
                    modified_str = superimpose_fun(base_struct, interact, copied_current_node, chain_str2_copy,
                                                   complex_id_copy, similar_seq, True, options)

                    # update the complex_id with clash if it couldn't add the chain
                    if modified_str == 'clash':
                        nodes.add_interaction("clash", interact)
                        modified_str = None

                    #  update the complex_id with a new node if it was able to add the chain
                    if modified_str:
                        if options.st and similar_seq[interact[0]] in stoichiometry_dict:
                            stoichiometry_dict[similar_seq[interact[0]]] -= 1

                        if len(complex_id_copy.get_nodes()) not in complex_id_dict:
                            complex_id_dict[len(complex_id_copy.get_nodes())] = []

                        # Recursively recall the function for the new node created
                        update_structure(base_struct, complex_id_copy, complex_id_dict, similar_seq, chains_str_dict,
                                         stoichiometry_dict, directory, options)
                        if options.verbose:
                            print('Popping')
                        if options.subunit_n:
                            options.subunit_n += 1
                        complex_id_copy.pop_structure(base_struct)
                else:
                    nodes.add_interaction("full", interact)

            # Same as before but in case the chains are not the same
            else:

                for i in interact:
                    if similar_seq[chains_str_dict[nodes.get_chain_type()]] == similar_seq[i]:
                        other_chain = [x for x in interact if x != i][0]
                        if not options.st or (similar_seq[other_chain] not in stoichiometry_dict or stoichiometry_dict[
                            similar_seq[other_chain]]):

                            modified_str = superimpose_fun(base_struct, interact, copied_current_node, i,
                                                           complex_id_copy, similar_seq, False, options)
                            if modified_str == 'clash':
                                nodes.add_interaction("clash", interact)
                                modified_str = None

                            if modified_str:
                                if options.st and similar_seq[other_chain] in stoichiometry_dict:
                                    stoichiometry_dict[similar_seq[other_chain]] -= 1

                                if len(complex_id_copy.get_nodes()) not in complex_id_dict:
                                    complex_id_dict[len(complex_id_copy.get_nodes())] = []

                                update_structure(base_struct, complex_id_copy, complex_id_dict, similar_seq,
                                                 chains_str_dict, stoichiometry_dict, directory, options)
                                if options.verbose:
                                    print('Popping')
                                if options.subunit_n:
                                    options.subunit_n += 1
                                complex_id_copy.pop_structure(base_struct)
                            break
                        else:
                            nodes.add_interaction("full", interact)

    # If the complex can't accept more chains print the pdb file resulting and go back a branch.
    for nodes in complex_id.get_nodes():
        verify = False
        if None not in nodes.interaction_dict.values():
            verify = True
            file_name = write_to_pdb(base_struct, directory)
            modeller_funcs(file_name, options)
            if options.subunit_n == 0:
                exit(0)

        if not options.intensive and verify:
            exit(0)

    branch_id.pop()
    if options.verbose:
        print('\nReturning to branch %s' % ".".join([str(x) for x in branch_id[:-1]]))


def macrocomplex_builder(id_dict, similar_seq, interaction_dict, seq_dict, directory, options):
    """
    This function rebuilds a complex with the interactions we obtained from the pdb files.
    :param str_dict: dictionary with all the interactions we want to build the complex with.
    :param id_dict: dictionary with all the chains with their specific key.
    :return: returns the built macro-complex/es
    """
    global branch_id
    # # returns a set with all the chains passed by the user
    # chains = set([item for sublist in [x for x in str_dict] for item in sublist])

    # unique structure selection
    unique_chains = set(similar_seq.values())

    # Create a dictionary with unique chains with chain name as key and str in the value
    chains_str_dict = {v: k for k, v in id_dict.items() if k in unique_chains}

    # initialize a complex id dictionary
    complex_id_dict = {}

    if not os.path.exists('result_' + directory):
        os.makedirs('result_' + directory)
    else:
        for the_file in os.listdir('result_' + directory):
            file_path = os.path.join('result_' + directory, the_file)
            if os.path.isfile(file_path):
                os.unlink(file_path)
    stoichiometry_dict = {}

    # Handle stoichiometry, let the user pass the stoichiometry
    if options.st:
        chain_set = set(similar_seq.values())

        print('\nWe have found %s different proteins in your input. Would you like to set sotoickiometry values for '
              'any of them?\n Enter \'q\' for skipping the process' % len(chain_set))
        reverse_similar_seq = reverse_dictionary(similar_seq)

        chain_counter = 1
        for chain in chain_set:
            print('\nChain %s:\n' % chain_counter)
            name_str = ', '.join(reverse_similar_seq[chain])
            print('\tNames: %s\n\tSequence:\n\t%s' % (name_str, seq_dict[chain]))
            copy_number = input('\tIntroduce number of copies:')
            if copy_number and copy_number != 'q':
                copy_number = int(copy_number)
                stoichiometry_dict[chain] = copy_number
            chain_counter += 1

    if options.subunit_n:
        options.subunit_n -= 1
    elif options.subunit_n == 0:
        exit(0)

    for chain in chains_str_dict:

        if options.verbose:
            print("\nStarting new Branch: %s" % ".".join([str(x) for x in branch_id]))

        # initialize an empty structure
        base_struct = Structure.Structure(directory)
        base_struct.add(Model.Model(0))
        chain_copied = copy_chain(chains_str_dict[chain], 1)

        # add chain to the new structure
        base_struct[0].add(chain_copied)
        similar_seq[chain_copied] = similar_seq[chains_str_dict[chain]]

        if options.st and similar_seq[chain_copied] in stoichiometry_dict:
            stoichiometry_dict[similar_seq[chain_copied]] -= 1

        complex_id = ComplexId(interaction_dict, id_dict, similar_seq, base_struct)
        complex_id_dict[len(complex_id.get_nodes())] = []
        update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict, stoichiometry_dict,
                         directory, options)
        branch_id[-1] += 1


def reverse_dictionary(dictionary):
    """
    Takes a dictionary, returns a dictionary with values as keys and arrays of keys as values.
    :param dictionary: dictionary we want to reverse.
    :return reverse_dict: dictionary with values as keys and keys as values.
    """

    reverse_dict = {}
    for key, value in dictionary.items():

        if value not in reverse_dict:
            reverse_dict[value] = set([])

        reverse_dict[value].add(key.get_id())

    return reverse_dict

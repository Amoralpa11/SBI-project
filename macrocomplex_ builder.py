from Bio.PDB import *
from Complex_breaker import *
from Complex_id import *
from ResidueDepth_copy import *
from Complex_breaker import trim_to_superimpose
from modeller_optimization import structure_optimization
import argparse

branch_id = [1]
pdb_counter = 1


def write_to_pdb(structure):
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
            model_counter +=1
            final_structure.add(Model.Model(model_counter))

        new_chain = chain.copy()
        new_chain.id = id_list.pop()

        final_structure[model_counter].add(new_chain)

    io = PDBIO()
    io.set_structure(final_structure)
    file_name = 'result/' + structure.id + str(pdb_counter) + '.pdb'
    io.save(file_name)
    pdb_counter += 1
    return file_name


def get_clash_chains(structure, chain, prev_chain):
    """
    This function recieves a hain and a structure and calculates if there is a clash between these.
    :param structure: structure where we want to add the chain
    :param chain: chain object we want to add to the structure
    :return: True or false, True if there is clash and false if there is no clash.
    """
    global options
    chain_atoms = [x for x in list(chain.get_atoms()) if x.get_id() == 'CA' or x.get_id() == 'P']
    atom_list = list(structure.get_atoms())
    ns = NeighborSearch(atom_list)

    clash_counter = 0
    for atom1 in chain_atoms:
        atom_produces_clash = False
        for atom in ns.search(atom1.get_coord(), 1.2, 'A'):
            print('%s(%s)-%s(%s)' % (atom1.get_id(),atom1.get_parent().get_id(), atom.get_id(),atom1.get_parent().get_id()))
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
                if options.verbose:
                    print('More than 5 clashes found')
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
    for atom in [atom for atom in ref_chain.get_atoms() if atom.get_id() == 'CA' or atom.get_id() == 'P']:  # For every alpha carbon in chain
        for atom2 in ns.search(atom.get_coord(), 8, level='A'):
            if atom2.get_id() == 'CA' or atom2.get_id() == 'P':  # for every alpha carbon at 8 armstrongs or less from atom
                chain2 = atom2.get_parent().get_parent()  # Getting to which chain it belongs
                if chain2 != ref_chain and chain2 not in neighbor_chains and chain2.get_id() != node.get_chain():
                    neighbor_chains.append(chain2)  # If it is not in the same chain and it is not already a
                    # key (we already have its interactions) we add the chain as a value
    if options.verbose:
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
    global options
    if options.verbose:
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

def update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict, stoichiometry_dict):
    global options
    global branch_id

    # TODO subunit limit
    if options.subunit_n:
        if options.subunit_n > 0:
            options.subunit_n -= 1
        else:
            file_name = write_to_pdb(base_struct)
            if options.optimize:
                structure_optimization(file_name)
    branch_id.append(0)

    for node in complex_id.get_nodes():
        if options.verbose:
            print('%s: %s' % (node.get_chain(),node))
        for interaction, value in node.get_interaction_dict().items():
            if options.verbose:
                print("%s: %s " % (interaction, value))

    for other_CI in [ident for ident in complex_id_dict[len(complex_id.get_nodes())]]:

        if complex_id.compare_with(other_CI, 4):
            if options.verbose:
                print('Repeated Complex id found')

            branch_id.pop()
            if options.verbose:
                print('\nReturning to branch %s' % ".".join([str(x) for x in branch_id[:-1]]))

            return
    complex_id_dict[len(complex_id.get_nodes())].append(complex_id)

    for nodes in complex_id.get_nodes():
        for interact in [ interaction[0] for interaction in nodes.get_interaction_dict().items() if interaction[1] is None ]:

            branch_id[-1] += 1

            print("\nStarting new Branch: %s" % ".".join([str(x) for x in branch_id]))

            if nodes.get_chain() == 17:
                print('stop')


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

            if similar_seq[interact[0]] == similar_seq[interact[1]]:
                if not options.st or (similar_seq[interact[0]] not in stoichiometry_dict or stoichiometry_dict[similar_seq[interact[0]]]):
                    # TODO len de complex id
                    chain_str2_copy = copy_chain(interact[0], len(complex_id.get_nodes()))
                    modified_str = superimpose_fun(base_struct, interact, copied_current_node, chain_str2_copy, complex_id_copy, similar_seq, True)

                    if modified_str == 'clash':
                        nodes.add_interaction("clash", interact)
                        modified_str = None

                    if modified_str:
                        if options.st and similar_seq[interact[0]] in stoichiometry_dict:
                            stoichiometry_dict[similar_seq[interact[0]]] -=1

                        if len(complex_id_copy.get_nodes()) not in complex_id_dict:
                            complex_id_dict[len(complex_id_copy.get_nodes())] = []

                        update_structure(base_struct, complex_id_copy, complex_id_dict, similar_seq, chains_str_dict, stoichiometry_dict)
                        if options.verbose:
                            print('Popping')
                        if options.subunit_n:
                            options.subunit_n += 1
                        complex_id_copy.pop_structure(base_struct)
                else:
                    nodes.add_interaction("full", interact)
            else:

                for i in interact:
                    if similar_seq[chains_str_dict[nodes.get_chain_type()]] == similar_seq[i]:
                        other_chain = [x for x in interact if x != i][0]
                        if not options.st or (similar_seq[other_chain] not in stoichiometry_dict or stoichiometry_dict[similar_seq[other_chain]]):

                            modified_str = superimpose_fun(base_struct, interact, copied_current_node, i, complex_id_copy, similar_seq, False)
                            if modified_str == 'clash':
                                nodes.add_interaction("clash", interact)
                                modified_str = None

                            if modified_str:
                                if options.st and similar_seq[other_chain] in stoichiometry_dict:
                                    stoichiometry_dict[similar_seq[other_chain]] -= 1

                                if len(complex_id_copy.get_nodes()) not in complex_id_dict:
                                    complex_id_dict[len(complex_id_copy.get_nodes())] = []

                                update_structure(base_struct, complex_id_copy, complex_id_dict, similar_seq, chains_str_dict, stoichiometry_dict)
                                if options.verbose:
                                    print('Popping')
                                if options.subunit_n:
                                    options.subunit_n += 1
                                complex_id_copy.pop_structure(base_struct)
                            break
                        else:
                            nodes.add_interaction("full", interact)

    for nodes in complex_id.get_nodes():
        verify = False
        if None not in nodes.interaction_dict.values():
            verify = True
            file_name = write_to_pdb(base_struct)

        if options.optimize and verify:
            structure_optimization(file_name)

        if not options.intensive and verify:
            exit(0)

    branch_id.pop()
    if options.verbose:
        print('\nReturning to branch %s' % ".".join([str(x) for x in branch_id[:-1]]))

def macrocomplex_builder(id_dict, similar_seq, interaction_dict, seq_dict, directory):
    """
    This function rebuilds a complex with the interactions we obtained from the pdb files.
    :param str_dict: dictionary with all the interactions we want to build the complex with.
    :param id_dict: dictionary with all the chains with their specific key.
    :return: returns a XXXXXXX with the macrocomplex built and... the middle steps??
    """

    global options
    global branch_id

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
    stoichiometry_dict = {}
    if options.st:

        chain_set = set(similar_seq.values())

        print('\nWe have found %s different proteins in your input. Would you like to set sotoickiometry values for any of them?\n Enter \'q\' for skipping the process' % len(chain_set))
        reverse_similar_seq = reverse_dictionary(similar_seq)

        chain_counter = 1
        for chain in chain_set:
            print('\nChain %s:\n' % chain_counter)
            name_str = ', '.join(reverse_similar_seq[chain])
            print('\tNames: %s\n\tSequence:\n\t%s' % (name_str, seq_dict[chain]))
            copy_number = input('\tIntroduce number of copies:' )
            if copy_number and copy_number != 'q':
                copy_number = int(copy_number)
                stoichiometry_dict[chain] = copy_number
            chain_counter += 1

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
        update_structure(base_struct, complex_id, complex_id_dict, similar_seq, chains_str_dict, stoichiometry_dict)
        branch_id[-1] += 1


def reverse_dictionary(dictionary):

    """
    Takes a dictionary, returns a dictionary with values as keys and arrays of keys as values
    :param dictionary:
    :return dictionary:
    """

    reverse_dict = {}
    for key, value in dictionary.items():
        if value not in reverse_dict:
            reverse_dict[value] = set([])
        reverse_dict[value].add(key.get_id())

    return reverse_dict

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
                        action="store_false")

    parser.add_argument('-int', "--intensive",
                        dest="intensive",
                        help="Perform an intensive search or just return the 1st structure found.",
                        action="store_true")

    parser.add_argument('-st', "--stoichiometry",
                        dest="st",
                        help="Allows the user to pass the stoichiometry once the chains have been processed.",
                        action='store_true')

    parser.add_argument('-br', '--break',
                        dest="break_complex",
                        action='store',
                        default=None,
                        help='Indicate if you want to obtain all the pairwise interactions, "all", or just one of '
                             'each type, "unique". '
                        )

    options = parser.parse_args()


    class WrongArgumentBreak(Exception):
        pass

    # get_all_interaction_pairs('')
    # print('starting')

    if not options.break_complex:
        result = get_interaction_pairs_from_input(options.infile)
        id_dict = result[1]
        interaction_dict = result[0]
        similar_sequences = result[2]
        seq_dict = result[3]
        macrocomplex_builder(id_dict, similar_sequences, interaction_dict, seq_dict, options.infile)

    elif options.break_complex == "all":
        get_all_interaction_pairs(options.infile)

    elif options.break_complex == "unique":
        get_interaction_pairs(options.infile)

    else:
        raise WrongArgumentBreak('%s is not an accepted argument for -br, please pass "all" or "unique" if you are '
                                 'passing the -br argument')
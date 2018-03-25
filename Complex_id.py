class Node(object):

    """This class is the main component of the ComplexId class. It represents a single chain in the structure and
    holds information about all the interactions it is having and to which chains (other nodes) it is interacting
    with """

    def __init__(self, chain, complex_id, chain_type = None):

        """ When a Node is created we have to link it to a chain in the structure and to its parent complex id"""

        if type(chain) == type(int()):
            self.chain = chain
        else:
            self.chain = chain.get_id()

        if not chain_type:

            self.chain_type = complex_id.id_dict[complex_id.similar_sequences[chain]]

        else:
            self.chain_type = chain_type

        self.interaction_dict = complex_id.get_all_interactions_of_chain(self.chain_type)  # From the parent complex id we get all the interactions that this chain cand have


    def add_interaction(self, node, interaction):
        self.interaction_dict[tuple(interaction)] = node

    def get_chain_type(self):
        return self.chain_type

    def get_interaction_dict(self):
        return self.interaction_dict

    def get_chain(self):
        return self.chain

    def __set_interaction__dict(self,interaction_dict):
        self.interaction_dict = interaction_dict

    def get_deep_interactions(self, cycles, prev_node = None):

        """This function returns a dictionary with interactions (tuple of two chains) as keys. As values it has
        another dictionary from the other node involved in the interaction. The number of times this pattern repeats
        is determined by the parameter cycle. The dictionaries do not include the interaction that is its key in the
        dictionary they belong to """

        deep_interactions_dict = {}
        if cycles > 0:
            for interaction in self.interaction_dict:
                next_node = self.interaction_dict[tuple(interaction)]
                if next_node not in [prev_node,None,'clash']: # for all the interactions stablished
                    deep_interactions_dict[tuple(interaction)] = self.interaction_dict[tuple(interaction)].get_deep_interactions(cycles-1,self)  # We create the key for the dictionary and start a recursive funcitión calling again get_deep_interactions with the other chain involved in the interaction, with a cycle less and the current node as previous node

            if len(deep_interactions_dict) > 0:
                return deep_interactions_dict



    def copy_interactions_from(self,old_node, nodes_dict):

        for interaction, value in old_node.get_interaction_dict().items():
            if type(value) == type(self):
                self.add_interaction(nodes_dict[value],interaction)
            else:
                self.add_interaction(value,interaction)



def get_nodes_from_structure(complex_id, structure):


    """This function takes a complex id and an structure and adds a non connected node for every chain in the structure"""

    for chain in structure.get_chains():
        complex_id.add_node(chain)



class ComplexId(object):


    """This class holds the information of protein interactions in a structure. Its main component are the nodes, stored in the nodes list. It also have the information about the interactions that can have every chain in the complex, a dictionary relating similar sequences in the user input"""

    def __init__(self, interaction_dict, id_dict, similar_sequences, structure=None):

        self.id_dict = id_dict
        self.interaction_dict = interaction_dict
        self.similar_sequences = similar_sequences
        self.nodes = []
        if structure:
            get_nodes_from_structure(self, structure)

    def get_all_interactions_of_chain(self, chain_type):


        """Takes an string that is the identifier of the main chain of a kind and extracts all the interactions from the interaction_dict where it is involved this group of chains"""


        interaction_list = {}

        for pair in [pair for pair in self.interaction_dict if chain_type in pair]:
            for interaction in self.interaction_dict[pair]:
                interaction_list[tuple(interaction)] = None

        return interaction_list

    def add_node(self, chain, node=None, interaction=None):


        """takes a chain, an existing node and an interaction and creates a new node with one interaction poionting to the existing node and adds an interaction to the existing node pointing to the new node"""

        new_node = Node(chain, self)


        if interaction and node:
            if interaction[0] != interaction[1] or interaction[::-1] not in node.get_interaction_dict():
                new_node.add_interaction(node, interaction)
                node.add_interaction(new_node, interaction)
            else:

                new_node.add_interaction(node, interaction)
                node.add_interaction(new_node, interaction[::-1])



        self.nodes.append(new_node)



    def get_nodes(self):
        return self.nodes

    def get_id_dict(self):
        return self.id_dict

    def get_interaction_dict(self):
        return self.interaction_dict

    def get_chain_type_list(self):
        return sorted([node.get_chain_type() for node in self.get_nodes()])

    def get_interaction_id(self,cycles):


        """Takes a number of cycles and calls the function get_deep_interactions for every node in the complex id. returns a list of deep_interaction (dicts of dicts) of 'cycle' levels"""

        deep_interaction_list = []
        for node in self.get_nodes():
            deep_interaction_list.append(node.get_deep_interactions(cycles))

        return deep_interaction_list


    def copy(self):

        """

        Returns a copy of the complex id keeping the chain objects in the interactions dicts of the nodes pointing to the same reference to avoid missusing memory
        :return:
        """

        new_complex_id = ComplexId(self.interaction_dict,self.id_dict,self.similar_sequences)
        nodes_dict = {}
        new_node_list = []
        for node in self.get_nodes():
            nodes_dict[node] = Node(node.get_chain(),new_complex_id,node.get_chain_type())
            new_node_list.append(nodes_dict[node])

        for old_node in self.get_nodes():
            nodes_dict[old_node].copy_interactions_from(old_node,nodes_dict)

        new_complex_id.__set_nodes__(new_node_list)

        return new_complex_id

    def compare_with(self,complex2,cycles):


        """
        Takes one complex and a number of cycles and compares self with the complex provided. Returns true if complex
        are equeal and and false if not.

        the comparison is done at three levels, First it checks if the number of nodes is the same, if it is,
        it checks if there is the same number of chain_types among all the nodes. If not it compares the
        deep_interaction lists increasing the level of deepnes up to cycles if they are the same.

        :param complex2: ComplexId
        :param cycles: int
        :return: Bool
        """
        if len(self.get_nodes()) != len(complex2.get_nodes()):
            return False
        if self.get_chain_type_list() != complex2.get_chain_type_list():
            return False
        for cycle in range(1, cycles):

            interaction_list2 = complex2.get_interaction_id(cycle)

            for interaction_dict in self.get_interaction_id(cycle):
                if interaction_dict not in interaction_list2:
                    return False
                else:
                    interaction_list2.remove(interaction_dict)

        return True

    def pop_structure(self,structure):

        """
        This function removes from structure the chain that has the same id than the later node added to this complex_id
        :param structure:
        :return: None
        """
        structure[0].detach_child(self.nodes[-1].get_chain())

    def __set_nodes__(self,nodes):
        self.nodes = nodes


    def __getitem__(self,chain_id):
        for node in self.get_nodes():
            if node.get_chain() == chain_id:
                return node

if __name__ == '__main__':
    pass

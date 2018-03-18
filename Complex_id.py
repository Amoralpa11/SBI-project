from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model


class Node(object):

    def __init__(self, chain_type, chain, complex_id, pos):
        self.chain = chain
        self.chain_type = chain_type
        self.interaction_dict = complex_id.get_all_interactions_of_chain(self.chain_type)
        self.pos = pos

    def add_interaction(self, node, interaction):
        self.interaction_dict[tuple(interaction)] = node

    def get_chain_type(self):
        return self.chain_type

    def get_interaction_list(self):
        return self.interaction_dict

    def get_deep_interactions(self, cycles, prev_node = None):
        deep_interactions_dict = {}
        if cycles > 0:
            for interaction in self.interaction_dict:
                if prev_node != self.interaction_dict[tuple(interaction)] != None:
                    deep_interactions_dict[tuple(interaction)] = self.interaction_dict[tuple(interaction)].get_deep_interactions(cycles-1,self)

            if len(deep_interactions_dict) > 0:
                return deep_interactions_dict


def get_nodes_from_structure(complex_id, structure):
    for chain in structure.get_chains():
        complex_id.add_node(chain)



def pop_structure(structure, complex_id):
    pass


class ComplexId(object):

    def __init__(self, structure, interaction_dict, id_dict, similar_sequences):
        self.id_dict = id_dict
        self.interaction_dict = interaction_dict
        self.similar_sequences = similar_sequences
        self.nodes = []
        get_nodes_from_structure(self, structure)

    def get_all_interactions_of_chain(self, chain_type):

        interaction_list = {}

        for pair in [pair for pair in self.interaction_dict if chain_type in pair]:
            for interaction in self.interaction_dict[pair]:
                interaction_list[tuple(interaction)] = None

        return interaction_list

    def add_node(self, chain, node=None, interaction=None):

        chain_type = self.id_dict[self.similar_sequences[chain]]

        new_node = Node(chain_type, chain, self, len(self.nodes))

        if interaction and node:
            new_node.add_interaction(node, interaction)
            node.add_interaction(new_node, interaction)

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
        deep_interaction_list = []
        for node in self.get_nodes():
            deep_interaction_list.append(node. get_deep_interactions(cycles))

        return deep_interaction_list


    def __copy__(self):
        pass

    def compare_with(self,complex2,cycles):
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


if __name__ == '__main__':
    pass

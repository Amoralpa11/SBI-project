
class Node(object):

    def __init__(self,chain,complex_id,interactions = None):
        self.chain = chain.get_id()
        self.interactions = complex_id

    def add_interaction(self, node,interaction):

             self.interactions.append((node,interaction))



def get_nodes_from_structure(structure):



class ComplexId(object):

    def __init__(self, structure,pdb_interact_dict):
        self.nodes = get_nodes_from_structure(structure)

if __name__ == '__main__':

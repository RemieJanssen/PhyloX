import networkx as nx

class DiNetwork(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        edges = kwargs.get("edges", [])
        super().__init__(edges, *args, **kwargs)
        self.add_nodes_from(kwargs.get("nodes", []))
        for label in kwargs.get("labels", []):
            self.nodes[label[0]]["label"] = label[1]

    @property
    def leaves(self):
        return set([node for node in self.nodes if self.is_leaf(node)])

    @property
    def roots(self):
        return set([node for node in self.nodes if self.is_root(node)])

    @property
    def reticulation_number(self):
        return sum([max(self.in_degree(node) - 1, 0) for node in self.nodes])
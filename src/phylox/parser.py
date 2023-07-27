import re
import json
from phylox import DiNetwork
from phylox.base import find_unused_node


def extended_newick_to_dinetwork(newick, internal_labels=False, network=None):
    """
    Converts a Newick string to a networkx DAG with leaf labels.
    The newick string may or may not have length:bootstrap:probability annotations.
    The newick string may or may not have internal node labels.
    The newick string may or may not have hybrid nodes.

    :param newick: a string in extended Newick format for phylogenetic networks.
    :param internal_labels: a boolean, indicating whether the internal nodes of the network are labeled.
    :return: a phylogenetic network, i.e., a networkx digraph with leaf labels represented by the `label' node attribute.

    :example:
    >>> newick = "(A:1.1,B:1.2,(C:1.3,D:1.4)E:1.6)F;"
    >>> network = extended_newick_to_dinetwork(newick)
    >>> set(network.leaves) == {'A', 'B', 'C', 'D'}
    True
    >>> p = network.parent('A')
    >>> network[p]['A']['length']
    1.1
    """

    network = network or DiNetwork()
    nested_list = newick_to_nested_list(newick)
    print(nested_list)




def newick_to_nested_list(newick):
    nested_list = [""]
    print(newick)
    while newick:
        character = newick[0]
        newick = newick[1:]
        print()
        print(character)
        print(nested_list)
        if character == "(":
            newick, child = newick_to_nested_list(newick)
            nested_list[-1] = child
            nested_list+=[""]
        elif character == ")":
            return newick, nested_list
        elif character == ",":
            nested_list+=[""]
        elif character == ";":
            pass
        else:
            nested_list[-1]+=character
        print(nested_list)
    return nested_list



def extended_newick_to_dinetwork2(newick, internal_labels=False):
    """
    Converts a Newick string to a networkx DAG with leaf labels.
    The newick string may or may not have length:bootstrap:probability annotations.
    The newick string may or may not have internal node labels.
    The newick string may or may not have hybrid nodes.

    :param newick: a string in extended Newick format for phylogenetic networks.
    :param internal_labels: a boolean, indicating whether the internal nodes of the network are labeled.
    :return: a phylogenetic network, i.e., a networkx digraph with leaf labels represented by the `label' node attribute.

    :example:
    >>> newick = "(A:1.1,B:1.2,(C:1.3,D:1.4)E:1.6)F;"
    >>> network = extended_newick_to_dinetwork(newick)
    >>> set(network.leaves) == {'A', 'B', 'C', 'D'}
    True
    >>> p = network.parent('A')
    >>> network[p]['A']['length']
    1.1
    """
    if newick[-1] != ";":
        raise ValueError("Newick string does not end with ;")
    newick = newick[:-1]

    # convert the newick string to json string
    # add ' or " to the names if necessary
    newick = re.sub(r"#H([d]+)", r"#R\1", newick)

    # replace [ or ( with {[ and ] or ) with ]}
    newick = re.sub(r"[\[\(]", r'{"children": [', newick)
    newick = re.sub(r"[\]\)]", r"]}", newick)

    # replace :a:b:c values with dictionary {"length":a,"bootstrap":b,"probability":c}
    # for internal nodes
    newick = re.sub(
        r"}([\da-zA-Z\_]+):([\d\.]+):([\d\.]+):([\d\.]+)",
        r', "label": "\1","length": \2,"bootstrap": \3,"probability": \4}',
        newick,
    )
    # for leaves
    newick = re.sub(
        r"([\da-zA-Z\_]+):([\d\.]+):([\d\.]+):([\d\.]+)",
        r'{"label": "\1", "length": \2,"bootstrap": \3,"probability": \4}',
        newick,
    )
    # replace :a values with dictionary {"length":a}
    # for internal nodes
    newick = re.sub(
        r"}([\da-zA-Z\_]+):([\d\.]+)", r', "label": "\1","length": \2}', newick
    )
    # for leaves
    newick = re.sub(
        r"([\da-zA-Z\_]+):([\d\.]+)", r'{"label": "\1", "length": \2}', newick
    )

    # add " around leaf names
    while newick != re.sub(r"([\[\,])([#a-zA-Z\d]+)([\]\,])", r'\1"\2"\3', newick):
        newick = re.sub(
            r"([\[\,])([#a-zA-Z\d]+)([\]\,])", r'\1{"label": "\2"}\3', newick
        )

    # replace internal node name with dictionary {"label": name}
    newick = re.sub(r"}([#a-zA-Z\d]+)", r', "label": "\1"}', newick)

    # convert the json string to a network
    newick_json = json.loads(newick)
    network = json_to_dinetwork(newick_json)
    if not internal_labels:
        network = remove_internal_labels(network)
    return network


def remove_internal_labels(network):
    """
    Removes the internal labels from a network.

    :param network: a phylogenetic network.
    :return: a phylogenetic network with the internal labels removed.
    """
    for node in network.nodes:
        if not network.is_leaf(node):
            network.nodes[node].pop("label", None)
    return network


def json_to_dinetwork(newick_json, network=None, root_node=None):
    """
    Converts a json string to a phylox DiNetwork

    :param newick_json: a string in json format for phylogenetic networks.
    """
    network = network or DiNetwork()

    node = newick_json.get("label") or root_node or find_unused_node(network)
    network.add_node(node)
    if label := newick_json.get("label"):
        network.nodes[node]["label"] = label
    for child_dict in newick_json.get("children", []):
        child_dict_without_label_and_children = {
            k: v for k, v in child_dict.items() if k not in ("label", "children")
        }
        child = child_dict.get("label", find_unused_node(network))
        network.add_edge(node, child, **child_dict_without_label_and_children)
        json_to_dinetwork(child_dict, network, root_node=child)
    return network


# def Newick_To_Network(newick):
#     newick=newick[:-1]

#     #remove internal labels
#     newick = re.sub(r"I([\d]+)", "", newick)

#     #remove lengths
#     newick = re.sub(r"([\d]+)\.([\d]+)", "", newick)
#     newick = re.sub(r"E-[\d]+", "", newick)
#     newick = re.sub(r":", "", newick)

#     #make into list format
#     newick = newick.replace("(","[")
#     newick = newick.replace(")","]")
#     newick = re.sub(r"\]\#H([\d]+)", r",#R\1]", newick)
#     newick = re.sub(r"#([RH])([\d]+)", r"'#\1\2'", newick)


#     #add "" if necessary
#     newick = re.sub(r"([ABCD])", r"'\1'", newick)
#     newick = re.sub(r" ", "", newick)


#     nestedtree = ast.literal_eval(newick)
#     edges, leaves, label_set, current_node = NestedList_To_Tree(nestedtree,1)
#     edges.append([0,1])
#     ret_labels = dict()
#     leaf_labels = dict()
#     for l in leaves:
#         if len(l)>2 and (l[:2]=="#H" or l[:2]=="#R"):
#             ret_labels[l[2:]]=[]
#         else:
#             leaf_labels[l]=[]
#     for l in label_set:
#         if len(l[0])>2 and (l[0][:2]=="#H" or l[0][:2]=="#R"):
#             if l[0][1]=='H':
#                 ret_labels[l[0][2:]]+=[l[1]]
#             else:
#                 ret_labels[l[0][2:]]=[l[1]]+ret_labels[l[0][2:]]
#         else:
#             leaf_labels[l[0]]+=[l[1]]
#     network = nx.DiGraph()
#     network.add_edges_from(edges)
#     for retic in ret_labels:
#         r = ret_labels[retic]
#         receiving = r[0]
#         parent_receiving = 0
#         for p in network.predecessors(receiving):
#             parent_receiving = p
#         network.remove_node(receiving)
#         for v in r[1:]:
#             network.add_edge(v,parent_receiving)
#             network = nx.contracted_edge(network,(v,parent_receiving))
#             network.remove_edge(v,v)
#             parent_receiving = v
#     leaves = set()
#     for l in leaf_labels:
#          leaf_labels[l]=leaf_labels[l][0]
#          leaves.add(l)
#     return network, leaves, leaf_labels


# def NestedList_To_Tree(nestedList,next_node):
#     edges = []
#     leaves = set()
#     labels = []
#     top_node = next_node
#     current_node = next_node+1
#     for t in nestedList:
#         edges.append((top_node,current_node))
#         if type(t)==list:
#             extra_edges, extra_leaves, extra_labels, current_node = NestedList_To_Tree(t,current_node)
#         else:
#             extra_edges = []
#             extra_leaves = set([str(t)])
#             extra_labels = [[str(t), current_node]]
#             current_node+=1
#         edges = edges + extra_edges
#         leaves = leaves.union(extra_leaves)
#         labels = labels + extra_labels
#     return edges, leaves, labels, current_node


# network,_,leaf_labels = Newick_To_Network(newick_string)
# reverse_labels = {x:y for (y,x) in leaf_labels.items()}


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                         I/O Functions                     #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# ########
# ######## Convert Newick to a networkx Digraph with labels (and branch lengths)
# ########
# # Write length newick: convert ":" to "," and then evaluate as list of lists using ast.literal_eval
# # Then, in each list, the node is followed by the length of the incoming arc.
# # This only works as long as each branch has length and all internal nodes are labeled.
# def Newick_To_Network(newick):
#     """
#     Converts a Newick string to a networkx DAG with leaf labels.

#     :param newick: a string in extended Newick format for phylogenetic networks.
#     :return: a phylogenetic network, i.e., a networkx digraph with leaf labels represented by the `label' node attribute.
#     """
#     # Ignore the ';'
#     newick = newick[:-1]
#     # If names are not in string format between ', add these.
#     if not "'" in newick and not '"' in newick:
#         newick = re.sub(r"\)#H([\d]+)", r",#R\1)", newick)
#         newick = re.sub(r"([,\(])([#a-zA-Z\d]+)", r"\1'\2", newick)
#         newick = re.sub(r"([#a-zA-Z\d])([,\(\)])", r"\1'\2", newick)
#     else:
#         newick = re.sub(r"\)#H([d]+)", r"'#R\1'\)", newick)
#     newick = newick.replace("(", "[")
#     newick = newick.replace(")", "]")
#     nestedtree = ast.literal_eval(newick)
#     edges, current_node = NestedList_To_Network(nestedtree, 0, 1)
#     network = nx.DiGraph()
#     network.add_edges_from(edges)
#     network = NetworkLeafToLabel(network)
#     return network


# # Returns a network in which the leaves have the original name as label, and all nodes have integer names.
# def NetworkLeafToLabel(network):
#     """
#     Renames the network nodes to integers, while storing the original node names in the `original' node attribute.

#     :param network: a phylogenetic network
#     :return: a phylogenetic network with original node names in the `original' node attribute.
#     """
#     for node in network.nodes():
#         if network.out_degree(node) == 0:
#             network.node[node]['label'] = node
#     return nx.convert_node_labels_to_integers(network, first_label=0, label_attribute='original')


# # Auxiliary function to convert list of lists to graph
# def NestedList_To_Network(nestedList, top_node, next_node):
#     """
#     Auxiliary function used to convert list of lists to graph.

#     :param nestedList: a nested list.
#     :param top_node: an integer, the node name of the top node of the network represented by the list.
#     :param next_node: an integer, the lowest integer not yet used as node name in the network.
#     :return: a set of edges of the network represented by the nested list, and an updated next_node.
#     """
#     edges = []
#     if type(nestedList) == list:
#         if type(nestedList[-1]) == str and len(nestedList[-1]) > 2 and nestedList[-1][:2] == '#R':
#             retic_node = '#H' + nestedList[-1][2:]
#             bottom_node = retic_node
#         else:
#             bottom_node = next_node
#             next_node += 1
#         edges.append((top_node, bottom_node))
#         for t in nestedList:
#             extra_edges, next_node = NestedList_To_Network(t, bottom_node, next_node)
#             edges = edges + extra_edges
#     else:
#         if not (len(nestedList) > 2 and nestedList[:2] == '#R'):
#             edges = [(top_node, nestedList)]
#     return edges, next_node


# # Sets the labels of the nodes of a network with a given label dictionary
# def Set_Labels(network, label_dict):
#     """
#     Sets the labels of the leaves of a network using a dictionary of labels.

#     :param network: a networkx digraph, a DAG.
#     :param label_dict: a dictionary, containing the labels (values) of the nodes of the network (keys).
#     :return: a phylogenetic network, obtained by labeling network with the labels.
#     """
#     for node, value in label_dict.items():
#         network.node[node]['label'] = value


# ###############################################################################
# ###############################################################################
# ###############################################################################
# #######                                                           #############
# #######                     AAE CutTree CLASS                     #############
# #######                                                           #############
# ###############################################################################
# ###############################################################################
# ###############################################################################


# #A class that represents a network as a tree where hybrid edges have been cut at the hybrid nodes.
# #Used as an intermediate to find the Newick string of a network.
# class CutTree:
#     def __init__(self, network = None, current_node = None, leaf_labels= dict()):
#          self.hybrid_nodes = dict()
#          self.no_of_hybrids = 0
#          self.root = None
#          self.nw = deepcopy(network)
#          self.current_node = current_node
#          self.leaf_labels = leaf_labels
#          if not self.current_node:
#              self.current_node = 2*len(self.nw)
#          if network:
#              self.Find_Root()
#              network_nodes = list(self.nw.nodes)
#              for node in network_nodes:
#                  if self.nw.in_degree(node)>1:
#                      self.no_of_hybrids+=1
#                      enumerated_parents = list(enumerate(self.nw.predecessors(node)))
#                      for i,parent in enumerated_parents:
#                          if i==0:
#                              self.hybrid_nodes[node]=self.no_of_hybrids
#                          else:
#                              self.nw.add_edges_from([(parent,self.current_node,self.nw[parent][node])])
#                              self.nw.remove_edge(parent,node)
#                              self.hybrid_nodes[self.current_node] = self.no_of_hybrids
#                              self.current_node+=1
# #             self.CheckLabelSet()

#     #Returns the root node of the tree
#     def Find_Root(self):
#         for node in self.nw.nodes:
#             if self.nw.in_degree(node)==0:
#                 self.root = node
#                 return node

#     #Returns a newick string for the tree
#     def Newick(self,probabilities = False):
#         return self.Newick_Recursive(self.root,probabilities = probabilities)+";"

#     #Returns the newick string for the subtree with given root
#     #does not append the; at the end, for the full newick string of the tree, use Newick()
#     # auxiliary function for finding the newick string for the tree
#     def Newick_Recursive(self,root,probabilities = False):
#         if self.nw.out_degree(root)==0:
#             if root in self.hybrid_nodes:
#                 return "#H"+str(self.hybrid_nodes[root])
#             elif root in self.leaf_labels:
#                 return self.leaf_labels[root]
#             return str(root)
#         Newick = ""
#         for v in self.nw.successors(root):
#             Newick+= self.Newick_Recursive(v,probabilities)+":"+str(self.nw[root][v]['length'])
#             if probabilities and v in self.hybrid_nodes:
#                 Newick+="::"+str(self.nw[root][v]['prob'])
#             Newick+= ","
#         Newick = "("+Newick[:-1]+")"
#         if root in self.hybrid_nodes:
#             Newick += "#H"+str(self.hybrid_nodes[root])
#         return Newick

#     '''
#     def CheckLabelSet(self):
#         for v in self.nw.nodes:
#             if self.nw.out_degree(v)==0:
#                 if v not in self.leaf_labels and v not in self.hybrid_nodes:
#                     print("non-labelled leaf!")
#                     return False
#         return True
#     '''

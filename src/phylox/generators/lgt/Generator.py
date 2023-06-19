# #By Pons et al.

# import matplotlib
# import networkx as nx
# import pylab as plt
# import random,numpy
# import itertools
# import os
# import numpy as np
# from networkx.drawing.nx_agraph import graphviz_layout


# def last_node(net):
#     return max(net.nodes())


# def speciate(net,leaf):
#     l = last_node(net)
#     net.add_edge(leaf,l+1)
#     net.add_edge(leaf,l+2)


# def lgt(net,leaf1,leaf2):
#     net.add_edge(leaf1,leaf2,secondary=True)
#     l = last_node(net)
#     net.add_edge(leaf1,l+1)
#     net.add_edge(leaf2,l+2)


# def draw(net):
#     positions = nx.drawing.nx_agraph.graphviz_layout(net, prog="dot")
#     # nx.draw(net, pos=graphviz_layout(net),prog='dot',args="-Grankdir=LR", with_labels=True)
#     nx.draw(net,positions,with_labels=True)
#     plt.show()


# def leaves(net):
#     return [u for u in net.nodes() if net.out_degree(u)==0]


# def non_trivial_blobs(net):
#     blobs = list(nx.biconnected_components(nx.Graph(net)))
#     return [bl for bl in blobs if len(bl) > 2]


# def internal_blobs(net):
#     internal_nodes = set([u for u in net.nodes() if net.out_degree(u)>0])
#     blobs = list(nx.biconnected_components(nx.Graph(net)))
#     blobs = [bl for bl in blobs if len(bl) > 2]
#     nodes_in_blobs = set().union(*blobs)
#     nodes_not_in_blobs = internal_nodes - nodes_in_blobs
#     blobs.extend([set([u]) for u in nodes_not_in_blobs])
#     return blobs


# def compute_hash(net):
#     mapping_blobs = {}
#     blobs = internal_blobs(net)
#     for blob in blobs:
#         for node in blob:
#             mapping_blobs[node] = blob

#     mapping = {}
#     for l in leaves(net):
#         parent = list(net.predecessors(l))[0]
#         mapping[l] = mapping_blobs[parent]
#     return mapping


# def internal_and_external_pairs(net):
#     lvs = leaves(net)
#     pairs = [(l1,l2) for l1 in lvs for l2 in lvs if l1 != l2]
#     mapping = compute_hash(net)
#     internal_pairs = []
#     external_pairs = []
#     for pair in pairs:
#         if mapping[pair[0]] == mapping[pair[1]]:
#             internal_pairs.append(pair)
#         else:
#             external_pairs.append(pair)
#     return internal_pairs, external_pairs


# def random_leaf(net):
#     return random.choice(leaves(net))


# def random_pair(net,wint,wext):
#     int_pairs, ext_pairs = internal_and_external_pairs(net)
#     return random.choices(int_pairs+ext_pairs, weights=[wint]*len(int_pairs)+[wext]*len(ext_pairs))[0]


# def simulation_1(num_steps,prob_lgt,wint,wext):
#     net = nx.DiGraph()
#     net.add_edge(1,2)
#     net.add_edge(1,3)
#     for i in range(num_steps):
#         event = random.choices(['spec','lgt'],[1-prob_lgt, prob_lgt])[0]
#         #event = numpy.random.choice(['spec','lgt'],p=[1-prob_lgt, prob_lgt])
#         if event == 'spec':
#             l = random.choice(leaves(net))
#             speciate(net,l)
#         else:
#             pair = random_pair(net,wint,wext)
#             lgt(net,pair[0],pair[1])
#     return net


# def simulation_2(leaves_goal,prob_lgt,wint,wext):
#     net = nx.DiGraph()
#     net.add_edge(1,2)
#     net.add_edge(1,3)
#     leaves_no = 2
#     while leaves_no < leaves_goal:
#         event = random.choices(['spec','lgt'],[1-prob_lgt, prob_lgt])[0]
#         #event = numpy.random.choice(['spec','lgt'],p=[1-prob_lgt, prob_lgt])
#         if event == 'spec':
#             l = random.choice(leaves(net))
#             speciate(net,l)
#             leaves_no+=1
#         else:
#             pair = random_pair(net,wint,wext)
#             lgt(net,pair[0],pair[1])
#     return net


# def simulation_3(leaves_goal,retics_goal,wint,wext):
#     #pick a number of extant lineages for each LGT event independently
#     retics_at_lineage = dict()
#     for r in range(retics_goal):
#         lin = random.choice(range(2,leaves_goal+1))
#         if lin in retics_at_lineage:
#             retics_at_lineage[lin]+=1
#         else:
#             retics_at_lineage[lin]=1
#     net = nx.DiGraph()
#     net.add_edge(1,2)
#     net.add_edge(1,3)
#     if 2 in retics_at_lineage:
#         for j in range(retics_at_lineage[2]):
#             pair = random_pair(net,wint,wext)
#             lgt(net,pair[0],pair[1])
#     for i in range(3,leaves_goal+1):
#         l = random.choice(leaves(net))
#         speciate(net,l)
#         if i in retics_at_lineage:
#             for j in range(retics_at_lineage[i]):
#                 pair = random_pair(net,wint,wext)
#                 lgt(net,pair[0],pair[1])
#     return net


# def reticulations(G):
#     return [v for v in G.nodes() if G.in_degree(v)==2]

# def local_level(G,bicc):
#     rets=list(set(reticulations(G)).intersection(bicc)) # reticulations present in the blob
#     if len(rets)==0:
#         return 0
#     else:
#         bicc_edges=[e for e in G.edges() if ((e[0] in bicc)&(e[1]in bicc))]
#         end_nodes=[e[1] for e in bicc_edges]
#         return len([ret for ret in rets if ret in end_nodes])


# def OutputNetwork(network,out_type,leaves = None):
#     output = ""
#     if out_type=="el":
#         for i,e in enumerate(network.edges):
#             if i!=0:
#                 output+="\r\n"
#             output += str(e[0])+" "+str(e[1])
#     elif out_type=="pl":
#         #find number of leaves if necessary
#         if leaves == None:
#             leaves = 0
#             for v in network.nodes:
#                 if network.out_degree(v)==0:
#                     leaves+=1

#         #relabel nodes
#         labelDict           = dict()
#         leaves_current      = 1
#         internal_current    = leaves+1

#         for v in network.nodes:
#             if network.out_degree(v)==0:
#                 labelDict[v]=leaves_current
#                 leaves_current+=1
#             else:
#                 labelDict[v]=internal_current
#                 internal_current+=1

#         #generate parent lists
#         for i,v in enumerate(network.nodes):
#             if i!=0:
#                 output+="\r\n"
#             output += str(v)+";"
#             parentList=""
#             for p in network.predecessors(v):
#                 parentList+=str(labelDict[p])+","
#             output+=parentList[:-1]
#     return output


# ##################################
# ########     ALL DATA     ########
# ##################################

# """
# number_of_experiments = 500
# values_of_n = [100]
# values_of_alpha = [0.0,0.1,0.2,0.3,0.4]
# values_of_beta = [0.01,0.1,1,10,100] #weight of outside blob pair, inside blob weight is 1
# foldername = "./Dataset_n_100/"


# for out_type in ["el","pl"]:
#     this_path = foldername+out_type
#     if not os.path.exists(this_path):
#         os.makedirs(this_path)

# for (n, alpha, beta) in itertools.product(values_of_n, values_of_alpha, values_of_beta):
#     print(alpha,beta)
#     for i in range(number_of_experiments):
#         resG = simulation_2(n, alpha, 1, beta)
#         for out_type in ["el","pl"]:
#             filename = foldername+out_type+"/"+str(alpha)+"_"+str(beta)+"_"+str(i)
#             f=open(filename, 'w')
#             f.write(OutputNetwork(resG,out_type))
#             f.close()
# """


# ##################################
# ########   FIXED RETICS   ########
# ##################################


# number_of_experiments = 500
# values_of_n     = [50]#[10,15,20,25] # [100]
# values_of_k     = [1,2,3,4,5,10,20,50,100,200,500]#[6,7,8,9,10]#[1,2,3,4,5]
# values_of_n_k   = list(itertools.product(values_of_n, values_of_k))
# values_of_alpha = [float(k)/(n+k+1) for n,k in values_of_n_k]#Making sure the desired number of retics is the mode of the binom distribution that determines it
# values_of_n_k_and_alpha     = zip(values_of_n_k,values_of_alpha)
# values_of_beta = [1.0,0.1,0.01]#[0.01,0.1,1,10,100] #weight of outside blob pair, inside blob weight is 1
# def foldername(n):
#     return "./Dataset_n_"+str(n)+"_fixed_k/"

# for n in values_of_n:
#     for out_type in ["el","pl"]:
#         this_path = foldername(n)+out_type
#         if not os.path.exists(this_path):
#             os.makedirs(this_path)

# for (n_k_and_alpha, beta) in itertools.product(values_of_n_k_and_alpha, values_of_beta):
#     n_k,alpha=n_k_and_alpha
#     n,k = n_k
#     print(n,k,beta)
#     for i in range(number_of_experiments):
# #        print(k,beta,i)
#         resG=nx.DiGraph()
#         while len(reticulations(resG))!=k:
# # naive method:
# #            resG = simulation_1(n+k, alpha, 1, beta)
# # fast method:
#             resG = simulation_3(n  ,     k, 1, beta)
#         for out_type in ["el","pl"]:
#             filename = foldername(n)+out_type+"/"+str(k)+"_"+str(beta)+"_"+str(i)
#             f=open(filename, 'w')
#             f.write(OutputNetwork(resG,out_type))
#             f.close()


# #Compare realistic data
# """
# number_of_experiments = 500
# values_of_n_k       = [(4,5),(5,5),(6,2),(6,3),(6,4),(6,6),(7,4),(7,32),(8,3),(8,5),(9,3),(9,4),(9,5),(11,5),(13,1),(13,2),(13,6),(15,6),(16,4),(17,2),(21,6),(22,9),(24,2),(25,2),(28,5),(29,2),(39,11)]
# values_of_alpha = [float(i[1])/(i[0]+i[1]+1) for i in values_of_n_k]#Making sure the desired number of retics is the mode of the binom distribution that determines it
# values_of_n_k_and_alpha     = zip(values_of_n_k,values_of_alpha)
# values_of_beta = [0.01,0.1,1] #weight of outside blob pair, inside blob weight is 1
# foldername = "./Dataset_Realistic_Comparison/"


# for out_type in ["el","pl"]:
#     this_path = foldername+out_type
#     if not os.path.exists(this_path):
#         os.makedirs(this_path)

# for (n_k_alpha, beta) in itertools.product(values_of_n_k_and_alpha, values_of_beta):
#     n_k,alpha=n_k_alpha
#     n,k=n_k
#     print(n,k,alpha,beta)
#     for i in range(number_of_experiments):
# #        print(k,beta,i)
#         resG=nx.DiGraph()
#         while len(reticulations(resG))!=k:
# # naive method:
# #            resG = simulation_1(n+k, alpha, 1, beta)
# # fast method:
#             resG = simulation_3(n  ,     k, 1, beta)
#         for out_type in ["el","pl"]:
#             filename = foldername+out_type+"/"+str(n)+"_"+str(k)+"_"+str(beta)+"_"+str(i)
#             f=open(filename, 'w')
#             f.write(OutputNetwork(resG,out_type))
#             f.close()
# """


# ##################################
# ########## PLOTTING STUFF ########
# ##################################


# # In[33]:

# """
# f = plt.figure()
# x_positions=range(len(values_of_beta))
# x_axis = values_of_beta
# symbols = "ox*+"
# i = 0
# dibs = [None]*4
# legends = [None]*4
# for n in values_of_n:
#     for alpha in values_of_alpha:
#         to_plot = [stats_numblobs[(n, alpha, beta)] for beta in values_of_beta]
#         dibs[i]=plt.scatter(x_positions, to_plot, marker = symbols[i], c='black')
#         legends[i] = "n=%d, alpha=%0.1f" % (n, alpha)
#         i += 1
# plt.xticks(x_positions, x_axis)
# plt.legend(dibs, legends)
# plt.xlabel('beta')
# plt.ylabel('avg number of nontrivial blobs')
# plt.show()
# f.savefig("fig10.pdf", bbox_inches='tight')


# # In[34]:


# f = plt.figure()
# x_positions=range(len(values_of_beta))
# x_axis = values_of_beta
# symbols = "ox*+"
# i = 0
# dibs = [None]*4
# legends = [None]*4
# for n in values_of_n:
#     for alpha in values_of_alpha:
#         to_plot = [stats_level[(n, alpha, beta)] for beta in values_of_beta]
#         dibs[i]=plt.scatter(x_positions, to_plot, marker = symbols[i], c='black')
#         legends[i] = "n=%d, alpha=%0.1f" % (n, alpha)
#         i += 1
# plt.xticks(x_positions, x_axis)
# plt.legend(dibs, legends)
# plt.xlabel('beta')
# plt.ylabel('avg level')
# plt.show()
# f.savefig("fig11.pdf", bbox_inches='tight')


# # In[35]:

# number_of_experiments2 = 500
# values_of_n2 = [10,20,30]
# values_of_alpha2 = [0,0.1,0.2,0.3,0.4]
# values_of_beta2 = [0.01, 1, 10]
# stats_level2 = {}
# stats_numblobs2 = {}

# for (n, alpha, beta) in itertools.product(values_of_n2, values_of_alpha2, values_of_beta2):
#     levels = []
#     numblobs = []
#     for _ in range(number_of_experiments):
#         resG = simulation(n, alpha, 1, beta)
#         bic_comp=list(nx.biconnected_components(nx.Graph(resG)))
#         rets_x_bicc=[local_level(resG,b) for b in bic_comp]
#         level=max(rets_x_bicc)
#         levels.append(level)
#         sparse=[1 for x in rets_x_bicc if x!=0]
#         numblobs.append(sum(sparse))
#     stats_level2[(n, alpha, beta)] = float(sum(levels))/len(levels)
#     stats_numblobs2[(n, alpha, beta)] = float(sum(numblobs))/len(numblobs)


# # In[ ]:


# f = plt.figure()
# x_positions=range(len(values_of_alpha2))
# x_axis = values_of_alpha2
# symbols = "ox*+"
# i = 0
# dibs = [None]*4
# legends = [None]*4
# for n in [10,30]:
#     for beta in [0.01,10]:
#         to_plot = [stats_level2[(n, alpha, beta)] for alpha in values_of_alpha2]
#         dibs[i]=plt.scatter(x_positions, to_plot, marker = symbols[i], c='black')
#         legends[i] = "n=%d, beta=%0.2f" % (n, beta)
#         i += 1
# plt.xticks(x_positions, x_axis)
# plt.legend(dibs, legends)
# plt.xlabel('alpha')
# plt.ylabel('avg level')
# plt.show()
# f.savefig("fig8.pdf", bbox_inches='tight')


# # In[ ]:


# f = plt.figure()
# x_positions=range(len(values_of_alpha2))
# x_axis = values_of_alpha2
# symbols = "ox*+"
# i = 0
# dibs = [None]*4
# legends = [None]*4
# for n in [10,30]:
#     for beta in [0.01,10]:
#         to_plot = [stats_numblobs2[(n, alpha, beta)] for alpha in values_of_alpha2]
#         dibs[i]=plt.scatter(x_positions, to_plot, marker = symbols[i], c='black')
#         legends[i] = "n=%d, beta=%0.2f" % (n, beta)
#         i += 1
# plt.xticks(x_positions, x_axis)
# plt.legend(dibs, legends)
# plt.xlabel('alpha')
# plt.ylabel('avg number of nontrivial blobs')
# plt.show()
# f.savefig("fig9.pdf", bbox_inches='tight')
# """

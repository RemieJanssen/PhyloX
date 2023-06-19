import ast
import itertools
import random
import time

import networkx as nx

import matplotlib.pyplot as plt

random.seed(a=1234321)


def find_cherry(N, x):
    lst = list()
    for p in N.predecessors(x):
        if N.in_degree(p) == 1:
            for pc in N.successors(p):
                if pc != x:
                    t = N.out_degree(pc)
                    if t == 0:
                        lst.append((pc, x))
                    if t == 1:
                        for pcc in N.successors(pc):
                            if N.out_degree(pcc) == 0:
                                lst.append((pcc, x))
    return lst


def find_ret_cherry(N, x):
    lst = list()
    for p in N.predecessors(x):
        if N.out_degree(p) == 1:
            for pp in N.predecessors(p):
                for ppc in N.successors(pp):
                    if ppc != p:
                        if N.out_degree(ppc) == 0:
                            lst.append((x, ppc))
    return lst


def check_cherry(N, x, y):
    if N.has_node(x):
        if N.has_node(y):
            for px in N.predecessors(x):
                for py in N.predecessors(y):
                    if px == py:
                        return 1
                    if N.out_degree(px) == 1:
                        if px in N.successors(py):
                            return 2
    return False


def reduce_pair(N, x, y):
    k = check_cherry(N, x, y)
    if k == 1:
        for px in N.predecessors(x):
            N.remove_node(x)
            for ppx in N.predecessors(px):
                N.remove_node(px)
                N.add_edge(ppx, y)
            return True
    if k == 2:
        for px in N.predecessors(x):
            for py in N.predecessors(y):
                N.remove_edge(py, px)
                if N.in_degree(px) == 1:
                    for ppx in N.predecessors(px):
                        N.add_edge(ppx, x)
                        N.remove_node(px)
                for ppy in N.predecessors(py):
                    N.add_edge(ppy, y)
                    N.remove_node(py)
                return True
    return False


def find_tcs(N):
    lst1 = list()
    for x in N.nodes():
        if N.out_degree(x) == 0:
            cherry1 = find_cherry(N, x)
            lst1.extend(cherry1)
    lst2 = list()
    while lst1:
        cherry = lst1.pop()
        k = check_cherry(N, *cherry)
        if (k == 1) or (k == 2):
            reduce_pair(N, *cherry)
            lst2.append(cherry)
            lst1.extend(find_cherry(N, cherry[1]))
            lst1.extend(find_ret_cherry(N, cherry[1]))
    return lst2


def cps_reduces_network(N, lst):
    for cherry in lst:
        reduce_pair(N, *cherry)
    if N.size() == 1:
        return True
    return False


def tcn_contains(N, M):
    return cps_reduces_network(M, find_tcs(N))


def tester(foldername):
    maxLeavesRetics = 1000
    stepLeavesRetics = 25
    reps = 4
    folder_name = foldername
    f = open("./" + folder_name + "/data.txt", "w+")
    f.write(
        "leaves;reticulations;reticulations_subnetwork;repetition;subnetwork;running_time\n"
    )
    f.close()

    allIndices = list(
        itertools.product(
            range(stepLeavesRetics, maxLeavesRetics + 1, stepLeavesRetics),
            range(stepLeavesRetics, maxLeavesRetics + 1, stepLeavesRetics),
            range(stepLeavesRetics, maxLeavesRetics + 1, stepLeavesRetics),
            range(reps),
        )
    )
    print(list(allIndices)[0])
    random.shuffle(allIndices)
    print(list(allIndices)[0])
    for i in allIndices:
        if i[2] <= i[1]:
            leaves, reticulations, reticulationsSubnetwork, rep = i
            index1 = "0000000" + str(leaves)
            index1 = index1[-4:]
            index2 = "0000000" + str(reticulations)
            index2 = index2[-4:]
            index3 = "0000000" + str(reticulationsSubnetwork)
            index3 = index3[-4:]
            index = (
                "n=" + index1 + "_k=" + index2 + "_kSub=" + index3 + "_rep=" + str(rep)
            )
            name = index + ".txt"

            test = open("./" + folder_name + "/" + name, "r")
            line1 = test.read()
            test.close()
            line1 = line1.split("\n")
            M = nx.DiGraph()
            N = nx.DiGraph()
            N.add_edges_from(ast.literal_eval(line1[0]))
            M.add_edges_from(ast.literal_eval(line1[1]))

            start = time.time()
            contains = tcn_contains(N, M)
            end = time.time()
            runningTime = end - start
            f = open("./" + folder_name + "/data.txt", "a+")
            f.write(
                index1
                + " ; "
                + index2
                + " ; "
                + index3
                + " ; "
                + str(rep)
                + " ; "
                + str(contains)
                + " ; "
                + str(runningTime)
                + "\n"
            )
            f.close()


# path of the folder with all input files
foldername = "tests"
# use network containment on all input files and store the data in data.txt in the same folder
tester(foldername)

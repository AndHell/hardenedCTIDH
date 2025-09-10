#!/usr/bin/env python3

import sys
import numpy as np
# import matplotlib.pyplot as plt
# import networkx as nx
from math import floor

from memoized import memoized
import costisog


bitlength = lambda x: len(bin(x)[2:])                       # number of bits
hamming_weight = lambda x: bin(x).count("1");               # hamming weight: number of bits equal 1

sign = lambda x: (1, -1)[x < 0]                             # Sign of an integer
isequal = { True : 1 , False : 0 }                          # Simulating constant-time integer comparison

SQR = 1
ADD = 0

measure = lambda x: (x[0] + SQR * x[1] + ADD * x[2])




def dac_search(target,r0,r1,r2,chain,chainlen,best,bestlen):
  if chainlen >= bestlen: 
    return best,bestlen
  if r2 > target: return best,bestlen
  if r2<<(bestlen-1-chainlen) < target: return best,bestlen
  if (r2 == target): 
    return chain,chainlen
  chain *= 2
  chainlen += 1
  best,bestlen = dac_search(target,r0,r2,r0+r2,chain+1,chainlen,best,bestlen)
  best,bestlen = dac_search(target,r1,r2,r1+r2,chain,chainlen,best,bestlen)
  return best,bestlen


@memoized
def daclen(target):
  best = None
  bestlen = -1
  while best == None:
    bestlen += 1
    best, bestlen = dac_search(target,1,2,3,0,0,best,bestlen)
  return bestlen


@memoized
def eval_cost(l):
    bg = bsgs(l)
    w0 = costisog.isog(l, 0,bg)
    w1 = costisog.isog(l, 1,bg)
    return w1 - w0

@memoized
def isog_cost(l):
    bg = bsgs(l)
    wo = costisog.isog(l, 0,bg)
    return wo


@memoized
def bsgs(prime):
  return costisog.optimize(prime,1)[1]

'''
    dynamic_programming_algorithm():
    inputs: the list of small odd primes to be processed and its length
    output: the optimal strategy and its cost of the input list of small odd primes
'''
def dynamic_programming_algorithm(L, C_dac, C_isog, C_eval):

    # If the approach uses dummy operations, to set DUMMY = 2.0;
    # otherwise, to set DUMMY = 1.0 (dummy free approach);

    daclens = [daclen(e) for e in L]

    n = len(L)

    S = { 1: {} }    # Initialization of each strategy
    C = { 1: {} }    # Initialization of the costs: 0.

    # cMUL  = lambda l: costisog.dac(daclen(l))
    # cEVAL = lambda l: eval_cost(l)
    # cISOG = lambda l: isog_cost(l)

    # cMUL  = lambda l: np.array( [ 4.0 * (daclens[L.index(l)] + 2), 2.0 * (daclens[L.index(l)] + 2), 6.0 * (daclens[L.index(l)] + 2) - 2.0] )
    # cEVAL = lambda l: np.array( [2.0*(l - 1.0), 2.0, (l + 1.0)] )
    # cISOG = lambda l: np.array( [(3.0*l + 2.0*hamming_weight(l) - 9.0 + isequal[l == 3]*4.0), (l + 2.0*bitlength(l) + 1.0 + isequal[l == 3]*2.0), (3.0*l - 7.0 + isequal[l == 3]*6.0)] ) 

    #C_xMUL  = list(map(cMUL,  L_dac))   # list of the costs of each [l]P
    #C_xEVAL = list(map(cEVAL, L))   # list of the costs of each degree-l isogeny evaluation
    #C_xISOG = list(map(cISOG, L))   # list of the costs of each degree-l isogeny construction


    C_xMUL = C_dac
    C_xISOG = C_isog
    C_xEVAL = C_eval
    # print(C_xMUL)
    # real cost:
    # C_xISOG = [154, 170, 178, 194, 218, 242, 250, 280, 296, 304, 536, 560, 568, 544, 552, 570, 594, 610, 634, 638, 824, 811, 786, 816, 840, 886, 872, 859, 875, 883, 1069, 1095, 1119, 1154, 1159, 1092, 1124, 1134, 1166, 1321, 1337, 1369, 1326, 1344, 1392, 1311, 1377, 1496, 1441, 1465, 1437, 1485, 1522, 1530, 1546, 1567, 1594, 1599, 1647, 1535, 1683, 1699, 1731, 1715, 1771, 1811, 1816, 1899, 1969, 2049, 2065, 2097, 2025, 2057, 2073, 2097, 2121, 2363, 2379, 2387, 2433, 2427, 2667, 2675, 2567,]
    # C_xISOG.reverse()
    # C_xEVAL = [74, 82, 86, 94, 106, 118, 122, 134, 142, 146, 262, 274, 278, 303, 307, 319, 331, 339, 351, 354, 482, 454, 454, 466, 478, 506, 494, 502, 510, 514, 626, 642, 654, 670, 674, 665, 681, 666, 682, 770, 778, 794, 788, 800, 824, 785, 821, 876, 853, 865, 851, 875, 908, 912, 920, 919, 944, 935, 959, 919, 1030, 1038, 1054, 1031, 1059, 1079, 1091, 1138, 1194, 1234, 1242, 1258, 1219, 1235, 1243, 1255, 1267, 1478, 1486, 1490, 1485, 1510, 1653, 1657, 1606,]
    # C_xEVAL.reverse()

    
    for i in range(n):
            S[1][tuple([L[i]])] = [];           # Strategy with a list with only one element (a small odd prime number l_i)
            C[1][tuple([L[i]])] = C_xISOG[i];   # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]

    for i in range(2, n+1):
        C[i] = { }
        S[i] = { }

    # Assuming #L = n, we proceed.
    get_neighboring_sets = lambda L, k: [ tuple(L[i:i+k]) for i in range(n-k+1)] # This function computes all the k-tuple: (l_1, l_2, ..., l_{k)),
                                                                                 # (l_2, l_3, ..., l_{k+1)), ..., (l_{n-k}, l_{n-k+1, ..., l_{n}).
                                                      
    for i in range(2, n+1):
        for Tuple in get_neighboring_sets(L, i):
            if C[i].get(Tuple) is None:
                alpha = [ (b, 
                    C[len(Tuple[:b])][Tuple[:b]] +       # Subtriangle on the right side with b leaves
                    C[len(Tuple[b:])][Tuple[b:]] +       # Subtriangle on the left side with (i - b) leaves
                    2.0*sum([ C_xMUL[L.index(t)] for t in Tuple[:b] ])   +   # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                    2.0*sum([ C_xEVAL[L.index(t)] for t in Tuple[b:] ])  +   # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                    1.0*sum([ C_xMUL[L.index(t)] for t in Tuple[b:] ])       # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                  # ^  2 = 1 (point on twist) + 1 (to account for dummies)
                    ) for b in range(1, i - 1) 
                    ] +\
                    [ (i - 1,
                    C[i - 1][Tuple[:(i - 1)]] +       # Subtriangle on the right side with (i - 1) leaves
                    C[1][Tuple[(i - 1):]]     +       # Subtriangle on the left side with 1 leaf (only one vertex)
                    1.0*sum([ C_xMUL[L.index(t)] for t in Tuple[:(i - 1)] ]) + # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with 1 leaf
                    2.0*C_xEVAL[L.index(Tuple[i - 1])] +                      # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                    1.0*C_xMUL[L.index(Tuple[i - 1])]                    # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                  # ^  2 = 1 (point on twist) + 1 (to account for dummies)
                    )
                    ]

                b, C[i][Tuple] = min(alpha,  key=lambda t:  t[1])     # We save the minimal cost corresponding to the triangle with leaves Tuple
                S[i][Tuple] = [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]    # We save the optimal strategy corresponding to the triangle with leaves Tuple
               
                
    return S[n][tuple(L)], C[n][tuple(L)]   # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to C_xISOG[global_L.index(L[0])].


def DRT(n):

	vertexes = dict()      # list of the position of each node
	vertex_colors = []    # color of each node: red for the leaves, otherwise color is set to white
	acc = 0

	# Different shape of the isogeny graph
	for i in range(n):
		for j in range(n - 1 - i):
			vertex_colors.append('black')
			vertexes[acc] = (i, -j)
			acc += 1

	
	return vertexes, vertex_colors, [], []



def strategy_evaluation(strategy, n):
    vertexes = dict()      # list of the position of each node
    edges    = []          # edges of the isogeny triangle

    edge_colors   = []    # color of each edge: blue for scalar multiplications, and orange for isogeny evalutions
    vertex_colors = []    # color of each node: red for the leaves, otherwise color is set to white

    vertexes[0] = (0.0,0.0)	# Root of the strategy
    ramifications = [] 	# nodes having two different types of edges
    moves = [0] 		# 
    k = 0       		# current element of the strategy
    t = 0				# current vertex

    # Strategy evaluation starts next
    ramifications.append([0, vertexes[0]])							# The root is always a ramification
    for i in range(len(strategy)):
        
        # Getting the new possible ramifications
        while sum(moves) < (n - 1 - i):
            vertex = ramifications[-1][1]
                
            moves.append(strategy[k])					# Increasing moves
            k += 1                                  	# Next element of the strategy
            t += 1										# Next vertex
            edges.append( (t - 1, t) )					# New edge to be added
            edge_colors.append('tab:blue')					# Color of this type of edge is always blue

            vertexes[t] = (i, vertex[1] - strategy[k-1])	# New vertex to be added
            ramifications.append([t, vertexes[t]])			# New ramification
            vertex_colors.append('black')
        
        # Removing the ramifications (it is not required more!)
        ramifications.pop()
        moves.pop()
        
        # Updating the ramifications
        for j in range(len(ramifications)):
            t += 1
            vertexes[t] = ( ramifications[j][1][0] + 1.0, ramifications[j][1][1] )
            edges.append( (ramifications[j][0], t))

            ramifications[j] = [t, vertexes[t]]
            edge_colors.append('tab:red')

    return vertexes, vertex_colors, edges, edge_colors
    
def plot_strat(S):
    pass
    # n = len(S) + 1

    # # Strategy written as a graph
    # vertexes, vertex_colors, edges, edge_colors = strategy_evaluation(S, n)

    # # Simba method written as a graph
    # # vertexes, vertex_colors, edges, edge_colors = simba(n, 3)

    # # All the Discrete Right Triangle
    # # vertexes, vertex_colors, edges, edge_colors = DRT(n)
    # G = nx.Graph()

    # # Adding nodes in specific positions
    # G.add_nodes_from(list(range(len(vertexes))))

    # nx.set_node_attributes(G, vertexes, 'pos')
    # # Adding edges with specific colors
    # for i in range(len(edges)):
    #     G.add_edge(edges[i][0], edges[i][1], color=edge_colors[i])

    # # Setting variables for a pretty plot of the graph
    # edges = G.edges()
    # edge_colors = [G[u][v]['color'] for u, v in edges]
    # weights = [6 for u, v in edges]
    # vertex_sizes = [24] * len(vertexes)

    # # Finally, the graph will be plotted
    # # plt.figure(1)
    # # nx.draw(
    # #     G,
    # #     vertexes,
    # #     node_color=['black'] * len(vertexes),
    # #     node_size=vertex_sizes,
    # #     edge_color=edge_colors,
    # #     width=weights,
    # # )

    # # # Saving the graph as a .PNG figure
    # # # file_name = 'out.png'
    # # # plt.savefig(file_name)
    # # # print("saving graph: " + file_name)
    # # plt.show()
    # # #plt.close()



# ells = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 
#         83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 
#         173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 
#         269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 
#         373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
#         467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 
#         593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 
#         691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 
#         821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 
#         937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033,
#         1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 
#         1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 
#         1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 
#         1327, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439,]

# ells_m7 = [ 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 7591, 8963, 10391,]

# batch_stop = [19, 39, 59, 76, 95, 111, 128, 148, 168, 184,]
# batch_keys = [10, 10, 10, 8, 9, 8, 8, 10, 9, 7,]

# L = []
# L_dac = []
# for b, k in zip(batch_stop, batch_keys):
#     batch = ells[b-k:b]
#     daccost = max([costisog.dac(daclen(ell)) for ell in batch])
#     L_dac += [daccost]*k
#     L += [ells[i] for i in range(b-k+1, b+1)]
          

# L_dac.reverse()
# L.reverse()
#print(L_dac)

# S, cost = dynamic_programming_algorithm(L, L_dac)
# print("Optimal strat:" ,int(cost))
# print(S)
# plot_strat(S)

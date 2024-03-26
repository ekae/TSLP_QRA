#!/usr/bin/env python
# coding: utf-8

import networkx as nx
import itertools
from graph_tools import *
from tslp import *
import time
import datetime


def run_simulation(n, path_k, lmax, padding=0, end_node_mode="convexhull", feed_value=25, filename="",filepath="", is_draw=True):
    print(f"date/time: {datetime.datetime.now():%y-%m-%d_%H-%M-%S}")
    
    if filename != "":
        print(f"-- sim {filename} graph [{end_node_mode}] -- ")
        RGG = read_graph_from_gml(filename)
        print(f"nodes : {len(RGG.nodes())}, lmax : {lmax}, padding : {padding}")
    else:
        print(f"-- sim gabriel graph [{end_node_mode}] -- ")
        print(f"nodes : {n}, lmax : {lmax}, padding : {padding}")
        
    is_graph_unsolvable = True
    try_loop = 0
    while (is_graph_unsolvable):

        try_loop += 1

        if (try_loop == 100):
            print("break run out of tries")
            break
            
        if filename == "":
            RGG = create_gabriel_graph(n, end_node_mode=end_node_mode, padding_percent=padding)

        if (RGG != None):

            GC = GraphContainer(RGG)
            if (check_neighbour_length(GC.graph, GC.end_nodes, lmax)):
                G = GC.graph.to_directed()

                end_nodes_pairs = list(itertools.combinations(GC.end_nodes, r=2)) #list of end nodes pairs
                end_nodes_pairs = sorted_end_node_pairs(RGG,end_nodes_pairs) #rearrange end nodes pair from highest distance to lowest
                chosen_edges_list = [] #a collection of edge of optimized path between a pair of nodes 
                var_node_list = [] #a collection of node from each iterations of paht optimization
                var_node_constraints_list = [] #a collection of node constraints from each iterations of paht optimization
                total_comp_time = 0 #compuation time for current iteration
                mcf_comp_time = 0 #path optimization time
                ncpc_comp_time = 0 #node optimization time
                rep_amt = 0 #amount of repeater selected
                sol_length = 0 #sum solution paht length from all subgraph

                for (i,j) in end_nodes_pairs:

                    splited_chained_nodes = GC.break_chained_list((i,j), GC.chained_end_nodes)
                    lp_mcf = LpFeedForwardMinCostFlowNoEndnodeAsRepeater(graph=G, end_nodes=GC.end_nodes, chained_end_nodes=GC.chained_end_nodes, 
                             end_nodes_pair=(i,j), splited_chained_nodes=splited_chained_nodes, lmax=lmax, path_k=path_k, previous_chosen_edges=chosen_edges_list, feed_value = feed_value)
                    lp_mcf.solve()

                    if lp_mcf.result:

                        chosen_edges_list.extend(lp_mcf.chosen_edges)
                        subG = nx.from_edgelist(lp_mcf.chosen_edges) #convert all nodes from multi MCF path to a graph

                        nc = NodeReachConstraints(GC, subG, (i,j))
                        nc.gen_node_cover_constraints(lmax)
                        var_node_list += nc.input_nodes
                        var_node_constraints_list += nc.repeator_node_constraints
                        is_graph_unsolvable = False
                    else:
                        is_graph_unsolvable = True
                        Print("infeasible")
                        break

                if not is_graph_unsolvable:

                    lp_ncpc = LpPathTraversal(var_node_list, var_node_constraints_list, GC.end_nodes)
                    vars_ncpc = lp_ncpc.exec(is_logging=False)

                    if lp_ncpc.result:

                        print("vars/cons", lp_ncpc.num_vars, lp_ncpc.num_cons)
                        print("avg node con", nx.average_node_connectivity(G))
                        print("repeater count", len(lp_ncpc.chosen_repeater_node))
                        print("solution length", sum(GC.graph[u][v]['length'] for (u, v) in chosen_edges_list))
                        if is_draw:
                            draw_graph(GC.graph, chosen_edges=chosen_edges_list, chosen_nodes=lp_ncpc.chosen_repeater_node)
                    else:
                        print("ncpc no solution")
                        continue
            else:
                print(f"tries[{try_loop}] no possible neighbour")


run_simulation(n=40, path_k=2, lmax=130, padding=1)


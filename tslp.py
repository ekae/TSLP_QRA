from pulp import *
import networkx as nx
import numpy
import time

class LPFormulation:
    """
    Base class for linear program (LP) containing basic LP operations.

    Parameters:
        lp_name (str): Name of the linear program.
        lp_var_name (str): Name of the linear program variables.

    Attributes:
        linear_program_name (str): Name of the linear program.
        linear_prog_varible_name (str): Name of the linear program variables.
        comp_time (float): Computation time for solving the LP.
        chosen_edges (list): List of chosen edges in the solution.
        chosen_nodes (list): List of chosen nodes in the solution.
        result (bool): Result of the LP formulation (True if successful, False otherwise).
        solutions (list): List of solutions obtained from solving the LP.
    """

    def __init__(self, lp_name, lp_var_name):
        self.linear_program_name = lp_name
        self.linear_prog_varible_name = lp_var_name
        self.comp_time = 0
        self.chosen_edges = []
        self.chosen_nodes = []
        self.result = False
        self.solutions = []

    def solve(self):
        """
        Set up variables and constraints for the linear program.
        """
        pass  # Placeholder for actual implementation
    
class LpMinCostFlow(LPFormulation):
    """
    Finding the best path through minimum-cost flow.

    Parameters:
        graph (networkx.Graph): The NetworkX directed graph representing the network (a set of nodes and edges).
        end_nodes_pair (list): A list of end nodes.
        path_k (int, optional): The number of shortest paths to find. Defaults to 2.
    """
    
    def __init__(self, graph, end_nodes_pair, path_k=2):
        """
        Initialize the LpMinCostFlow class.

        Parameters:
            graph (networkx.Graph): The NetworkX directed graph representing the network (a set of nodes and edges).
            end_nodes_pair (list): A list of end nodes.
            path_k (int, optional): The number of shortest paths to find. Defaults to 2.
        """
        super().__init__(lp_name=self.__class__.__name__, lp_var_name="Flow")
        
        self.G = graph
        self.end_nodes_pair = end_nodes_pair
        self.number_of_path_k = path_k
        
    def solve(self, is_logging=False):
        """
        Solve the linear program to find k shortest paths from the given subgraph.

        Parameters:
            logging (bool, optional): Whether to enable logging. Defaults to False.

        Returns:
            list: The selected edges by MCF.
        """
        
        #Prepare vars for LP
        nodes = list(self.G.nodes)
        edges = list(self.G.edges)
        
        #Create the problem
        prob = LpProblem(self.linear_program_name, LpMinimize)
        vars = LpVariable.dicts(self.linear_prog_varible_name,edges,0,None,LpInteger)

        #add objective function
        prob += lpSum((vars[(i,j)]*(self.G[i][j]["length"]) for (i,j) in edges))
        
        #add constraints: flow cap all equal 1
        for (i,j) in edges:
            prob += (vars[(i,j)] <= 1)
            
       # Add constraints: required flow (s,t) and flow conservation constraints (u != (s,t))
        for n in nodes:
            if n == self.end_nodes_pair[0]:  # if node is source
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) - lpSum(
                    [vars[(i, j)] for (i, j) in self.G.out_edges(n)]) == self.number_of_path_k * -1
            elif n == self.end_nodes_pair[1]:  # if node is destination
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) - lpSum(
                    [vars[(i, j)] for (i, j) in self.G.out_edges(n)]) == self.number_of_path_k
            else:
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) - lpSum(
                    [vars[(i, j)] for (i, j) in self.G.out_edges(n)]) == 0


        #add constraints:incoming flow less than one for each edges except endnodes pair
        for n in nodes:
            if n not in self.end_nodes_pair:
                 prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) <= 1
        
        if is_logging: print(prob)
        
        #set timer and solve the LP
        start_time = time.time()
        prob.solve(CPLEX_PY(msg=0))
        self.comp_time = time.time() - start_time
        if is_logging: print("MCF Status:", LpStatus[prob.status])
        if prob.status == 1:   
            self.result = True
            for v in vars.items():
                if v[1].varValue == 1:
                    self.chosen_edges.append(v[0])
        
        return self.chosen_edges

class FeedForwardEdge:
    """
    Class representing a feedforward edge in a graph.
    """
    
    def __init__(self):
        self.data = {}

    def __getitem__(self, i):
        if i not in self.data:
            self.data[i] = {}
        return self.data[i]

    def __setitem__(self, i):
        if i not in self.data:
            self.data[i] = 0
        self.data[i] = value
    
    
    def __contains__(self, keys):
        if not isinstance(keys, tuple) or len(keys) != 2:
            raise ValueError("Expected a tuple of two keys.")
        
        i, j = keys
        return i in self.data and j in self.data[i]

class LpFeedForwardMinCostFlow(LPFormulation):
    """
    Class representing a linear program formulation for finding the best path through minimum-cost flow,
    enchancing with Feed forward cost reduction.

    Attributes:
        G (networkx.Graph): The NetworkX directed graph representing the network (a set of nodes and edges).
        end_nodes_pair (list): A list of end nodes.
        number_of_path_k (int): The number of shortest paths to find.
        pre_chosen_edges (list): A list of previously chosen edges.
        feed_value (int): The feed value.
    """
    
    def __init__(self, graph, end_nodes_pair, path_k = 2, previous_chosen_edges = [], feed_value = 0):
        super().__init__(lp_name = self.__class__.__name__, lp_var_name = "Flow")
        """
        Initialize the LpFeedForwardMinCostFlow object.
    
        Parameters:
            graph (networkx.Graph): The NetworkX directed graph representing the network (a set of nodes and edges).
            end_nodes_pair (list): A list of end nodes.
            path_k (int, optional): The number of shortest paths to find. Defaults to 2.
            previous_chosen_edges (list, optional): A list of previously chosen edges. Defaults to [].
            feed_value (int, optional): The feed value. Defaults to 0.
        """
        
        self.G = graph
        self.end_nodes_pair = end_nodes_pair
        self.number_of_path_k = path_k
        self.pre_chosen_edges = previous_chosen_edges
        self.feed_value = feed_value
        
    def solve(self, is_logging=False):
        """
        Find k shortest paths from the given subgraph.

        Parameters:
            is_logging (bool, optional): Flag to indicate whether logging information should be printed. Defaults to False.

        Returns:
            list: List of selected edges by MCF.
        """
        # Prepare vars for LP
        nodes = list(self.G.nodes)
        edges = list(self.G.edges)

        # Feed forward discount
        FeedFwdEdges = FeedForwardEdge()
        for (i, j) in edges:
            if (i, j) in self.pre_chosen_edges:
                FeedFwdEdges[i][j] = self.feed_value
            else:
                FeedFwdEdges[i][j] = 0

        # Create the LP problem
        prob = LpProblem(self.linear_program_name, LpMinimize)
        vars = LpVariable.dicts(self.linear_prog_varible_name, edges, 0, None, LpInteger)

        # Add objective function
        if self.pre_chosen_edges == []:
            prob += lpSum((vars[(i, j)] * (self.G[i][j]["length"]) for (i, j) in edges))
        else:
            prob += lpSum((vars[(i, j)] * (self.G[i][j]["length"] * (100 - FeedFwdEdges[i][j]) / 100) for (i, j) in edges))

        # Add constraints: flow cap all equal 1
        for (i, j) in edges:
            prob += (vars[(i, j)] <= 1)

        # Add constraints: required flow (s,t) and flow conservation constraints (u != (s,t))
        for n in nodes:
            if n == self.end_nodes_pair[0]:  # If node is source
                if is_logging:
                    print(f"Check end node {n}")
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) - \
                        lpSum([vars[(i, j)] for (i, j) in self.G.out_edges(n)]) == self.number_of_path_k * -1
            elif n == self.end_nodes_pair[1]:  # If node is destination
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) - \
                        lpSum([vars[(i, j)] for (i, j) in self.G.out_edges(n)]) == self.number_of_path_k
            else:
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) - \
                        lpSum([vars[(i, j)] for (i, j) in self.G.out_edges(n)]) == 0

        # Add constraints: incoming flow less than one for each edges except endnodes pair
        for n in nodes:
            if n not in self.end_nodes_pair:
                prob += lpSum([vars[(i, j)] for (i, j) in self.G.in_edges(n)]) <= 1

        if is_logging:
            print(prob)  # or prob.writeLP("Min_Cost_flow.lp")

        # Set timer and solve the LP
        start_time = time.time()
        prob.solve(CPLEX_PY(msg=0))
        self.comp_time = time.time() - start_time
        if is_logging:
            print("MCF Status:", LpStatus[prob.status])
        if prob.status == 1:
            self.result = True
            for v in vars.items():
                if v[1].varValue == 1:
                    self.chosen_edges.append(v[0])

        return self.chosen_edges

class LpFeedForwardMinCostFlowNoEndnodeAsRepeater(LPFormulation):
    """
    Class representing a linear program formulation for finding the best path through minimum-cost flow without end nodes acting as repeaters.

    Attributes:
        G (networkx.Graph): The NetworkX directed graph representing the network (a set of nodes and edges).
        end_nodes_list (list): A list of end nodes.
        splited_chained_nodes (list): A list of split chained nodes.
        lmax (int): The maximum distance between nodes.
        end_nodes_pair (list): A list of end nodes pairs.
        chained_end_nodes (list): A list of chained end nodes.
        number_of_path_k (int): The number of shortest paths to find.
        pre_chosen_edges (list): A list of previously chosen edges.
        feed_value (int): The feed value.
        prob (LpProblem): The linear program problem object.
    """

    def __init__(self, graph, end_nodes, chained_end_nodes, end_nodes_pair, splited_chained_nodes, lmax, path_k=2, previous_chosen_edges=[], feed_value=0):
        """
        Initialize the LpFeedForwardMinCostFlowNoEndnodeAsRepeater object.

        Parameters:
            graph (networkx.Graph): The NetworkX directed graph representing the network (a set of nodes and edges).
            end_nodes (list): A list of end nodes.
            chained_end_nodes (list): A list of chained end nodes.
            end_nodes_pair (list): A list of end nodes pairs.
            splited_chained_nodes (list): A list of split chained nodes.
            lmax (int): The maximum distance between nodes.
            path_k (int): The number of shortest paths to find. Defaults to 2.
            previous_chosen_edges (list, optional): A list of previously chosen edges. Defaults to [].
            feed_value (int): The feed value. Defaults to 0.
        """
        super().__init__(lp_name = self.__class__.__name__, lp_var_name = "Flow")
        
        self.G = graph
        self.end_nodes_list = end_nodes
        self.splited_chained_nodes = splited_chained_nodes
        self.lmax = lmax
        self.end_nodes_pair = end_nodes_pair
        self.chained_end_nodes = chained_end_nodes
        self.number_of_path_k = path_k
        self.pre_chosen_edges = previous_chosen_edges
        self.feed_value = feed_value
        self.prob = LpProblem(self.linear_program_name, LpMinimize)
        
    def solve(self, is_logging=False):
        """
        Find k shortest paths from the given subgraph and return the selected edges by MCF.

        Parameters:
            is_logging (bool, optional): Flag indicating whether to print debug logs. Defaults to False.

        Returns:
            list: List of chosen edges by MCF.
        """
        #Prepare vars for LP
        nodes = list(self.G.nodes)
        edges = list(self.G.edges)
        
        if is_logging: print( self.end_nodes_pair)
        
        #####feed fwd discount########
        FeedFwdEdges = FeedForwardEdge()
        for (i,j) in edges:
            if (i,j) in self.pre_chosen_edges:
                FeedFwdEdges[i][j] = self.feed_value
            else:
                FeedFwdEdges[i][j] = 0
            
        #Create the 'prob'
        vars = LpVariable.dicts(self.linear_prog_varible_name,edges,0,None,LpInteger)

        #add objective function
        if self.pre_chosen_edges == []:
            self.prob += lpSum((vars[(i,j)]*(self.G[i][j]["length"]) for (i,j) in edges))
        else:
            self.prob += lpSum((vars[(i,j)]*(self.G[i][j]["length"]*(100-FeedFwdEdges[i][j])/100) for (i,j) in edges))
        
        #add constraints: flow cap all equal 1
        for (i,j) in edges:
            self.prob += (vars[(i,j)] <= 1)
            
        #add constraints: required flow (s,t) and flow convervation constraints (u != (s,t))
        for n in nodes:
            
            if n == self.end_nodes_pair[0]: #if node is source
                # self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == self.number_of_path_k
                self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) - lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == self.number_of_path_k * -1
            elif n == self.end_nodes_pair[1]: #if node is destination
                # self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) == self.number_of_path_k
                self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) - lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == self.number_of_path_k#check if the total length of chosen edge in and out exeed Lmax
            else:
                self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) - lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == 0

        #add constraints: selecting edge of an end node must not exceed Lmax
        for n in self.end_nodes_list:
            if n not in self.end_nodes_pair:
                self.prob += lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.in_edges(n)) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.out_edges(n)) <= self.lmax
        
        #add constraints: selecting edge of chain end nodes must not exceed Lmax
        # head of chain[in/out edges] - intermediate chain [brigde edges] + end of chain[in/out edges]
        for chain in self.splited_chained_nodes:
            intermediate_chain = []
            for idx_node, node in enumerate(chain):
                if idx_node != len(chain)-1:
                    intermediate_chain.append((chain[idx_node],chain[idx_node+1]))
                    intermediate_chain.append((chain[idx_node+1],chain[idx_node]))
            self.prob += lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.in_edges(chain[0])) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.out_edges(chain[0])) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.in_edges(chain[len(chain)-1])) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.out_edges(chain[len(chain)-1])) - lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in intermediate_chain) <= self.lmax
            for (i,j) in intermediate_chain:
                lenght = self.G[i][j]["length"]
                if is_logging: print(f"chain{chain} >> {lenght}*{(i,j)}")
            

        #add constraints:incoming flow less than one for each edges except endnodes pair
        for n in nodes:
            if n not in self.end_nodes_pair:
                 self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) <= 1
        
        if is_logging: print(self.prob) # or prob.writeLP("Min_Cost_flow.lp")   
        
        #set timer and solve the LP
        start_time = time.time()
        self.prob.solve(CPLEX_PY(msg=0)) #msg=0
        self.comp_time = time.time() - start_time
        if is_logging: print("MCF Status:", LpStatus[self.prob.status])
        if self.prob.status == 1:   
            self.result = True
            for v in vars.items():
                if is_logging: print(f"{v[0]} = {v[1].varValue}")
                if v[1].varValue > 0:
                    self.chosen_edges.append(v[0])
                if not v[1].varValue.is_integer():
                    print("Non integer solution detected: ", v[1], "=", v[1].varValue)
        
        return self.chosen_edges    

class NodeReachConstraints:
    """
    Class that holds the related data for calculating the NCPC constraints.
    """

    def __init__(self, graph_container, subgraph, end_nodes_pair):
        """
        Initialize the object.

        Parameters:
            graph_container (GraphContainer): Graph container object containing the network graph and relevant properties.
            subgraph (nx.Graph): Subgraph representing a portion of the network.
            end_nodes_pair (list): A list containing the pair of end nodes.

        Attributes:
            G (nx.Graph): The network graph.
            end_nodes (list): List of end nodes in the graph.
            end_nodes_pair (list): Pair of end nodes.
            input_paths (list): List to store input paths.
            input_nodes (list): List to store input nodes.
            left_end_nodes_coverage (list): List to store left end nodes coverage.
            right_end_nodes_coverage (list): List to store right end nodes coverage.
            repeater_node_constraints (list): List to hold repeater node constraints.
        """
        self.G = graph_container.graph
        self.end_nodes = graph_container.end_nodes
        self.end_nodes_pair = end_nodes_pair
        self.input_paths = []
        self.input_nodes = []
        self.left_end_nodes_coverage = []
        self.right_end_nodes_coverage = []
        self.repeator_node_constraints = []            
        self.input_paths = self.gen_input_path_list(subgraph)
        self.input_nodes = self.gen_input_nodes_list()
        self.repeator_node_constraints = self.gen_rep_con_holder()
        
    #select node with linklenght for each Zi    
    def gen_node_cover_constraints(self, Lmax=130.0, is_logging=False):
        """
        Generate node coverage constraints.

        Parameters:
            Lmax (float): Maximum link length.
            is_logging (bool): Flag indicating whether to log intermediate steps.

        Returns:
            None
        """
        # Iterate over paths
        for idx_path, path in enumerate(self.input_paths): # [[a,b,c],[d,e,f]] >> [[ab,bc,],[de,ef]]
            #finding constraint for each-path
            if is_logging: print("Path: {} {}" .format(idx_path, path))

            #For loop of  element Zi
            for idx_node, node in enumerate(self.input_nodes[idx_path]):
                #prepare coverage of end nodes in the path
                if node in self.end_nodes:
                    if is_logging: print("end node: ",node)
                    self.node_cover(idx_path, idx_node, Lmax, is_logging)
            for idx_node, node in enumerate(self.input_nodes[idx_path]):
                #finding the remaining coverage for candinate repeater node
                if node not in self.end_nodes:
                    if is_logging: print("repeator node: ",node)
                    self.node_cover(idx_path, idx_node, Lmax, is_logging)
    
    def node_cover(self, idx_path, idx_node, Lmax, is_logging):
        """
        Calculate node coverage for a given path and node index.

        Parameters:
            idx_path (int): Index of the path.
            idx_node (int): Index of the node.
            Lmax (float): Maximum link length.
            is_logging (bool): Flag indicating whether to log intermediate steps.

        Returns:
            None
        """
        #For loop each next node from Zi     
            #check if right-node within linklenght
                #yes append Zi right list
        distance = 0.0
        idx_i = idx_node
        idx_j = idx_node + 1 #start from the next edge         
        repeator_node_right_constraints = []


        while distance <= Lmax and idx_j < len(self.input_nodes[idx_path]): #traverse next node

            distance += self.G[self.input_nodes[idx_path][idx_i]][self.input_nodes[idx_path][idx_j]]["length"]        
            if is_logging: print("[Checking]{}-{} dist={}".format(self.input_nodes[idx_path][idx_i],self.input_nodes[idx_path][idx_j],distance))
            if distance <= Lmax:
                if idx_node == 0: #self.input_nodes[idx_path][idx_node] in self.end_nodes:
                    self.right_end_nodes_coverage.append(self.input_nodes[idx_path][idx_j])                    
                    if is_logging: print("%s[%d] is in right-cover of %s at %f" %(self.input_nodes[idx_path][idx_j],idx_j,self.input_nodes[idx_path][idx_node], distance))
                idx_i += 1
                idx_j += 1


        #For loop start end node from Zi
                #check if left-node within linklenght
                    #yes append Zi left list
        distance = 0.0
        idx_j = idx_node
        idx_i = idx_node - 1 #start from the previous edge     
        repeator_node_left_constraints = []

        while distance <= Lmax and idx_i > 0: #and idx_node >= 0: #traverse previous node

            distance += self.G[self.input_nodes[idx_path][idx_i]][self.input_nodes[idx_path][idx_j]]["length"]
            if is_logging: print("[Checking]{}-{} dist={}".format(self.input_nodes[idx_path][idx_i],self.input_nodes[idx_path][idx_j],distance))
            if distance <= Lmax:
                if self.input_nodes[idx_path][idx_node] == self.end_nodes:
                    self.left_end_nodes_coverage.append(self.input_nodes[idx_path][idx_i])
                    if is_logging: print("%s[%d] is in left-cover of %s at %f" %(self.input_nodes[idx_path][idx_j],idx_j,self.input_nodes[idx_path][idx_node], distance))
                elif self.input_nodes[idx_path][idx_node] not in self.right_end_nodes_coverage:
                    # if self.input_nodes[idx_path][idx_i] not in self.end_nodes:
                    repeator_node_left_constraints.append(self.input_nodes[idx_path][idx_i])
                    if is_logging: print("%s[%d] is in left-range of %s at %f" %(self.input_nodes[idx_path][idx_i],idx_i,self.input_nodes[idx_path][idx_node], distance))
                    # else:
                    #     if is_logging: print("%s[%d] is in left-range of %s at %f ,as an endnode putting defualt 1 value" %(self.input_nodes[idx_path][idx_i],idx_i,self.input_nodes[idx_path][idx_node], distance))
                    #     repeator_node_left_constraints.append(1)
                    #     #if is_logging: print("%s[%d] is in left-range of %s at %f **Skipped**" %(self.input_nodes[idx_path][idx_i],idx_i,self.input_nodes[idx_path][idx_node], distance))
                idx_i -= 1
                idx_j -= 1



        #a Left constraint for Node_i is collect, Let put it in Node_i constraint
        self.repeator_node_constraints[idx_path][idx_node].append(repeator_node_left_constraints)
        #a Right constraint for Node_i is collect, Let put it in Node_i constraint
        self.repeator_node_constraints[idx_path][idx_node].append(repeator_node_right_constraints)
                    
    #select node within linklenght for each Zi
    #left only and no endnodes as repeater
    def gen_node_cover_constraints_no_endnode_repeater(self, Lmax=136.0, is_logging=False):
        """
        Generate node coverage constraints without end nodes as repeaters.

        Parameters:
            Lmax (float): Maximum link length.
            is_logging (bool): Flag indicating whether to log intermediate steps.

        Returns:
            None
        """
        #For Loop path
        for idx_path, path in enumerate(self.input_paths): # [[a,b,c],[d,e,f]] >> [[ab,bc,],[de,ef]]
            #finding constraint for each-path
            if is_logging: print("Path: {} {}" .format(idx_path, path))

            #For loop of  element Zi
            for idx_node, node in enumerate(self.input_nodes[idx_path]):
                #prepare coverage of end nodes in the path
                if (node in self.end_nodes):
                    if is_logging: print("end node: ",node)
                    self.node_cover_no_endnode_repeater(idx_path, idx_node, Lmax, is_logging)
            for idx_node, node in enumerate(self.input_nodes[idx_path]):
                #finding the remaining coverage for candinate repeater node
                if node not in self.end_nodes:
                    if is_logging: print("repeator node: ",node)
                    self.node_cover_no_endnode_repeater(idx_path, idx_node, Lmax, is_logging)
                    
    def node_cover_no_endnode_repeater(self, idx_path, idx_node, Lmax, is_logging):
        """
        Calculate node coverage for a given path and node index without end nodes as repeaters.

        Parameters:
            idx_path (int): Index of the path.
            idx_node (int): Index of the node.
            Lmax (float): Maximum link length.
            is_logging (bool): Flag indicating whether to log intermediate steps.

        Returns:
            None
        """
        #For loop each next node from Zi     
            #check if right-node within linklenght
                #yes append Zi right list
        distance = 0.0
        idx_i = idx_node
        idx_j = idx_node + 1 #start from the next edge         
        repeator_node_right_constraints = []


        while distance <= Lmax and idx_j < len(self.input_nodes[idx_path]): #traverse next node

            distance += self.G[self.input_nodes[idx_path][idx_i]][self.input_nodes[idx_path][idx_j]]["length"]        
            if is_logging: print("[Checking]{}-{} dist={}".format(self.input_nodes[idx_path][idx_i],self.input_nodes[idx_path][idx_j],distance))
            if distance <= Lmax:
                if self.input_nodes[idx_path][idx_node] in self.end_nodes:
                    if idx_node == 0:
                        self.right_end_nodes_coverage.append(self.input_nodes[idx_path][idx_j])                    
                        if is_logging: print("%s[%d] is in right-cover of %s at %f" %(self.input_nodes[idx_path][idx_j],idx_j,self.input_nodes[idx_path][idx_node], distance))
                idx_i += 1
                idx_j += 1
                #old
                # if self.input_nodes[idx_path][idx_node] in self.end_nodes:
                #     if idx_node != 0 and idx_node != len(self.input_nodes[idx_path])-1:
                #         self.right_end_nodes_coverage.append(self.input_nodes[idx_path][idx_j])                    
                #         if is_logging: print("%s[%d] is in right-cover of %s at %f" %(self.input_nodes[idx_path][idx_j],idx_j,self.input_nodes[idx_path][idx_node], distance))
                # elif self.input_nodes[idx_path][idx_node] not in self.left_end_nodes_coverage:
                #     if self.input_nodes[idx_path][idx_j] not in self.end_nodes:
                #         repeator_node_right_constraints.append(self.input_nodes[idx_path][idx_j])
                #         if is_logging: print("%s[%d] is in right-range of %s at %f" %(self.input_nodes[idx_path][idx_j],idx_j,self.input_nodes[idx_path][idx_node], distance))
                # idx_i += 1
                # idx_j += 1


        #For loop each previous node from Zi
                #check if left-node within linklenght
                    #yes append Zi left list
        distance = 0.0
        idx_j = idx_node
        idx_i = idx_node - 1 #start from the next edge     
        repeator_node_left_constraints = []

        while distance <= Lmax and idx_i > 0: #and idx_node >= 0: #traverse previous node

            distance += self.G[self.input_nodes[idx_path][idx_i]][self.input_nodes[idx_path][idx_j]]["length"]
            if is_logging: print("[Checking]{}-{} dist={}".format(self.input_nodes[idx_path][idx_i],self.input_nodes[idx_path][idx_j],distance))
            if distance <= Lmax:
                if self.input_nodes[idx_path][idx_node] == self.end_nodes:
                    if idx_node == len(self.input_nodes[idx_path])-1:
                        self.left_end_nodes_coverage.append(self.input_nodes[idx_path][idx_i])
                        if is_logging: print("%s[%d] is in left-cover of %s at %f" %(self.input_nodes[idx_path][idx_j],idx_j,self.input_nodes[idx_path][idx_node], distance))
                elif self.input_nodes[idx_path][idx_node] not in self.right_end_nodes_coverage:
                    if self.input_nodes[idx_path][idx_i] not in self.end_nodes:
                        repeator_node_left_constraints.append(self.input_nodes[idx_path][idx_i])
                        if is_logging: print("%s[%d] is in left-range of %s at %f" %(self.input_nodes[idx_path][idx_i],idx_i,self.input_nodes[idx_path][idx_node], distance))
                idx_i -= 1
                idx_j -= 1



        #a Left constraint for Node_i is collect, Let put it in Node_i constraint
        self.repeator_node_constraints[idx_path][idx_node].append(repeator_node_left_constraints)
        #a Right constraint for Node_i is collect, Let put it in Node_i constraint
        self.repeator_node_constraints[idx_path][idx_node].append(repeator_node_right_constraints)
    
    
    def gen_input_path_list(self,subG):
        """
        Generate a list of input paths from a subgraph of optimal route.

        Parameters:
            subG (nx.Graph): The subgraph representing the optimal route.

        Returns:
            list: A list of input paths.
        """        
        input_paths_list = []
        for path in sorted(nx.all_simple_edge_paths(subG, self.end_nodes_pair[0], self.end_nodes_pair[1])):
            input_paths_list.append(path)
        # print(input_paths_list)
        return input_paths_list
    
    def gen_input_nodes_list(self):
        """
        Generate a list of input nodes from a subgraph of optimal route.
        """
        input_nodes_list = []
        for i in range(len(self.input_paths)):
            input_nodes_list.append([])
        for idx_path, path in enumerate(self.input_paths):
            for idx_edge, edge in enumerate(path):        
                if idx_edge == 0:
                    input_nodes_list[idx_path].append(edge[0])
                input_nodes_list[idx_path].append(edge[1])
        # print(input_nodes_list)
        return input_nodes_list
                
    def gen_rep_con_holder(self):
        """
        Generate a list of holder for constraints.
        """
        holder = []
        
        for i in range(len(self.input_nodes)):
            holder.append([])
            for j in range(len(self.input_nodes[i])):
                # print(i,"-",j)
                holder[i].append([])
        return holder
   
class LpPathTraversal:
    """
    Class for solving the path traversal (PT) problem using linear programming.
    """
    def __init__(self, input_nodes_list ,repeator_node_constraints, end_nodes):
        """
        Initialize the LpPathTraversal object.

        Parameters:
            input_nodes_list (list): A list containing input nodes.
            repeater_node_constraints (list): A list containing repeater node constraints.
            end_nodes (list): A list of end nodes.

        Attributes:
            lp_name (str): The name of the linear programming problem.
            lp_var_name (str): The name of the LP variables.
            input_nodes (list): A list containing input nodes.
            end_nodes (list): A list of end nodes.
            repeator_node_cons (list): A list containing repeater node constraints.
            comp_time (float): The computation time.
            chosen_repeater_node (list): The chosen repeater nodes.
            result (bool): The result of the LP problem.
            num_vars (int): The number of variables.
            num_cons (int): The number of constraints.
        """
        self.lp_name = "nodecover_pathcover_Problem"
        self.lp_var_name = "Node"
        self.input_nodes = input_nodes_list
        self.end_nodes = end_nodes
        self.repeator_node_cons = repeator_node_constraints
        self.comp_time = 0
        self.chosen_repeater_node = []
        self.result = False
        self.num_vars = 0
        self.num_cons = 0    
     
    def exec(self, is_logging=False):
        """
        Find the least repeater location covering all the paths.

        Parameters:
            is_logging (bool): Whether to enable logging.

        Returns:
            dict: Dictionary containing the LP variables.
        """
        #Prepare vars for LP
        var_nodes = set()
        var_repeater_nodes = set()
        
        #add only repeaters which appear in solution node
        for path_num in range(len(self.input_nodes)):
            for node in self.input_nodes[path_num]:
                var_nodes.add(node)
        for node in var_nodes:
            if node not in self.end_nodes:
                var_repeater_nodes.add(node)
            
        if is_logging: print(var_nodes)
        
        #Create linear problem
        prob = LpProblem(self.lp_name, LpMinimize)
        
        #Create variable dict
        vars = LpVariable.dicts(self.lp_var_name,var_nodes,0,None,LpInteger)
        for node in self.end_nodes: #set end nodes value to 1
            vars[node].setInitialValue(1)
        
        #First, add obj function that has only repeater node
        prob += lpSum((vars[node]) for node in var_repeater_nodes)
        
        #add repeater location constraints
        for path_num in range(len(self.input_nodes)):
            for node_cons_leftrigth in self.repeator_node_cons[path_num]:
                if len(node_cons_leftrigth[0]) > 0:
                    prob += lpSum((vars[node]) for node in node_cons_leftrigth[0]) >= 1
                if len(node_cons_leftrigth[1]) > 0:
                    prob += lpSum((vars[node]) for node in node_cons_leftrigth[1]) >= 1
        if is_logging: print(prob)

        # prob.writeLP("nodecover_pathcover_Problem.lp")
        start_time = time.time()
        prob.solve(CPLEX_PY(msg=0))
        self.comp_time = time.time() - start_time
        if is_logging: print("NCPC Status:", LpStatus[prob.status])
        if prob.status == 1:   
            self.result = True
            self.num_vars = len(var_repeater_nodes)
            self.num_cons = prob.numConstraints()
            if is_logging: print("Number of var-con:", self.num_vars, ' ', self.num_cons)
            obj_value = 0   
            for v in vars.items():
                if v[0] in var_repeater_nodes and v[1].varValue > 0:
                    self.chosen_repeater_node.append(v[0])
                    if is_logging: 
                        obj_value += v[1].varValue
                        if not v[1].varValue.is_integer():
                            print("Non integer solution detected: ", v[1], "=", v[1].varValue)
                        else:
                            print(v[1], "=", v[1].varValue)
            if is_logging: print("objective value = ",obj_value)
                    
        return vars
    
class LpFeedForwardMinCostFlowNoEndnodeAsRepeater33(LPFormulation):
    """
    Finding best path through mininum-cost flow.

    graph : networkx.Graph
        The NetworkX directed graph representing the network (a set of nodes and edges).
    end_nodes_pair : list
        a list of end nodes
    """
    
    def __init__(self, graph, end_nodes, chained_end_nodes, end_nodes_pair, splited_chained_nodes, lmax, path_k = 2, previous_chosen_edges = [], feed_value = 0):
        super().__init__(lp_name = self.__class__.__name__, lp_var_name = "Flow")
        
        self.G = graph
        self.end_nodes_list = end_nodes
        self.splited_chained_nodes = splited_chained_nodes
        self.lmax = lmax
        self.end_nodes_pair = end_nodes_pair
        self.chained_end_nodes = chained_end_nodes
        self.number_of_path_k = path_k
        self.pre_chosen_edges = previous_chosen_edges
        self.feed_value = feed_value
        self.prob = LpProblem(self.linear_program_name, LpMinimize)
        
    def solve(self, is_logging=False):
        
        """
        find k shortest path from given subgraph
        
        return chosen_edges: the selected edges by MCF
        """
        
        #Prepare vars for LP
        nodes = list(self.G.nodes)
        edges = list(self.G.edges)
        
        if is_logging: print( self.end_nodes_pair)
        
        #####feed fwd discount########
        FeedFwdEdges = FeedForwardEdge()
        for (i,j) in edges:
            if (i,j) in self.pre_chosen_edges:
                FeedFwdEdges[i][j] = self.feed_value
            else:
                FeedFwdEdges[i][j] = 0
            
        #Create the 'prob'
        vars = LpVariable.dicts(self.linear_prog_varible_name,edges,0,None,LpInteger)

        #add objective function
        if self.pre_chosen_edges == []:
            self.prob += lpSum((vars[(i,j)]*(self.G[i][j]["length"]) for (i,j) in edges))
        else:
            self.prob += lpSum((vars[(i,j)]*(self.G[i][j]["length"]*(100-FeedFwdEdges[i][j])/100) for (i,j) in edges))
        
        #add constraints: flow cap all equal 1
        for (i,j) in edges:
            self.prob += (vars[(i,j)] <= 1)
            
        #add constraints: required flow (s,t) and flow convervation constraints (u != (s,t))
        for n in nodes:
            
            if n == self.end_nodes_pair[0]: #if node is source
                # self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == self.number_of_path_k
                self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) - lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == self.number_of_path_k * -1
            elif n == self.end_nodes_pair[1]: #if node is destination
                # self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) == self.number_of_path_k
                self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) - lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == self.number_of_path_k#check if the total length of chosen edge in and out exeed Lmax
            else:
                self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) - lpSum([vars[(i,j)] for (i,j) in self.G.out_edges(n)]) == 0

        #add constraints: selecting edge of an end node must not exceed Lmax
        for n in self.end_nodes_list:
            if n not in self.end_nodes_pair:
                self.prob += lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.in_edges(n)) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.out_edges(n)) <= self.lmax
        
        #add constraints: selecting edge of chain end nodes must not exceed Lmax
        # head of chain[in/out edges] - intermediate chain [brigde edges] + end of chain[in/out edges]
        for chain in self.splited_chained_nodes:
            intermediate_chain = []
            for idx_node, node in enumerate(chain):
                if idx_node != len(chain)-1:
                    intermediate_chain.append((chain[idx_node],chain[idx_node+1]))
                    intermediate_chain.append((chain[idx_node+1],chain[idx_node]))
            self.prob += lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.in_edges(chain[0])) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.out_edges(chain[0])) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.in_edges(chain[len(chain)-1])) + lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in self.G.out_edges(chain[len(chain)-1])) - lpSum(vars[(i,j)]*self.G[i][j]["length"] for (i,j) in intermediate_chain) <= self.lmax
            for (i,j) in intermediate_chain:
                lenght = self.G[i][j]["length"]
                if logging: print(f"chain{chain} >> {lenght}*{(i,j)}")
            

        #add constraints:incoming flow less than one for each edges except endnodes pair
        for n in nodes:
            if n not in self.end_nodes_pair:
                 self.prob += lpSum([vars[(i,j)] for (i,j) in self.G.in_edges(n)]) <= 1
        
        if logging: print(self.prob) # or prob.writeLP("Min_Cost_flow.lp")   
        
        #set timer and solve the LP
        start_time = time.time()
        self.prob.solve(CPLEX_PY(msg=0)) #msg=0
        self.comp_time = time.time() - start_time
        if logging: print("MCF Status:", LpStatus[self.prob.status])
        if self.prob.status == 1:   
            self.result = True
            for v in vars.items():
                if logging: print(f"{v[0]} = {v[1].varValue}")
                if v[1].varValue > 0:
                    self.chosen_edges.append(v[0])
                if not v[1].varValue.is_integer():
                    print("Non integer solution detected: ", v[1], "=", v[1].varValue)
        
        return self.chosen_edges    
import numpy as np
import math
from scipy.stats import norm, rankdata
from scipy.spatial import ConvexHull, distance, Delaunay, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import ast
import os
import random
import heapq


class GraphContainer:
    """
    Class that holds the graph and some relevant properties to solve the repeater allocation problem.

    Parameters
    ----------
    graph (Networkx Graph): NetworkX undirected graph.

    Attributes
    ----------
    graph : networkx.Graph
        Graph object representing the network.
    end_nodes : list
        List of end nodes in the graph, consisted of node in set T.
    possible_rep_nodes : list
        List of possible repeater node locations, consisted of node in set R.
    unique_end_node_pairs : list
        List of all unique source-destination pairs.
    num_end_nodes : int
        Total number of end nodes.
    num_repeater_nodes : int
        Total number of repeater nodes.
    num_unique_pairs : int
        Total number of unique source-destination pairs.
    chained_end_nodes : list
        List of chained end nodes for omitting them in the LP.
    """
    def __init__(self, graph):
        self.graph = graph
        self.end_nodes = []
        self.possible_rep_nodes = []
        
        for node, nodedata in graph.nodes.items():
            if nodedata["type"] == 'end_node':
                self.end_nodes.append(node)
            else:
                self.possible_rep_nodes.append(node)
        
        self.num_end_nodes = len(self.end_nodes)
        self.num_repeater_nodes = len(self.possible_rep_nodes)
        if self.num_repeater_nodes == 0:
            # Trivial graph
            return
        self.unique_end_node_pairs = list(itertools.combinations(self.end_nodes, r=2))
        self.num_unique_pairs = len(list(self.unique_end_node_pairs))
        self.chained_end_nodes = self._find_chained_end_node()
        
        # Add length parameter to edges if this is not defined yet
        for i, j in graph.edges():
            if 'length' not in graph[i][j]:
                if 'Longitude' in graph.nodes[self.possible_rep_nodes[0]]:
                    self._compute_dist_lat_lon(graph)
                else:
                    self._compute_dist_cartesian(graph)
                break
    
    def _compute_dist_lat_lon(self, graph):
        """Compute the distance in km between two points based on their latitude and longitude.
        Assumes both are given in radians.

        Parameters:
            graph (NetworkX graph): Graph containing nodes with latitude and longitude attributes.

        Returns:
            None
        """
        R = 6371  # Radius of the earth in km
        for edge in graph.edges():
            node1, node2 = edge
            lon1 = np.radians(graph.nodes[node1]['Longitude'])
            lon2 = np.radians(graph.nodes[node2]['Longitude'])
            lat1 = np.radians(graph.nodes[node1]['Latitude'])
            lat2 = np.radians(graph.nodes[node2]['Latitude'])
            delta_lat = lat2 - lat1
            delta_lon = lon2 - lon1
            a = np.sin(delta_lat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin(delta_lon / 2) ** 2)
            c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
            dist = np.round(R * c, 5)
            graph.edges[node1, node2]['length'] = dist

    def _compute_dist_cartesian(self, graph):
        """
        Compute the distance in kilometers between two points based on their Cartesian coordinates.

        Parameters:
        graph (NetworkX graph): Graph containing nodes with Cartesian coordinates.

        Returns:
        None
        """
        for edge in graph.edges():
            node1, node2 = edge
            # dx = np.abs(graph.nodes[node1]['xcoord'] - graph.nodes[node2]['xcoord'])
            # dy = np.abs(graph.nodes[node1]['ycoord'] - graph.nodes[node2]['ycoord'])
            dx = np.abs(graph.nodes[node1]['pos'][0] - graph.nodes[node2]['pos'][0])
            dy = np.abs(graph.nodes[node1]['pos'][1] - graph.nodes[node2]['pos'][1])
            dist = np.round(np.sqrt(np.square(dx) + np.square(dy)), 5)
            graph.edges[node1, node2]['length'] = dist

    def _print_graph_data(self):
        """
        Print various statistics about the graph.

        Parameters:
            None

        Returns:
            None
        """
        total_length, min_length, max_length = 0, float('inf'), 0
        num_edges = len(self.graph.edges())

        for i, j in self.graph.edges():
            length = self.graph[i][j]['length']
            total_length += length
            if length > max_length:
                max_length = length
            if length < min_length:
                min_length = length

        print("Total number of nodes:", len(self.graph.nodes()))
        print("Number of end nodes:", self.num_end_nodes)
        print("Number of possible repeater nodes:", self.num_repeater_nodes)
        print("Total number of edges:", num_edges)
        print("Average length is", total_length / num_edges)
        print("Maximum edge length is", max_length)
        print("Minimum edge length is", min_length)
    
    def _find_chained_end_node(self, is_logging=False):
        """
        Find chained end nodes in the graph.

        Parameters:
            is_logging (bool): Whether logging is enabled.

        Returns:
            list: List of lists containing chained end nodes.
        """
        chained_end_node = []
        
        for end_node in self.end_nodes:
            result = []
            result.append(end_node)
            non_visited_end_node = self.end_nodes
            visited_node = None
            there_is_more = True
            while (there_is_more):
                non_visited_end_node = [node for node in non_visited_end_node if node not in result]
                if is_logging: print(f"non_visit_end_node: {non_visited_end_node}")
                result = self._deep_find_chained_end_node(self.graph, end_node, non_visited_end_node, visited_node, is_logging)
                if result:
                    chained_end_node.append(result)
                else:
                    there_is_more = False
        #clean the output
        chained_end_node = self._remove_reversed_list(chained_end_node)
        
        return chained_end_node
    
     #a recursive of find chained end node
    def _deep_find_chained_end_node(self, G, start_node, end_nodes, visited_node=None, is_logging=False):
        """
        Recursively find chained end nodes starting from a given node.

        Parameters:
            G (NetworkX graph): Input graph.
            start_node (str): The node to start the search from.
            end_nodes (list): List of end nodes.
            visited_node (set): Set of visited nodes.
            is_logging (bool): Whether logging is enabled.

        Returns:
            list: List containing chained end nodes.
        """
        if visited_node is None:
            visited_node = set()

        visited_node.add(start_node)
        if is_logging:
            print(f"start node: {start_node}, visited: {visited_node}")
        for n_node in list(G.neighbors(start_node)):
            if is_logging:
                print(f"check on {n_node}")
            if n_node in end_nodes and n_node not in visited_node:
                if is_logging:
                    print(f"{start_node}---", end="")
                path = self._deep_find_chained_end_node(G, n_node, end_nodes, visited_node, is_logging)
                if path:
                    if is_logging:
                        print(f"return path: {path}")
                    return [start_node] + path

        if start_node in end_nodes:
            return [start_node]
        else:
            return None
    
    # Remove reversed lists
    def _remove_reversed_list(self, mix_reverse_list):
        """
        Remove reversed sublists from a list of sublists.

        Parameters:
            mix_reverse_list (list): List of sublists.

        Returns:
            list: List of unique sublists without their reversed counterparts.
        """
        # Convert each sublist to a tuple to make them hashable
        lists_tuples = [tuple(sublist) for sublist in mix_reverse_list]
        unique_lists_tuples = []

        for sublist_tuple in lists_tuples:
            reverse_sublist_tuple = tuple(reversed(sublist_tuple))
            if reverse_sublist_tuple not in unique_lists_tuples:
                unique_lists_tuples.append(sublist_tuple)
        
        return unique_lists_tuples
    
    @staticmethod
    def break_chained_list(end_node_pair, chained_end_nodes, is_logging=False):
        """
        Break chained end nodes list into separate chains.

        Parameters:
            end_node_pair (list): List of end nodes.
            chained_end_nodes (list): List of chained end nodes.
            is_logging (bool): Whether logging is enabled.

        Returns:
            list: List of updated chained end nodes.
        """
        
        updated_chained_end_nodes = []
            
        for chain_tuple in chained_end_nodes:
            chain = list(chain_tuple)
            if is_logging: print(f"\nchain = {chain}")
            is_end_node_found = False
            for end_node in end_node_pair:
                try:
                    index = chain.index(end_node)
                    is_end_node_found = True
                    if index == 0 or index == (len(chain)-1):
                        if is_logging: print(f"--{end_node} is founded {chain} >>", end="")
                        chain.remove(end_node)
                        if is_logging: print(chain) 
                        if chain not in updated_chained_end_nodes:
                            updated_chained_end_nodes.append(chain)
                            if is_logging: print(f"--updated_chained_end_nodes = {updated_chained_end_nodes}")
                    else:
                        if is_logging: print(f"--{end_node} is founed in between..")
                        if (chain[:index] not in updated_chained_end_nodes):
                            updated_chained_end_nodes.append(chain[:index])
                            if is_logging: print(f"--updated_chained_end_nodes = {updated_chained_end_nodes}")
                        chain = chain[index+1:]  #remove the first part(chain[:index])
                        if (chain not in updated_chained_end_nodes):
                            updated_chained_end_nodes.append(chain)
                            if is_logging: print(f"--updated_chained_end_nodes = {updated_chained_end_nodes}")
                        chain.remove(end_node)
                except ValueError:
                    if is_logging: print(f"--{end_node} not found in chain {chain}")
            if len(chain) > 0 and is_end_node_found == False and chain not in updated_chained_end_nodes:
                updated_chained_end_nodes.append(chain)
                if is_logging: print(f"--updated_chained_end_nodes = {updated_chained_end_nodes}")

        updated_chained_end_nodes = [chain for chain in updated_chained_end_nodes if chain]
        updated_chained_end_nodes = [chain for chain in updated_chained_end_nodes if len(chain) != 1]

        return updated_chained_end_nodes

def graph_to_gml(graph, foldername, filename):
    """
    Save a graph in GML format.

    Parameters:
        graph (NetworkX graph): The graph to be saved.
        foldername (str): Name of the folder where the file will be saved.
        filename (str): Name of the file to be saved.
    """
    os.makedirs(foldername, exist_ok=True)
    filepath = os.path.join(foldername, filename)
    nx.write_gml(graph, filepath)
       
    

def read_graph_from_gml(file, end_node_pair=[], is_draw = False):
    """
    Read a graph from a GML file.

    Parameters:
        file (str): Path to the GML file.
        end_node_pair (list): List of end nodes pairs.
        is_draw (bool): Whether to draw the graph.

    Returns:
        NetworkX graph: The read graph.
    """
    G = nx.read_gml(file)
    file_name = file[:-4]
    
            
    # SurfnetCore nodes list
    #  "A'dam 1",A'dam 2, 'Almere', 'Amersfoort', 'Apeldoorn', 'Arnhem', 'Breda 1', 
    #  'Breda 2', 'Delft 1', 'Delft 2', 'Den Bosch', 'Den Haag 1', 'Den Haag 2', 
    #  'Deventer', 'Dordrecht', 'Dwingeloo', 'Eindhoven 1', 'Eindhoven 2', 'Enschede 1', 
    #  'Enschede 2', 'Gorinchem', 'Groningen 1', 'Groningen 2', 'Harlingen', 'Heerlen', 
    #  'Hilversum', 'Hoorn', 'Leeuwarden', 'Leiden 1', 'Leiden 2', 'Leiden 3', 'Lelystad 1', 
    #  'Lelystad 2', 'Maasbracht', 'Maastricht', 'Meppel', 'Middenmeer', 'Nieuwegein', 
    #  'Nijmegen 1', 'Nijmegen 2', "R'dam 1", "R'dam 2", "R'dam 3", 'Tilburg 1', 'Tilburg 2', 
    #  'Utrecht 1', 'Utrecht 2', 'Venlo', 'Wageningen 1', 'Wageningen 2', 'Westerbork', 
    #  'Zutphen', 'Zwolle 1', 'Zwolle 2'    

    if file_name == "SurfnetCore":
        if end_node_pair:
            print("Manual endnodes input")
            end_node_list = end_node_pair
        else:
            end_node_list = ["Delft 1", "Groningen 1", "Maastricht", "Enschede 2"]
    else:
        raise NotImplementedError(f"Unknown Dataset: {file_name}")

    pos = {}
    for node, nodedata in G.nodes.items():
        if "position" in nodedata:
            pos[node] = ast.literal_eval(nodedata["position"])
        elif "Longitude" in nodedata and "Latitude" in nodedata:
            pos[node] = [nodedata['Longitude'], nodedata['Latitude']]
        else:
            raise ValueError("Cannot determine node position.")
        nodedata['type'] = 'end_node' if node in end_node_list else 'repeater_node'
        
    nx.set_node_attributes(G, pos, name='pos')
    
    if is_draw:
        draw_graph(G)
    return G

def find_two_largest_multipliers(n):
    """
    Find the two largest multipliers of a given number.

    Parameters:
        n (int): The number to find the multipliers for.

    Returns:
        tuple: A tuple containing the two largest multipliers.
    """
    largest = 1
    for i in range(int(math.sqrt(n)), 0, -1):
        if n % i == 0:
            largest = i
            break
    second_largest = n // largest
    return largest, second_largest

def generate_grid(square_area_side_length, grid_number):
    """
    Generate a grid of rectangles within a square area.

    Parameters:
        square_area_side_length (int): The side length of the square area.
        grid_number (int): The number of rectangles to generate.

    Returns:
        tuple: A tuple containing the width and height of each rectangle,
               and a list of coordinates representing the grid.
    """
    cols, rows = find_two_largest_multipliers(grid_number)

    # Calculate the width and height of each rectangle.
    rect_width = square_area_side_length / cols
    rect_height = square_area_side_length / rows

    # Create the grid by placing rectangles at their respective positions.
    grid = []
    for row in range(rows):
        for col in range(cols):
            x = col * rect_width
            y = row * rect_height
            grid.append((x, y))

    return rect_width, rect_height, grid


def ccw(A,B,C):
    """
    Determines if three points are oriented in a counterclockwise direction.

    Parameters:
        A (tuple): Coordinates of point A (x, y).
        B (tuple): Coordinates of point B (x, y).
        C (tuple): Coordinates of point C (x, y).

    Returns:
        bool: True if the points are oriented counterclockwise, False otherwise.
    """
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def DelaunayGraph(points):
    tri = Delaunay(points)
    G = nx.Graph()
    
    # Iterate over triangles and add edges to the graph
    for simplex in tri.simplices:
        for i in range(len(simplex)):
            for j in range(i + 1, len(simplex)):
                G.add_edge(simplex[i], simplex[j])
                
    return G

def create_gabriel_graph(n, end_node_mode=None, padding_percent=0, node_number_step = 20):
    """
    Generates a Gabriel graph based on a set of points in a grid.

    Parameters:
        n (int): The number of nodes in the graph.
        end_node_mode (str): Mode for adding end nodes.
        padding_percent (float): Percentage of padding to apply to the grid.
        node_number_step (int) : number of node increase in each step, using 
                                 for proportional expand sim area.

    Returns:
        NetworkX Graph: The generated Gabriel graph.
    """
    x_padding = 0
    y_padding = 0
    
    graph_side_length = math.sqrt(n/node_number_step)*101.82
    x_length, y_length, grid = generate_grid(graph_side_length, n)
    
    if padding_percent > 0:
        x_padding = x_length*padding_percent/100
        y_padding = y_length*padding_percent/100
    
    x_coordinates = np.random.uniform(0+x_padding, x_length-x_padding, n)
    y_coordinates = np.random.uniform(0+y_padding, y_length-y_padding, n)
    
    #shift into grid
    for i, (x, y) in enumerate(grid):
        x_coordinates[i] += x
        y_coordinates[i] += y
    points = np.asarray([[x_coordinates[i], y_coordinates[i]] for i in range(len(x_coordinates))])
    
    pos_dict = dict(zip(list(range(n)),points.tolist())) #create a dict of node coordinates
    
    delaunayGraph = DelaunayGraph(points)
    voronoiDiagram = Voronoi(points)
    voronoiVertices = voronoiDiagram.vertices
    voronoiCenter = voronoiDiagram.points.mean(axis=0)
    ptp_bound = voronoiDiagram.points.ptp(axis=0)
    
    gabrielGraph = nx.Graph()
    for node in range(n):
        gabrielGraph.add_node(node)
    nx.set_node_attributes(gabrielGraph, pos_dict, name='pos')
    
    for u, v in delaunayGraph.edges():
        uRegion = set(voronoiDiagram.regions[voronoiDiagram.point_region[u]])
        vRegion = set(voronoiDiagram.regions[voronoiDiagram.point_region[v]])
        boundary = sorted(list(uRegion.intersection(vRegion)))[-2:]
        boundaryVertices = [None, voronoiVertices[boundary[1]]]
        
        if (boundary[0] == -1):
            tangent = points[u] - points[v]
            tangent /= np.linalg.norm(tangent)
            normal = np.array([-tangent[1], tangent[0]])
            
            midPoint = 0.5 * (points[u] + points[v])
            direction = np.sign(np.dot(midPoint - voronoiCenter, normal)) * normal
            farPoint = voronoiVertices[boundary[1]] + direction * ptp_bound.max()
            boundaryVertices[0] = farPoint
        else:
            boundaryVertices[0] = voronoiVertices[boundary[0]]
        
        if intersect(points[u], points[v], boundaryVertices[0], boundaryVertices[1]): 
            gabrielGraph.add_edge(u, v)
    

    gabrielGraph = finalize_graph_creation(gabrielGraph, end_node_mode, graph_side_length)
    
    return gabrielGraph

def create_waxman_graph(num_nodes, beta=0.4, alpha=0.1, L=None, end_node_mode=None, draw=False):
    """
    Creates a Waxman random graph.

    Parameters:
        num_nodes (int): The number of nodes in the graph.
        beta (float): Controls link density (0 < beta <= 1).
        alpha (float): Sensitivity of link formation to distance (0 < alpha <= 1).
        L (float or None): Maximum distance between nodes.
        end_node_mode (str or None): Mode for determining end nodes.
        draw (bool): Whether to draw the graph.

    Returns:
        networkx.Graph: The generated Waxman random graph.
    """
    G = nx.waxman_graph(n=num_nodes, beta=beta, alpha=alpha, L=L, domain=(0, 0, 1, 1), metric=None, seed=None)
    
    if (has_isolated_elements(G)):
        return None

    G = finalize_graph_creation(G, end_node_mode)

    if draw:
        draw_graph(G)
    return G

def has_isolated_elements(G):
    """
    Checks for isolated nodes and isolated graphs in the given graph.

    Parameters:
        G (networkx.Graph): The graph to check.

    Returns:
        bool: True if isolated nodes or isolated graphs are found, False otherwise.
    """
    
    isolates = list(nx.isolates(G))
    if len(isolates) > 0:
        return True
    # Check for isolated graphs
    if nx.number_connected_components(G) >1:
        return True
    
    return False #no isolated node and graph

def find_convex_hull(G, graph_side_length):
    """
    Finds the nodes forming the convex hull of the graph G.

    Parameters:
        G (NetworkX graph): Input graph.
        graph_side_length (float): Length of the graph side.

    Returns:
        list: List of nodes forming the convex hull.
    """
    end_nodes = []
    pos = nx.get_node_attributes(G, 'pos')
    
    hull = ConvexHull(np.array(list(pos.values()))) 
    for node in hull.vertices:
        end_nodes.append(str(node))
    
    return end_nodes


def find_4corners_node(G, graph_side_length):
    """
    Finds the four corner nodes of the graph G.

    Parameters:
        G (NetworkX graph): Input graph.
        graph_side_length (float): Length of the graph side.

    Returns:
        list: List of four corner nodes.
    """
    end_nodes = []
    second_end_nodes = []
    
    # LB
    corner_point = 2 * graph_side_length
    corner_node = ''
    second_corner_node = ''
    second_corner_point = 2 * graph_side_length
    
    for n, d in G.nodes.items():
        current_point = G.nodes[n]['pos'][0] + G.nodes[n]['pos'][1]
        if corner_point > current_point:
            if corner_node != '':
                second_corner_node = corner_node
                second_corner_point = corner_point
            corner_point = current_point
            corner_node = n
        elif corner_node != '' and second_corner_point > current_point:
                second_corner_point = current_point
                second_corner_node = n
            
    end_nodes.append(corner_node)
    second_end_nodes.append(second_corner_node)
    
    # LT
    corner_point = -1 * graph_side_length
    corner_node = ''
    second_corner_node = ''
    second_corner_point = -1 * graph_side_length
    
    for n, d in G.nodes.items():
        current_point = - G.nodes[n]['pos'][0] + G.nodes[n]['pos'][1]
        if corner_point < current_point:
            if corner_node != '':
                second_corner_node = corner_node
                second_corner_point = corner_point
            corner_point = current_point
            corner_node = n
        elif corner_node != '' and second_corner_point < current_point:
                second_corner_point = current_point
                second_corner_node = n
                
    end_nodes.append(corner_node)
    second_end_nodes.append(second_corner_node)
    
    # RT
    corner_point = 0
    corner_node = ''
    second_corner_node = ''
    second_corner_point = 0
    
    for n, d in G.nodes.items():
        current_point = G.nodes[n]['pos'][0] + G.nodes[n]['pos'][1]
        
        if corner_point < current_point:
            if corner_node != '':
                second_corner_node = corner_node
                second_corner_point = corner_point
            corner_point = current_point
            corner_node = n
        elif corner_node != '' and second_corner_point < current_point:
                second_corner_point = current_point
                second_corner_node = n
                
    end_nodes.append(corner_node)
    second_end_nodes.append(second_corner_node)
    
    # RB
    corner_point = -1 * graph_side_length
    corner_node = ''
    second_corner_node = ''
    second_corner_point = -1 * graph_side_length
    
    for n, d in G.nodes.items():
        current_point = G.nodes[n]['pos'][0] - G.nodes[n]['pos'][1]        
        if corner_point < current_point:
            if corner_node != '':
                second_corner_node = corner_node
                second_corner_point = corner_point
            corner_point = current_point
            corner_node = n
        elif corner_node != '' and second_corner_point < current_point:
                second_corner_point = current_point
                second_corner_node = n
                
    end_nodes.append(corner_node)
    second_end_nodes.append(second_corner_node)
    
    # Check if same corner node was selected, if so, step back to the 2nd candidate corner node
    for i, n_i in enumerate(end_nodes):
        for j, n_j in enumerate(end_nodes):
            if i != j and n_i == n_j:
                end_nodes[j] = second_end_nodes[j]
    
    return end_nodes

def finalize_graph_creation(G, end_node_mode, graph_side_length):
    """
    Assigns type and coordinates for nodes and changes node labels to strings.

    Parameters:
        G (NetworkX graph): Input graph.
        end_node_mode (str): Mode for selecting end nodes ("4corners" or "convexhull").
        graph_side_length (float): Length of the graph side.

    Returns:
        NetworkX graph: Updated graph with assigned types and coordinates.
    """
    # Convert node labels to strings
    label_remapping = {key: str(key) for key in G.nodes() if type(key) is not str}
    G = nx.relabel_nodes(G, label_remapping)
    
    # Assign type for all nodes
    for node in G.nodes():
        G.nodes[node]['type'] = 'repeater_node'
    
    # Determine end nodes based on selected mode
    if end_node_mode == "4corners":
        for node in find_4corners_node(G, graph_side_length):
            G.nodes[node]['type'] = 'end_node'
    elif end_node_mode == "convexhull":
        for node in find_convex_hull(G, graph_side_length):
            G.nodes[node]['type'] = 'end_node'
    else:
        end_nodes = []
    
    # Assign coordinates for all nodes
    for node in G.nodes():
        G.nodes[node]['xcoord'] = G.nodes[node]['pos'][0]
        G.nodes[node]['ycoord'] = G.nodes[node]['pos'][1]
    
    return G

def draw_graph(G, title="", chosen_edges=[], updated_edges=[], chosen_nodes=[], show_node_label = True, show_edge_label=True):
    """
    Draws the input graph.

    Parameters:
        G (NetworkX graph): Input graph.
        title (str): Title of the plot.
        edge_label (bool): Whether to display edge labels.

    Returns:
        None
    """
    node_pos = nx.get_node_attributes(G, 'pos')
    edge_pos = node_pos
    label_pos = pos_nudge(node_pos, 1, 1)
    length = nx.get_edge_attributes(G, 'length')
    
    none_type_node = []
    repeater_nodes = []
    end_nodes = []
    
    node_labels = {}
    end_node_labels = {}
    repeater_node_labels = {}
    chosen_node_labels = {}
    
    for node, nodedata in G.nodes.items():
        if not nodedata['type']:
            none_type_node.append(node)
            if show_node_label:
                node_labels[node] = node
        else:
            if nodedata['type'] == 'repeater_node':
                repeater_nodes.append(node)
                if show_node_label:
                    repeater_node_labels[node] = node
            elif nodedata['type'] == 'end_node':
                end_nodes.append(node)
                if show_node_label:
                    end_node_labels[node] = node
            else:
                NotImplementedError(f"Unknowed node type: {nodedata['type']}")
        
        if node in chosen_nodes:
            chosen_node_labels[node] = node
    
    fig, ax = plt.subplots(figsize=(15, 15))
    
    #draw nodes
    if end_nodes:
        end_nodes = nx.draw_networkx_nodes(G=G, pos=node_pos, nodelist=end_nodes, node_shape='s', node_size=750,
                                       node_color=[[1.0, 120 / 255, 0.]], label="End Node", linewidths=3)
        end_nodes.set_edgecolor(["k"])
        if end_node_labels:
            nx.draw_networkx_labels(G=G, pos=label_pos, labels=end_node_labels, font_size=20, 
                                font_weight="bold", font_color="b", font_family='serif')
    if repeater_nodes:
        rep_nodes = nx.draw_networkx_nodes(G=G, pos=node_pos, nodelist=repeater_nodes, node_size=750,
                                       node_color=[[1, 1, 1]], label="Repeater Node")
        rep_nodes.set_edgecolor(["k"])
        if repeater_node_labels:
            nx.draw_networkx_labels(G=G, pos=label_pos, labels=repeater_node_labels, font_size=12, font_weight="bold")
        
    if chosen_nodes:
        chosen_nodes = nx.draw_networkx_nodes(G=G, pos=node_pos, nodelist=chosen_nodes, node_size=750,
                                       node_color="#2CBAFA", label="Chosen Node")
        chosen_nodes.set_edgecolor(["k"])
        if chosen_node_labels and show_node_label:
            nx.draw_networkx_labels(G=G, pos=label_pos, labels=chosen_node_labels, font_size=12, 
                                    font_weight="bold", font_color="k", font_family='serif')
    
    #draw edges
    
    nx.draw_networkx_edges(G=G, pos=node_pos, width=1, node_size=750)
    if chosen_edges:
        nx.draw_networkx_edges(G=G, pos=node_pos, edgelist=chosen_edges, edge_color="#2CBAFA", width=3, node_size=750)
    if show_edge_label:
        nx.draw_networkx_edge_labels(G=G, pos=node_pos, edge_labels=length)
        
    ax.set_title(title, fontsize=20)   
    ax.set_aspect(1)
    ax.axis('equal') 
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    plt.axis('on')
    fig.tight_layout()
    plt.show()
    
def pos_nudge(pos, x_shift, y_shift):
    """
    Adjusts node positions by a given x and y shift.

    Parameters:
        pos (dict): Dictionary containing node positions.
        x_shift (float): Amount to shift node positions along the x-axis.
        y_shift (float): Amount to shift node positions along the y-axis.

    Returns:
        dict: Updated node positions.
    """
    if (x_shift == 1 and y_shift == 1):
        return pos
    else:
        return {n:(x*x_shift, y*y_shift) for n,(x,y) in pos.items()}    
    
def prune_edge_exceeding_lmax(G, lmax):
    """
    Removes edges from the graph exceeding the given maximum length.

    Parameters:
        G (NetworkX graph): Input graph.
        lmax (float): Maximum allowed edge length.

    Returns:
        NetworkX graph: Graph with edges pruned.
    """
    edges_to_remove = []
    for (i,j) in G.edges:
        if G[i][j]["length"] > lmax:
            edges_to_remove.append((i,j))
    for (i,j) in edges_to_remove:
        G.remove_edge(i, j)
    return G
    
def convert_bidirectional(G, is_logging=False):
    """
    Convert the graph to a bidirectional graph.

    Parameters:
        G (NetworkX graph): Input graph.
        is_logging (bool): Whether to log the conversion process.

    Returns:
        NetworkX graph: Bidirectional graph.
    """
    G = G.to_directed()
    for edge, edge_data in G.edges.items():
        chk_edge_start = edge[0]
        chk_edge_end = edge[1]
    
        for node, node_data in G.nodes.items():
            if node == chk_edge_start:
                chk_edge_start = node_data["old_label"]
            elif node == chk_edge_end:
                chk_edge_end = node_data["old_label"]  
        if edge_data['start'] == chk_edge_end and edge_data['end'] == chk_edge_start:
            if is_logging: 
                print(edge, "start:", edge_data['start'], " end:", edge_data['end'], " >> start:", chk_edge_start, " end:", chk_edge_end) 
            edge_data['start'] = chk_edge_start
            edge_data['end'] = chk_edge_end
    return G

def check_neighbour_length(G, end_nodes, lmax, is_logging=False):
    """
    Checks if all neighbors of end nodes have lengths less than lmax.

    Parameters:
        G (NetworkX graph): Input graph.
        end_nodes (list): List of end nodes.
        lmax (float): Maximum allowed edge length.
        is_logging (bool): Whether to log the checking process.

    Returns:
        bool: True if all neighbor lengths are less than lmax, False otherwise.
    """
    
    # Initialize the maximum length to zero
    max_length = 0
    num_possible_neighbour = 2

    # Iterate over each node in the graph
    for e_node in end_nodes:
        possible_neighbour = 0
        # Get the neighbors of the current node
        neighbors = G.neighbors(e_node)
        if is_logging: print(f"node:{e_node} ", end=">>")

        # Iterate over each neighbor
        for neighbor in neighbors:
            # Get the shortest path length between the current node and its neighbor
            length = G[e_node][neighbor]["length"]
            if is_logging: print(f"{neighbor}[{length}]",end=' ')

            # Return False if neighbors length more than lmax
            if length < lmax:
                possible_neighbour += 1
            
        if is_logging: print("")
        if possible_neighbour < num_possible_neighbour:
            return False
        
    return True

def sorted_end_node_pairs(G, end_nodes_pair):
    
    dist_end_nodes_pair = []
    
    for (u,v) in end_nodes_pair:
        dx = np.abs(G.nodes[u]['pos'][0] - G.nodes[v]['pos'][0])
        dy = np.abs(G.nodes[u]['pos'][1] - G.nodes[v]['pos'][1])
        dist = np.round(np.sqrt(np.square(dx) + np.square(dy)), 5)
        dist_end_nodes_pair.append(dist)
    
    ranks = rankdata(dist_end_nodes_pair)
    sorted_data = sorted(zip(end_nodes_pair, ranks), key=lambda x: x[1], reverse=True)
    
    return [x[0] for x in sorted_data]
    
def shortest_path_repeaters(graph, start, rep_node, lmax, is_logging=False):
    """
    Checks if all neighbors of end nodes have lengths less than lmax.

    Parameters:
        G (NetworkX graph): Input graph.
        end_nodes (list): List of end nodes.
        lmax (float): Maximum allowed edge length.
        is_log (bool): Whether to log the checking process.

    Returns:
        bool: True if all neighbor lengths are less than lmax, False otherwise.
    """
    priority_queue = [(0, start)] # Priority queue to store (distance, node) pairs
    edges_to_repeater = []
    edges_to_repeater_distance = {}
    founded_repeater = []
    distances = {node: float('infinity') for node in graph.nodes} # Dictionary to store the shortest distances from the start node
    predecessors = {node: None for node in graph.nodes}
    distances[start] = 0
    
    if is_logging: print(f"targeted: {start} to {rep_node}")
    while priority_queue:
        current_distance, current_node = heapq.heappop(priority_queue)
        if is_logging: print(f"pop {current_node}-{current_distance}")

        if current_distance > lmax: # Stop if the node type is rep_node or the distance exceeds lmax
            continue
        elif current_node in rep_node and current_node != start and current_node not in founded_repeater:
            nodes_to_repeater = []
            founded_repeater.append(current_node)
            traced_back_node = current_node
            edges_to_repeater_distance[current_node] = current_distance
            if is_logging: print(f"  current_node is {current_node}")
            # while traced_back_node is not None:# and traced_back_node:
            #     nodes_to_repeater.append(traced_back_node)
            #     edges_to_repeater.append((predecessors[traced_back_node], traced_back_node))
            #     traced_back_node = predecessors[traced_back_node]
            if is_logging: print(f"nodes to repeater {nodes_to_repeater}")
            continue
        else:        
            # Explore neighbors
            for neighbor in graph.neighbors(current_node):
                length = graph[current_node][neighbor].get('length')
                distance = current_distance + length

                # Update the distance if a shorter path is found
                if distance < distances[neighbor]:

                    if is_logging: print(f"  push {current_node}->{distance}->{neighbor}")
                    distances[neighbor] = distance
                    predecessors[neighbor] = current_node
                    heapq.heappush(priority_queue, (distance, neighbor))
    
    return edges_to_repeater_distance
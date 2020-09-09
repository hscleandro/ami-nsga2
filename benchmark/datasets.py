#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    module creating and initializing various datasets
    
    the created graphs are networkx graphs with the following additional characteristics:
    
    each node can have a weight which is specified with the attribute 'weight'
    each graph can be associated with ground truth information that can be
    a clustering or a specific module. The attribute storing the ground truth is specified
    with the graph attribute 'groups' 
    
"""

import pandas as pd
import networkx as nx
import os
from scipy.stats import truncnorm
from zipfile import ZipFile
from random import Random
import sys


from typing import Dict, Set
try:
    import toml
    TOML_INSTALLED = True
except ImportError:
    TOML_INSTALLED = False


# the directory where all the datasets are located
DATADIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")

class Datasets:

    @staticmethod
    def init_graph(G: nx.Graph, name="", node_weight='unchanged', edge_weight='unchanged', default_groups=None):
        """
        Initialize graph, nodes and edges attributes

        :param G: a graph
        :type G: networkx.Graph
        :param name: the name of the graph
        :type name: String
        :param node_weight: specify the weight of nodes
        :type node_weight: 'unchanged' or None for nodes with uniform weights equals to 1
        :param edge_weight: specify the weight of edges
        :type edge_weight: 'unchanged' or None for edges with uniform weights equals to 1
        :param default_groups: the node attribute representing the known group(s)
        :type default_groups: a string or None if the default_groups are unknown
        """

        mapping_needed = False
        # node labels must consist on consecutive integers starting at 0
        # if it is not the case, a mapping is necessary

        # convert nodes label to int if all labels are numbers and are
        # represented as strings

        # convert string labels to int values if possible
        try:
            if all([x.isdigit() for x in G.nodes]):
                mapping = {x: int(x) for x in G.nodes()}
                G = nx.relabel_nodes(G, mapping, copy=False)
        except AttributeError:
            pass

        mapping_needed = False

        if all(isinstance(n, int) for n in G.nodes()):
            if sorted([x for x in G.nodes()]) == list(range(G.number_of_nodes())):
                pass  # perfect
            else:
                mapping_needed = True
        else:
            mapping_needed = True

        if mapping_needed:
            nodes_mapping = dict(enumerate(sorted(G.nodes)))
            mapping = {y: x for x, y in nodes_mapping.items()}
            G = nx.relabel_nodes(G, mapping)
            for node in G.nodes:
                G.nodes[node]['label'] = nodes_mapping[node]
        for node in G.nodes:
            if 'weight' not in G.nodes[node] or node_weight == None:
                G.nodes[node]['weight'] = 1
        for edge in G.edges:
            if 'weight' not in G.edges[edge] or edge_weight == None:
                G.edges[edge]['weight'] = 1

        G.graph['default_groups'] = default_groups
        G.graph['nb_nodes'] = G.number_of_nodes()
        G.graph['nb_edges'] = G.number_of_edges()
        G.graph['name'] = name
        return G

    @staticmethod
    def get_groups(G: nx.Graph, groups_attribute=None) -> Dict[int, Set[int]]:
        """
        Return the groups identified in the graph.

        :param G: a graph
        :type G: networkx.Graph
        :param groups_attribute: the node attribute representing the known group(s)
        :type groups_attribute: a string of None to use default group
        """
        if not groups_attribute:
            groups_attribute = G.graph['default_groups']
        if not groups_attribute:
            return {}
        groups = {}
        for n in G.nodes:
            group_nb = int(G.nodes[n][groups_attribute])
            if group_nb == 0:
                continue
            try:
                groups[group_nb].add(n)
            except KeyError:
                groups[group_nb] = set([n])
        return groups

    @staticmethod
    def get_guyon_graph(nb: int) -> nx.Graph:
        zip_file = ZipFile(os.path.join(DATADIR, "guyon.zip"))
        df = pd.read_csv(zip_file.open(
            "pp1_" + str(nb) + ".csv"), sep=',', header=0)
        df = df.drop(df.columns[0], axis=1)
        A = df.values
        G = nx.from_numpy_matrix(A)
        df = pd.read_csv(zip_file.open("pvalues.csv"), sep=',', header=0)
        df = df.drop(df.columns[0], axis=1)
        for ctr in range(G.number_of_nodes()):
            G.nodes[ctr]['weight'] = 1 - df.iloc[ctr, nb - 1]

        df = pd.read_csv(zip_file.open("truehits.csv"), sep=',', header=0)
        df = df.drop(df.columns[0], axis=1)
        for ctr in range(G.number_of_nodes()):
            G.nodes[ctr]['truehit'] = df.iloc[ctr, nb - 1]

        return Datasets.init_graph(G, edge_weight=None, default_groups='truehit')

    @staticmethod
    def get_scale_free_graph(nb_nodes: int, nb_initial_nodes: int, nb_seeds: int, nb_to_select: int, p_prob: float,
                             q_prob: float, rng: Random = None) -> nx.Graph:
        """
        nb_nodes = 1000 # number of nodes
        nb_initial_nodes = 1 # number of initial nodes
        nb_seeds = 3 # umber of seed nodes
        nb_to_select = 10 # number of selected nodes in a group
        p_prob = 0.09 # probability to add an edge between existing nodes
        q_prob = 0.04 # Probability value of rewiring of existing edges
        """

        MODULE_SIZE = nb_to_select
        min_espace_between_cluster = 1

        G = nx.extended_barabasi_albert_graph(nb_nodes, nb_initial_nodes, p_prob, q_prob, rng)
        G = Datasets.connect_graph(G, rng)

        selected = list(G.nodes())
        no_select = []
        # n: is the seed vertex. If we have more than one group (nb_seeds > 1) we use most_distance() selection.
        # cluster: is all the vertices selected at the end of the execution
        # group: is each group selected in each interaction (each module).
        cluster = set()
        for i in range(nb_seeds):
            if i > 0:
                # most_distance() tries to identify a distant seed based on the other vertices already selected
                n = Datasets.most_distance(G, selected, cluster)
            else:
                n = rng.sample(selected, 1)[0]
            group = set()
            group.add(n)
            print("seed:" + str(n))

            while group.__len__() < MODULE_SIZE:
                # get_order() returns the next neighborhood level after the vertices selected in the module
                # if the neighborhood level 'nb_to_select' still contains vertices already selected previously
                # the function returns -1
                order = Datasets.get_order(G, n, nb_to_select, no_select)

                if order < 0:
                    print("ERRO [1]! The density and size of the network prevent you from creating a network with "
                          + str(nb_seeds) + " modules of size " + str(MODULE_SIZE) + " each.\n")
                    print("The modules need to be far away on the network!\n")
                    sys.exit(1)

                rnd = rng.random()
                # checks to neighborhood order according with the function 1/pow(10, dist)

                if 1 >= rnd > 0.1:
                    order = order + 1
                elif 0.01 < rnd <= 0.1:
                    order = order + 2
                elif 0.001 < rnd <= 0.01:
                    order = order + 3
                elif 0.0001 < rnd <= 0.001:
                    order = order + 4

                # neighbors_order returns all neighbors of the node n of a specific order
                population = Datasets.neighbors_order(G, n, order)

                # Removes neighbors and nodes from other groups already identified
                # besides the nodes of this same group that have already been chosen.
                population = population - set(no_select)

                ## Erro! When a neighborhood level 'order' is composed only of vertices already selected.
                if population.__len__() == 0:
                    print("ERRO [2]! The density and size of the network prevent you from creating a network with "
                          + str(nb_seeds) + " modules of size " + str(MODULE_SIZE) + " each.\n")
                    print("The modules need to be far away on the network!\n")
                    sys.exit(1)
                    # order = Datasets.get_order(G, n, nb_to_select, no_select)
                    # population = Datasets.neighbors_order(G, n, order)

                else:
                    # one node is chosen randomly in the population
                    sel = rng.sample(list(population), 1)[0]
                    group.add(sel)
                    no_select = list(set(no_select).union(group))
                    #nb_selected += 1

            print("group " + str(i + 1) + ":")
            print(group)

            cluster.update(group)
            no_select = set()
            for j in cluster:
                no_select = no_select.union(Datasets.neighborhood(G, j, min_espace_between_cluster))
            selected = list(set(selected) - set(no_select))

        # assign a uniform distribution for all nodes
        low, up = 0, 1
        for node in G.nodes:
            G.nodes[node]['weight'] = rng.uniform(low, up)
            G.nodes[node]['truehit'] = 0

        # assign a truncated normal distribution for each node in cluster
        mean, std = 0, 0.05
        a, b = (low - mean) / std, (up - mean) / std
        for x in cluster:
            G.nodes[x]['weight'] = 1 - truncnorm.rvs(
                loc=mean, scale=std, a=a, b=b)
            G.nodes[x]['truehit'] = 1

        return Datasets.init_graph(G, edge_weight=None, default_groups='truehit')

    ## Remove x% of edges
    @staticmethod
    def get_scale_free_graph_edge(nb_nodes: int, nb_initial_nodes: int, nb_seeds: int, nb_to_select: int, p_prob: float,
                                  q_prob: float, exc_edge, rng: Random = None) -> nx.Graph:
        """
        nb_nodes = 1000 # number of nodes
        nb_initial_nodes = 1 # number of initial nodes
        nb_seeds = 3 # umber of seed nodes
        nb_to_select = 10 # number of selected nodes in a group
        p_prob = 0.09 # probability to add an edge between existing nodes
        q_prob = 0.04 # Probability value of rewiring of existing edges
        exc_edge = 0.2 # Percent of excluded edges
        """

        G = Datasets.get_scale_free_graph(nb_nodes, nb_initial_nodes, nb_seeds, nb_to_select, p_prob, q_prob, rng)
        if exc_edge > 0:
            list_edges = list(G.edges)
            size = len(list_edges)
            x1 = size
            exc = round(size * exc_edge)
            for i in range(exc):
                exclude_edge = rng.sample(range(size), 1)[0]
                e = list_edges[exclude_edge]
                G.remove_edge(*e)
                list_edges = list(G.edges)
                size = len(list_edges)
            x2 = len(list(G.edges))

            G.nodes[0]['e'] = x1
            G.nodes[1]['e'] = x2

            G = Datasets.connect_graph(G, rng)

        return Datasets.init_graph(G, edge_weight=None, default_groups='truehit')


    ## Suport Functions

    def connect_graph(G: nx.Graph, rng: Random = None) -> nx.Graph:
        """
        Input: Graph containing some disconnected vertices
        Output: Graph connected in a single component

        """
        # Identifying vertices
        list_deg = list(G.nodes())

        # Obtaining vertices that compose the max connected coponent
        candidates = max(nx.connected_components(G), key=len)
        #candidates = max(Gc)

        # Identifying vertices that wil be reconnected in the network
        list_null = list(set(list_deg) - set(candidates))

        ## Add the isolate vertex in the max connected component
        for i in list_null:
            j = rng.sample(list(candidates), 1)[0]
            G.add_edge(i, j)

        return Datasets.init_graph(G)

    def neighbors_order(G, start, k):
        """
        Input: G: Network
               start: seed vertex
               k: Neigboors level from the seed
        Output: Neighbors level k from the seed vertex

        """
        nbrs = set([start])
        for l in range(k):
            nbrs = set((nbr for n in nbrs for nbr in G[n]))
        return nbrs

    def neighborhood(G, start, k):
        """
        Input: G: Network
               start: seed vertex
               k: Neigboors level from the seed
        Output: Neighborhood k from the seed vertex.
                For example: neighborhood(v,2) = {neighbors_leve1(v) U neighbors_leve2(v)}
        """
        nbrs = set([start])
        neighbors = set([start])
        for l in range(k):
            nbrs = set((nbr for n in nbrs for nbr in G[n]))
            neighbors.update(nbrs)
        return neighbors

    def most_distance(G, selected, cluster):
        """
        Input: G: Network
               selected: vertices that can be chosen as seed
               cluster: vertices already selected for other groups
        Output: Vertex contained in "selected" furthest from vertices contained in "cluster"
        """
        min_dist = 0
        choice_vertex = selected[0]
        for sel in selected:
            test_dist = 1000
            for cl in cluster:
                dist = nx.shortest_path_length(G, source=sel, target=cl)
                if dist < test_dist:
                    test_dist = dist
            if test_dist > min_dist:
                min_dist = test_dist
                choice_vertex = sel

        return choice_vertex

    def get_order(G, n, max_order, no_select):
        """
        Input: G: Network
               n: seed vertex
               max_order: max number of vertex in a module
               no_select: vertices that cannot be chosen because they have already been chosen or
                          are close to vertices already chosen
        Output: minimum order of neighborhood so that a vertex can be chosen.
                when the network contains vertices with low density and large groups,
                the neighborhood orders are being occupied very quickly.
        """
        order = 1
        population = Datasets.neighbors_order(G, n, order)
        population = population - set(no_select)
        while population.__len__() == 0 and order <= max_order:
            order += 1
            population = Datasets.neighbors_order(G, n, order)
            population = population - set(no_select)
        if order == max_order:
            return -1
        else:
            return order - 1
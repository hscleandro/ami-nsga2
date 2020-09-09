#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    various scoring methods
    
    A first group of functions computes the score of a set of nodes on a graph,
    optionnaly specifying the attribute that is used to store the nodes' weights
    
    A second group of functions measures the accuracy of a prediction given
    a predicted set of nodes and a set of nodes considered as ground truth 
"""

from sklearn import metrics
import networkx as nx
from typing import Set
import json
import time
import os
import itertools
import re

# get the hiperparameters for each tool
json_load = os.getcwd() + '/set_variables.json'
f = open(json_load)
data = json.load(f)

class Scores:

    @staticmethod
    def average_shortest_path(G: nx.Graph, cluster: Set):
        """compute average shortest path in modules composed of vertices belonging to cluster"""
        sum_shortest_path_length = 0
        for n1 in cluster:
            for n2 in cluster:
                if n1 == n2:
                    continue
                sum_shortest_path_length += nx.shortest_path_length(G,source=n1,target=n2)
        nb_nodes = len(cluster)
        return sum_shortest_path_length/(nb_nodes*(nb_nodes-1))
    
    @staticmethod
    def measure_f1(G: nx.Graph, real_cluster: Set, pred_cluster: Set) -> float:
        prediction = [0] * G.graph['nb_nodes']
        for x in pred_cluster:
            prediction[x] = 1
        real = [0] * G.graph['nb_nodes']
        for x in real_cluster:
            real[x] = 1
        return metrics.f1_score(real, prediction)

    @staticmethod
    def measure_precision(G: nx.Graph, real_cluster: Set, pred_cluster: Set) -> float:
        prediction = [0] * G.graph['nb_nodes']
        for x in pred_cluster:
            prediction[x] = 1
        real = [0] * G.graph['nb_nodes']
        for x in real_cluster:
            real[x] = 1
        return metrics.precision_score(real, prediction)

    @staticmethod
    def measure_recall(G: nx.Graph, real_cluster: Set, pred_cluster: Set) -> float:
        prediction = [0] * G.graph['nb_nodes']
        for x in pred_cluster:
            prediction[x] = 1
        real = [0] * G.graph['nb_nodes']
        for x in real_cluster:
            real[x] = 1
        return metrics.recall_score(real, prediction)

    @staticmethod
    def bionet_result(network_path: str, weight_path: str, fdrs: str, output: str):
        # Execute the bionet tool
        start_time = time.time()
        bionet_command = "Rscript " + os.getcwd() + "/tools/BIONET_test.R" + \
                         " -n " + network_path + \
                         " -w " + weight_path + \
                         " --fdrs " + fdrs + \
                         " -o " + output

        os.system(bionet_command)

        # get the bionet results
        bionet_json = output + 'bionet_result.json'
        f = open(bionet_json)
        end_time = time.time()
        final_time = end_time - start_time
        return json.load(f), final_time

    @staticmethod
    def cosine_result(network_path: str, weight_path: str, ninter: str,
                      popSize: str, minsize: str, output: str):
        start_time = time.time()
        # Execute the bionet tool
        cosine_command = "Rscript " + os.getcwd() + "/tools/COSINE_test.R" + \
                         " -n " + network_path + \
                         " -w " + weight_path + \
                         " --popSize " + popSize + \
                         " --min " + minsize + \
                         " --int " + ninter + \
                         " -o " + output

        os.system(cosine_command)

        # get the COSINE results
        cosine_json = output + 'cosine_result.json'
        f = open(cosine_json)
        end_time = time.time()
        final_time = end_time - start_time

        return json.load(f), final_time

    @staticmethod
    def giga_result(network_path: str, weight_path: str, minsize: str, output: str):
        start_time = time.time()
        # Execute the bionet tool
        giga_command = "Rscript " + os.getcwd() + "/tools/GIGA_test.R" + \
                         " -n " + network_path + \
                         " -w " + weight_path + \
                         " --th " + minsize + \
                         " -o " + output

        os.system(giga_command)

        # get the COSINE results
        giga_json = output + 'giga_result.json'
        f = open(giga_json)
        end_time = time.time()
        final_time = end_time - start_time

        return json.load(f), final_time

    @staticmethod
    def knode_result(network_path: str, weight_path: str, minsize: str, output: str):
        start_time = time.time()
        # Execute the bionet tool
        knode_command = "Rscript " + os.getcwd() + "/tools/KNODE_test.R" + \
                       " -n " + network_path + \
                       " -w " + weight_path + \
                       " --th " + minsize + \
                       " -o " + output

        os.system(knode_command)

        # get the COSINE results
        knode_json = output + 'giga_result.json'
        f = open(knode_json)

        end_time = time.time()
        final_time = end_time - start_time

        return json.load(f), final_time

    @staticmethod
    def aminsga2_result(network_path: str, weight_path: str, ngeneration: str, psize: str, minsize: str, maxsize: str,
                        offsize: str, chmut: str, delmut: str, tourn: str, dconnect: str, output: str):
        start_time = time.time()

        aminsga2_command = "Rscript " + re.sub('benchmark','',os.getcwd()) + "aminsga2.R" + \
                        " --net " + network_path + \
                        " --profile " + weight_path + \
                        " --ngeneration " + ngeneration + \
                        " --popsize " + psize + \
                        " --minsize " + minsize + \
                        " --maxsize " + maxsize + \
                        " --offspingsize " + offsize + \
                        " --changemutationtax " + chmut + \
                        " --delmutationtax " + delmut + \
                        " --tournament " + tourn + \
                        " --dconnection " + dconnect + \
                        " --output " + output

        os.system(aminsga2_command)

        # get the COSINE results
        aminsga2_json = output + 'aminsga2_result.json'
        f = open(aminsga2_json)

        end_time = time.time()
        final_time = end_time - start_time

        return json.load(f), final_time

    @staticmethod
    def baseline(G: nx.Graph, minsize: int):
        start_time = time.time()
        list_node = {}
        # Execute the bionet tool
        for node in G.nodes:
            d = {node:G.nodes[node]['weight']}
            list_node.update(d)

        new_list = {k: v for k, v in sorted(list_node.items(), key=lambda item: item[1])}
        end_time = time.time()
        final_time = end_time - start_time

        return list(dict(itertools.islice(new_list.items(), minsize)).keys()), final_time

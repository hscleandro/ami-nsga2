
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Script used for comparison (benchmark) between the tools for identification of active modules. 

    The structure of the artificial networks used in this benchmark is made using as  
    reference the study of (Barabási, 2000).  
    The parameters bellow Can be defined by the user.

    Density [-i];
    Network size[-s]; 
    Active module size [-m], 
    Number of executions [-r];
    Number of random removed edges (for null test estimation) [-d];
    Probability to add an edge between existing nodes [-p]; and 
    Probability value of rewiring of existing edges [-q]

    Use the parameter m = 0 to execute a network structure based in real biological STRING dabase 
    (combined_score > 700 and cooexpression > 0) from a approximation using linear regression.

    The signal variation of the simulated genes was calculated using the same strategy as the study
    made by (Robinson, 2017).

    Four other state-of-the-art strategies (Bionet, GIGA, COSINE, Knode), as well as the baseline 
    containing only the best values obtained by the signal variation, are analyzed. The 
    hyperparameters and the selection of each tool are defined in the set_variables.json file.

    Albert,  Réka  and  Albert-László  Barabási  (2000).  
    “Topology  of  evolving  networks:local events and universality”. 
    In:Physical review letters85.24, p. 5234.

    Robinson, Sean et al. (2017). 
    “Incorporating interaction networks into the determina-tion of 
     functionally related hit genes in genomic experiments with Markov ran-dom fields.” 
    In:Bioinformatics (Oxford, England)33.14, pp. i170–i179.ISSN: 1367-4811
    
"""

from datasets import Datasets
from scores import Scores
import argparse
import statistics
from random import Random
from sklearn.linear_model import LinearRegression
import numpy as np
from datetime import datetime
import networkx as nx
import os
import json

import networkx as nx

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--nbruns', dest='number_of_runs', type=int, required=False, default=1000, help="specifies the number of runs (default=1000)")
#   parser.add_argument('-n', '--nbmodules', dest='nb_modules', type=int, required=False, default=1, help="specifies the number of modules to generate (default=1)")
    parser.add_argument('-m', '--modulesize', dest='module_size', type=int, required=False, default=10, help="specifies the size of each modules (default=10)")
    parser.add_argument('-i', '--initialmod', dest='initial_module', type=int, required=False, default=1, help="specifies the number of initial nodes for the preferential attachment (default=1)")
    parser.add_argument('-s', '--networksize', dest='network_size', type=int, required=False, default=1000, help="specifies the number of vertices in graph (default=1000)")
    parser.add_argument('-p', '--probp', dest='prob_p', type=float, required=False, default=0.0, help=" Probability value for adding an edge between existing nodes  p + q < 1 (default = 0")
    parser.add_argument('-q', '--probq', dest='prob_q', type=float, required=False, default=0.0, help=" Probability value of rewiring of existing edges  p + q < 1  (default = 0")
    parser.add_argument('-d', '--deledges', dest='removed_edges', type=float, required=False, default=0.0, help="proportion of edges to remove (default = 0.0")
    parser.add_argument('-v', '--verbose', dest='verbose', required=False, action="store_true", default=False, help="displays results on screen")
    parser.add_argument('-o', '--outfile', dest='outfile', required=False, default=None, help="name of the output file (default=no output)")
    parser.add_argument('-t', '--outtype', dest='output_type', type=int, required=False, default=None, help=" [1] to salve sim networks in adjacency matrix; [2] to save sim network in edge list (default=no output)")
    parser.add_argument('-g', '--graphgen', dest='graph_generation', required=False, choices=["guyondata"], help="""Specify the data to use.
                        If 'guyondata' is chosen the predictions are done on the dataset specified in Guyon's paper and all other graph parameters are ignored.""")
    return parser.parse_args()

def write_network_edge_list(G, output, file_name="adj_teste"):

    network_name = output + file_name + '_connection_list.edgelist'
    weight_list = output + file_name + '_weight.csv'

    fh = open(network_name, 'wb')
    nx.write_edgelist(G, fh)

    with open(weight_list, 'w') as output_finalw:
        for node in G.nodes:
            outline = []
            outline.append(str(node)+","+str(float(G.nodes[node]['weight'])))
            output_finalw.write("".join(outline) + "\n")

    with open(output + file_name + '_hittrue.csv', 'w') as output_finalt:
        k = 0
        for node in G.nodes:
            outline = []
            if G.nodes[node]['truehit'] == 1:
                outline.append(str(k))
                output_finalt.write("".join(outline) + "\n")
            k +=1
    return network_name, weight_list

if __name__ == '__main__':
    arg = parse_arguments()
    rng = Random()
    nb_modules = 1
    number_activegenes = nb_modules * arg.module_size

    json_load = os.getcwd() + '/set_variables.json'
    f = open(json_load)
    data = json.load(f)

    if arg.initial_module == 0:
        ## Simulation based in the extraction of modules in the real STRING network
        # Number of vertex extracted from STRING
        x = np.array([1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]).reshape((-1, 1))
        # Rate of conncetion in each subnetwork size indicated above
        y = np.array([2.376260, 3.597987, 4.849775, 5.982308, 7.258233, 8.374246, 9.597230, 10.786962, 11.996532])
        # Linear function to predict the rate of connection in each network size (based in STRING)
        predit_model = LinearRegression()
        predit_model.fit(x, y)
        y_pred = predit_model.predict(np.array([arg.network_size]).reshape((-1, 1)))
    else:
        initial_module = arg.initial_module

    f1_scores = []

    if arg.outfile:
        str_date = datetime.today().isoformat('_')
        outfile = open(str(arg.outfile + "result_" + str_date + ".csv"), 'w')

        initial_head = "#graph, real_size, active_module_size, rate_connection, average_shortest_paths"
        if 'bionet' in data['tools']:
            initial_head = initial_head + ", th_bionet, size_result_bionet, f1_score_bionet, recall_bionet, runtime_bionet"
        if 'cosine' in data['tools']:
            initial_head = initial_head + ", th_cosine, size_result_cosine, f1_score_cosine, recall_cosine, runtime_cosine"
        if 'giga' in data['tools']:
            initial_head = initial_head + ", th_giga, size_result_giga, f1_score_giga, recall_giga, runtime_giga"
        if 'knode' in data['tools']:
            initial_head = initial_head + ", th_knode, size_result_knode, f1_score_knode, recall_knode, runtime_knode"
        if 'aminsga2' in data['tools']:
            initial_head = initial_head + ", th_aminsga2, size_result_aminsga2, f1_score_aminsga2, recall_aminsga2, runtime_aminsga2"

        initial_head = initial_head + ", th_baseline, f1_score_baseline, recall_baseline, runtime_baseline \n"

        outfile.write(initial_head)

    for ctr in range(arg.number_of_runs):
        if arg.graph_generation == "guyondata":
            G = Datasets.get_guyon_graph(ctr + 1)
        else:
            G = Datasets.get_scale_free_graph_edge(arg.network_size,
                                                   initial_module,
                                                   nb_modules,
                                                   arg.module_size,
                                                   arg.prob_p,
                                                   arg.prob_q,
                                                   arg.removed_edges,
                                                   rng)


        rate_conection = len(G.edges)/len(G.nodes)

        average_shortest_paths = []
        for _, cluster in Datasets.get_groups(G).items():
            nodes = list(cluster)
            average_shortest_paths.append(Scores.average_shortest_path(G, nodes))
        result = str(ctr) +","+ str(arg.network_size)+","+ str(arg.module_size) +","+ str(rate_conection) +","+ str(average_shortest_paths[0])
        # save the network and the simulation transcripts in the file
        network_path, weight_path = write_network_edge_list(G, arg.outfile, "sim" + str(ctr + 1))
        truehits = Datasets.get_groups(G)
        th = set([j for i in truehits.values() for j in i])

        # geting the result of others tools (verify set_variables.json file)
        if 'bionet' in data['tools']:
            print("\nExecuting Bionet\n")
            bionet_output, time = Scores.bionet_result(network_path,
                                                       weight_path,
                                                       str(data['bionet_fdrs']),
                                                       arg.outfile)

            th_bionet = len(th.intersection(bionet_output))
            len_bionet = len(bionet_output)
            f1_score_bionet = Scores.measure_f1(G, th, bionet_output)
            recall_bionet = Scores.measure_recall(G, th, bionet_output)
            runtime_bionet = time

            result = result +","+ str(th_bionet) +","+ str(len_bionet) +","+ str(f1_score_bionet) +","+ str(recall_bionet) +","+\
                     str(runtime_bionet)

        if 'cosine' in data['tools']:
            print("\nExecuting COSINE\n")
            cosine_output, time = Scores.cosine_result(network_path,
                                                       weight_path,
                                                       str(data['cosine_n_iter']),
                                                       str(data['cosine_popSize']),
                                                       str(number_activegenes),
                                                       arg.outfile)

            th_cosine = len(th.intersection(cosine_output))
            len_cosine = len(cosine_output)
            f1_score_cosine = Scores.measure_f1(G, th, set(cosine_output))
            recall_cosine = Scores.measure_recall(G, th, cosine_output)
            runtime_cosine = time

            result = result +","+ str(th_cosine) +","+ str(len_cosine) +","+str(f1_score_cosine) +","+ str(recall_cosine) +","+\
                     str(runtime_cosine)

        if 'giga' in data['tools']:
            print("\nExecuting GIGA\n")
            giga_output, time = Scores.giga_result(network_path,
                                                   weight_path,
                                                   str(number_activegenes),
                                                   arg.outfile)

            th_giga = len(th.intersection(giga_output))
            len_giga = len(giga_output)
            f1_score_giga = Scores.measure_f1(G, th, giga_output)
            recall_giga = Scores.measure_recall(G, th, giga_output)
            runtime_giga = time

            result = result +","+ str(th_giga) +","+ str(len_giga) +","+ str(f1_score_giga) +","+ str(recall_giga) +","+\
                str(runtime_giga)

        if 'knode' in data['tools']:
            print("\nExecuting Knode\n")
            knode_output, time = Scores.knode_result(network_path,
                                                     weight_path,
                                                     str(number_activegenes),
                                                     arg.outfile)

            th_knode = len(th.intersection(knode_output))
            len_knode = len(knode_output)
            f1_score_knode = Scores.measure_f1(G, th, knode_output)
            recall_knode = Scores.measure_recall(G, th, knode_output)
            runtime_knode = time

            result = result +","+ str(th_knode) +","+ str(len_knode) +","+ str(f1_score_knode) +","+ str(recall_knode) +","+\
                str(runtime_knode)

        if 'aminsga2' in data['tools']:
            print("\nExecuting AMINSGA2\n")

            aminsga2_output, time = Scores.aminsga2_result(network_path,
                                                           weight_path,
                                                           str(data['aminsga2_number_of_generations']),
                                                           str(data['aminsga2_population_size']),
                                                           str(number_activegenes),
                                                           str(data['aminsga2_upper_range']),
                                                           str(data['aminsga2_offspring_size']),
                                                           str(data['aminsga2_change_mutation_tax']),
                                                           str(data['aminsga2_deletion_mutation_tax']),
                                                           str(data['aminsga2_tournament']),
                                                           str(data['aminsga2_d_connection']),
                                                           arg.outfile)

            th_aminsga2 = len(th.intersection(aminsga2_output))
            len_aminsga2 = len(aminsga2_output)
            f1_score_aminsga2 = Scores.measure_f1(G, th, aminsga2_output)
            recall_aminsga2 = Scores.measure_recall(G, th, aminsga2_output)
            runtime_aminsga2 = time

            result = result +","+ str(th_aminsga2) +","+ str(len_aminsga2) +","+ str(f1_score_aminsga2) +","+ str(recall_aminsga2) +","+\
                str(runtime_aminsga2)

        pvalue_baseline, time = Scores.baseline(G,number_activegenes)

        th_baseline = len(th.intersection(pvalue_baseline))
        f1_score_baseline = Scores.measure_f1(G, th, pvalue_baseline)
        recall_baseline = Scores.measure_recall(G, th, pvalue_baseline)
        runtime_baseline = time

        result = result +","+ str(th_baseline) +","+ str(f1_score_baseline) +","+ str(recall_baseline) +","+\
            str(runtime_baseline)

        if arg.outfile:
            outfile.write("".join(result) + "\n")
            outfile.flush()
#
    if arg.outfile:
        outfile.close()

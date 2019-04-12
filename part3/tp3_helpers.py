import networkx as nx
import numpy as np
import random
from collections import defaultdict

def generate_random_graph(cluster_sizes, p = 0.3, q = 0.1):
    graph = nx.Graph()
    num_cluster = len(cluster_sizes)
    num_node = sum(cluster_sizes)
    start_idx = 0
    label = 0
    for size in cluster_sizes:
        graph.add_nodes_from(list(np.arange(start_idx, start_idx+size)), community=label)
        label += 1
        start_idx += size
    for i in range(num_node):
        for j in range(i+1, num_node):
            x = np.random.uniform()
            if graph.nodes[i]["community"] == graph.nodes[j]["community"] and x<p:
                graph.add_edge(i,j, type="in")
            elif graph.nodes[i]["community"] != graph.nodes[j]["community"] and x<q:
                graph.add_edge(i,j, type="out")  
    return graph

def label_propagation(graph, max_iter=50, num_community=4):
    for idx, node in enumerate(graph.nodes):
        graph.nodes[node]["label"] = idx
    
    for i in range(max_iter):
        nodes = [node for node in graph.nodes]
        random.seed(123)
        random.shuffle(nodes)
        for node in nodes:
            neighbors = graph.neighbors(node)
            scores = {}
            for neighbor in neighbors:
                if graph.nodes[neighbor]["label"] in scores:
                    scores[graph.nodes[neighbor]["label"]] = scores[graph.nodes[neighbor]["label"]] + 1
                else:
                    scores[graph.nodes[neighbor]["label"]] = 1
            top = [key for key,val in scores.items() if val == max(scores.values())]
            new_label = random.sample(top,1)[0]
            graph.nodes[node]["label"] = new_label
        labels = [graph.nodes[node]["label"] for node in graph.nodes]
        curr_label_count = len(set(labels))
        if curr_label_count == num_community:
            break
    return graph

def new_algorithm(graph, max_size=2000, num_community=4):
    edge_list = []
    for i in range(1):
        edge_list += [edge for edge in graph.edges]
    random.seed(123)
    random.shuffle(edge_list)
    
    deg = defaultdict(int)
    com_size = defaultdict(int)
    labels = defaultdict(int)
    label = 1
    for i in range(1):
        random.shuffle(edge_list)
        for n1, n2 in edge_list:
            if labels[n1] == 0:
                labels[n1] = label
                label += 1
            if labels[n2] == 0:
                labels[n2] = label
                label += 1

            deg[n1] += 1
            deg[n2] += 1
            com_size[labels[n1]] += 1
            com_size[labels[n2]] += 1

            if com_size[labels[n1]] <= max_size and com_size[labels[n2]] <= max_size:
                if com_size[labels[n1]] <= com_size[labels[n2]]:
                    com_size[labels[n2]] += deg[n1]
                    com_size[labels[n1]] -= deg[n1]
                    labels[n1] = labels[n2]
                else:
                    com_size[labels[n1]] += deg[n2]
                    com_size[labels[n2]] -= deg[n2]
                    labels[n2] = labels[n1]
        if i>=1 and len(set(labels.values())) <= num_community:
            break
    for node in graph.nodes:
        graph.nodes[node]["label"] = labels[node]
    return graph
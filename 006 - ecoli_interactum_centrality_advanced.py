import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import igraph as ig
from scipy import stats
from random import sample
import numpy as np
from scipy.stats import ttest_ind

# v1 = np.random.normal(size=100)
# v2 = np.random.normal(size=100)
#
#
# res = ttest_ind(v1, v2).pvalue
#
# print(res)

def pairwiseDisconnectivityIndex(g, v=0):
    """
    g: Graph
    v: Vertex, default = 0
    return: The pairwise disconnectivity of vertex v
    """
    N0 = 0
    vertices = g.vcount()
    for i in range(vertices):  # count number of ordered pairs of vertices
        for j in range(vertices):
            if (i != j):
                if (g.vertex_disjoint_paths(i, j,
                                            neighbors="ignore") != 0):  # If there is a path between vertex i and j
                    N0 += 1

    aux = g.copy()
    aux.delete_vertices([v])
    Nv = 0
    nvertices = aux.vcount()
    for i in range(nvertices):  # count number of ordered pairs of vertices
        for j in range(nvertices):
            if (i != j):
                if (aux.vertex_disjoint_paths(i, j,
                                              neighbors="ignore") != 0):  # If there is a path between vertex i and j
                    Nv += 1

    return (N0 - Nv) / N0


try:

    # Read the csv data
    df_all = pd.read_csv('data/511145.protein.links.full.v11.5.txt', sep=' ')
    print("Row count of raw data:  {}".format(df_all.shape[0]))

    df_deg_data = pd.read_csv('data/511145.protein.aliases.v11.5.txt', sep='\t', usecols=[0])
    print("Row count of raw data:  {}".format(df_deg_data.shape[0]))
    df_deg_data = df_deg_data.drop_duplicates()
    print("Row count of non duplicate raw data:  {}".format(df_deg_data.shape[0]))

    # Print first 10 rows
    print("Display first 10 rows of raw data")
    print(df_all.head(10))
    # Print all column names
    print("Display all column names")
    print(df_all.columns)

    print("Display first 10 rows of raw data")
    print(df_deg_data.head(10))
    print("Display all column names")
    print(df_deg_data.columns)

    # filter dataset to include only data with > 700 combined_score
    df_filtered = df_all[df_all['combined_score'] > 700]

    # row count of filtered rows
    print("Row count of filtered data combined_score > 700:  {}".format(df_filtered.shape[0]))

    df_working_data = df_filtered[['protein1', 'protein2', 'combined_score']]
    print("Display working data column names")
    print(df_working_data.columns)

    G = nx.from_pandas_edgelist(df_working_data, "protein1", "protein2", ["combined_score"])

    labels = df_deg_data['#string_protein_id'].values.tolist()
    print(labels)
    nx.set_node_attributes(G, labels, "NODETYPE")
    labels.append("ESSENTIAL")

    print(G.nodes['511145.b0001']["NODETYPE"])
    selected_nodes = [n for n, v in G.nodes(data=True) if v['NODETYPE'] == 'ESSENTIAL']
    print(selected_nodes)

    # Computing centrality
    print('Computing degree of centrality')
    degCent = nx.degree_centrality(G)
    # print('Degree Centrality is {}'.format(degCent))

    # Descending order sorting centrality
    # degCent_sorted = dict(sorted(degCent.items(), key=lambda item: item[1], reverse=True))
    # print(degCent_sorted)

    # degCent_sorted_top_10_items = itertools.islice(degCent_sorted.items(), 0, 10)
    #
    # for key, value in degCent_sorted_top_10_items:
    #     print(key, value)

    # # Computing closeness
    # print('Computing closeness centrality')
    # closeCent = nx.closeness_centrality(G)
    # print('Closeness Centrality is {}'.format(closeCent))
    # #
    # # # Computing betweeness
    # print('Computing betweenness centrality')
    # betCent = nx.betweenness_centrality(G)
    # print('Betweenness Centrality is {}'.format(betCent))
    # #
    # # compute eigenvector centrality
    # print('Computing eigenvector centrality')
    # eigenCent = nx.eigenvector_centrality(G)
    # print(eigenCent)

    # compute assortavity index

    # remove randomly selected nodes (to make example fast)
    # num_to_remove = int(len(G) / 1.5)
    # nodes = sample(list(G.nodes), num_to_remove)
    # G.remove_nodes_from(nodes)
    #
    # # remove low-degree nodes
    # low_degree = [n for n, d in G.degree() if d < 10]
    # G.remove_nodes_from(low_degree)

    # largest connected component
    # components = nx.connected_components(G)
    # largest_component = max(components, key=len)
    # H = G.subgraph(largest_component)

    # convert to igraph
    # h = ig.Graph.from_networkx(G)
    # test_value = pairwiseDisconnectivityIndex(h)
    #
    # print('pairwisedisconnectivityindex is {}'.format(test_value))

    # # compute centrality
    # centrality = nx.betweenness_centrality(H, k=10, endpoints=True)
    #
    # # compute community structure
    # lpc = nx.community.label_propagation_communities(H)
    # community_index = {n: i for i, com in enumerate(lpc) for n in com}
    #
    # #### draw graph ####
    # fig, ax = plt.subplots(figsize=(20, 15))
    # pos = nx.spring_layout(H, k=0.15, seed=4572321)
    # node_color = [community_index[n] for n in H]
    # node_size = [v * 20000 for v in centrality.values()]
    # nx.draw_networkx(
    #     H,
    #     pos=pos,
    #     with_labels=False,
    #     node_color=node_color,
    #     node_size=node_size,
    #     edge_color="gainsboro",
    #     alpha=0.4,
    # )
    #
    # # Title/legend
    # font = {"color": "k", "fontweight": "bold", "fontsize": 20}
    # ax.set_title("E. coli protein interactome", font)
    # # Change font color for legend
    # font["color"] = "r"
    #
    # ax.text(
    #     0.80,
    #     0.10,
    #     "node color = community structure",
    #     horizontalalignment="center",
    #     transform=ax.transAxes,
    #     fontdict=font,
    # )
    # ax.text(
    #     0.80,
    #     0.06,
    #     "node size = betweeness centrality",
    #     horizontalalignment="center",
    #     transform=ax.transAxes,
    #     fontdict=font,
    # )
    #
    # r = nx.degree_assortativity_coefficient(G)
    # print(f"{r:3.1f}")
    #
    # # Resize figure for label readibility
    # ax.margins(0.1, 0.05)
    # fig.tight_layout()
    # plt.axis("off")
    # plt.show()

    print('Program exited successfully')
    # Descending order sorting betweeness
    # betCent_sorted = dict(sorted(betCent.items(), key=lambda item: item[1], reverse=True))
    #
    # betCent_sorted_top_10_items = itertools.islice(betCent_sorted.items(), 0, 10)
    #
    # for key, value in betCent_sorted_top_10_items:
    #     print(key, value)

#     # G = nx.Graph([(1, 2), (1, 3), (1, 5), (3, 4), (4, 5)])
#     # G = nx.read_edgelist("data/edge_list.txt", nodetype=str)
#     #
#     num_nodes = G.number_of_nodes()
#     print('Number of Nodes {}'.format(num_nodes))
#     num_edges = G.number_of_edges()
#     print('Number of Edges {}'.format(num_edges))
#     #
#     # print(nx.info(G))
#     # subax1 = plt.subplot(121)
#     # nx.draw(G, with_labels=True, font_weight='bold')
#     # fig, ax = plt.subplots(1, 1, figsize=(8, 6));
#     # nx.draw_networkx(G, ax=ax)
#     # #
#     # Computing centrality
#     degCent = nx.degree_centrality(G)
#     # Descending order sorting centrality
#     degCent_sorted = dict(sorted(degCent.items(), key=lambda item: item[1], reverse=True))
#
#     # Computing betweeness
#     # betCent = nx.betweenness_centrality(G, normalized=True, endpoints=True)
#     #
#     # # Descending order sorting betweeness
#     # betCent_sorted = dict(sorted(betCent.items(), key=lambda item: item[1], reverse=True))
#
#     # Color for regular nodes
#     N_nodes = 50
#     color_list = N_nodes * ['lightsteelblue']
#
#     # Getting indices on top 10 nodes for each measure
#     N_top = 10
#     colors_top_10 = ['tab:orange', 'tab:blue', 'tab:green', 'lightsteelblue']
#     keys_deg_top = list(degCent_sorted)[0:N_top]
#     # keys_bet_top = list(betCent_sorted)[0:N_top]
#
#     # Computing centrality and betweeness intersection
#     # inter_list = list(set(keys_deg_top) & set(keys_bet_top))
#     # inter_list = list(set(keys_deg_top))
#     #
#     # # Setting up color for nodes
#     # for i in inter_list:
#     #     color_list[i] = colors_top_10[1]
#     #
#     # for i in range(N_top):
#     #     if keys_deg_top[i] not in list(set(keys_deg_top)):
#     #         color_list[keys_deg_top[i]] = colors_top_10[0]
#         # if keys_bet_top[i] not in inter_list:
#         #     color_list[keys_bet_top[i]] = colors_top_10[1]
#
#     # Draw graph
#     pos = nx.circular_layout(G)
#     nx.draw(G, pos, with_labels=True)
#
#     # Setting up legend
#     labels = ['Top 10 deg cent', 'Top 10 bet cent', 'Top 10 deg and bet cent', 'no top 10']
#     # for i in range(len(labels)):
#     #     plt.scatter([], [], label=labels[i])
#     plt.legend(loc='center')

# selected_nodes = [n for n,v in G.nodes(data=True) if v['since'] == 'December 2008']
# print (selected_nodes)
#
# selected_edges = [(u,v) for u,v,e in G.edges(data=True) if e['since'] == 'December 2008']
# print (selected_edges)

except Exception as e:
    print(e)
    print('Exception occurred')
#
# plt.show()

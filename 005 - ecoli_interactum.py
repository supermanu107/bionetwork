import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import itertools



try:

    # Read the csv data
    df_all = pd.read_csv('data/EColiProteinDataset.csv')
    print("Row count of raw data:  {}".format(df_all.shape[0]))

    # Print first 10 rows
    print("Display first 10 rows of raw data")
    print(df_all.head(10))

    # Print all column names
    print("Display all column names")
    print(df_all.columns)

    # filter dataset to include only data with > 700 combined_score
    df_filtered = df_all.loc[df_all['combined_score'] > 700]

    # row count of filtered rows
    print("Row count of filtered data combined_score > 700:  {}".format(df_filtered.shape[0]))

    df_working_data = df_filtered[['protein1', 'protein2', 'combined_score']]
    print("Display working data column names")
    print(df_working_data.columns)

    G = nx.from_pandas_edgelist(df_working_data, "protein1", "protein2", ["combined_score"])

    # Computing centrality
    print('Computing degree of centrality')
    degCent = nx.degree_centrality(G)

    # Descending order sorting centrality
    degCent_sorted = dict(sorted(degCent.items(), key=lambda item: item[1], reverse=True))
    # print(degCent_sorted)

    degCent_sorted_top_10_items = itertools.islice(degCent_sorted.items(), 0, 10)

    for key, value in degCent_sorted_top_10_items:
        print(key, value)

    # Computing closeness
    print('Computing closeness centrality')
    closeCent = nx.closeness_centrality(G)
    print('Closeness of Centrality is {}'.format(closeCent))

    # Computing betweeness
    print('Computing betweenness centrality')
    betCent = nx.betweenness_centrality(G)
    print(betCent)

    # compute eigenvector centrality
    print('Computing eigenvector centrality')
    eigenCent = nx.eigenvector_centrality(G)
    print(eigenCent)

    # compute assortavity index


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
except:
    print('Exception occurred')
#
# plt.show()

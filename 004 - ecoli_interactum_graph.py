import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
# Import library
from d3graph import d3graph, vec2adjmat
import numpy as np

# print(pd.options.display.max_rows)

G = nx.karate_club_graph()
adjmat = nx.adjacency_matrix(G).todense()
adjmat=pd.DataFrame(index=range(0,adjmat.shape[0]), data=adjmat, columns=range(0,adjmat.shape[0]))
adjmat.columns=adjmat.columns.astype(str)
adjmat.index=adjmat.index.astype(str)

# Make the dataframe
df = pd.DataFrame(index=adjmat.index)

# Add some columns. Note that columns that start with: node_ are removed from the information.

df['degree']=np.array([*G.degree()])[:,1]
df['other info']=np.array([*G.degree()])[:,1]

node_color = []
for i in range(0,len(G.nodes)):
    node_color.append(G.nodes[i]['club'])
    node_name=node_color
df['name']=node_name

node_size = df.degree.values*2

# Make some graphs
out = d3graph(adjmat, df=df, node_color=node_size, node_size=node_size)

out.show()


#
#
# # Set source and target nodes
# # source = ['node A','node F','node B','node B','node B','node A','node C','node Z']
# # target = ['node F','node B','node J','node F','node F','node M','node M','node A']
# # weight = [5.56, 0.5, 0.64, 0.23, 0.9, 3.28, 0.5, 0.45]
# #
# # # Create adjacency matrix
# # adjmat = vec2adjmat(source, target, weight=weight)
#
# # target  node A  node B  node F  node J  node M  node C  node Z
# # source
# # node A    0.00     0.0    5.56    0.00    3.28     0.0     0.0
# # node B    0.00     0.0    1.13    0.64    0.00     0.0     0.0
# # node F    0.00     0.5    0.00    0.00    0.00     0.0     0.0
# # node J    0.00     0.0    0.00    0.00    0.00     0.0     0.0
# # node M    0.00     0.0    0.00    0.00    0.00     0.0     0.0
# # node C    0.00     0.0    0.00    0.00    0.50     0.0     0.0
# # node Z    0.45     0.0    0.00    0.00    0.00     0.0     0.0
#
#
# #
# try:
#     df = pd.read_csv('data/EColiProteinDataset.csv', index_col=0)
#     print(df.head(10))
#     print(df.columns)
# #     # print(df.keys())
# #     # print(df.to_string())
# #
# #     # print(df.head(10))
# #     # print(df.shape[0])
# #     # print(df.iloc[:, 15])
# #
#     df2 = df.loc[df['combined_score'] > 700]
#     print(df2.shape[0])
#
#     df3 = df2[['protein2', 'neighborhood', 'combined_score']]
#
#     df4 = pd.concat([df3])
#
#
#     # print(df2['combined_score'])
#
#     d3 = d3graph(charge=1000)
#
#     # # print(df2['protein2'].to_numpy())
#     #
#     # # Create adjacency matrix
#     # adjmat = vec2adjmat(df3['protein2'].to_numpy(), df3['neighborhood'].to_numpy(), weight=df3['combined_score'].to_numpy())
#     #
#     # # Build force-directed graph with default settings
#     # d3.graph(adjmat)
#     # d3.show()
# #
#     G = nx.from_pandas_edgelist(df2, "protein2", "neighborhood", ["combined_score"])
#     adjmat = nx.adjacency_matrix(G).todense()
#     d3.graph(adjmat)
#     # d3.set_node_properties(size=df3['combined_score'].to_numpy())
#     d3.show(showfig=True)
#
# #
# #     # G = nx.Graph([(1, 2), (1, 3), (1, 5), (3, 4), (4, 5)])
# #     # G = nx.read_edgelist("data/edge_list.txt", nodetype=str)
# #     #
# #     num_nodes = G.number_of_nodes()
# #     print('Number of Nodes {}'.format(num_nodes))
# #     num_edges = G.number_of_edges()
# #     print('Number of Edges {}'.format(num_edges))
# #     #
# #     # print(nx.info(G))
# #     # subax1 = plt.subplot(121)
# #     # nx.draw(G, with_labels=True, font_weight='bold')
# #     # fig, ax = plt.subplots(1, 1, figsize=(8, 6));
# #     # nx.draw_networkx(G, ax=ax)
# #     # #
# #     # Computing centrality
# #     degCent = nx.degree_centrality(G)
# #     # Descending order sorting centrality
# #     degCent_sorted = dict(sorted(degCent.items(), key=lambda item: item[1], reverse=True))
# #
# #     # Computing betweeness
# #     # betCent = nx.betweenness_centrality(G, normalized=True, endpoints=True)
# #     #
# #     # # Descending order sorting betweeness
# #     # betCent_sorted = dict(sorted(betCent.items(), key=lambda item: item[1], reverse=True))
# #
# #     # Color for regular nodes
# #     N_nodes = 50
# #     color_list = N_nodes * ['lightsteelblue']
# #
# #     # Getting indices on top 10 nodes for each measure
# #     N_top = 10
# #     colors_top_10 = ['tab:orange', 'tab:blue', 'tab:green', 'lightsteelblue']
# #     keys_deg_top = list(degCent_sorted)[0:N_top]
# #     # keys_bet_top = list(betCent_sorted)[0:N_top]
# #
# #     # Computing centrality and betweeness intersection
# #     # inter_list = list(set(keys_deg_top) & set(keys_bet_top))
# #     # inter_list = list(set(keys_deg_top))
# #     #
# #     # # Setting up color for nodes
# #     # for i in inter_list:
# #     #     color_list[i] = colors_top_10[1]
# #     #
# #     # for i in range(N_top):
# #     #     if keys_deg_top[i] not in list(set(keys_deg_top)):
# #     #         color_list[keys_deg_top[i]] = colors_top_10[0]
# #         # if keys_bet_top[i] not in inter_list:
# #         #     color_list[keys_bet_top[i]] = colors_top_10[1]
# #
# #     # Draw graph
# #     pos = nx.circular_layout(G)
# #     nx.draw(G, pos, with_labels=True)
# #
# #     # Setting up legend
# #     labels = ['Top 10 deg cent', 'Top 10 bet cent', 'Top 10 deg and bet cent', 'no top 10']
# #     # for i in range(len(labels)):
# #     #     plt.scatter([], [], label=labels[i])
# #     plt.legend(loc='center')
# except:
#     print('Exception occured')
# #
# # plt.show()
# # Initialize
#
#
# #
# # import networkx as nx
# # import pandas as pd
# # from d3graph import d3graph
# #
# # G = nx.karate_club_graph()
# # adjmat = nx.adjacency_matrix(G).todense()
# # adjmat=pd.DataFrame(index=range(0,adjmat.shape[0]), data=adjmat, columns=range(0,adjmat.shape[0]))
# # adjmat.columns=adjmat.columns.astype(str)
# # adjmat.index=adjmat.index.astype(str)
# # adjmat.iloc[3,4]=5
# # adjmat.iloc[4,5]=6
# # adjmat.iloc[5,6]=7
# #
# # from tabulate import tabulate
# # print(tabulate(adjmat.head(), tablefmt="grid", headers="keys"))
# #
# # df = pd.DataFrame(index=adjmat.index)
# # df['degree']=np.array([*G.degree()])[:,1]
# # df['other info']=np.array([*G.degree()])[:,1]
# # node_size=df.degree.values*2
# # node_color=[]
# # for i in range(0,len(G.nodes)):
# #     node_color.append(G.nodes[i]['club'])
# #     node_name=node_color
# #
# # # Make some graphs
# # d3 = d3graph()
# #
# # d3.graph(adjmat)
# # d3.set_node_properties(color=node_color, cmap='Set1')
# # d3.show()
# #
# # d3.set_node_properties(label=node_name, color=node_color, cmap='Set1')
# # d3.show()
# #
# # d3.set_node_properties(adjmat, size=node_size)
# # d3.show()
# #
# # d3.set_node_properties(color=node_size, size=node_size)
# # d3.show()
# #
# # d3.set_edge_properties(edge_distance=100)
# # d3.set_node_properties(color=node_size, size=node_size)
# # d3.show()
# #
# # d3 = d3graph(charge=1000)
# # d3.graph(adjmat)
# # d3.set_node_properties(color=node_size, size=node_size)
# # d3.show()
# #
# # d3 = d3graph(collision=1, charge=250)
# # d3.graph(adjmat)
# # d3.set_node_properties(color=node_name, size=node_size, edge_size=node_size, cmap='Set1')
# # d3.show()
# #
# # d3 = d3graph(collision=1, charge=250)
# # d3.graph(adjmat)
# # d3.set_node_properties(color=node_name, size=node_size, edge_size=node_size, edge_color='#00FFFF', cmap='Set1')
# # d3.show()

import networkx as nx
import matplotlib.pyplot as plt

G = nx.Graph([(1, 2), (1, 3), (1, 5), (3, 4), (4, 5)])

num_nodes = G.number_of_nodes()
print('Number of Nodes {}'.format(num_nodes))
num_edges = G.number_of_edges()
print('Number of Edges {}'.format(num_edges))

subax1 = plt.subplot(121)
nx.draw(G, with_labels=True, font_weight='bold')

degree_c = nx.degree_centrality(G)
print('Degrees of Centrality is {}'.format(degree_c))

closeness_c = nx.closeness_centrality(G)
print('Closeness of Centrality is {}'.format(closeness_c))

# compute centrality
betweenness_centrality = nx.betweenness_centrality(G)
print('Betweenness of Centrality is {}'.format(betweenness_centrality))

plt.show()

import networkx as nx
import matplotlib.pyplot as plt

G = nx.Graph([(1, 2), (1, 3), (1, 5), (3, 4), (4, 5)])

num_nodes = G.number_of_nodes()
print('Number of Nodes {}'.format(num_nodes))
num_edges = G.number_of_edges()
print('Number of Edges {}'.format(num_edges))

# subax1 = plt.subplot(121)
# nx.draw(G, with_labels=True, font_weight='bold')

degree_c = nx.degree_centrality(G)
print('Degrees of Centrality is {}'.format(degree_c))

closeness_c = nx.closeness_centrality(G)
print('Closeness of Centrality is {}'.format(closeness_c))

# largest connected component
components = nx.connected_components(G)
largest_component = max(components, key=len)
H = G.subgraph(largest_component)

# compute centrality
centrality = nx.betweenness_centrality(G)
print('Betweenness of Centrality is {}'.format(centrality))

# compute community structure
lpc = nx.community.label_propagation_communities(H)
community_index = {n: i for i, com in enumerate(lpc) for n in com}

#### draw graph ####
fig, ax = plt.subplots(figsize=(20, 15))
pos = nx.spring_layout(H, k=0.15, seed=4572321)
node_color = [community_index[n] for n in H]
node_size = [v * 20000 for v in centrality.values()]
nx.draw_networkx(
    H,
    pos=pos,
    with_labels=False,
    node_color=node_color,
    node_size=node_size,
    edge_color="gainsboro",
    alpha=0.4,
)

# Title/legend
font = {"color": "k", "fontweight": "bold", "fontsize": 20}
ax.set_title("Network Diagram with Weighted Centrality", font)
# Change font color for legend
font["color"] = "r"

ax.text(
    0.80,
    0.06,
    "node size = betweeness centrality",
    horizontalalignment="center",
    transform=ax.transAxes,
    fontdict=font,
)

# Resize figure for label readibility
ax.margins(0.1, 0.05)
# fig.tight_layout()
# plt.axis("off")
plt.show()
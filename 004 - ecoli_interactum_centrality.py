import networkx as nx
import pandas as pd

try:
    df_all = pd.read_csv('data/511145.protein.links.full.v11.5.txt', sep=' ')
    print("Total Rows:  {}".format(df_all.shape[0]))
    print(df_all.head(10))
    print(df_all.columns)

    # filter  > 700 combined_score
    df_filtered = df_all[df_all['combined_score'] > 700]
    print("Total rows after filter:  {}".format(df_filtered.shape[0]))

    df_data = df_filtered[['protein1', 'protein2', 'combined_score']]
    print(df_data.columns)

    # df_alias_all = pd.read_csv('data/511145.protein.aliases.v11.5.txt', sep='\t', usecols=[0])
    # print("Row count of alis raw data:  {}".format(df_alias_all.shape[0]))

    G = nx.from_pandas_edgelist(df_data, "protein1", "protein2", ["combined_score"])

    print('Nodes {}  Edges {}'.format(G.number_of_nodes(), G.number_of_edges()))

    degCent = nx.degree_centrality(G)
    print('Degree Centrality is {}'.format(degCent))

    closeCent = nx.closeness_centrality(G)
    print('Closeness Centrality is {}'.format(closeCent))

    betCent = nx.betweenness_centrality(G)
    print('Betweenness Centrality is {}'.format(betCent))

    assortativity = nx.degree_assortativity_coefficient(G)
    print('assortavity is {}'.format(assortativity))

except Exception as e:
    print(e)

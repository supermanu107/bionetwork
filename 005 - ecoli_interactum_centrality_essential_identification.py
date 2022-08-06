import networkx as nx
import pandas as pd
import numpy as np

try:
    # Read the E-Coli data from STRING Database
    df_STRING_all = pd.read_csv('data/511145.protein.links.full.v11.5.txt', sep=' ')
    print("Row count of STRING raw data:  {}".format(df_STRING_all.shape[0]))

    # Print first 10 rows
    print("Display first 10 rows of STRING raw data")
    print(df_STRING_all.head(10))

    # Print all column names
    print("Display all STRING column names")
    print(df_STRING_all.columns)

    # filter dataset to include only data with > 700 combined_score
    df_STRING_filtered = df_STRING_all[df_STRING_all['combined_score'] > 700]
    # row count of filtered rows
    print("Row count of STRING filtered data combined_score > 700:  {}".format(df_STRING_filtered.shape[0]))

    df_STRING_data = df_STRING_filtered[['protein1', 'protein2', 'combined_score']]
    print("Display df_STRING working data column names")
    print(df_STRING_data.columns)

    # Read the E-Coli data from DEG Database
    df_DEG_all = pd.read_csv('data/511145.protein.aliases.v11.5.txt', sep='\t', usecols=[0])
    print("Row count of DEG raw data:  {}".format(df_DEG_all.shape[0]))

    # Remove the duplicate values from DEG
    df_DEG_data = df_DEG_all.drop_duplicates()
    print("Row count of non duplicate DEG raw data:  {}".format(df_DEG_data.shape[0]))

    G = nx.from_pandas_edgelist(df_STRING_data, "protein1", "protein2", ["combined_score"])

    list_essential_nodes = df_DEG_data['#string_protein_id'].values.tolist()
    # print(list_essential_nodes)

    # ID node as essential or not essential and add node attribute
    print('number of nodes {} number of edges {}'.format(G.number_of_nodes(), G.number_of_edges()))
    for node in G.nodes():
        if node not in list_essential_nodes:
            print("{} is not essential".format(node))
        # else:
        #     print("{} is not essential".format(node))

    # Computing centrality
    # print('Computing degree of centrality')
    # degCent = nx.degree_centrality(G)
    # print('Degree Centrality is {}'.format(degCent))
    #
    # # Computing closeness
    # print('Computing closeness centrality')
    # closeCent = nx.closeness_centrality(G)
    # print('Closeness Centrality is {}'.format(closeCent))
    #
    # # Computing betweeness
    # print('Computing betweenness centrality')
    # betCent = nx.betweenness_centrality(G)
    # print('Betweenness Centrality is {}'.format(betCent))

    print('Program exited successfully')

except Exception as e:
    print('Exception occurred')
    print(e)

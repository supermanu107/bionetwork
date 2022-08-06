import networkx as nx
import pandas as pd

try:
    # Read the txt data
    df_all = pd.read_csv('data/511145.protein.links.full.v11.5.txt', sep=' ')
    print("Row count of raw data:  {}".format(df_all.shape[0]))

    # Print first 10 rows
    print("Display first 10 rows of raw data")
    print(df_all.head(10))

    # Print all column names
    print("Display all column names")
    print(df_all.columns)

    # filter dataset to include only data with > 700 combined_score
    df_filtered = df_all[df_all['combined_score'] > 700]
    # row count of filtered rows
    print("Row count of filtered data combined_score > 700:  {}".format(df_filtered.shape[0]))

    df_working_data = df_filtered[['protein1', 'protein2', 'combined_score']]
    print("Display working data column names")
    print(df_working_data.columns)

    G = nx.from_pandas_edgelist(df_working_data, "protein1", "protein2", ["combined_score"])

    # Computing centrality
    print('Computing degree of centrality')
    degCent = nx.degree_centrality(G)
    print('Degree Centrality is {}'.format(degCent))

    # Computing closeness
    print('Computing closeness centrality')
    closeCent = nx.closeness_centrality(G)
    print('Closeness Centrality is {}'.format(closeCent))

    # Computing betweeness
    print('Computing betweenness centrality')
    betCent = nx.betweenness_centrality(G)
    print('Betweenness Centrality is {}'.format(betCent))

    print('Program exited successfully')

except Exception as e:
    print('Exception occurred')
    print(e)

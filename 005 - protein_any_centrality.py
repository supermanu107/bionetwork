import urllib
import networkx as nx
import pandas as pd
import gzip
import requests
import shutil
import re
import os
from urllib import request


# Download protein links from STRING website
def download_protein_url(url, protein_link_gz_file_name):
    print("\nDownloading proteins from STRING database", url)
    request.urlretrieve(url=url, filename=protein_link_gz_file_name)
    protein_link_file_name = re.split(pattern=r'\.gz', string=protein_link_gz_file_name)[0]

    # Unzip ".gz" file and save in data folder
    with gzip.open(protein_link_gz_file_name, 'rb') as f_in:
        with open("data/" + protein_link_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # delete the *.gz file
    if os.path.exists(protein_link_gz_file_name):
        os.remove(protein_link_gz_file_name)


# Download essential genes from netgenes website
def download_netgenes_url(url, organism_name):
    print("Downloading essential genes from Netgenes database: ", url)
    with open("data/" + organism_name + '.csv', "wb") as file:
        response = requests.get(url)
        file.write(response.content)


try:
    taxonomyID = int(input("Please enter NCBI Taxonomy ID: "))

    # Download all protein connection links from STRING database
    protein_link_file_name = str(taxonomyID) + ".protein.links.v11.5.txt.gz";
    protein_url = "https://stringdb-static.org/download/protein.links.v11.5/" + protein_link_file_name
    download_protein_url(protein_url, protein_link_file_name)

    # Select organism name matching taxonomy id
    df_netgenes_species = pd.read_csv("data/netgenes_species_list.csv", sep=',')
    df_netgenes_selected = df_netgenes_species.loc[df_netgenes_species['NCBI Taxonomy ID'] == taxonomyID]
    selected_organism_name = df_netgenes_selected['Organism Name'].values[0]

    df_all = pd.read_csv('data/' + str(taxonomyID) + '.protein.links.v11.5.txt', sep=' ')
    # filter  > 700 combined_score
    df_filtered = df_all[df_all['combined_score'] >= 700]
    df_data = df_filtered[['protein1', 'protein2', 'combined_score']]
    G1 = nx.from_pandas_edgelist(df_data, "protein1", "protein2", ["combined_score"])
    num_of_nodes = G1.number_of_nodes()
    num_of_edges = G1.number_of_edges()

    # read essential gens from NetGenes Data
    netgenes_url = "https://rbc-dsai-iitm.github.io/NetGenes/CSV/" + urllib.parse.quote(selected_organism_name) + ".csv"
    download_netgenes_url(netgenes_url, selected_organism_name)
    df_essential_genes_all = pd.read_csv('data/' + selected_organism_name + '.csv', header=0,
                                         names=['node', 'alias', 'notes', 'score'], sep=',')
    num_of_essential_nodes = df_essential_genes_all.shape[0]


    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 200)

    essential_nodes_ratio = "{:.1f} %".format((num_of_essential_nodes/num_of_nodes)*100)
    table1_dict = {
        "Organism (NCBI taxonomy ID)": [selected_organism_name + " (" + str(taxonomyID) + ")"],
        "Nodes (proteins)" : [num_of_nodes],
        "Essential nodes": [num_of_essential_nodes],
        "(%)": [essential_nodes_ratio],
        "Edges (interactions)": [num_of_edges]
            }

    print("\nTable 1 - Summary of the networks considered in this study")
    print("----------------------------------------------------------")
    df_table1 = pd.DataFrame(table1_dict)
    print(df_table1)

    df_essential_gene_data_list = df_essential_genes_all['node'].values.tolist()

    # match protein names against essential genes and find match
    # Essential to Essential (EE)
    df_ee_data = df_data.loc[df_data['protein1'].isin(df_essential_gene_data_list) & df_data['protein2'].isin(df_essential_gene_data_list)]
    Gee = nx.from_pandas_edgelist(df_ee_data, "protein1", "protein2", ["combined_score"])
    num_ee = Gee.number_of_edges()
    num_ee_ratio = "{:.3f}".format((num_ee/num_of_edges))

    # Non-Essential to Essential (NE)
    df_ne_data = df_data.loc[-df_data['protein1'].isin(df_essential_gene_data_list) & df_data['protein2'].isin(df_essential_gene_data_list)]
    Gne = nx.from_pandas_edgelist(df_ne_data, "protein1", "protein2", ["combined_score"])
    num_ne = Gne.number_of_edges()
    num_ne_ratio = "{:.3f}".format((num_ne/num_of_edges))

    # Non-Essential to Non-Essential (NN)
    df_nn_data = df_data.loc[~df_data['protein1'].isin(df_essential_gene_data_list) & ~df_data['protein2'].isin(df_essential_gene_data_list)]
    Gnn = nx.from_pandas_edgelist(df_nn_data, "protein1", "protein2", ["combined_score"])
    num_nn = Gnn.number_of_edges()
    num_nn_ratio = "{:.3f}".format((num_nn/num_of_edges))

    assortativity_coefficient = nx.degree_assortativity_coefficient(G1)

    print("\nTable 2 - Assortativity coefficients for the networks considered in this study")
    print("--------------------------------------------------------------------------------")
    table2_dict = {
        "Organism": [selected_organism_name],
        "Nodes (proteins)": [num_of_nodes],
        "EE" : [num_ee_ratio],
        "NE": [num_ne_ratio],
        "NN": [num_nn_ratio],
        "Total edges": [num_of_edges],
        "Assortavity (r)": assortativity_coefficient
            }
    df_table2 = pd.DataFrame(table2_dict)
    print(df_table2)

    # Centrality measures for EE data
    G = nx.from_pandas_edgelist(df_ee_data, "protein1", "protein2", ["combined_score"])

    # "https://stackoverflow.com/questions/50243215/networkx-how-to-write-to-csv-multiple-centrality-metrics"
    df_centrality_measures = pd.DataFrame(dict(
        DEGREE_CENTRALITY=nx.degree_centrality(G),
        EIGENVECTOR=nx.eigenvector_centrality(G),
        CLOSENESS_CENTRALITY=nx.closeness_centrality(G),
        BETWEENNESS_CENTRALITY=nx.betweenness_centrality(G)
    ))
    
    #
    output_file_name = str(taxonomyID) + '-centrality.csv'
    print("\nSaving all centrality values in {}".format(output_file_name))
    df_centrality_measures.to_csv(output_file_name, index_label='Essential Node')

    print('\nDone')

except Exception as e:
    print(e)


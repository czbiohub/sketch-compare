import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster
import scipy.cluster.hierarchy as shc
# from scipy.cluster.hierarchy import dendrogram, linkage

from sklearn.metrics import silhouette_score


def clean_data(distances_tsv, tailing, metadata):
    '''
    clean data and remove water samples
    '''
    
    # read in data tables with sample distances computed with ska (split kmer analysis)
    dd_tsv = pd.read_csv(distances_tsv, sep='\t')
    dd_tsv.columns = dd_tsv.columns.str.replace(" ", "_")
    
    # Add prefix columns and size
    dd_tsv["Sample_1"] = dd_tsv.apply(lambda row: row["Sample_1"].split(tailing)[0], axis=1)
    dd_tsv["Sample_2"] = dd_tsv.apply(lambda row: row["Sample_2"].split(tailing)[0], axis=1)

    # re-organize columns
    dd_tsv = dd_tsv[["Sample_1", 
                     "Sample_2", 
                     "Matches", 
                     "Mismatches", 
                     "Jaccard_Index",
                     "Mash-like_distance",
                     "SNPs",
                     "SNP_distance",
                    ]]
    
    # fix samples that need their file names changed
    # replace row values
    dd_tsv.replace({
        'CMS_0015_RNA_A_S13': 'CMS_015_RNA_A_S13',
        'CMS_001_16_S5_L001': 'CMS_006_RNA_A_S5',  # wierd case
        'CMS_001_17_S6_L001': 'CMS_017_RNA_A_S6',
        'CMS_001_32_S7_L001': 'CMS_032_RNA_A_S7',
        'CMS_037_RNA_A_S21': 'CMS_037_RNA_A_S21',
        'CMS_058_RNA_A_S9':  'CMS_058_RNA_A_S9'
    }, inplace=True)

    water_samples = metadata[
        (metadata.visual_species.isnull()) 
        & 
        (metadata.sourmash_species.isnull())
    ].index.to_list()
    
    # remove data with no species label (water samples)
    filtered_out_water = dd_tsv[~(
        (dd_tsv.Sample_1.isin(water_samples)) 
        | 
        (dd_tsv.Sample_2.isin(water_samples))
    )]

    return filtered_out_water


def create_pivot(df, values="SNP_distance", diagonal=0):    
    ''' create pivot tables with sample 1 and sample 2 correlations'''
    
    dd_tsv_pivot = df.pivot(
        columns="Sample_1", 
        index="Sample_2", 
        values=values
    )
    dd_tsv_pivot = dd_tsv_pivot.loc[dd_tsv_pivot.columns, :]
    np.fill_diagonal(dd_tsv_pivot.values, diagonal)
    dd_tsv_pivot = dd_tsv_pivot.fillna(0) + dd_tsv_pivot.T.fillna(0) - np.diag(dd_tsv_pivot.values.diagonal())
    return dd_tsv_pivot


def add_metadata_to_pivot(df, metadata_fields, metadata):
    """join metadata onto pivot table for axis labeling"""
    
    dd_tsv_pivot_w_labels = pd.merge(
        df.reset_index(), 
        metadata[metadata_fields], 
        how='left',
        left_on="Sample_1", 
        right_index=True, 
    ).reset_index().set_index(
        ["Sample_1"] + metadata_fields
    ).drop(["index"], axis=1)
    
    return dd_tsv_pivot_w_labels


def get_linkage(pivot_df):
    '''clustering linkage'''
    
    cdist = scipy.spatial.distance.squareform(pivot_df)
    return scipy.cluster.hierarchy.linkage(cdist, method="ward")


def get_cluster_map(pivot_df, linkage, num_clusters, extra=True):
    """cluster distance matrix and re-assign species as cluster species mode"""
    
    clusters = fcluster(linkage, num_clusters, criterion='maxclust')
    
    cluster_map = pd.DataFrame()
    
    if extra:
        cluster_map['visual_genus'] = pivot_df.index.get_level_values(1)
        cluster_map['visual_species'] = pivot_df.index.get_level_values(2)
        cluster_map['sourmash_genus'] = pivot_df.index.get_level_values(3)
        cluster_map['sourmash_species'] = pivot_df.index.get_level_values(4)
        
    cluster_map["Sample_1"] = pivot_df.columns
    cluster_map['cluster'] = clusters
    cluster_map.set_index("Sample_1", inplace=True)

    cluster_chunks = []
    for clust in cluster_map.cluster.unique():
        subset = cluster_map[cluster_map.cluster == clust]
        species_mode = subset.visual_species.mode()[0]
        genus_mode = subset[subset["visual_species"] == species_mode].visual_genus.values[0]
        subset["ska_species"] = species_mode
        subset["ska_genus"] = genus_mode
        cluster_chunks.append(subset)

    cluster_map_w_ska = pd.concat(cluster_chunks)
#     cluster_map_w_ska = cluster_map_w_ska[[
#         "genus", 
#         "species", 
#         "sourmash_genus",
#         "sourmash_species", 
#         "ska_genus",
#         "ska_species"
#     ]]
    return cluster_map_w_ska

    
def correlation_matrix(pivot_w_labels, linkage, cluster_map, figsize=(60,60)):
    '''correlatian matrix comparing visual, sourmash, and ska species assignment'''
    
    original_species_coloring = dict(
        zip(
            pivot_w_labels.index.get_level_values(2).unique(), 
            sns.color_palette("bright", 12)
        )
    )
    
    sourmash_species_coloring = dict(
        zip(
            pivot_w_labels.index.get_level_values(4).unique(), 
            sns.color_palette("deep", 12)
        )
    )
    
    ska_species_coloring = dict(
        zip(
            pivot_w_labels.index.get_level_values(6).unique(), 
            sns.color_palette("dark", 13)
        )
    )   
    
    pivot_index = pivot_w_labels.index
    
    key_colors_original = pivot_index.get_level_values(2).map(original_species_coloring)
    key_colors_sourmash = pivot_index.get_level_values(4).map(sourmash_species_coloring)
    key_colors_ska = pivot_index.get_level_values(6).map(ska_species_coloring)
    
    key_colors_labels = [
        key_colors_original, 
        key_colors_sourmash, 
        key_colors_ska
    ]

    return sns.clustermap(
            pivot_w_labels.reset_index().drop([
                "visual_genus", 
                "visual_species", 
                "sourmash_genus",
                "sourmash_species",
                "ska_genus",
                "ska_species"
            ], axis=1).set_index("Sample_1"),
            metric="correlation", 
            row_colors=key_colors_labels, 
            col_colors=key_colors_labels,
            figsize=figsize,
            col_linkage=linkage,
            row_linkage=linkage,
            cbar_kws={'label': 'distance'},
            xticklabels=pivot_w_labels.index.get_level_values(4),
            yticklabels=pivot_w_labels.index.get_level_values(6)
    )


def join_on_ska_labels(fp, metadata, tailing, k, values="SNP_distance"):
    '''read in fp and create pivot distance matrix with old and new species labels'''
    
    df_cleaned = clean_data(
        fp, 
        metadata=metadata,
        tailing=tailing
    )

    dd_pivot = create_pivot(
        df_cleaned, 
        values=values, 
        diagonal=0
    )

    metadata_fields = [
        "visual_genus", 
        "visual_species", 
        "sourmash_genus",
        "sourmash_species"
    ]

    dd_pivot_w_labels = add_metadata_to_pivot(
        dd_pivot, 
        metadata_fields, 
        metadata=metadata
    )

    linkage = get_linkage(dd_pivot)

    cluster_map = get_cluster_map(
        dd_pivot_w_labels, 
        linkage, 
        k, 
        extra=True
    )

    ska_metadata_fields = metadata_fields + ["ska_genus", "ska_species"]

    final_pivot_w_labels = add_metadata_to_pivot(
        dd_pivot, 
        ska_metadata_fields, 
        metadata=cluster_map
    )
    
    return final_pivot_w_labels, cluster_map, linkage


def wrap_clustermap_and_mismatches(fp, metadata, tailing, k=10, figsize=(60,60), values="SNP_distance", diagonal=0):
    '''create clustermap with new and old species assigned labels'''
    
    dd_pivot_w_ska_labels, cluster_map, linkage = join_on_ska_labels(
        fp, 
        metadata, 
        tailing,
        k,
        values=values
    )
    
    c = correlation_matrix(
        dd_pivot_w_ska_labels, 
        linkage, 
        cluster_map, 
        figsize=(60,60)
    )
    
    #return dd_pivot_w_ska_labels
    return dd_pivot_w_ska_labels.query(
        'visual_species != sourmash_species | ' + 
        'visual_species != ska_species | ' + 
        'sourmash_species != ska_species'
    ), c


def hierarchical_clustering(fp, metadata, tailing, k, values="SNP_distance"):

    dd_pivot_w_ska_labels, *_ = join_on_ska_labels(
        fp, 
        metadata, 
        tailing, 
        k,
        values="SNP_distance"
    )
    
    cdist = scipy.spatial.distance.squareform(dd_pivot_w_ska_labels)
    Z = shc.linkage(cdist, method="ward")
    fig = plt.figure(figsize=(25, 10))

    ax = plt.gca()
#     xlbls = ax.get_ymajorticklabels()
#     num=-1
#     for lbl in xlbls:
#         num+=1
#         val=species_coloring[lbl]
#         lbl.set_color(my_palette(val))
    
    dn = shc.dendrogram(
        Z, 
        leaf_font_size=8, 
        labels=dd_pivot_w_ska_labels.index.get_level_values(6), 
        color_threshold=0.008
    )
    
    return dn


def get_silhouette_score(fp, metadata, tailing, values="SNP_distance", correlation=False):
    '''use silhouette metric to see how compact and distinct clusters are '''
    
    clean_df = clean_data(
        fp, 
        metadata=metadata,
        tailing=tailing
    )

    pivot_df = create_pivot(
        clean_df, 
        values=values, 
        diagonal=0
    )
    
    cdist = scipy.spatial.distance.squareform(pivot_df)  
    linkage = scipy.cluster.hierarchy.linkage(cdist, method="ward")
    clusters = fcluster(linkage, 10, criterion='maxclust')
    
    cluster_map = pd.DataFrame()
    
    cluster_map["Sample_1"] = pivot_df.columns
    cluster_map['cluster'] = clusters
    cluster_map.set_index("Sample_1", inplace=True)
    
    return silhouette_score(pivot_df, metric="precomputed", labels=cluster_map["cluster"])

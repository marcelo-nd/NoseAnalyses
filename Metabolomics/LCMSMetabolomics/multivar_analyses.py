#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 09:56:33 2024

@author: marcelo
"""
import numpy as np
import pandas as pd
from scipy.spatial import distance
import skbio
import plotly.express as px
from skbio.stats.distance import permdisp
import plotly.io as pio
from PIL import Image
import io
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.cluster import KMeans
from yellowbrick.cluster import KElbowVisualizer
from sklearn.metrics import silhouette_score
from PyComplexHeatmap import *

def pcoa_metabolomics(cleaned_data, metadata):
    # In order to perform a PCoA as described below,
    # it is important that the filenames in our metadata are identical as well as in the same order
    # as the filenames in our feature table.
    DataCheck = np.isin(metadata.index, cleaned_data.index)

    # Check if all elements are True
    if DataCheck.all():
        print("All elements meet the condition")
        print(np.bincount(DataCheck.astype(int)))
    else:
        print("Not all elements meet the condition.")
    
    print('Dimension: ', cleaned_data.shape)
    print(cleaned_data.head(n=2))
    
    # Then we will calculate pairwise distances across all samples in our data 
    # using the Euclidean distance metric. When we use the Euclidean distance,
    # it's referred to as Principal Component Analysis (PCA). For other distance metrics,
    # it corresponds to a Principal Coordinates Analysis (PCoA).
    raw_pcoa = distance.squareform(distance.pdist(cleaned_data, metric="euclidean"))
    distm = skbio.stats.distance.DistanceMatrix(raw_pcoa, ids=cleaned_data.index)
    
    # All pairwise Euclidean distances are now stored within our distance (dissimilarity) matrix (distm).
    # The distance matrix is then used as input for the PCoA.
    #computing multi-dimensional scaling on distance matrix to 10 components
    PCoA = skbio.stats.ordination.pcoa(distm, number_of_dimensions=10) # 10 principal coordinates

    return PCoA

def pca_plot(pcoa_obj, metadata, attribute, attribute_type='categorical'):
    pcoa_df= pcoa_obj.samples
    explained_var= np.round(pcoa_obj.proportion_explained*100, decimals =2)
    
    title = "PCoA plot representing the sample distribution across different groups"
    df = pd.merge(pcoa_df[['PC1', 'PC2']], metadata[attribute].apply(str if attribute_type == 'categorical'
                                                                      else float),
                  left_index=True,
                  right_index=True)

    if attribute_type == 'continuous':
        fig = px.scatter(df, x='PC1', y='PC2', template='plotly_white', width=600, height=400,
                         color=attribute,
                         color_continuous_scale='Viridis')
    else:
        fig = px.scatter(df, x='PC1', y='PC2', template='plotly_white', width=600, height=400,
                         color=attribute)

        fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={"text": title, 'x': 0.5, "font_color": "#3E3D53", "xanchor": "center"},  # Center title
        xaxis_title=f'PCoA1 {explained_var[0]}%',
        yaxis_title=f'PCoA2 {explained_var[1]}%',
        legend_title_text='',
        legend=dict(
            title=dict(font=dict(size=18, family='bold')),
            itemclick='toggleothers',
            itemdoubleclick='toggle'
        )
    )

    fig.update_traces(marker=dict(size=6))

    return fig

def homoscedasticity_check(cleaned_data, metadata, attribute_permanova):
    # Before conducting PERMANOVA, it's essential to check the homogeneity of group dispersions
    # (also known as 'Homoscedasticity'), as PERMANOVA assumes that group dispersions are equal across groups.
    # Violation of this assumption can increase the Type I error rate.
    group = metadata.loc[:, attribute_permanova]
    group.unique()
    
    raw_pcoa = distance.squareform(distance.pdist(cleaned_data, metric="euclidean")) #compute distance
    distm = skbio.stats.distance.DistanceMatrix(raw_pcoa, ids=cleaned_data.index)
    
    np.random.seed(0)
    print(permdisp(distm, group))
    return None

def permanova_metab(cleaned_data, metadata, attribute_permanova, distmetric="euclidean"):
    # Permutational multivariate analysis of variance (PERMANOVA) is a non-parametric method
    # for multivariate ANOVA, where P-values are obtained using permutations.
    
    # Before conducting PERMANOVA, it's essential to check the homogeneity of group dispersions
    # (also known as 'Homoscedasticity'), as PERMANOVA assumes that group dispersions are equal across groups.
    # Violation of this assumption can increase the Type I error rate.
    
    distm = skbio.stats.distance.DistanceMatrix(distance.squareform(distance.pdist(cleaned_data, metric=distmetric)), ids=cleaned_data.index)
    np.random.seed(0)
    print("\n")
    print("Homoscedasticity check!\n")
    print(permdisp(distm, metadata[attribute_permanova]))
    print("\n")
    
    # Perform PERMANOVA test on the chosen attribute
    permanova = skbio.stats.distance.permanova(distm, metadata[attribute_permanova])
    print("Permanova results!\n")
    print(permanova)
    print("\n")
    
custom_palette = [
    "orange", "green", "red", "blue", "black", "slategray", "purple", "cyan",
    "maroon", "yellow", "magenta", "teal", "pink", "lavender", "olive",
    "turquoise", "brown", "navy", "goldenrod", "indigo", "darkorange"
]

def pcoa_w_metrics(data, meta, distmetric, attribute,
             col_attribute, mdtype='categorical', cols=custom_palette,
             title='Principal coordinates plot', plot=True, print_perm=True):
     
    metrices = ['euclidean','cityblock','canberra','braycurtis','jaccard','minkowski']
    
    if not distmetric in metrices:
      print('Wrong metric. Please choose one out of the following list: ', metrices)

    # Create the distance matrix from the original data
    distance_matrix = skbio.stats.distance.DistanceMatrix(distance.squareform(distance.pdist(data, metric=distmetric)),
                                                         ids=data.index)

    np.random.seed(0)
    # perform PERMDISP test
    permdisp = skbio.stats.distance.permdisp(distance_matrix, meta[attribute])
    if print_perm:
        print("\n")
        print("Homoscedasticity check!\n")
        print(permdisp)
        print("\n")

    # perform PERMANOVA test
    permanova = skbio.stats.distance.permanova(distance_matrix, meta[attribute])

    # SumOfSquares_ratio is the ratio of sum of squares between groups (SSB) and the sum of squares within groups (SSW)
    SumOfSquares_ratio = (permanova['test statistic'] * (permanova['number of groups'] - 1)) / (permanova['sample size'] - permanova['number of groups'])
    permanova['R2'] = SumOfSquares_ratio / (SumOfSquares_ratio + 1)
    if print_perm:
      print("Permanova results!\n")
      print(permanova)
      print("\n")

    # perfom PCoA
    pcoa = skbio.stats.ordination.pcoa(distance_matrix)
    df = pcoa.samples[['PC1', 'PC2']]
    df = df.set_index(meta.index)
    df = pd.merge(df[['PC1', 'PC2']], meta[attribute].apply(str), left_index=True, right_index=True)

    title_text = (
    f'{title}<br>'  # First line with existing title
    f'PERMISP p-value: {permdisp["p-value"]}<br>'  # Second line with PERMISP p-value
    f'PERMANOVA p-value: {permanova["p-value"]}  R2: {permanova["R2"]:.4f}<br>'  # Third line with PERMANOVA p-value and R2
    )

    if mdtype == 'categorical':
      pcoa_w_metrics_plot = px.scatter(df, x='PC1', y='PC2', template='plotly_white', width=600, height=400, color=col_attribute, color_discrete_sequence=cols)
    elif mdtype == 'continuous':
      pcoa_w_metrics_plot = px.scatter(df, x='PC1', y='PC2', template='plotly_white', width=600, height=400, color=col_attribute, color_continuous_scale=cols)
    else:
      print('Wrong mdtype parameter. Please choose from categorical or continuous')
    pcoa_w_metrics_plot.update_layout(font={"color":"grey", "size":12, "family":"Sans"},
                      title={"text": title_text, 'x':0.18, "font_color":"#3E3D53"},
                      xaxis_title=f'PC1 {round(pcoa.proportion_explained[0]*100, 1)}%',
                      yaxis_title=f'PC2 {round(pcoa.proportion_explained[1]*100, 1)}%')
    if plot:
      pcoa_w_metrics_plot.show()
    
    return pcoa_w_metrics_plot
    
def pcoa_explore(cleaned_data_choice, category_permanova, category_type, category_colors, pcoa_titles, distance_metric = "euclidean"):
    

    # Initialize a list to store the PCoA plots
    pcoa_plots = []

# Specifing a higher scale factor for better image quality
    scale_factor = 20  # You can adjust this value as needed
    
    # Iterate through the variables and create the plots
    for i in range(5):
        print(i)
        #plot_data = globals()[cleaned_data_choice.loc[i, 'Variable_name']]
        plot_data = cleaned_data_choice.loc[i, 'Variable_name']
        #pcoa_metadata  = globals()[cleaned_data_choice.loc[i, 'metadata_name']]
        pcoa_metadata  = cleaned_data_choice.loc[i, 'metadata_name']
        #print(np.array_equal(pcoa_metadata['filename'],plot_data.index)) #should return True now
        #print(pcoa_metadata.head())
        #print(plot_data.head())
        if category_type == 'categorical':
            pcoa_metadata[category_permanova] = pcoa_metadata[category_permanova].astype(str)
        
        print(pcoa_metadata.head(2))
        # Generate the PCoA plot
        pcoa_plot = pcoa_w_metrics(data=plot_data, 
                             meta=pcoa_metadata, 
                             distmetric=distance_metric,
                             attribute=category_permanova, 
                             col_attribute=category_colors,
                             mdtype=category_type, 
                             title=pcoa_titles[i], plot=False, print_perm=False)
    
        # Convert Plotly figure to a static image (PNG format) with the specified scale
        img_bytes = pio.to_image(pcoa_plot, format="png", scale=scale_factor)
        img = Image.open(io.BytesIO(img_bytes))
        pcoa_plots.append(img)
    
    # Create subplots for two scatter plots side by side
    fig, axs = plt.subplots(2, 3, figsize=(12, 6))  # 2 rows, 3 columns for six plots
    
    # Display the images in the subplots
    axs[0, 0].imshow(pcoa_plots[0])
    axs[0, 1].imshow(pcoa_plots[1])
    axs[0, 2].imshow(pcoa_plots[2])
    axs[1, 0].imshow(pcoa_plots[3])
    axs[1, 1].imshow(pcoa_plots[4])
    
    # Remove axis labels and ticks
    for row in axs:
        for ax in row:
            ax.axis('off')
    
    # Adjust spacing between subplots
    plt.tight_layout()
    
    return plt

def h_cluster(cleaned_data, cluster_num_method = None):
    np.random.seed(1234) # Setting a seed for reproducing the same result
    linkage_data = linkage(cleaned_data, method='complete', metric='euclidean')
    
    if(cluster_num_method == None):
        dendo_obj = dendrogram(linkage_data, labels=cleaned_data.index, no_plot=False)
        plt.grid(False)
        plt.show()
        return None
    elif(cluster_num_method == "elbow"):
        # Calculating elbow
        # Instantiate the clustering model and visualizer
        model = KMeans()
        visualizer = KElbowVisualizer(model, k=(1,10), timings=False)
        clusters_value = int(visualizer.fit(cleaned_data).elbow_value_)
        #print(elbow_value)
    elif(cluster_num_method == "silhouette"):
        silhouette_scores = []
        for i in range(2,11):
            km = KMeans(n_clusters=i, random_state=42)
            km.fit_predict(cleaned_data)
            silhouette_scores.append(silhouette_score(cleaned_data, km.labels_))
        #print(max(range(len(silhouette_scores)), key=silhouette_scores.__getitem__)+2)
        clusters_value = (max(range(len(silhouette_scores)), key=silhouette_scores.__getitem__)+2)
        
    
    clustering = dendrogram(linkage_data, truncate_mode='lastp', p=clusters_value, show_contracted=True, no_plot=True)
    cluster_labels = fcluster(linkage_data, t=clusters_value, criterion='maxclust')
    # Pair these cluster labels with original labels (which are in cleaned_data.index)
    label_cluster_pairs = list(zip(cleaned_data.index, cluster_labels))
    cluster_results = pd.DataFrame(label_cluster_pairs, columns=['Sample Name', 'Assigned Cluster'])
    # Get the y-axis level corresponding to cluters_value clusters
    clust_array = np.array(clustering['dcoord'])
    clust_threshold = np.floor(np.min(clust_array[np.nonzero(clust_array)]))
    # Plot the dendrogram with the 4 clusters coloured
    dendo_obj = dendrogram(linkage_data, color_threshold=clust_threshold, labels=cleaned_data.index, no_plot=False)
    plt.grid(False)
    plt.show()
        
def metabo_heatmap(cleaned_data, input_list, ins_lev, meta):
    heatmap_attributes = list(ins_lev.iloc[input_list]['ATTRIBUTES'])
    ann = meta.loc[:,heatmap_attributes]
    ann.columns = ann.columns.str.replace('ATTRIBUTE_','') #removing the ATTRBUTE suffix from column names for easier heatmap visualization
    heatmap_attributes = [attr.replace('ATTRIBUTE_','') for attr in heatmap_attributes]
    
    print(heatmap_attributes)
    
    def generate_colors(attributes):
      colors = {}
      for attribute in attributes:
        unique_levels = list(set(ann.loc[:,attribute]))
        n_colors = len(unique_levels)
        selected_colors = custom_palette[:n_colors]
        colors_levels = {unique_levels[i]: selected_colors[i] for i in range(n_colors)}
        colors[attribute] = colors_levels
      return(colors)
    
    colors = generate_colors(heatmap_attributes)
    print(colors)
    
    col_ha = HeatmapAnnotation(df=ann, colors=colors, verbose=0)
    
    cm = ClusterMapPlotter(data=cleaned_data.T, #heatmap_data.T,
                       top_annotation=col_ha,
                       col_dendrogram=True,
                       col_cluster_method='complete',
                       col_cluster_metric='euclidean', # you can change the distance here
                       row_dendrogram=True,
                       row_cluster_method='complete',
                       row_cluster_metric='euclidean', # you can change the distance here
                       cmap='seismic',
                       vmin=-1,  # Set minimum value for colormap
                       vmax=1,  # Set maximum value for colormap
                       verbose=0)
    
    return None

def random_forest():
    return None
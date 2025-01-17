#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:37:16 2025

@author: marcelo
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os

#from sklearn.metrics import silhouette_score
#from sklearn.cluster import KMeans
#from yellowbrick.cluster import KElbowVisualizer

#### From scripts
import sys
sys.path.append('/mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/Metabolomics/LCMSMetabolomics/')
from data_prep import InsideLevels, MergingAnnotationsFT, tidyTables, ft_md_merging, blank_removal, imputation, normalize_ft, scale_ft
from multivar_analyses import pcoa_metabolomics, pca_plot, permanova_metab, pcoa_w_metrics, custom_palette, pcoa_explore, h_cluster, metabo_heatmap, random_forest
from univar_analyses import norm_test, gen_anova_data, anova_vis, anova_boxplots, tukey_post_hoc_test, volcano, tukey_boxplots, gen_ttest_data
from univar_analyses import ttest_volcano, gen_kruskal_data, kruskal_viz, kruskal_boxplots, dunn_post_hoc_test, dunn_volcano, dunn_boxplots

### Get files

taskID = "1bf212a7af8c40e198afabd9c4fcf538"

ft_url = os.path.join('https://gnps2.org/resultfile?task='+taskID+'&file=nf_output/clustering/featuretable_reformated.csv')
md_url = os.path.join('https://gnps2.org/resultfile?task='+taskID+'&file=nf_output/metadata/merged_metadata.tsv')
an_gnps_url = os.path.join('https://gnps2.org/resultfile?task='+taskID+'&file=nf_output/library/merged_results_with_gnps.tsv')

ft = pd.read_csv(ft_url)
md = pd.read_csv(md_url, sep = "\t")
an_gnps = pd.read_csv(an_gnps_url, sep = "\t")


ft = pd.read_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/FBMN/170125_all_samples_blks_ctrl_removed/featuretable_reformated.csv")
md = pd.read_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/FBMN/170125_all_samples_blks_ctrl_removed/merged_metadata.tsv", sep = "\t")
an_gnps = pd.read_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/FBMN/170125_all_samples_blks_ctrl_removed/merged_results_with_gnps.tsv", sep = "\t")
an_analog = an_gnps.copy() # since GNPS2 offers a single merged file containing both actual annotations and analogs

##### Explore original tables
print('Dimension: ', ft.shape) #gets the dimension (number of rows and columns) of ft
print(ft.head()) # gets the first 5 rows of ft

print('Dimension: ', md.shape)
print(md.head())

print(an_gnps.head(n=2))

print(an_analog.head(n=2))

print('Dimension: ', an_gnps.shape)
print('Dimension: ', an_analog.shape)

print(InsideLevels(md))

##### Step 1. Merging annotations from GNPS

ft_w_an = MergingAnnotationsFT(ft, an_gnps, an_analog)

print('Dimension: ', ft_w_an.shape)
ft_w_an.head(n=2)

##### Step 2. Tidy tables
new_tables = tidyTables(ft, md, ft_w_an)

new_ft_tidy = new_tables[0]

new_md_tidy = new_tables[1]

#new_ft_tidy.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/FBMN/170125_all_samples_blks_ctrl_removed/annotated_QuantTable.csv") # save to file

# Data-cleanup

# Merge the metadata and feature table (transposed) together.
ft_merged_with_md = ft_md_merging(new_ft_tidy, new_md_tidy) # This file can be used for batch correction


# Examine Metadata Attributes
ins_lvls = InsideLevels(new_md_tidy)
print(ins_lvls) #skipping the 0th column (filename) and looking at the summary table

#####

##### Step 3. Blank removal. This step creates a new ft table and a new metadata table.
blk_tables = blank_removal(md_df=new_md_tidy, ft_t_df=new_ft_tidy, cutoff = 0.3, sample_type_index=2, blank_level_index=[1], sample_level_index=[2])

blk_rem = blk_tables[0]

md_Samples = blk_tables[1]

print('Dimension: ',blk_rem.shape)
blk_rem.head(n=3)

# metadata without the blanks info:
print('Dimension: ', md_Samples.shape)
md_Samples.head(n=3)

##### Step 4. Imputation (dealing with missing values)
imp_ft = imputation(blk_rem)

##### Step 5. Either do Normalization or Scaling.
### Normalization
tic_norm_ft = normalize_ft(imp_ft)

### Scaling
scaled_ft = scale_ft(imp_ft)

#scaled_ft.to_csv("/mnt/d/2_OtherProjects/SY_bioreactors/SY_bioreactors/annotated_QuantTable_scaled.csv") # save to file

# Check if the feature table and metadata have the same sample names!
if (md_Samples.index == scaled_ft.index).all():
    print("pass")
else:
    print("WARNING: Sample names in feature and metadata table are NOT the same!")

#Setting 'filename' column as 'index' for md_Samples
#md_Samples.set_index('filename', inplace=True)
#md_Samples.index.name = None

# put the rows in the feature table and metadata in the same order
#imp_ft.sort_index(inplace=True)
#md_Samples.sort_index(inplace=True)


########################################################### Multivariate analyses
######## Principal coordinates analysis (PCoA)
pcoa = pcoa_metabolomics(cleaned_data = scaled_ft, metadata = md_Samples)

pca_plot_fig = pca_plot(pcoa_obj = pcoa, metadata = md_Samples, attribute = 'ATTRIBUTE_Sample')

pca_plot_fig.show(renderer="png")

pca_plot_fig.write_image("/mnt/d/2_OtherProjects/SY_bioreactors/Figures/pca_plot.svg")

# Scree plot?

######## Permanova
permanova_metab(cleaned_data = scaled_ft, metadata=md_Samples, distmetric = "euclidean", attribute_permanova="ATTRIBUTE_Volunteer")

################


####### Permanova + PCoA

pcoa_w_metrics_plot = pcoa_w_metrics(data = scaled_ft, meta = md_Samples, distmetric = "braycurtis",
                                     attribute = "ATTRIBUTE_Sample", col_attribute = "ATTRIBUTE_Sample",
                                     mdtype="categorical", cols=custom_palette,
                                     title="Principal coordinates plot", plot=True, print_perm=True)

pcoa_w_metrics_plot.show(renderer="png")

sys.path.append('/mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/Metabolomics/LCMSMetabolomics/')
from multivar_analyses import pcoa_metabolomics, pca_plot, permanova_metab, pcoa_w_metrics, custom_palette, pcoa_explore, h_cluster, metabo_heatmap, random_forest

pcoa_w_metrics_plot = pcoa_w_metrics(data = scaled_ft, meta = md_Samples, distmetric = "braycurtis",
                                     attribute = "ATTRIBUTE_Sample", col_attribute = "ATTRIBUTE_Sample",
                                     mdtype="categorical", title="Principal coordinates plot",
                                     plot=True, print_perm=True, pWidth= 1200, pHeight=850, dot_size=12)

pcoa_w_metrics_plot.show(renderer="png")

pcoa_w_metrics_plot.write_image("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/pcoa_w_metrics.pdf", format="pdf")

####### Permanova + PCoA (Attribute and Data cleanup strategies exploration)

# Define the choice of distance metric, category for permanova calculation,
# category for coloring the PCoA plot and the type of the chosen category in the cell below.

distance_metric = "euclidean"
category_permanova = "ATTRIBUTE_Volunteer"
category_colors = "ATTRIBUTE_Volunteer"
category_type = 'categorical'

pcoa_titles = ["a) Before Data cleanup", "b) After Blank removal", "c) After Imputation", "d) After TIC Normalization", "e) After Scaling"]

# Define the dataframe of data choices
cleaned_data_choice = pd.DataFrame({
    'INDEX': range(1, 6),
    'Description': ["Raw data", "Blank removed data", "Imputed data", "TIC Normalized", "Scaled"],
    'Variable_name': [new_ft_tidy, blk_rem, imp_ft, tic_norm_ft, scaled_ft],
    'metadata_name': [new_md_tidy, md_Samples, md_Samples, md_Samples, md_Samples]
})

cleaned_data = scaled_ft

pcoa_explore_plot = pcoa_explore(cleaned_data_choice = cleaned_data_choice, category_permanova = category_permanova,
                                 category_type = category_type, category_colors = category_colors,
                                 pcoa_titles = pcoa_titles,
                                 distance_metric = distance_metric)

pcoa_explore_plot.savefig("/mnt/d/2_OtherProjects/SY_bioreactors/Figures/figtest3.svg")
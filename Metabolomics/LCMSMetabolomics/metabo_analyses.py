#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:27:01 2024

@author: marcelo
"""

import pandas as pd
import matplotlib.pyplot as plt
#import numpy as np


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

ft = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/SD_BeachSurvey_GapFilled_quant.csv")

md = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/20221125_Metadata_SD_Beaches_with_injection_order.txt", sep = "\t")

an_gnps = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/GNPS_result_FBMN.tsv.txt", sep = "\t")

an_analog = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/GNPS_analog_result_FBMN.tsv.txt", sep = "\t")

### Soyoungs results

ft = pd.read_csv("/mnt/c/Users/marce/Desktop/Sy_bioreactors/Sy_bioreactors_FBMN/09f095a559eb4d909fc032685b6290c7-featuretable_reformated.csv")

md = pd.read_csv("/mnt/c/Users/marce/Desktop/Sy_bioreactors/Sy_bioreactors_FBMN/09f095a559eb4d909fc032685b6290c7-merged_metadata.tsv", sep = "\t")

an_gnps = pd.read_csv("/mnt/c/Users/marce/Desktop/Sy_bioreactors/Sy_bioreactors_FBMN/09f095a559eb4d909fc032685b6290c7-merged_results_with_gnps.tsv", sep = "\t")

an_analog = an_gnps.copy()

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

new_ft_tidy.to_csv("/mnt/c/Users/marce/Desktop/Sy_bioreactors/annotated_QuantTable.csv") # save to file

# Data-cleanup

# Merge the metadata and feature table (transposed) together.
ft_merged_with_md = ft_md_merging(new_ft_tidy, new_md_tidy) # This file can be used for batch correction


# Examine Metadata Attributes
ins_lvls = InsideLevels(new_md_tidy)
print(ins_lvls) #skipping the 0th column (filename) and looking at the summary table

##### Step 3. Blank removal. This step creates a new ft table and a new metadata table.
blk_tables = blank_removal(md_df=new_md_tidy, ft_t_df=new_ft_tidy, cutoff = 0.3, sample_type_index=1, blank_level_index=[1], sample_level_index=[2])

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

scaled_ft.to_csv("/mnt/c/Users/marce/Desktop/Sy_bioreactors/annotated_QuantTable_scaled.csv") # save to file

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

pca_plot_fig = pca_plot(pcoa_obj = pcoa, metadata = md_Samples, attribute = 'ATTRIBUTE_operation_day')

pca_plot_fig.show(renderer="png")

pca_plot_fig.write_image("/mnt/c/Users/marce/Desktop/figtest.svg")

# Scree plot?

######## Permanova
permanova_metab(cleaned_data = scaled_ft, metadata=md_Samples, distmetric = "euclidean", attribute_permanova="ATTRIBUTE_Month")

################ this gives FALSE


####### Permanova + PCoA

pcoa_w_metrics_plot = pcoa_w_metrics(data = scaled_ft, meta = md_Samples, distmetric = "euclidean",
                                     attribute = "ATTRIBUTE_Month", col_attribute = "ATTRIBUTE_Month",
                                     mdtype="categorical", cols=custom_palette,
                                     title="Principal coordinates plot", plot=True, print_perm=True)

pcoa_w_metrics_plot.show(renderer="png")

pcoa_w_metrics_plot.write_image("/mnt/c/Users/marce/Desktop/figtest2.svg")

####### Permanova + PCoA (Attribute and Data cleanup strategies exploration)

# Define the choice of distance metric, category for permanova calculation,
# category for coloring the PCoA plot and the type of the chosen category in the cell below.

distance_metric = "euclidean"
category_permanova = "ATTRIBUTE_Month"
category_colors = "ATTRIBUTE_Month"
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

pcoa_explore_plot.savefig("/mnt/c/Users/marce/Desktop/figtest3.svg")

####### Hierarchial Clustering Algorithm
hc_plot = h_cluster(cleaned_data=new_ft_tidy, cluster_num_method=None)

#hc_plot = h_cluster(cleaned_data=scaled_ft, cluster_num_method="elbow")

#hc_plot = h_cluster(cleaned_data=scaled_ft, cluster_num_method="silhouette")

plt.savefig("/mnt/c/Users/marce/Desktop/colored_dendrogram.svg", format="svg")
####### Heatmaps
print(InsideLevels(new_md_tidy.iloc[:, 1:]))
heatmap_attributes = [2, 3, 9]

metabo_heatmap(cleaned_data=scaled_ft, meta=md_Samples, input_list= heatmap_attributes, ins_lev=ins_lvls)

###### Supervised learning with Random Forest



########################################################### Univariate analyses

# Normality Test
# In order to decide whether to go for parametric or non-parametric tests, we test for normality. 
norm_test_df = norm_test(imp_ft, new_md_tidy)

# ANOVA
anova_attribute = 'ATTRIBUTE_Sample_Area'

anova_results = gen_anova_data(cleaned_data= imp_ft, metadata=new_md_tidy, groups_col = anova_attribute)

results_df_anova = anova_results[0]

results_df_anova

significant_features_anova = anova_results[1]

significant_features_anova

### Step 63: Visualize ANOVA Results
fig = anova_vis(results_df_anova)

fig.show(renderer="png")
# save fig as pdf
fig.write_image("/mnt/c/Users/marce/Desktop/plot_ANOVA.pdf", scale=3) # you can use different file types here e.g. svg, png

# Boxplots for the first 5 
anova_boxplots(imp_ft, new_md_tidy, anova_attribute=anova_attribute, anova_results=results_df_anova, save_to=None)
#fig.show()

# Tukey's post-hoc test

# Applying the tukey_post_hoc_test function
tukey = tukey_post_hoc_test(cleaned_data=scaled_ft, metadata=new_md_tidy,
                            anova_attribute = anova_attribute,
                            contrasts = [("Mission_Bay", "La_Jolla Reefs")],
                            metabolites= significant_features_anova)

tukey.head()

# Visualize Results with a Volcano Plot
volcano_tukey_fig = volcano(sig_results=tukey, anova=results_df_anova)

volcano_tukey_fig.show(renderer="png")
# save image as pdf
volcano_tukey_fig.write_image("/mnt/c/Users/marce/Desktop/TukeyHSD.pdf", scale=3)


tukey_boxplots(cleaned_data=scaled_ft, metadata=new_md_tidy, anova_attribute = anova_attribute,tukey=tukey)

# T-tests
ttest_attribute = 'ATTRIBUTE_Month'
target_group = 'Jan'

ttest = gen_ttest_data(cleaned_data = scaled_ft, metadata = new_md_tidy, ttest_attribute=ttest_attribute, target_group=target_group)
ttest.head(5)

ttest_volcano_plot = ttest_volcano(ttest=ttest, target_group=target_group)

ttest_volcano_plot.show(renderer="png")

# Kruskall wallis (Non-parametric anova)
# select an attribute to perform Kruskal Wallis
kruskal_attribute = 'ATTRIBUTE_Sample_Area'

kruskal_results = gen_kruskal_data(cleaned_data = scaled_ft, metadata = new_md_tidy, groups_col = kruskal_attribute)

kruskal = kruskal_results[0]
print(kruskal.head())

kruskal_signif = kruskal_results[1]
print(kruskal_signif)

kruskal_fig = kruskal_viz(kruskal = kruskal)

kruskal_fig.show(renderer="png")

kruskal_boxplots(cleaned_data=scaled_ft, metadata=new_md_tidy, kruskal = kruskal, kruskal_attribute = kruskal_attribute, features = 5)

# Dunn's post hoc test
dunn_group1 = "La_Jolla Reefs"
dunn_group2 = "Mission_Bay"

dunn_contrasts = [(dunn_group1, dunn_group2)]

# Applying the dunn_post_hoc_test function
dunn = dunn_post_hoc_test(cleaned_data=scaled_ft, metadata=new_md_tidy,
                          kruskal_attribute = kruskal_attribute,
                          contrasts = dunn_contrasts,
                          metabolites= kruskal_signif)

dunn.head()

dunn_volc_fig = dunn_volcano(dunn, kruskal)

dunn_volc_fig.show(renderer="png")

dunn_boxplots(cleaned_data=scaled_ft, metadata=new_md_tidy, dunn = dunn, kruskal_attribute = kruskal_attribute, features = 3)

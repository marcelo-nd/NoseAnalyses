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

#### From scripts
import sys
sys.path.append('/mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/Metabolomics/LCMSMetabolomics/')
from data_prep import InsideLevels, MergingAnnotationsFT, tidyTables, ft_md_merging, blank_removal, imputation, normalize_ft, scale_ft
from multivar_analyses import pcoa_metabolomics, pca_plot, permanova_metab, pcoa_w_metrics, custom_palette, pcoa_explore, h_cluster, metabo_heatmap, random_forest
from univar_analyses import norm_test, gen_anova_data, anova_vis, anova_boxplots, tukey_post_hoc_test, volcano, tukey_boxplots, gen_ttest_data
from univar_analyses import ttest_volcano, gen_kruskal_data, kruskal_viz, kruskal_boxplots, dunn_post_hoc_test, dunn_volcano, dunn_boxplots

### Get files
# no analog search
taskID = "115e1b26a65b436e9057a52305f5a315"

# with analog search
taskID = "4bf27f73090641c4ad8acea79edc8855"

ft_url = os.path.join('https://gnps2.org/resultfile?task='+taskID+'&file=nf_output/clustering/featuretable_reformated.csv')
md_url = os.path.join('https://gnps2.org/resultfile?task='+taskID+'&file=nf_output/metadata/merged_metadata.tsv')
an_gnps_url = os.path.join('https://gnps2.org/resultfile?task='+taskID+'&file=nf_output/library/merged_results_with_gnps.tsv')

ft = pd.read_csv(ft_url)
md = pd.read_csv(md_url, sep = "\t")
an_gnps = pd.read_csv(an_gnps_url, sep = "\t")

# No analog search
ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/FBMN/featuretable_reformated.csv")
md.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/FBMN/merged_metadata.tsv")
an_gnps.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/FBMN/merged_results_with_gnps.tsv")

# With analog search
ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/FBMN/featuretable_reformated.csv")
md.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/FBMN/merged_metadata.tsv")
an_gnps.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/FBMN/merged_results_with_gnps.tsv")


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

# Merge the metadata and feature table (transposed) together. FIX!!!!!!!!!!!!!!!!!

#ft_merged_with_md = ft_md_merging(new_ft_tidy, new_md_tidy) # This file can be used for batch correction

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

# Check if the feature table and metadata have the same sample names!
if (md_Samples.index == scaled_ft.index).all():
    print("pass")
else:
    print("WARNING: Sample names in feature and metadata table are NOT the same!")

# Write the tables to disk
# no analog search
imp_ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_quantTable_imp.csv")
tic_norm_ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_quantTable_ticNorm.csv")
scaled_ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/annotated_QuantTable_scaled.csv") # save to file

# analog search
imp_ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_quantTable_imp.csv", sep="*")
tic_norm_ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_quantTable_ticNorm.csv", sep="*")
scaled_ft.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/annotated_QuantTable_scaled.csv", sep="*") # save to file


# Read tables



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

pca_plot_fig.write_image("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/Figures/pca_plot.pdf", width=1200, height=850, format="pdf")

# Scree plot?

######## Permanova
permanova_metab(cleaned_data = scaled_ft, metadata=md_Samples, distmetric = "euclidean", attribute_permanova="ATTRIBUTE_Volunteer")

####### Permanova + PCoA
pcoa_w_metrics_plot = pcoa_w_metrics(data = scaled_ft, meta = md_Samples, distmetric = "euclidean",
                                     attribute = "ATTRIBUTE_Sample", attribute2 = "ATTRIBUTE_Time", col_attribute = "ATTRIBUTE_Sample",
                                     mdtype="categorical", title="Principal coordinates plot",
                                     plot=True, print_perm=True, pWidth= 1200, pHeight=850, dot_size=12)

pcoa_w_metrics_plot.show(renderer="png", pheight = 800, pwidth = 1200)

pcoa_w_metrics_plot.write_image("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/Figures/pcoa_w_metrics_euclidean.pdf", format="pdf")

# bray curtis
pcoa_w_metrics_plot = pcoa_w_metrics(data = scaled_ft, meta = md_Samples, distmetric = "braycurtis",
                                     attribute = "ATTRIBUTE_Sample", attribute2 = "ATTRIBUTE_Time", col_attribute = "ATTRIBUTE_Sample",
                                     mdtype="categorical", title="Principal coordinates plot",
                                     plot=True, print_perm=True, pWidth= 1200, pHeight=850, dot_size=12)

pcoa_w_metrics_plot.show(renderer="png", pheight = 800, pwidth = 1200)

pcoa_w_metrics_plot.write_image("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/Figures/pcoa_w_metrics_euclidean_braycurtis.pdf", format="pdf")



###### Supervised learning with Random Forest
rf_results = random_forest(feature_table = imp_ft, metadata_table = md_Samples, attribute = "ATTRIBUTE_Sample")

# Order random forest results by features importance
rf_results.sort_values(by="Importance", ascending=False, inplace = True)

# Write random forest to csv
rf_results.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_results.csv")
# with Analog search
rf_results.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_results.csv")

# read rf results
rf_results = pd.read_csv('/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_results.csv', delimiter=',', index_col=0)
# with analog search
rf_results = pd.read_csv('/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_results.csv', delimiter=',', index_col=0)

# Make a list (Series) of the features (this is already ordered by importance)
rf_features_list = rf_results["Feature"]

# Filter feature table to contain only the 100 most important features according to random forest classification.
rf_feature_table = imp_ft[rf_features_list[1:101]]

# Write random forest filtered feature table
rf_feature_table.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_imp.csv")
rf_feature_table.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_imp.csv", sep="*")


# Scaled filtered feature table.
rf_ft_scaled = scale_ft(rf_feature_table)
#rf_feature_table_scaled = scaled_ft[rf_features_list[1:101]]

# Write rf filtered feature table to csv
rf_ft_scaled.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_scaled.csv")
# With analog search
rf_ft_scaled.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_scaled.csv", sep="*")

# PCoA from rf filtered scaled feature table.
#pcoa = pcoa_metabolomics(cleaned_data = rf_ft_scaled, metadata = md_Samples)
#pca_plot_fig = pca_plot(pcoa_obj = pcoa, metadata = md_Samples, attribute = 'ATTRIBUTE_Sample')

pca_plot_fig = pcoa_w_metrics(data = rf_ft_scaled, meta = md_Samples, distmetric = "euclidean",
                              attribute = 'ATTRIBUTE_Sample', col_attribute = 'ATTRIBUTE_Sample',
                              attribute2='ATTRIBUTE_Time', plot=True, print_perm=True, pWidth= 1200, pHeight=850, dot_size=12)

pca_plot_fig.show(renderer="png")

pca_plot_fig.write_image("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/Figures/pcoa_rf_euclidean.pdf", format = "pdf") # save figure

# bray curtis
pca_plot_fig = pcoa_w_metrics(data = rf_ft_scaled, meta = md_Samples, distmetric = "braycurtis",
                              attribute = 'ATTRIBUTE_Sample', col_attribute = 'ATTRIBUTE_Sample',
                              attribute2='ATTRIBUTE_Time', plot=True, print_perm=True, pWidth= 1200, pHeight=850, dot_size=12)

pca_plot_fig.show(renderer="png")
pca_plot_fig.write_image("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/Figures/pcoa_rf_braycurtis.pdf", format = "pdf") # save figure

print(pcoa.eigvals)

# Heatmap between samples and features.
heatmap_attributes = [0]
metabo_heatmap(cleaned_data=rf_ft_scaled, meta=md_Samples, input_list= heatmap_attributes, ins_lev=ins_lvls)

metabo_heatmap(cleaned_data=scaled_ft, meta=md_Samples, input_list= heatmap_attributes, ins_lev=ins_lvls)




##### Filter unannotated features. Only run to keep annotated features!
rf_results_an = rf_results[~rf_results['Feature'].str.contains('nan', na=False)]

# Write filtered random forest results to csv
rf_results_an.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_results_an.csv")
# with analog search
rf_results.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_results_an.csv")

rf_feature_table.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_imp_an.csv")
# With analog search
rf_feature_table.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_imp_an.csv")
# Write rf filtered feature table to csv
#rf_feature_table_scaled.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/3_200125_No_QCs_noSinStrs/DA/rf_featuretable_scaled_an.csv")
# With analog search
#rf_feature_table_scaled.to_csv("/mnt/c/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Metabolomics/UT_LCMS/SC100/Results/4_no.qcs_no.sin.strs_an.search/DA/rf_featuretable_scaled_an.csv")





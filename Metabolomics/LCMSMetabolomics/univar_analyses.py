#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:19:36 2024

@author: marcelo
"""
import pandas as pd
from scipy.stats import shapiro
import statsmodels.api as sm
import pingouin as pg

def norm_test(cleaned_data, metadata):
    #getting 'cleaned_data_with_md' and selecting all columns that starts with X into the 'norm_test'
    # merging metadata with cleaned data
    cleaned_data_with_md = metadata.merge(cleaned_data, left_index=True, right_index=True, how="inner")
    norm_test = cleaned_data_with_md.loc[:, cleaned_data_with_md.columns.str.startswith('X')]
    #norm_test.iloc[:, 0]
    
    uni_data= norm_test.copy()
    op_shapiro = pd.DataFrame(uni_data.columns) # Get all column to a dataframe called "op_shapiro"
    op_shapiro.rename(columns={ op_shapiro.columns[0]: "Metabolites" }, inplace = True) # Convert row "Metabolites" to header and remove the row
    #op_shapiro
    pd.set_option('display.max_colwidth', None)
    
    # Compute Shapiro-Wilk test for each column
    p_values = []
    for col in uni_data.columns:
        _, p_value = shapiro(uni_data[col])
        p_values.append(p_value)

    p_adjusted = sm.stats.fdrcorrection(p_values)[1]
    op_shapiro['p_adj'] = p_adjusted #adding a column "p_adj" with corrected p-values. The method used is FDR
    op_shapiro['p_adj'] = op_shapiro['p_adj'].astype(float)
    #op_shapiro
    # Classify distributions as normal or non-normal based on adjusted p-values
    op_shapiro['distribution'] = ['Normal' if p_adj >= 0.05 else 'Non-normal' for p_adj in op_shapiro['p_adj']]
    # Compute and print number of features with normal and non-normal distribution
    num_normal = sum(op_shapiro['distribution'] == 'Normal')
    num_non_normal = sum(op_shapiro['distribution'] == 'Non-normal')
    print("No. of features with normal distribution:", num_normal)
    print("No. of features with non-normal distribution:", num_non_normal)
    return op_shapiro

def gen_anova_data(cleaned_data, metadata, groups_col):
    
    cleaned_data_with_md = metadata.merge(cleaned_data, left_index=True, right_index=True, how="inner")
    
    norm_test = cleaned_data_with_md.loc[:, cleaned_data_with_md.columns.str.startswith('X')]
    
    results = []

    for col in norm_test.columns:
        result = pg.anova(data=cleaned_data_with_md, dv=col, between=groups_col, detailed=True).set_index('Source')
        p = result.loc[groups_col, 'p-unc']
        f = result.loc[groups_col, 'F']
        meansq = result.loc[groups_col, 'MS']
        significance = 'Significant' if p < 0.05 else 'Non-significant'

        results.append([col, p, f, meansq, significance])
        anova_columns = ['metabolite', 'p', 'F', 'MS', 'Significance']
        # Add Benjamini-Hochberg corrected p-values for multiple testing correction
        
    results = pd.DataFrame(results, columns=anova_columns)
    
    if 'p_fdr_bh' not in results.columns:
        results.insert(2, 'p_fdr_bh', pg.multicomp(results['p'], method='fdr_bh')[1])

    # Add significance based on the corrected p-values
    if 'significant' not in results.columns:
        results.insert(3, 'bh_significant', results['p_fdr_bh'] < 0.05)

    # Sort by p-value
    results.sort_values('p_fdr_bh', inplace=True)


    return results


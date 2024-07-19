#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:19:36 2024

@author: marcelo
"""
import pandas as pd
import numpy as np
from scipy.stats import shapiro
import statsmodels.api as sm
import pingouin as pg
import plotly.express as px
import plotly.graph_objects as go

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
    
    print("Total no.of features on which ANOVA test was performed:", len(results))
    print("No.of features that showed significant difference:", len(results[results["bh_significant"]==True]))
    print("No.of features that did not show significant difference:", len(results[results["bh_significant"]==False]))
    
    # collecting the significant features in an array
    anova_significant_metabolites = results[results["bh_significant"]==True]["metabolite"].to_numpy()

    return [results, anova_significant_metabolites]

def anova_vis(anova_results):
    # first plot insignificant features
    fig = px.scatter(x=anova_results[anova_results['bh_significant'] == False]['F'].apply(np.log),
                y=anova_results[anova_results['bh_significant'] == False]['p_fdr_bh'].apply(lambda x: -np.log(x)),
                template='plotly_white', width=600, height=600)
    fig.update_traces(marker_color="#696880")

    # plot significant features
    fig.add_scatter(x=anova_results[anova_results['bh_significant']]['F'].apply(np.log),
                y=anova_results[anova_results['bh_significant']]['p_fdr_bh'].apply(lambda x: -np.log(x)),
                mode='markers+text',
                text=anova_results['metabolite'].iloc[:4],
                textposition='top left', textfont=dict(color='#ef553b', size=7), name='bh_significant')

    fig.update_layout(font={"color":"grey", "size":12, "family":"Sans"},
                  title={"text":"ANOVA - FEATURE SIGNIFICANCE", 'x':0.5, "font_color":"#3E3D53"},
                  xaxis_title="log(F)", yaxis_title="-log(p)", showlegend=False)

    return fig

def anova_boxplots(cleaned_data, metadata, anova_attribute, anova_results, save_to=None):
    # boxplots with top 5 metabolites from ANOVA
    cleaned_data_with_md = metadata.merge(cleaned_data, left_index=True, right_index=True, how="inner")
    for metabolite in anova_results.sort_values('p_fdr_bh').iloc[:5, 0]:
        fig = px.box(cleaned_data_with_md, x=anova_attribute, y=metabolite, color=anova_attribute)
        fig.update_layout(showlegend=False, title=metabolite, xaxis_title="", yaxis_title="intensity", template="plotly_white", width=500)
        fig.show(renderer="png")

def tukey_post_hoc_test(cleaned_data, metadata, anova_attribute, contrasts, metabolites):
    cleaned_data_with_md = metadata.merge(cleaned_data, left_index=True, right_index=True, how="inner")
    print(cleaned_data_with_md.head)
    # Ensure metabolites is a list
    metabolites = [metabolites] if isinstance(metabolites, str) else metabolites

    results = []
    for metabolite in metabolites:
        for contrast in contrasts:
            filtered_df = cleaned_data_with_md[cleaned_data_with_md[anova_attribute].isin(contrast)][[metabolite, anova_attribute]]
            pairwise_tukey = pg.pairwise_tukey(filtered_df, dv=metabolite, between=anova_attribute)

            # Getting the order from Tukey result
            group_A = pairwise_tukey['A'].iloc[0]
            group_B = pairwise_tukey['B'].iloc[0]

            diff_value = pairwise_tukey['diff'].iloc[0]
            #print(diff_value)
            p_tukey_value = pairwise_tukey['p-tukey'].iloc[0]

            results.append([f'{group_A}-{group_B}', metabolite, int(metabolite.split('_')[0].replace('X', '')), diff_value, p_tukey_value])

    tukey = pd.DataFrame(results, columns=['contrast', 'stats_metabolite', 'stats_ID', 'stats_diff', 'stats_p'])

    # add Benjamini-Hochberg corrected p-values
    tukey['stats_p_bh'] = pg.multicomp(tukey['stats_p'], method='fdr_bh')[1]

    # add significance
    tukey['stats_significant'] = tukey['stats_p_bh'] < 0.05

    # sort by p-value
    tukey.sort_values('stats_p_bh', inplace=True)
    
    print("Total no.of features on which ANOVA test was performed:", len(tukey))
    print("No.of features that showed significant difference:", len(tukey[tukey["stats_significant"]==True]))
    print("No.of features that did not show significant difference:", len(tukey[tukey["stats_significant"]==False]))

    return tukey

def volcano(sig_results, anova):
    fig = px.scatter(template='plotly_white', width=1000, height=600)

    # plot insignificant values
    fig.add_trace(go.Scatter(x=sig_results[sig_results['stats_significant'] == False]['stats_diff'],
                             y=sig_results[sig_results['stats_significant'] == False]['stats_p'].apply(lambda x: -np.log10(x)),
                             mode='markers', marker_color='#696880', name='insignificant'))
    
    # plot significant values
    fig.add_trace(go.Scatter(x=sig_results[sig_results['stats_significant']]['stats_diff'],
                             y=sig_results[sig_results['stats_significant']]['stats_p'].apply(lambda x: -np.log10(x)),
                             mode='markers+text', text=anova['metabolite'].iloc[:5], textposition='top left',
                             textfont=dict(color='#ef553b', size=8), marker_color='#ef553b', name='significant'))
    
    # Grabbing the groups from the first contrast
    group_A = sig_results['contrast'].str.split('-').str[0].iloc[0]
    group_B = sig_results['contrast'].str.split('-').str[1].iloc[0]
    
    # Add annotations to the figure to indicate upregulated and downregulated groups
    fig.add_annotation(
        x=0.85, # Position on x-axis (between 0 and 1)
        y=1.02, # Position on y-axis (between 0 and 1)
        xref="paper",
        yref="paper",
        text=f"{group_A} upregulated",
        showarrow=False
    )
    
    fig.add_annotation(
        x=0.15, # Position on x-axis (between 0 and 1)
        y=1.02, # Position on y-axis (between 0 and 1)
        xref="paper",
        yref="paper",
        text=f"{group_B} upregulated",
        showarrow=False
    )
    
    fig.update_layout(font={"color":"grey", "size":12, "family":"Sans"},
                      title={"text":"TUKEY PLOT", 'x':0.5, "font_color":"#3E3D53"},
                      xaxis_title="stats_diff", yaxis_title="-log(p)")
    
    return fig

def tukey_boxplots(cleaned_data, metadata, anova_attribute, tukey):
    cleaned_data_with_md = metadata.merge(cleaned_data, left_index=True, right_index=True, how="inner")
    rightside_top_metabolites = tukey[tukey['stats_significant']].sort_values('stats_diff', ascending=False)['stats_metabolite'].iloc[:3]
    for stats_metabolite in rightside_top_metabolites:
        fig = px.box(cleaned_data_with_md, x=anova_attribute, y=stats_metabolite, color=anova_attribute)
        fig.update_layout(showlegend=False, title=stats_metabolite, xaxis_title="", yaxis_title="intensity", template="plotly_white", width=500)
        fig.show(renderer="png")
    leftside_top_metabolites = tukey[tukey['stats_significant']].sort_values('stats_diff')['stats_metabolite'].iloc[:3]
    for stats_metabolite in leftside_top_metabolites:
        fig = px.box(cleaned_data_with_md, x=anova_attribute, y=stats_metabolite, color=anova_attribute)
        fig.update_layout(showlegend=False, title=stats_metabolite, xaxis_title="", yaxis_title="intensity", template="plotly_white", width=500)
        fig.show(renderer="png")
        
def gen_ttest_data(cleaned_data, metadata, ttest_attribute, target_group):
    cleaned_data_with_md = metadata.merge(cleaned_data, left_index=True, right_index=True, how="inner")
    uni_data = cleaned_data_with_md.loc[:, cleaned_data_with_md.columns.str.startswith('X')]
    columns = uni_data.columns
    #norm_test.iloc[:, 0]

    ttest = []
    for col in columns:
        group1 = cleaned_data_with_md[col][cleaned_data_with_md[ttest_attribute]==target_group]
        group2 = cleaned_data_with_md[col][cleaned_data_with_md[ttest_attribute]!=target_group]
        result = pg.ttest(group1, group2)
        result['Metabolite'] = col

        ttest.append(result)

    ttest = pd.concat(ttest).set_index('Metabolite')

    ttest.insert(8, 'p-bh', pg.multicomp(ttest['p-val'], method='fdr_bh')[1])
    # add significance
    ttest.insert(9, 'Significance', ttest['p-bh'] < 0.05)
    ttest.sort_values('p-bh', inplace=True)

    return ttest

def ttest_volcano(ttest, target_group):
    # To avoid taking -log10(0)
    ttest['log_p'] = ttest['p-bh'].apply(lambda x: -np.log10(x + 1e-50))
    
    fig = px.scatter(template='plotly_white', width=1000, height=600)
    
    # Plot insignificant values
    fig.add_trace(go.Scatter(x=ttest[ttest['Significance'] == False]['T'],
                             y=ttest[ttest['Significance'] == False]['log_p'],
                             mode='markers', marker_color='#696880', name='insignificant'))
    
    # Plot significant values
    fig.add_trace(go.Scatter(x=ttest[ttest['Significance']]['T'],
                             y=ttest[ttest['Significance']]['log_p'],
                             mode='markers+text', text=ttest.index[:4], textposition='top left',
                             textfont=dict(color='#ef553b', size=8), marker_color='#ef553b', name='significant'))
    
    fig.update_layout(font={"color":"grey", "size":12, "family":"Sans"},
                      title={"text":"T-TEST VOLCANO PLOT", 'x':0.5, "font_color":"#3E3D53"},
                      xaxis_title="T-statistic", yaxis_title="-log(p)")
    
    # Adding annotations
    fig.add_annotation(
        x=0.85,  # Position on x-axis (between 0 and 1)
        y=1.02,  # Position on y-axis (between 0 and 1)
        xref="paper",
        yref="paper",
        text=f"{target_group} upregulated",
        showarrow=False,
        font=dict(size=12, color="black")
    )
    
    fig.add_annotation(
        x=0.15,  # Position on x-axis (between 0 and 1)
        y=1.02,  # Position on y-axis (between 0 and 1)
        xref="paper",
        yref="paper",
        text=f"{target_group} downregulated",
        showarrow=False,
        font=dict(size=12, color="black")
    )
    
    # save image as pdf
    #fig.write_image(os.path.join(data_dir, "T-Test_Volcano_Plot.pdf"), scale=3)
    return fig
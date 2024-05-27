#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:58:36 2024

@author: marcelo
"""
import pandas as pd
import copy
import numpy as np

def InsideLevels(df):
    # get all the columns (equals all attributes) -> will be number of rows
    levels = []
    types = []
    count = []
    for col in df.columns:
        types.append(type(df[col][0]))
        levels.append(sorted(set(df[col].dropna())))
        tmp = df[col].value_counts()
        count.append([tmp[levels[-1][i]] for i in range(len(levels[-1]))])
    return pd.DataFrame({"ATTRIBUTES": df.columns, "LEVELS": levels, "COUNT":count, "TYPES": types}, index=range(1, len(levels)+1))


### Merging annotations with feature table

# Consolidate multiple annotations for a single '#Scan#' into one combined name
def combine_annotation_names(row):
    if row['Compound_Name'] == row['Analog_Compound_Name']:
        return row['Compound_Name']
    return ';'.join([str(row['Compound_Name']), str(row['Analog_Compound_Name'])])

def MergingAnnotationsFT(ft_df, an_df, an_ana_df):
    # deep copies to prevent modification of orignal DFs
    ft_df2 = copy.deepcopy(ft_df)
    an_df2 = copy.deepcopy(an_df)
    an_ana_df2 = copy.deepcopy(an_ana_df)
    
    # Checking if both "#Scan#" columns are of similar class to "row ID" column in FT dataframe
    if(ft_df["row ID"].dtype == an_df2["#Scan#"].dtype and ft_df2["row ID"].dtype == an_ana_df2["#Scan#"].dtype):
        
        # Rename the columns of 'an_analog' with a prefix 'Analog' (excluding the '#Scan#' column)
        an_ana_df2.columns = ['Analog_' + col if col != '#Scan#' else col for col in an_ana_df2.columns]
        
        # Merge 'an_analog' with 'an_gnps' using a full join on the '#Scan#' column
        an_final = pd.merge(an_df2, an_ana_df2, on='#Scan#', how='outer')
        
        # Group by Scan column
        an_final_single = an_final.groupby("#Scan#").apply(lambda group: pd.Series({"Combined_Name": combine_annotation_names(group.iloc[0])})).reset_index()
        
        # To get the DataFrame with that exact column name (without automatic renaming)
        an_final_single.columns = an_final_single.columns.str.replace('.', '_')
        
        # This annotation information is then merged with the feature table.
        ft_an = ft_df2.merge(an_final_single, left_on= "row ID",  how='left', right_on= "#Scan#", sort=True)
        # To do Sirius annotations merging
    return(ft_an)

def tidyTables(ft_p, md_p, ft_an_p):
    new_ft = copy.deepcopy(ft_p)
    new_md = copy.deepcopy(md_p)
    new_ft_an = copy.deepcopy(ft_an_p)
    
    # Cleaning the files
    new_ft.columns = new_ft.columns.str.replace(' Peak area', '') # Removing " Peak area" extensions from the column names of new_ft
    new_ft = new_ft.sort_values(by='row ID') # Arranging the rows of new_ft by ascending order of "row ID"
    
    new_ft = new_ft.loc[:, new_ft.notna().sum() > 0] # Removing columns in new_ft where all values are NaN
    new_md = new_md.loc[:, new_md.notna().sum() > 0] # Removing columns in new_md where all values are NaN
    
    new_md = new_md[new_md.apply(lambda row: all(item != "" for item in row), axis=1)] # Remove rows where all the elements are empty strings
    new_md = new_md.applymap(lambda x: x.strip() if isinstance(x, str) else x) # Remove leading and trailing spaces from each column of new_md
    
    # Update the row names of feature table
    # Compare "row ID" columns from new_ft_an with 
    comparison_result=(new_ft_an['row ID'].values == new_ft['row ID'].values).all()
    print("Comparison result: ", comparison_result)
    if(comparison_result):
        print("Should print True if you have an annotation file")
        # Changing the index (row names) of new_ft into the combined name as "XID_mz_RT":
        new_name = 'X' + new_ft['row ID'].astype(str) + '_' + new_ft['row m/z'].round(3).astype(str) + '_' + new_ft['row retention time'].round(3).astype(str)
        new_name_values = new_name.values
        combined_name_ft = new_ft_an['Combined_Name'].astype(str).values
        underscore_added = ['_' + x for x in combined_name_ft] #add a underscore prefix
        new_name_values = np.core.defchararray.add(new_name_values.astype(str), underscore_added)
        # Set the new index and remove trailing underscore if present
        new_ft.index = new_name_values
        # After adding the annotation information, the feature table "new_ft" now includes the annotation in the row names of the features.
        
        #### Selecting Relevant Columns
        # Selecting only the columns with names containing 'mzXML' or 'mzML'
        new_ft = new_ft.loc[:, new_ft.columns.str.contains('\\.mzXML$|\\.mzML$')]
        
        new_ft.shape

    # If either .mzXML or .mzML files are present, print a message
        if new_ft.columns.str.contains('\\.mzXML$').any() and new_ft.columns.str.contains('\\.mzML$').any():
            #print("Both .mzXML and .mzML file types are present in the data")
            # Checking the files again to see if the above changes have been made:
            #print('Dimension: ', new_ft.shape)
            #print(new_ft.head(n=2))
            #print('Dimension: ', new_md.shape)
            #new_md.head(n=2)
            ### Step 13: Verifying File Consistency
            new_ft = new_ft.reindex(columns=sorted(new_ft.columns)) # Ordering the columns of 'new_ft' by their names
            new_md = new_md.sort_values(by='filename').reset_index(drop=True) #ordering the md by the 1st column filename
            # how many files in the metadata are also present in the feature table
            if (not new_md['filename'].isin(new_ft.columns).value_counts()):  #if this returns False it means some files are missing
                print
    # check if new_ft column names and md row names are the same
    if sorted(new_ft.columns) == sorted(new_md['filename']):
        print(f"All {len(new_ft.columns)} files are present in both new_md & new_ft.")
    else:
        print("Not all files are present in both new_md & new_ft.\n")
        # print the md rows / ft column which are not in ft columns / md rows and remove them
        ft_cols_not_in_md = [col for col in new_ft.columns if col not in new_md['filename']]
        print(f"These {len(ft_cols_not_in_md)} columns of feature table are not present in metadata table and will be removed:\n{', '.join(ft_cols_not_in_md)}\n")
        new_ft.drop(columns=ft_cols_not_in_md, inplace=True)
        md_rows_not_in_ft = [row for row in new_md['filename'] if row not in new_ft.columns]
        print(f"These {len(md_rows_not_in_ft)} rows of metadata table are not present in feature table and will be removed:\n{', '.join(md_rows_not_in_ft)}\n")
        new_md.drop(md_rows_not_in_ft, inplace=True)
    
    return (new_ft, new_md)

# Data-cleanup

def ft_md_merging(new_ft_p, new_md_p):
    # As a first step of data-cleanup step, lets merge the metadata and feature table (transposed) together.
    # This can be used later for other purposes such as batch correction.
    ft_trans = pd.DataFrame(new_ft_p).T #transposing the ft
    ft_trans = ft_trans.apply(pd.to_numeric) #converting all values to numeric
    np.array_equal(new_md_p['filename'],ft_trans.index) #should return True now
    
    print(ft_trans.head(n=3))
    
    #merging metadata (new_md) and transposed feature table based on the sample names
    ft_merged_with_md = pd.merge(new_md_p, ft_trans.reset_index(), left_on='filename', right_on='index', how='left')
    #ft_merged_with_md.head(n=2)
    
    # option to export?
    #ft_merged_with_md.to_csv(os.path.join(data_dir,"Ft_md_merged.csv")) # This file can be used for batch correction
    return(ft_merged_with_md)
    


def batch_correction():
    return ""


def blank_removal(md_df, ft_t_df, sample_index_input, blank_num_input, sample_num_input):
    sample_attribute = md_df.columns[sample_index_input]
    
    df = pd.DataFrame({"LEVELS": InsideLevels(md_df.iloc[:, 1:]).iloc[sample_index_input-1]["LEVELS"]})
    df.index = [*range(1, len(df)+1)]
    print(df.head())
    
    blank_num = InsideLevels(md_df.iloc[:,1:]).iloc[sample_index_input-1]['LEVELS'][(blank_num_input-1)]
    sample_num = InsideLevels(md_df.iloc[:,1:]).iloc[sample_index_input-1]['LEVELS'][(sample_num_input-1)]
    
    print(blank_num)
    print(sample_num)
    
    md_Blank = md_df[md_df[sample_attribute] == blank_num]
    print(md_Blank.head())
    print(md_Blank.shape)
    
    Blank = ft_t_df[ft_t_df.index.isin(md_Blank['filename'])]
    print(Blank.head())
    print(Blank.shape)

    md_Samples  = md_df[md_df[sample_attribute] == sample_num]
    print(md_Samples.head())
    print(md_Samples.shape)
    
    Samples = ft_t_df[ft_t_df.index.isin(md_Samples['filename'])]
    print(Samples.head())
    print(Samples.shape)
    
    print(ft_t_df.head())
    return(Blank, Samples)




    
    

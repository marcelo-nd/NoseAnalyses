#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:58:36 2024

@author: marcelo
"""

import pandas as pd
import copy
import numpy as np
from sklearn.preprocessing import StandardScaler

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
    """"Merge the annotations from GNPS: Actual library hits and analogs. This is given in 'an_final' dataframe.
        The compound names of 'an_gnps' and 'an_analog' are combined by a separator ';'. As a result, for each
        #Scan#, we have a single associated compound name. This result is given in 'an_final_single' dataframe."""

    # deep copies to prevent modification of orignal DFs
    ft_df2 = copy.deepcopy(ft_df)
    an_df2 = copy.deepcopy(an_df)
    an_ana_df2 = copy.deepcopy(an_ana_df)
    
    # Checking if both "#Scan#" columns are of similar class to "row ID" column in FT dataframe
    if(ft_df2["row ID"].dtype == an_df2["#Scan#"].dtype and ft_df2["row ID"].dtype == an_ana_df2["#Scan#"].dtype):
        
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
    """ In the next cells, we bring feature table and metadata in the correct format such that the rownames of the
        metadata and column names of the feature table are the same. Filenames and the order of files need to correspond
        in both tables, as we will match metadata attributes to the feature table. In that way, both metadata and feature
        table, can easily be filtered. """
    
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
    #print("Comparison result: ", comparison_result)
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
        # After adding the annotation information, the feature table "new_ft" now
        # includes the annotation in the row names of the features.
        
        #### Step 16. Selecting Relevant Columns
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
            ### Step 17: Verifying File Consistency
            new_ft = new_ft.reindex(columns=sorted(new_ft.columns)) # Ordering the columns of 'new_ft' by their names
            new_md = new_md.sort_values(by='filename').reset_index(drop=True) #ordering the md by the 1st column filename
            # how many files in the metadata are also present in the feature table
            if (not new_md['filename'].isin(new_ft.columns).value_counts()):  #if this returns False it means some files are missing
                #print(not new_md['filename'].isin(new_ft.columns).value_counts())
                print()
    # check if new_ft column names and md row names are the same
    if sorted(new_ft.columns) == sorted(new_md['filename']):
        print(f"\nAll {len(new_ft.columns)} files are present in both new_md & new_ft.")
    else:
        print("\nNot all files are present in both new_md & new_ft.\n")
        # print the md rows / ft column which are not in ft columns / md rows and remove them
        ft_cols_not_in_md = [col for col in new_ft.columns if col not in new_md['filename']]
        print(f"\nThese {len(ft_cols_not_in_md)} columns of feature table are not present in metadata table and will be removed:\n{', '.join(ft_cols_not_in_md)}\n")
        new_ft.drop(columns=ft_cols_not_in_md, inplace=True)
        md_rows_not_in_ft = [row for row in new_md['filename'] if row not in new_ft.columns]
        print(f"\nThese {len(md_rows_not_in_ft)} rows of metadata table are not present in feature table and will be removed:\n{', '.join(md_rows_not_in_ft)}\n")
        new_md.drop(md_rows_not_in_ft, inplace=True)
        print("\nWill stop execution")
        return
    
    
    # ft_trans is the transposed features table
    ft_trans = pd.DataFrame(new_ft).T #transposing the ft
    ft_trans = ft_trans.apply(pd.to_numeric) #converting all values to numeric
    
    ft_trans = ft_trans.sort_index()
    ft_trans.index.name = None
    #print(ft_trans.head(n=3))
    
    
    # Sort new metadata table.
    new_md = new_md.sort_values("filename")
    # Set index of metadata table with the values of filenames.
    new_md.set_index('filename', inplace=True, drop = True)
    new_md.index.name = None
    #print(new_md.head(n=3))
    
    #print(np.array_equal(new_md['filename'],ft_trans.index)) #should return True now
    if(np.array_equal(new_md.index,ft_trans.index)):
        print("All files in MD and FT are the same.")
    else:
        print("Not all files in MD and FT are the same.")
    
    #checking the dimensions of our new ft and md:
    print("\nThe number of rows and columns in our original ft is:", ft_p.shape,"\n")
    #print("The number of rows and columns in our new transposed ft is:", new_ft.shape,"\n")
    print("The number of rows and columns in our new transposed ft is:", ft_trans.shape,"\n")
    print("The number of rows and columns in our new md is:", new_md.shape)
    
    #print(np.array_equal(new_md.index,ft_trans.index)) #should return True now
    if(np.array_equal(new_md.index,ft_trans.index)):
        print("\nIndices of md and ft tables are the same!")
    
    return (ft_trans, new_md)

# Data-cleanup

def ft_md_merging(new_ft_p, new_md_p):
    # As a first step of data-cleanup step, lets merge the metadata and feature table (transposed) together.
    # This can be used later for other purposes such as batch correction.
       
    #merging metadata (new_md) and transposed feature table based on the sample names
    ft_merged_with_md = pd.merge(new_md_p, new_ft_p.reset_index(), left_on='filename', right_on='index', how='left')
    #ft_merged_with_md.head(n=2)
    
    # option to export?
    #ft_merged_with_md.to_csv(os.path.join(data_dir,"Ft_md_merged.csv")) # This file can be used for batch correction
    return(ft_merged_with_md)

def batch_correction():
    return ""

def blank_removal(md_df, ft_t_df, cutoff, sample_type_index, blank_level_index, sample_level_index):
    # This function returns a new ft table and a new metadata table.
    sample_attribute = md_df.columns[sample_type_index-1]
    print(sample_attribute)
    
    df = pd.DataFrame({"LEVELS": InsideLevels(md_df).iloc[sample_type_index-1]["LEVELS"]})
    df.index = [*range(1, len(df)+1)]
    print(df.head())
    
    print(InsideLevels(md_df).iloc[sample_type_index-1]['LEVELS'])
    
    # Input for blanks, allowing multiple indices
    #blank_indices = [int(index.strip()) - 1 for index in blank_level_index]
    blank_indices= blank_level_index
    #blank_nums = [InsideLevels(md_df).iloc[sample_type_index-1]['LEVELS'][index] for index in blank_indices]
    blank_nums = [InsideLevels(md_df).iloc[sample_type_index-1]['LEVELS'][index-1] for index in blank_indices]
    print(blank_indices)
    print(blank_nums)
    
    # Input for samples, allowing multiple indices
    #sample_indices = [int(index.strip()) - 1 for index in sample_level_index]
    sample_indices = sample_level_index
    #sample_nums = [InsideLevels(md_df).iloc[sample_type_index-1]['LEVELS'][index] for index in sample_indices]
    sample_nums = [InsideLevels(md_df).iloc[sample_type_index-1]['LEVELS'][index-1] for index in sample_indices]
    print(sample_indices)
    print(sample_nums)

    md_Blank = md_df[md_df[sample_attribute].isin(blank_nums)]
    
    #Blank = ft_t_df[ft_t_df.index.isin(md_Blank['filename'])]
    Blank = ft_t_df[ft_t_df.index.isin(md_Blank.index)]

    md_Samples  = md_df[md_df[sample_attribute].isin(sample_nums)]
    
    Samples = ft_t_df[ft_t_df.index.isin(md_Samples.index)]
    
    #return(Blank, Samples)
    # Getting mean for every feature in blank and Samples in a DataFrame named 'Avg_ft'
    Avg_ft = pd.DataFrame({'Avg_blank': Blank.mean(axis=0, skipna=False)}) # Set skipna=False to check if there are NA values
    Avg_ft['Avg_samples'] = Samples.mean(axis=0, skipna=False) # Adding another column 'Avg_samples' for feature means of samples
    
    # Getting the ratio of blank vs Sample
    Avg_ft['Ratio_blank_Sample'] = (Avg_ft['Avg_blank'] + 1) / (Avg_ft['Avg_samples'] + 1)

    # Creating a bin with 1s when the ratio > Cutoff, else put 0s
    Avg_ft['Bg_bin'] = (Avg_ft['Ratio_blank_Sample'] > cutoff).astype(int)

    # Calculating the number of background features and features present
    print("Total no.of features:", Avg_ft.shape[0])
    print("No.of Background or noise features:", Avg_ft['Bg_bin'].sum())
    print("No.of features after excluding noise:", (Samples.shape[1] - Avg_ft['Bg_bin'].sum()))

    # Merging Samples with Avg_ft and selecting only the required rows and columns
    blk_rem = pd.concat([Samples.T, Avg_ft], axis=1, join='inner')
    blk_rem = blk_rem[blk_rem['Bg_bin'] == 0]  # Picking only the features
    blk_rem = blk_rem.drop(columns=['Avg_blank', 'Avg_samples', 'Ratio_blank_Sample', 'Bg_bin'])  # Removing the last 4 columns
    blk_rem = blk_rem.T
    
    return(blk_rem, md_Samples)

def imputation(ft_p):
    # get the lowest intensity (that is not zero) as a cutoff LOD value
    cutoff_LOD = round(ft_p.replace(0, np.nan).min(numeric_only=True).min())
    # Step 25: Random Value Generation & Zero Replacement
    
    # Set the seed for random number generation
    np.random.seed(141222)

    imp_ft = ft_p.copy()
    imp_ft = imp_ft.applymap(lambda x: np.random.randint(1, cutoff_LOD) if x == 0 else x)
    
    print('Dimension: ', imp_ft.shape)
    print(imp_ft.head(n=3))
    
    (imp_ft == 0).sum().sum() # checking if there are any zeros in our imputed table
    
    return(imp_ft)
    
def normalize_ft(p_ft):
    # Step 27: Total Ion Current (TIC) or sample-centric normalization
    normalized = p_ft.copy()
    
    # Dividing each element of a particular row (as each row is the sample) with its row sum
    norm_TIC = normalized.apply(lambda x: x/np.sum(x), axis=1)
    print('No NA values in Normalized data:',norm_TIC.isnull().values.any() == False)
    print('Dimension: ', norm_TIC.shape)
    print(norm_TIC.head(n=3))
    
    return(norm_TIC)

def scale_ft(ft_p):
    # center and scale filtered data
    ft_scaled = pd.DataFrame(StandardScaler().fit_transform(ft_p),
                              index=ft_p.index, columns=ft_p.columns)
    print(ft_scaled.head(n=2))
    return(ft_scaled)


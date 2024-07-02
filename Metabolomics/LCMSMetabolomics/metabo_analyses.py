#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:27:01 2024

@author: marcelo
"""

import pandas as pd
import numpy as np

import glob
import os
import re
import itertools
import io
import session_info
import shutil
import warnings
import matplotlib.pyplot as plt
import platform
import math
import scipy.stats as stats

import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.io as pio

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import distance
from scipy.stats import lognorm
from scipy.stats import shapiro
from scipy.stats import lognorm
from scipy.stats import spearmanr

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.preprocessing import OrdinalEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import class_weight
from sklearn.metrics import classification_report
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans

import pingouin as pg
import scikit_posthocs as sp
from PyComplexHeatmap import *
from yellowbrick.cluster import KElbowVisualizer

from PIL import Image

import statsmodels.api as sm

import time
from sklearn.inspection import permutation_importance

#### Not imported anymore?
from skbio.stats.distance import permdisp
import skbio


#### From scripts
import sys, os
sys.path.append('/mnt/c/Users/marce/Documents/GitHub/NoseAnalyses/Metabolomics/LCMSMetabolomics/')
from helper_functions import InsideLevels, combine_annotation_names, MergingAnnotationsFT, tidyTables, ft_md_merging, blank_removal


ft = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/SD_BeachSurvey_GapFilled_quant.csv")

md = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/20221125_Metadata_SD_Beaches_with_injection_order.txt", sep = "\t")

an_gnps = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/GNPS_result_FBMN.tsv.txt", sep = "\t")

an_analog = pd.read_csv("/mnt/d/1_NoseSynComProject/Metabolomics Data/metaboData/GNPS_analog_result_FBMN.tsv.txt", sep = "\t")

print('Dimension: ', ft.shape) #gets the dimension (number of rows and columns) of ft
print(ft.head()) # gets the first 5 rows of ft

print('Dimension: ', md.shape)
print(md.head())

print(an_gnps.head(n=2))

print(an_analog.head(n=2))

print('Dimension: ', an_gnps.shape)
print('Dimension: ', an_analog.shape)

print(InsideLevels(md))

ft_w_an = MergingAnnotationsFT(ft, an_gnps, an_analog)

new_tables = tidyTables(ft, md, ft_w_an)

new_ft_tidy = new_tables[0]

new_md_tidy = new_tables[1]


#checking the dimensions of our new ft and md:
print("The number of rows and columns in our original ft is:", ft.shape,"\n")
print("The number of rows and columns in our new ft is:", new_ft_tidy.shape,"\n")
print("The number of rows and columns in our new md is:", new_md_tidy.shape)

# Data-cleanup

# As a first step of data-cleanup step, lets merge the metadata and feature table (transposed) together.
ft_t = ft_md_merging(new_ft_tidy, new_md_tidy)


# Examine Metadata Attributes
ins_lvls = InsideLevels(new_md_tidy.iloc[:, 1:])
print(ins_lvls) #skipping the 0th column (filename) and looking at the summary table

blank_tables = blank_removal(md_df=new_md_tidy, ft_t_df=ft_t, sample_index_input=1, blank_num_input=1, sample_num_input=2)

blanks_table = blank_tables[0]

samples_table = blank_tables[1]

# Display the chosen blanks
print('Dimension: ',blanks_table.shape)
blanks_table.head(2)

# Display the chosen samples
print('Dimension: ',samples_table.shape)
samples_table.head(2)


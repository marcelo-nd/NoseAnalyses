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
from skbio.stats.distance import permdisp

#### From scripts
from helper_functions import InsideLevels, MergingAnnotationsFT, tidyTables


ft = pd.read_csv("/mnt/c/metaboData/SD_BeachSurvey_GapFilled_quant.csv")

md = pd.read_csv("/mnt/c/metaboData/20221125_Metadata_SD_Beaches_with_injection_order.txt", sep = "\t")

an_gnps = pd.read_csv("/mnt/c/metaboData/GNPS_result_FBMN.tsv.txt", sep = "\t")

an_analog = pd.read_csv("/mnt/c/metaboData/GNPS_analog_result_FBMN.tsv.txt", sep = "\t")

#print('Dimension: ', ft.shape) #gets the dimension (number of rows and columns) of ft
#print(ft.head()) # gets the first 5 rows of ft

#print('Dimension: ', md.shape)
#print(md.head())

#print(an_gnps.head(n=2))

#print(an_analog.head(n=2))

#print('Dimension: ', an_gnps.shape)
#print('Dimension: ', an_analog.shape)

#print(InsideLevels(md))

ft_w_an = MergingAnnotationsFT(ft, an_gnps, an_analog)

print('Dimension: ', ft_w_an.shape)
print(ft_w_an.head())

new_tables = tidyTables(ft, md, ft_w_an)

new_ft = new_tables[0]

new_md = new_tables[1]

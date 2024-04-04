# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 19:15:22 2024

@author: marce
"""

import random

import pandas as pd

# Set a fixed random seed for reproducibility
random.seed(42)

# Define the number of species and strains per species
num_species = 10
num_strains_per_species = 3

# Assign unique IDs to each strain
strains = [{'species_id': i, 'strain_id': j, 'unique_id': (i - 1) * num_strains_per_species + j} for i in range(1, num_species + 1) for j in range(1, num_strains_per_species + 1)]

# Generate 100 combinations
combinations = []
for _ in range(100):
    # Choose a random size for the combination (between 5 and 10)
    combination_size = random.randint(5, 10)
    
    # Shuffle the list of strains for randomness
    random.shuffle(strains)

    # Choose randomly one strain per species
    chosen_strains = []
    species_ids_chosen = set()
    for strain in strains:
        if strain['species_id'] not in species_ids_chosen:
            chosen_strains.append(strain)
            species_ids_chosen.add(strain['species_id'])
            if len(chosen_strains) == combination_size:  # Stop when the desired combination size is reached
                break

    # Return the list of strains using the unique id number
    combination = sorted([strain['unique_id'] for strain in chosen_strains])  # Sort the IDs in numerical order
    combinations.append(combination)

# Sort combinations in numerical order
combinations.sort()

#print(combinations)

# Print the combinations
#for idx, combination in enumerate(combinations, start=1):
#    print(f'Combination {idx}: {combination}')

new_combinations = []

for idx, combination in enumerate(combinations, start=1):
    new_combination = []
    #print(combination)
    for idx2 in range(1, (num_species*num_strains_per_species)+1):
        if idx2 in combination:
            new_combination.append(idx2)
        else:
            new_combination.append(0)
    new_combinations.append(new_combination)

#print(new_combinations)


# Assuming your list of lists is named 'list_of_lists'
# and all sublists have the same size
#new_combinations = [ [1, 2, 3, ...], [4, 5, 6, ...], ... ]  # Your data here

# Creating a DataFrame
df = pd.DataFrame(new_combinations)

# Setting the index from 1 to 100
df.index = range(1, 101)

df = df.transpose()

#new_rownames = ['Strain_{}'.format(i) for i in range(1, 31)]

new_rownames = ["Anaerococcus octavius 133",
                "Anaerococcus octavius 211",
                "Anaerococcus octavius 259",
                "Corynebacterium accolens 99",
                "Corynebacterium accolens 157",
                "Corynebacterium accolens 184",
                "Corynebacterium propinquum 16",
                "Corynebacterium propinquum 70",
                "Corynebacterium propinquum 265",
                "Corynebacterium pseudodiphtheriticum DSM44287",
                "Corynebacterium pseudodiphtheriticum 242",
                "Corynebacterium pseudodiphtheriticum 244",
                "Corynebacterium tuberculostearicum DSM44922",
                "Corynebacterium tuberculostearicum 102",
                "Corynebacterium tuberculostearicum 223",
                "Cutibacterium acnes 33",
                "Cutibacterium acnes 86",
                "Cutibacterium acnes 149",
                "Cutibacterium avidum 32",
                "Cutibacterium avidum 181",
                "Cutibacterium avidum 208",
                "Dolosigranulum pigrum 21",
                "Dolosigranulum pigrum 61",
                "Dolosigranulum pigrum 245",
                "Staphylococcus epidermidis 28",
                "Staphylococcus epidermidis 231",
                "Staphylococcus epidermidis 251",
                "Staphylococcus lugdunensis 81",
                "Staphylococcus lugdunensis 115",
                "Staphylococcus lugdunensis 239"
                ]

df.rename(index=dict(zip(df.index, new_rownames)), inplace=True)

print(df)

# Replace all values above 0 with 1

df2 = df.applymap(lambda x: 1 if x > 0 else x)

print(df)
print(df2)

df2.to_excel('C:/Users/marce/Desktop/output.xlsx', index=True)  # Set index=True if you want to include the DataFrame index in the Excel file

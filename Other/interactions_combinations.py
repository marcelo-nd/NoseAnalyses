# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:15:19 2023

@author: Marcelo
"""

set1 = ["1.1", "1.2", "1.3"]

set2 = ["2.1", "2.2", "2.3"]

set3 = ["3.1", "3.2", "3.3"]

set4 = ["4.1", "4.2", "4.3"]

set5 = ["5.1", "5.2", "5.3"]

set6 = ["6.1", "6.2", "6.3"]

set7 = ["7.1", "7.2", "7.3"]

set8 = ["8.1", "8.2", "8.3"]

set9 = ["9.1", "9.2", "9.3"]

set10 = ["10.1", "10.2", "10.3"]


sets_list = [set1, set2, set3, set4, set5, set6, set7, set8, set9, set10]

pairs = []

for set_origin in range(0, 10):
    destination_sets = list(range(0 , set_origin)) + list(range(set_origin+1, 4))
    #print("dest_sets: ", destination_sets)
    for set_destination in destination_sets:
        for element1 in range(0, 3):
            for element2 in range(0, 3):
                #print(sets_list[set_origin][element1], set_destination[set_destination][element2])
                #print("set_or: ", set_origin)
                #print("set_dest: ", set_destination)
                #print("dest element: ", sets_list[set_destination])
                #print({sets_list[set_origin][element1], sets_list[set_destination][element2]})
                current_pair = {sets_list[set_origin][element1], sets_list[set_destination][element2]}
                #print(current_pair)
                #print(pairs)
                if current_pair not in pairs:
                    pairs.append(current_pair)
                    #print(current_pair)

print(pairs)

print(len(pairs))

print(len(pairs)/9)

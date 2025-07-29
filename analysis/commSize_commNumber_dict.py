#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:33:24 2024

@author: magdalena
"""

# import random
# import matplotlib.pyplot as plt
# import numpy as np
# import cobra
from os import listdir
from os.path import isfile, join

# import xlrd
import misosoup_analysis_module as mym

# from os.path import isdir
# import pickle
# import scipy.cluster.hierarchy as sch
import yaml
from cobra.io import read_sbml_model

cs_sorted_fba_mc = [
    "pep",
    "icit",
    "glu__D",
    "3pg",
]  #'cyst__L',

cs_sorted_fba_msc = [
    "r5p",
    "met__L",
    "akg",
    "oaa",
    "f6p",
    "phe__L",
    "trp__L",
    "rmn",
    "for",
    "ile__L",
    "rib__D",
    "ser__D",
    "acgam",
    "cit",
    "ac",
    "pyr",
    "glc__D",
    "pro__L",
]  #'cys__L'#carbon sources sorted according to increasing number of strains that are viable in isolation (obtained with FBA)

cs_sorted_fba = cs_sorted_fba_mc + cs_sorted_fba_msc

all_species = [
    "E3R11",
    "C1R06",
    "m_3D05",
    "E3R18",
    "F3R11",
    "m_4F10",
    "C1M14",
    "A2M03",
    "m_5C01",
    "D3R19",
    "m_6D02",
    "B3R18",
    "B3M02",
    "C3R17",
    "m_6D03",
    "D2R18",
    "G2R10",
    "B3R10",
    "D2M02",
    "A3R16",
    "m_4A10",
    "D2R05",
    "C3M10",
    "D2R04",
    "A3R04",
    "A3M17",
    "m_3C02",
    "D2M19",
    "A1R12",
    "C2M11",
    "C2R02",
    "E3M17",
    "I3M07",
    "B2M17",
    "G3M19",
    "m_4B03",
    "E3M18",
    "m_1A01",
    "F3M07",
    "m_3B05",
    "B2M06",
    "E3M09",
    "C3R15",
    "I2M19",
    "C3R12",
    "C3M06",
    "m_4E07",
    "G2R14",
    "G2M07",
    "F3M07_",
    "A3R12",
    "E2M01",
    "F2R14",
    "I2R16",
    "A3R06",
    "E3R09",
    "F3R08",
    "E3R11_",
    "C2R09",
    "E3R01",
    "B2M13",
    "m_6C06",
]

path_results_script = "out/analysis/"


###
# load results needed
tolerance = 7

with open(
    path_results_script
    + "d_community_members_tolerance_1e-"
    + str(tolerance)
    + ".yaml",
    "r",
) as file:
    d_community_members = yaml.load(file, Loader=yaml.SafeLoader)


###
# for each cs and species, store number of communities of each size
d_comm_size_and_number = {}

for cs in cs_sorted_fba_mc:
    if len(d_community_members[cs]) == 0:
        d_comm_size_and_number[cs] = []
    else:
        d_comm_size_and_number[cs] = {}
        for comm in d_community_members[cs]:
            comm_size = len(comm)
            if comm_size not in d_comm_size_and_number[cs]:
                d_comm_size_and_number[cs][comm_size] = 1
            else:
                d_comm_size_and_number[cs][comm_size] += 1

for cs in cs_sorted_fba_msc:
    d_comm_size_and_number[cs] = {}
    for s in all_species:
        if len(d_community_members[cs][s]) == 0:
            d_comm_size_and_number[cs][s] = []
        else:
            d_comm_size_and_number[cs][s] = {}
            for comm in d_community_members[cs][s]:
                comm_size = len(comm)
                if comm_size not in d_comm_size_and_number[cs][s]:
                    d_comm_size_and_number[cs][s][comm_size] = 1
                else:
                    d_comm_size_and_number[cs][s][comm_size] += 1

###
# print (non trivial) results
for cs in d_comm_size_and_number:
    if cs in cs_sorted_fba_msc:
        for s in d_comm_size_and_number[cs]:
            if (
                len(d_comm_size_and_number[cs][s]) != 0
                and 1 not in d_comm_size_and_number[cs][s]
            ):
                print(cs, s, d_comm_size_and_number[cs][s])

    else:
        if len(d_comm_size_and_number[cs]) != 0:
            print(cs, d_comm_size_and_number[cs])

# total number of minimal communities
mc_2sp = 0
mc_3sp = 0
for cs in cs_sorted_fba_mc:
    if len(d_comm_size_and_number[cs]) != 0:
        for k in d_comm_size_and_number[cs]:
            if k == 2:
                mc_2sp += d_comm_size_and_number[cs][k]
            elif k == 3:
                mc_3sp += d_comm_size_and_number[cs][k]
            else:
                print("Comm size not taken into account: ", k)

# total number of minimal supplying communities
msc_1sp = 0
msc_2sp = 0
msc_3sp = 0
for cs in cs_sorted_fba_msc:
    if len(d_comm_size_and_number[cs]) != 0:
        for s in d_comm_size_and_number[cs]:
            if len(d_comm_size_and_number[cs][s]) != 0:
                for k in d_comm_size_and_number[cs][s]:
                    if k == 1:
                        msc_1sp += d_comm_size_and_number[cs][s][k]
                    elif k == 2:
                        msc_2sp += d_comm_size_and_number[cs][s][k]
                    elif k == 3:
                        msc_3sp += d_comm_size_and_number[cs][s][k]
                    else:
                        print("Comm size not taken into account: ", k)


### WITH TOLERANCE=7: 2210 minimal (supplying) communities

# 15 minimal communities:
# mc_2sp=15
# mc_3sp=0

# 2195 minimal supplying communities with 2 or 3 species:
# msc_1sp=214
# msc_2sp=1912
# msc_3sp=283
###

###
# save dictionary
with open(
    path_results_script
    + "d_cs_commSize_commNumber_tolerance_1e-"
    + str(tolerance)
    + ".yaml",
    "w",
) as file:
    documents = yaml.dump(d_comm_size_and_number, file)


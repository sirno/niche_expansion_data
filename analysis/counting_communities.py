#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 11:10:26 2025

@author: magdalena
"""

import itertools
import math
import random
from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt

# import xlrd
import misosoup_analysis_module as mym
import numpy as np

# from os.path import isdir
# import pickle
# import scipy.cluster.hierarchy as sch
import yaml

# import cobra
from cobra.io import read_sbml_model
from matplotlib import gridspec
from scipy import stats

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
path_results_misosoup = "out/misosoup/"


###
# load results needed

# for tolerance in [7,8,9]:
tolerance = 7
print("Tolerance: ", tolerance)

with open(
    path_results_script
    + "d_community_members_tolerance_1e-"
    + str(tolerance)
    + ".yaml",
    "r",
) as file:
    d_community_members = yaml.load(file, Loader=yaml.SafeLoader)


Comms = []
for c in d_community_members:
    if c in cs_sorted_fba_mc:
        for comm in d_community_members[c]:
            if ("y_F3M07_" not in comm) and ("y_E3R11_" not in comm):
                Comms.append(comm)
    elif c in cs_sorted_fba_msc:
        for sp in d_community_members[c]:
            for comm in d_community_members[c][sp]:
                if ("y_F3M07_" not in comm) and ("y_E3R11_" not in comm):
                    Comms.append(comm)


Comms.sort()
print("Total number of misosoup solutions: " + str(len(Comms)))  # 2318

###
# number of strains in communities
str1_comm = 0
str2_comm = 0
str3_comm = 0

for comm in Comms:
    if len(comm) == 1:
        str1_comm += 1
    elif len(comm) == 2:
        str2_comm += 1
    elif len(comm) == 3:
        str3_comm += 1

print(str(str1_comm) + " 1-strain community")  # 214
print(str(str2_comm) + " 2-strain community")  # 1853
print(str(str3_comm) + " 3-strain community")  # 251


###
# unique comms
Comms_unique = list(
    Comms for Comms, _ in itertools.groupby(Comms)
)  # 735 unique communities
print(str(len(Comms_unique)) + " unique communities: ")

str1_comm = 0
str2_comm = 0
str3_comm = 0

for comm in Comms_unique:
    if len(comm) == 1:
        str1_comm += 1
    elif len(comm) == 2:
        str2_comm += 1
    elif len(comm) == 3:
        str3_comm += 1

print(str(str1_comm) + " 1-strain community")  # 30
print(str(str2_comm) + " 2-strain community")  # 482
print(str(str3_comm) + " 3-strain community")  # 223


###
# number of complete minimal supplying communities
for c in cs_sorted_fba_mc:
    n_comm = 0
    for comm in d_community_members[c]:
        if ("y_F3M07_" not in comm) and ("y_E3R11_" not in comm):
            n_comm += 1
    print(c + ": " + str(n_comm) + " complete minimal supplying communities")

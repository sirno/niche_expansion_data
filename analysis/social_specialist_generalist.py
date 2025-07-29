#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:39:32 2024

@author: magdalena
"""

# import random
import math
from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt

# import xlrd
import misosoup_analysis_module as mym
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns

# from os.path import isdir
# import pickle
# import scipy.cluster.hierarchy as sch
import yaml

# import cobra
from cobra.io import read_sbml_model
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
    path_results_script + "d_misosoup_growth_tolerance_1e-" + str(tolerance) + ".yaml",
    "r",
) as file:
    d_growth_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

with open(
    path_results_script
    + "d_community_members_tolerance_1e-"
    + str(tolerance)
    + ".yaml",
    "r",
) as file:
    d_community_members = yaml.load(file, Loader=yaml.SafeLoader)


###
# every community where a species participates
d_communities = {}
for s in all_species:
    d_communities[s] = []


# first check species in mc
for cs in cs_sorted_fba_mc:
    if len(d_community_members[cs]) != 0:
        for comm in d_community_members[cs]:
            for y_s in comm:
                s = y_s.split("y_")[1]
                d_communities[s].append(comm)


# now species in msc
for cs in cs_sorted_fba_msc:
    for s in d_community_members[cs]:
        if len(d_community_members[cs][s]) != 0:
            for comm in d_community_members[cs][s]:
                if len(comm) > 1:
                    for y_s in comm:
                        s = y_s.split("y_")[1]
                        d_communities[s].append(comm)


###
# distance between communities
def jaccard_similarity(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    # intersection of two sets
    intersection = len(set1.intersection(set2))
    # Unions of two sets
    union = len(set1.union(set2))

    return intersection / union


total_communities = []
jaccard = []

for s in all_species:
    tc = len(d_communities[s])
    total_communities.append(tc)
    if tc != 0:
        jaccard_s = []
        for i1 in range(len(d_communities[s]) - 1):
            comm1 = d_communities[s][i1]
            for i2 in range(1, len(d_communities[s])):
                comm2 = d_communities[s][i2]
                jaccard_s.append(jaccard_similarity(comm1, comm2))

        jaccard.append(np.mean(jaccard_s))
    else:
        jaccard.append(0)


###
# order species according to jaccard distance
jaccard_ordered = [x for x, y in sorted(zip(jaccard, all_species))]
all_species_ordered = [y for x, y in sorted(zip(jaccard, all_species))]
total_communities_ordered = [y for x, y in sorted(zip(jaccard, total_communities))]

###
# plot
plt.plot(range(len(all_species_ordered)), total_communities_ordered, "ko")
plt.xlabel("Species", size=12, family="Arial")
plt.ylabel("Communities", size=12, family="Arial")

plt.show()

###
# plot jaccard similarity
plt.plot(range(len(all_species_ordered)), jaccard_ordered, "ko")
plt.xlabel("Species", size=12, family="Arial")
plt.ylabel("Jaccard similarity", size=12, family="Arial")

plt.show()

###
# histogram haccard similarity
f = plt.figure()
plt.hist(jaccard, 50, color="black", ec="white")
plt.xlim((0.3, 0.85))
plt.xlabel("Jaccard similarity", size=12, family="Arial")
plt.show()
f.savefig(
    path_results_script + "Jaccard_similarity_hist_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# jaccard vs number of communities
plt.plot(total_communities_ordered, jaccard_ordered, "ko")
plt.xlabel("Communities", size=12, family="Arial")
plt.ylabel("Jaccard similarity", size=12, family="Arial")
plt.show()

###
# jaccard vs fundamental niche

inisolation = [0] * len(all_species)
incomm_or_isolation = [0] * len(all_species)

for si, s in enumerate(all_species):
    for cs in cs_sorted_fba:
        if type(d_growth_misosoup[cs][s]) == float:
            inisolation[si] = inisolation[si] + 1
            incomm_or_isolation[si] = incomm_or_isolation[si] + 1
        elif d_growth_misosoup[cs][s] == "in community":  # find comm size
            incomm_or_isolation[si] = incomm_or_isolation[si] + 1

# 100*fn_ordered_rn[si]/len(cs_sorted_fba)

f = plt.figure()
plt.plot([100 * iv / len(cs_sorted_fba) for iv in inisolation], jaccard, "ko")
plt.xlabel("Fundamental niche size", size=12, family="Arial")
plt.ylabel("Jaccard similarity", size=12, family="Arial")
plt.ylim((0.3, 0.85))
plt.xlim((-5, 105))

plt.show()
f.savefig(
    path_results_script + "Jaccard_similarity_vs_fn_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)

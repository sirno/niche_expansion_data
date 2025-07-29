#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:59:28 2025

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

###
path_results_script = "out/analysis/"

###
# load data
# species niche size
df_fn = pd.read_csv("out/analysis/Fundamental_niche.csv")

# alternative communities
df_altcomm = pd.read_csv("out/analysis/Alternatice_communities_species_cs.csv")
df_altcomm.index = df_altcomm["Unnamed: 0"]
df_altcomm = df_altcomm.drop("Unnamed: 0", axis=1)

###
# get total comms for each carbon source
altcomms_for_each_cs = []
f_grow_for_each_cs = []
for cs in df_altcomm.index:
    altcomm_here = sum(df_altcomm.loc[cs].dropna())
    altcomms_for_each_cs.append(altcomm_here)

    f_grow = sum(df_fn[cs]) / len(df_fn[cs])
    f_grow_for_each_cs.append(f_grow)

###
# plot
fig, ax = plt.subplots()

plt.plot(f_grow_for_each_cs, altcomms_for_each_cs, "ko")

plt.ylabel("Alternative communities", size=12, family="Arial")
plt.xlabel("Fraction of species that thrive", size=12, family="Arial")

# Remove axes splines
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

# Add padding between axes and labels
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=10)

# Setting the number of ticks
plt.locator_params(axis="x", nbins=4)
plt.locator_params(axis="y", nbins=4)

plt.show()
fig.savefig(
    "out/analysis/alternative_comms_vs_fraction_that_thrive.pdf", bbox_inches="tight"
)

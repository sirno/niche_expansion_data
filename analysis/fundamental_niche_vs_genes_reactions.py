#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:11:17 2024

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
import yaml

# import cobra
from cobra.io import read_sbml_model

# from os.path import isdir
# import pickle
# import scipy.cluster.hierarchy as sch
from scipy.stats import pearsonr

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

all_species = [s for s in all_species if ((s != "F3M07_") and (s != "E3R11_"))]

print("########################################")
print("## F3M07_ and E3R11_ not included ##")
print("## WORKING WITH " + str(len(all_species)) + " STRAINS ##")
print("########################################")

path_results_script = "out/analysis/"


###
# load results needed
with open(path_results_script + "d_species_info_extended.yaml", "r") as file:
    d_species_info = yaml.load(file, Loader=yaml.SafeLoader)

with open(path_results_script + "d_misosoup_growth_tolerance_1e-7.yaml", "r") as file:
    d_growth_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

# species niche size
df_niche = pd.read_csv("out/analysis/Species_fundamental_realized_niche_size.csv")
df_niche = df_niche.set_index(df_niche["Unnamed: 0"])

###
# fundamental niche
fundamental_niche = df_niche["Fundamental niche size (%)"].loc[all_species]

###
# realized niche
realized_niche = df_niche["Realized niche size (%)"].loc[all_species]

###
# gene/reaction number
sp_reactions = []
sp_genes = []

for s in all_species:
    if s.startswith("m_"):
        s = s.split("m_")[1]

    sp_reactions.append(len(d_species_info[s]["reactions"]))
    sp_genes.append(len(d_species_info[s]["genes"]))


###
# plot for fundamental niche

# reactions
f = plt.figure()
plt.plot(fundamental_niche.values, sp_reactions, "ko")
plt.xlabel("Fundamental niche size (%)", size=12, family="Arial")
plt.ylabel("Reactions", size=12, family="Arial")
plt.show()
f.savefig(path_results_script + "fn_vs_reactions.pdf", bbox_inches="tight")

# genes
f = plt.figure()
plt.plot(fundamental_niche.values, sp_genes, "ko")
plt.xlabel("Fundamental niche size (%)", size=12, family="Arial")
plt.ylabel("Genes", size=12, family="Arial")
plt.show()
f.savefig(path_results_script + "fn_vs_genes.pdf", bbox_inches="tight")


###
# apply the pearsonr()

# reactions
corr, _ = pearsonr(fundamental_niche.values, sp_reactions)
print("Pearsons correlation fn-reactions: %.3f" % corr)  # 0.59

# genes
corr, _ = pearsonr(fundamental_niche.values, sp_genes)
print("Pearsons correlation fn-rgenes: %.3f" % corr)  # 0.46

###
# plot for realized niche

# reactions
f = plt.figure()
plt.plot(realized_niche.values, sp_reactions, "ko")
plt.xlabel("Realized niche size (%)", size=12, family="Arial")
plt.ylabel("Reactions", size=12, family="Arial")
plt.show()
f.savefig(path_results_script + "rn_vs_reactions.pdf", bbox_inches="tight")

# genes
f = plt.figure()
plt.plot(fundamental_niche.values, sp_genes, "ko")
plt.xlabel("Realized niche size (%)", size=12, family="Arial")
plt.ylabel("Genes", size=12, family="Arial")
plt.show()
f.savefig(path_results_script + "rn_vs_genes.pdf", bbox_inches="tight")


###
# apply the pearsonr()

# reactions
corr, _ = pearsonr(realized_niche.values, sp_reactions)
print("Pearsons correlation rn-reactions: %.3f" % corr)  # 0.24

# genes
corr, _ = pearsonr(realized_niche.values, sp_genes)
print("Pearsons correlation rn-rgenes: %.3f" % corr)  # 0.07

###
# plot for niche expansion
niche_expansion_degree = [
    realized_niche.iloc[i] - fundamental_niche.iloc[i] for i in range(len(all_species))
]

# reactions
f = plt.figure()
plt.plot(sp_reactions, niche_expansion_degree, "ko")
plt.ylabel("Niche expansion degree", size=12, family="Arial")
plt.xlabel("Reactions", size=12, family="Arial")
plt.show()
f.savefig(path_results_script + "reactions_vs_niche_expansion.pdf", bbox_inches="tight")

# genes
f = plt.figure()
plt.plot(sp_genes, niche_expansion_degree, "ko")
plt.ylabel("Niche expansion degree", size=12, family="Arial")
plt.xlabel("Genes", size=12, family="Arial")
plt.show()
f.savefig(path_results_script + "genes_vs_niche_expansion.pdf", bbox_inches="tight")


###
# apply the pearsonr()

# reactions
corr, _ = pearsonr(sp_reactions, niche_expansion_degree)
print("Pearsons correlation rn-reactions: %.3f" % corr)  # -0.086

# genes
corr, _ = pearsonr(sp_genes, niche_expansion_degree)
print("Pearsons correlation rn-rgenes: %.3f" % corr)  # -0.223


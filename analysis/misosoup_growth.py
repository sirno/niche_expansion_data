#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 14:35:46 2022

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
path_results_misosoup = "out/misosoup/"


d_growth_misosoup = {}  # for each cs and species, store: 'no growth'/growth value in isolation/'in community'

###
# load misosoup results
tolerance = 7

with open(
    path_results_misosoup + "tolerance_1e-" + str(tolerance) + ".yaml", "r"
) as file:
    d_misosoup = yaml.load(file, Loader=yaml.SafeLoader)


###
# find if species grow in isolation or in communities

# growth in isolation/in community for species in MC
for c in cs_sorted_fba_mc:
    d_growth_misosoup[c] = {}
    species_grow_in_MC = mym.mc_growth_incommunity(d_misosoup, c)
    for s in all_species:
        if s in species_grow_in_MC:
            d_growth_misosoup[c][s] = "in community"
        else:
            d_growth_misosoup[c][s] = "no growth"

# growth in isolation/in community for MSC
for c in cs_sorted_fba_msc:
    if c not in d_growth_misosoup:
        d_growth_misosoup[c] = {}
    for s in d_misosoup[c].keys():
        d_growth_misosoup[c][s] = {}

        does_it_grow = mym.msc_nogrowth_inisolation_incommunity(d_misosoup, c, s)
        if does_it_grow == "in isolation":
            d_growth_misosoup[c][s] = mym.growth_value_in_inisolation(d_misosoup, c, s)
        else:
            d_growth_misosoup[c][s] = does_it_grow

# save dictionary
with open(
    path_results_script + "d_misosoup_growth_tolerance_1e-" + str(tolerance) + ".yaml",
    "w",
) as file:
    documents = yaml.dump(d_growth_misosoup, file)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:36:03 2025

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
import pandas as pd

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

with open(path_results_script + "d_cobra_growth.yaml", "r") as file:
    d_cobra_growth = yaml.load(file, Loader=yaml.Loader)  # Loader=yaml.SafeLoader)

###
# create dictionary
d_species_fundamental_niche = {}
for s in d_cobra_growth:
    d_species_fundamental_niche[s] = []
    for c in cs_sorted_fba:
        if d_cobra_growth[s][c] != 0 and d_cobra_growth[s][c].objective_value > 0.00001:
            d_species_fundamental_niche[s].append(1)
        else:
            d_species_fundamental_niche[s].append(0)

# change it to dataframe
pd_fn = pd.DataFrame.from_dict(
    d_species_fundamental_niche, orient="index", columns=cs_sorted_fba
)

###
# save as csv
pd_fn.to_csv(path_results_script + "Fundamental_niche.csv")


###
# realized niche
tol = 7
with open(
    path_results_script + "/d_misosoup_growth_tolerance_1e-" + str(tol) + ".yaml", "r"
) as file:
    d_growth_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

with open(
    path_results_script + "d_cs_commSize_commNumber_tolerance_1e-" + str(tol) + ".yaml",
    "r",
) as file:
    d_comm_size_and_number = yaml.load(file, Loader=yaml.SafeLoader)

# transform to simpler 0/1 dictionary
d_realized_niche_2stcomm = {}
d_realized_niche_3stcomm = {}
for cs in d_growth_misosoup.keys():
    d_realized_niche_2stcomm[cs] = {}
    d_realized_niche_3stcomm[cs] = {}

    # minimal supplying communities
    for st in d_growth_misosoup[cs].keys():
        if st.startswith("m_"):
            st_name = st.split("m_")[1]
        else:
            st_name = st

        # set to 0
        d_realized_niche_2stcomm[cs][st_name] = 0
        d_realized_niche_3stcomm[cs][st_name] = 0

        if d_growth_misosoup[cs][st] == "in community":
            # minimal communities
            if cs in cs_sorted_fba_mc:
                if 2 in d_comm_size_and_number[cs]:
                    d_realized_niche_2stcomm[cs][st_name] = 1

                if 3 in d_comm_size_and_number[cs]:
                    d_realized_niche_3stcomm[cs][st_name] = 1
            else:
                if 2 in d_comm_size_and_number[cs][st]:
                    d_realized_niche_2stcomm[cs][st_name] = 1

                if 3 in d_comm_size_and_number[cs][st]:
                    d_realized_niche_3stcomm[cs][st_name] = 1

# Convert to DataFrame - 2 strains comm realized niche
pd_rn = pd.DataFrame.from_dict(d_realized_niche_2stcomm, orient="columns")
# order rows as in fn dataframe
pd_rn = pd_rn.reindex(pd_fn.index)
# order columns as in fn dataframe
pd_rn = pd_rn.loc[:, pd_fn.columns]

# save as csv
pd_rn.to_csv(path_results_script + "Realized_niche_2st_comms.csv")


# Convert to DataFrame - 3 strains comm realized niche
pd_rn = pd.DataFrame.from_dict(d_realized_niche_3stcomm, orient="columns")
# order rows as in fn dataframe
pd_rn = pd_rn.reindex(pd_fn.index)
# order columns as in fn dataframe
pd_rn = pd_rn.loc[:, pd_fn.columns]

# save as csv
pd_rn.to_csv(path_results_script + "Realized_niche_3st_comms.csv")


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:18:26 2024

@author: magdalena
"""

# import random
import json
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

df_carbon_sources = pd.read_csv("../data/Carbon_sources_used.csv")
cs_sorted_fba_mc = df_carbon_sources[~df_carbon_sources["Carbon sources CMSC"].isna()][
    "Carbon sources CMSC"
]
cs_sorted_fba_msc = df_carbon_sources[~df_carbon_sources["Carbon sources MSC"].isna()][
    "Carbon sources MSC"
]

cs_sorted_fba = [v for v in cs_sorted_fba_mc] + [v for v in cs_sorted_fba_msc]

# strains
with open("../data/Strains_used", "r") as fp:
    all_species = json.load(fp)


path_results_script = "out/analysis/"
path_results_misosoup = "out/misosoup/"


d_commSize_alternatives = {}  # for each cs and species that grows in community, store the number of alternative communities found of each size
d_tolerances = {}

###
# load results needed

color_green = "#b9d7b1"

for tolerance in [7]:  # ,8,9]:
    # with open(path_results_misosoup+'tolerance_1e-'+str(tolerance)+'.yaml', 'r') as file:
    #    d_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

    with open(
        path_results_script
        + "d_misosoup_growth_tolerance_1e-"
        + str(tolerance)
        + ".yaml",
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

    with open(
        path_results_script
        + "d_cs_commSize_commNumber_tolerance_1e-"
        + str(tolerance)
        + ".yaml",
        "r",
    ) as file:
        d_comm_size_and_number = yaml.load(file, Loader=yaml.SafeLoader)

    ###
    # find maximum community size (number of comm members)
    max_commSize = 0
    for cs in cs_sorted_fba:
        if cs in cs_sorted_fba_mc.values:
            if len(d_comm_size_and_number[cs]) != 0:
                comSize = max(d_comm_size_and_number[cs].keys())
                if comSize > max_commSize:
                    max_commSize = comSize
        elif cs in cs_sorted_fba_msc.values:
            for s in all_species:
                if len(d_comm_size_and_number[cs][s]) != 0:
                    comSize = max(d_comm_size_and_number[cs][s].keys())
                    if comSize > max_commSize:
                        max_commSize = comSize

    ###
    #
    d_commNumb_commSize_cs = {}
    for commsize in range(1, max_commSize + 1):
        d_commNumb_commSize_cs[commsize] = [0] * len(cs_sorted_fba)

    for ics, cs in enumerate(cs_sorted_fba):
        for commsize in range(1, max_commSize + 1):
            if cs in cs_sorted_fba_mc.values:
                if (
                    len(d_comm_size_and_number[cs]) != 0
                    and commsize in d_comm_size_and_number[cs]
                ):
                    d_commNumb_commSize_cs[commsize][ics] = (
                        d_commNumb_commSize_cs[commsize][ics]
                        + d_comm_size_and_number[cs][commsize]
                    )

            if cs in cs_sorted_fba_msc.values:
                for s in all_species:
                    if (
                        len(d_comm_size_and_number[cs][s]) != 0
                        and commsize in d_comm_size_and_number[cs][s]
                    ):
                        d_commNumb_commSize_cs[commsize][ics] = (
                            d_commNumb_commSize_cs[commsize][ics]
                            + d_comm_size_and_number[cs][s][commsize]
                        )

    d_tolerances[tolerance] = d_commNumb_commSize_cs

    ###
    # plot

    f1a = plt.figure()

    p1 = plt.plot(
        range(len(cs_sorted_fba)),
        [
            d_commNumb_commSize_cs[2][i] + d_commNumb_commSize_cs[3][i]
            for i in range(len(d_commNumb_commSize_cs[2]))
        ],
        marker="o",
        color="#000000",
    )
    p2 = plt.plot(
        range(len(cs_sorted_fba)),
        d_commNumb_commSize_cs[2],
        marker="o",
        color=color_green,
    )
    p3 = plt.plot(
        range(len(cs_sorted_fba)),
        d_commNumb_commSize_cs[3],
        marker="o",
        color="black",
        alpha=0.4,
    )  # , color='#909090', edgecolor='white')

    plt.legend(
        (p1[0], p2[0], p3[0]),
        ("total", "2-species communities", "3-species communities"),
    )
    plt.xlabel("Carbon sources", size=12, family="Times New Roman")
    plt.ylabel("Alternative communities", size=12, family="Arial")
    plt.xticks(
        range(len(cs_sorted_fba)),
        labels=cs_sorted_fba,
        rotation="vertical",
        size=12,
        family="Arial",
    )
    # plt.ylim([0,0.7])
    plt.yticks(size=12, family="Arial")
    plt.show()

    f1a.savefig(
        path_results_script + "alternativeComms_vs_cs_" + str(tolerance) + ".pdf",
        bbox_inches="tight",
    )

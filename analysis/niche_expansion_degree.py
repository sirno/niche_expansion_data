#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:18:26 2024

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
for tolerance in [7]:  # ,8,9]:
    print("Tolerance: ", tolerance)

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
    #

    # check if comms found on MC are allways of 2 species
    print("CHECK MCs ARE ALWAYS MADE UP OF 2 SPECIES:")
    for cs in cs_sorted_fba_mc:
        print(d_comm_size_and_number[cs])

    print("----------------")

    inisolation = [0] * len(all_species)
    incomm_size2 = [0] * len(all_species)
    incomm_size3 = [0] * len(all_species)
    incomm_sizeany = [0] * len(all_species)

    for si, s in enumerate(all_species):
        for cs in cs_sorted_fba:
            if type(d_growth_misosoup[cs][s]) == float:
                inisolation[si] = inisolation[si] + 1
            elif d_growth_misosoup[cs][s] == "in community":  # find comm size
                incomm_sizeany[si] = incomm_sizeany[si] + 1

                if cs in cs_sorted_fba_msc:
                    if 2 in d_comm_size_and_number[cs][s]:
                        incomm_size2[si] = incomm_size2[si] + 1

                    elif 3 in d_comm_size_and_number[cs][s]:
                        incomm_size3[si] = incomm_size3[si] + 1

                elif cs in cs_sorted_fba_mc:  ## OJO CON ESTO QUE NO DEJO MANUAL
                    incomm_size2[si] = incomm_size2[si] + 1

    color_green = "#b9d7b1"
    color_orange = "#eca04f"

    ###
    # plot degree of niche expansion
    f = plt.figure()
    bins = np.linspace(0, 100, 22)
    _, _, bars = plt.hist(
        [x * 100 / len(cs_sorted_fba) for ix, x in enumerate(incomm_size2)],
        bins=bins,
        color=color_green,
        ec="white",
        alpha=1,
        label="2 species",
    )
    plt.axvline(
        x=np.mean([x * 100 / len(cs_sorted_fba) for ix, x in enumerate(incomm_size2)]),
        linestyle="--",
        linewidth=1.5,
        color=color_green,
    )

    _, _, bars = plt.hist(
        [x * 100 / len(cs_sorted_fba) for ix, x in enumerate(incomm_sizeany)],
        bins=bins,
        color="black",
        ec="white",
        alpha=0.4,
        label="up to 3 species",
    )

    average_niche_expansion_degree = np.mean(
        [x * 100 / len(cs_sorted_fba) for ix, x in enumerate(incomm_sizeany)]
    )
    plt.axvline(
        x=average_niche_expansion_degree,
        linestyle="--",
        linewidth=1.5,
        color="black",
        alpha=0.4,
    )

    plt.xticks(
        list(range(0, 120, 20)), list(range(0, 120, 20)), size=12, family="Arial"
    )
    plt.yticks([0, 15], [0, 15], size=12, family="Arial")

    plt.xlabel("Niche expansion degree", size=12, family="Arial")
    plt.legend()
    plt.show()
    f.savefig(
        path_results_script
        + "niche_expansion_degree_hist_comm2or3_tol_"
        + str(tolerance)
        + ".pdf",
        bbox_inches="tight",
    )

    ###
    # plot scatter with hists
    def round_to_next5(n):
        return n + (5 - n) % 5

    def scatter_hist(x, y, ax, ax_histx, ax_histy):
        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        # the scatter plot:
        ax.scatter(x, y, c="black")

        # now determine nice limits by hand:
        binwidth = 5
        xymax = 100
        lim = (int(xymax / binwidth) + 1) * binwidth

        bins = np.arange(0, lim, binwidth)
        histx_vals = ax_histx.hist(x, bins=bins, color="black", ec="white")
        histy_vals = ax_histy.hist(
            y, bins=bins, orientation="horizontal", color="black", ec="white"
        )

        ylim = max(max(histx_vals[0]), max(histy_vals[0]))
        ylim = round_to_next5(ylim)
        ax_histx.set_ylim((0, ylim))
        ax_histy.set_xlim((0, ylim))

    # Start with a square Figure.
    fig = plt.figure(figsize=(6, 6))
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal axes and the main axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        left=0.1,
        right=0.9,
        bottom=0.1,
        top=0.9,
        wspace=0.05,
        hspace=0.05,
    )
    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx.set_yticks([0, 20], [0, 20], size=12, family="Arial")
    ax_histy.set_xticks([0, 20], [0, 20], size=12, family="Arial")

    # Draw the scatter plot and marginals.
    x_data = [x * 100 / len(cs_sorted_fba) for x in inisolation]
    iso_and_comm = [
        (inisolation[i] + incomm_sizeany[i]) * 100 / len(cs_sorted_fba)
        for i in range(len(inisolation))
    ]
    y_data = iso_and_comm
    scatter_hist(x_data, y_data, ax, ax_histx, ax_histy)

    # Add diagonal
    ax.plot([0, 100], [0, 100], ":", color="grey")

    ax.set_xticks(
        list(range(0, 120, 20)), list(range(0, 120, 20)), size=12, family="Arial"
    )
    ax.set_yticks(
        list(range(0, 120, 20)), list(range(0, 120, 20)), size=12, family="Arial"
    )

    ax.set_xlabel(
        "Fundamental niche size \n (percentage of env. colonized in isolation)",
        size=12,
        family="Arial",
    )
    ax.set_ylabel(
        "Realized niche size \n (percentage of env. colonized either in isolation or communities)",
        size=12,
        family="Arial",
    )
    plt.show()

    fig.savefig(
        path_results_script + "niche_expansion_degree_tol_" + str(tolerance) + ".pdf",
        bbox_inches="tight",
    )

    ###
    # stats
    print(
        "Average fundamental niche size: "
        + str(np.mean(x_data))
        + "+/-"
        + str(np.std(x_data))
    )
    print(
        str(len([i for i in range(len(all_species)) if x_data[i] == 0]))
        + "species have a fundamental niche size of 0"
    )
    max_fundamental_niche_size = max(x_data)
    print("Largest fundamental niche size: ", max_fundamental_niche_size)
    print(
        str(
            len(
                [
                    i
                    for i in range(len(all_species))
                    if x_data[i] == max_fundamental_niche_size
                ]
            )
        )
        + "species have the largest fundamental niche size"
    )
    print("Average niche expansion degree: ", average_niche_expansion_degree)
    species_no_niche_expansion = [
        i for i in range(len(all_species)) if x_data[i] == y_data[i]
    ]
    species_no_niche_expansion_number = len(species_no_niche_expansion)
    species_no_niche_expansion_percentage = len(species_no_niche_expansion) / len(
        all_species
    )
    print(
        "Species without niche expansion: "
        + str(species_no_niche_expansion_number)
        + " ("
        + str(species_no_niche_expansion_percentage)
        + "%)"
    )
    print(
        "Number of species that dont grow in any carbon source, not even when in communities: ",
        len([i for i in range(len(all_species)) if y_data[i] == 0]),
    )

"""
Tolerance:  7
CHECK MCs ARE ALWAYS MADE UP OF 2 SPECIES:
[]
[]
{2: 2}
{2: 13}
----------------
Average fundamental niche size: 15.689149560117304+/-19.08807092705026
32species have a fundamental niche size of 0
Largest fundamental niche size:  54.54545454545455
2species have the largest fundamental niche size
Average niche expansion degree:  30.645161290322584
Species without niche expansion: 21 (0.3387096774193548%)
Number of species that dont grow in any carbon source, not even when in communities:  20
Tolerance:  8
CHECK MCs ARE ALWAYS MADE UP OF 2 SPECIES:
[]
{2: 4}
[]
{2: 6}
----------------
Average fundamental niche size: 15.689149560117304+/-19.08807092705026
32species have a fundamental niche size of 0
Largest fundamental niche size:  54.54545454545455
2species have the largest fundamental niche size
Average niche expansion degree:  32.47800586510264
Species without niche expansion: 22 (0.3548387096774194%)
Number of species that dont grow in any carbon source, not even when in communities:  21
Tolerance:  9
CHECK MCs ARE ALWAYS MADE UP OF 2 SPECIES:
[]
{2: 2}
{2: 2}
[]
"""


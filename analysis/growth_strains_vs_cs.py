#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 07:22:14 2023

@author: magdalena
"""

# import random
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
# load misosoup results
for tol in [7, 8, 9]:
    print("Tolerance: ", tol)

    with open(
        path_results_script + "/d_misosoup_growth_tolerance_1e-" + str(tol) + ".yaml",
        "r",
    ) as file:
        d_growth_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

    with open(
        path_results_script
        + "d_cs_commSize_commNumber_tolerance_1e-"
        + str(tol)
        + ".yaml",
        "r",
    ) as file:
        d_comm_size_and_number = yaml.load(file, Loader=yaml.SafeLoader)

    ###
    # find species that do not grow in isolation nor in communities
    strains_nogrowth_allcss = []
    for si, species in enumerate(all_species):
        count_cs = 0
        for ic, carbon_source in enumerate(cs_sorted_fba):
            if d_growth_misosoup[carbon_source][species] == "no growth":
                count_cs += 1
        if count_cs == len(cs_sorted_fba):
            strains_nogrowth_allcss.append(species)

    print(
        "Strains that do not grow in isolation nor in communities: ",
        len(strains_nogrowth_allcss),
    )  # 20
    print(
        strains_nogrowth_allcss
    )  # ['m_4F10', 'A2M03', 'm_5C01', 'B3R18', 'D2R18', 'm_4A10', 'C3M10', 'D2R04', 'A3M17', 'B2M17', 'E3M09', 'C3R15', 'I2M19', 'm_4E07', 'G2R14', 'G2M07', 'F2R14', 'A3R06', 'C2R09', 'B2M13']
    strains_do_grow_number = len(all_species) - len(strains_nogrowth_allcss)
    strains_do_grow = [s for s in all_species if s not in strains_nogrowth_allcss]
    print(
        "Strains that grow in isolation or communities: ", strains_do_grow_number
    )  # 42

    #######
    # heatmap
    #######
    M = np.zeros((len(strains_do_grow), len(cs_sorted_fba)))
    for si, s in enumerate(strains_do_grow):
        for ci, c in enumerate(cs_sorted_fba):
            if d_growth_misosoup[c][s] != "no growth":
                M[si, ci] = 1

    f0, ax = plt.subplots(figsize=(30, 25))
    im = ax.imshow(M)
    # add ticks...
    ax.set_yticks(np.arange(len(strains_do_grow)))
    ax.set_xticks(np.arange(len(cs_sorted_fba)))
    # ... and label them with the respective list entries
    ax.set_yticklabels(strains_do_grow)
    ax.set_xticklabels(cs_sorted_fba)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    plt.show()
    f0.savefig(
        path_results_script + "heatmap_growth" + str(tol) + ".pdf", bbox_inches="tight"
    )

    #######
    # plt growth in isolation and in community
    #######

    nogrowth_css = [0] * len(cs_sorted_fba)
    inisolation_css = [0] * len(cs_sorted_fba)
    incommunities_css = [0] * len(cs_sorted_fba)
    incommunities_size2_css = [0] * len(cs_sorted_fba)
    incommunities_size3_css = [0] * len(cs_sorted_fba)

    for ic, carbon_source in enumerate(cs_sorted_fba):
        for species in all_species:
            if (
                d_growth_misosoup[carbon_source][species] == "no growth"
            ):  # it does not grow
                nogrowth_css[ic] += 1
            elif (
                type(d_growth_misosoup[carbon_source][species]) == float
            ):  # grows in isolation
                inisolation_css[ic] += 1
            else:  # grows in community
                incommunities_css[ic] += 1
                # check comm size
                if carbon_source in cs_sorted_fba_mc:  # ALL MC ARE OF SIZE 2
                    incommunities_size2_css[ic] += 1
                else:
                    if 2 in d_comm_size_and_number[carbon_source][species]:
                        incommunities_size2_css[ic] += 1
                    elif 3 in d_comm_size_and_number[carbon_source][species]:
                        incommunities_size3_css[ic] += 1

    ###
    # some stats

    # Average number of strains growing in isolation
    print(
        "Average number of strains growing in isolation: "
        + str(np.mean(inisolation_css))
        + "+/-"
        + str(np.std(inisolation_css))
    )  # 9.333333333333334+/-8.844332774281067
    print(
        "Average number of strains growing in communities: "
        + str(np.mean(incommunities_css))
        + "+/-"
        + str(np.std(incommunities_css))
    )  # 20.75+/-8.926785535678562
    inisolation_or_incommunities_css = [
        inisolation_css[i] + incommunities_css[i] for i in range(len(cs_sorted_fba))
    ]
    print(
        "Average number of strains growing in isolation or communities: "
        + str(np.mean(inisolation_or_incommunities_css))
        + "+/-"
        + str(np.std(inisolation_or_incommunities_css))
    )  # 30.083333333333332+/-9.924366758416154

    # how many strains grow in each cs?
    print("Species growth in isolation/ in communities / in isolation or community")
    for ic, carbon_source in enumerate(cs_sorted_fba):
        print(
            carbon_source,
            inisolation_css[ic],
            incommunities_css[ic],
            inisolation_css[ic] + incommunities_css[ic],
        )

    ###
    # plot fraction considering all species
    f_inisolation_css = [v / len(all_species) for iv, v in enumerate(inisolation_css)]
    f_incommunity_css = [v / len(all_species) for iv, v in enumerate(incommunities_css)]
    f_incommunitysize2_css = [
        v / len(all_species) for iv, v in enumerate(incommunities_size2_css)
    ]
    f_incommunitysize3_css = [
        v / len(all_species) for iv, v in enumerate(incommunities_size3_css)
    ]
    f_inisolation_or_incommunities = [
        v / len(all_species) for iv, v in enumerate(inisolation_or_incommunities_css)
    ]
    inisolation_or_incommunities_css = [
        inisolation_css[i] + incommunities_size2_css[i]
        for i in range(len(cs_sorted_fba))
    ]
    f_inisolation_or_incommunitiessize2 = [
        v / len(all_species) for iv, v in enumerate(inisolation_or_incommunities_css)
    ]

    color_green = "#b9d7b1"
    color_orange = "#eca04fff"

    f1a = plt.figure()
    p1 = plt.bar(
        range(len(cs_sorted_fba)),
        f_inisolation_css,
        color="#000000",
        edgecolor="white",
        alpha=0.7,
    )
    p2 = plt.bar(
        range(len(cs_sorted_fba)),
        f_incommunitysize2_css,
        bottom=f_inisolation_css,
        color=color_green,
        edgecolor="white",
        alpha=0.7,
    )
    p3 = plt.bar(
        range(len(cs_sorted_fba)),
        f_incommunitysize3_css,
        bottom=[
            f_inisolation_css[i] + f_incommunitysize2_css[i]
            for i in range(len(f_inisolation_css))
        ],
        color=color_orange,
        edgecolor="white",
        alpha=0.7,
    )
    # p3=plt.bar(range(len(cs_sorted_fba)), f_incommunity_css, bottom=f_inisolation_css, color='black', edgecolor='white', alpha=0.4)

    plt.legend(
        (p1[0], p2[0], p3[0]),
        ("isolation", "2-species communities", "3-species communities"),
    )

    plt.plot(
        range(len(cs_sorted_fba)),
        [np.mean(f_inisolation_css)] * len(cs_sorted_fba),
        ":",
        color="#000000",
    )
    plt.plot(
        range(len(cs_sorted_fba)),
        [np.mean(f_inisolation_or_incommunitiessize2)] * len(cs_sorted_fba),
        ":",
        color=color_green,
    )
    plt.plot(
        range(len(cs_sorted_fba)),
        [np.mean(f_inisolation_or_incommunities)] * len(cs_sorted_fba),
        ":",
        color="black",
        alpha=0.4,
    )
    plt.xlabel("Carbon sources", size=12, family="Arial")

    plt.ylabel("Fraction of species", size=12, family="Arial")
    plt.xticks(
        range(len(cs_sorted_fba)),
        labels=cs_sorted_fba,
        rotation="vertical",
        size=12,
        family="Arial",
    )
    plt.ylim([0, 0.7])
    plt.yticks([0, 0.35, 0.7], [0, 0.35, 0.7], size=12, family="Arial")
    plt.show()
    f1a.savefig(
        path_results_script + "growth_strains_vs_cs" + str(tol) + ".pdf",
        bbox_inches="tight",
    )

"""
Tolerance:  7
Strains that do not grow in isolation nor in communities:  20
['m_4F10', 'A2M03', 'm_5C01', 'B3R18', 'D2R18', 'm_4A10', 'C3M10', 'D2R04', 'A3M17', 'B2M17', 'E3M09', 'C3R15', 'I2M19', 'm_4E07', 'G2R14', 'G2M07', 'F2R14', 'A3R06', 'C2R09', 'B2M13']
Strains that grow in isolation or communities:  42
Average number of strains growing in isolation: 9.727272727272727+/-9.011012546181368
Average number of strains growing in communities: 18.818181818181817+/-9.537659154570388
Average number of strains growing in isolation or communities: 28.545454545454547+/-12.44525200711746
Species growth in isolation/ in communities / in isolation or community
pep 0 0 0
icit 0 0 0
glu__D 0 2 2
3pg 0 11 11
r5p 1 30 31
met__L 1 27 28
akg 2 29 31
oaa 3 24 27
f6p 4 28 32
phe__L 4 30 34
trp__L 6 28 34
rmn 9 27 36
for 12 22 34
ile__L 12 24 36
rib__D 13 25 38
ser__D 14 18 32
acgam 15 20 35
cit 19 20 39
ac 22 15 37
pyr 24 13 37
glc__D 26 11 37
pro__L 27 10 37
Tolerance:  8
Strains that do not grow in isolation nor in communities:  21
['m_4F10', 'A2M03', 'm_5C01', 'B3R18', 'D2R18', 'D2M02', 'm_4A10', 'C3M10', 'D2R04', 'A3M17', 'B2M17', 'E3M09', 'C3R15', 'I2M19', 'm_4E07', 'G2R14', 'G2M07', 'F2R14', 'A3R06', 'C2R09', 'B2M13']
Strains that grow in isolation or communities:  41
Average number of strains growing in isolation: 9.727272727272727+/-9.011012546181368
Average number of strains growing in communities: 20.136363636363637+/-11.132652442846814
Average number of strains growing in isolation or communities: 29.863636363636363+/-13.133217616360346
Species growth in isolation/ in communities / in isolation or community
pep 0 0 0
icit 0 4 4
glu__D 0 0 0
3pg 0 5 5
r5p 1 35 36
met__L 1 33 34
akg 2 34 36
oaa 3 31 34
f6p 4 31 35
phe__L 4 29 33
trp__L 6 31 37
rmn 9 29 38
for 12 22 34
ile__L 12 26 38
rib__D 13 25 38
ser__D 14 20 34
acgam 15 22 37
cit 19 17 36
ac 22 14 36
pyr 24 14 38
glc__D 26 11 37
pro__L 27 10 37
Tolerance:  9
Strains that do not grow in isolation nor in communities:  20
['m_4F10', 'A2M03', 'm_5C01', 'B3R18', 'D2R18', 'm_4A10', 'C3M10', 'D2R04', 'A3M17', 'B2M17', 'E3M09', 'C3R15', 'I2M19', 'm_4E07', 'G2R14', 'G2M07', 'F2R14', 'A3R06', 'C2R09', 'B2M13']
Strains that grow in isolation or communities:  42
Average number of strains growing in isolation: 9.727272727272727+/-9.011012546181368
Average number of strains growing in communities: 20.454545454545453+/-11.857974964045667
Average number of strains growing in isolation or communities: 30.181818181818183+/-13.71010365200644
Species growth in isolation/ in communities / in isolation or community
pep 0 0 0
icit 0 3 3
glu__D 0 2 2
3pg 0 0 0
r5p 1 35 36
met__L 1 34 35
akg 2 35 37
oaa 3 35 38
f6p 4 31 35
phe__L 4 32 36
trp__L 6 32 38
rmn 9 27 36
for 12 22 34
ile__L 12 25 37
rib__D 13 25 38
ser__D 14 21 35
acgam 15 24 39
cit 19 20 39
ac 22 14 36
pyr 24 12 36
glc__D 26 12 38
pro__L 27 9 36

"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:01:49 2025

@author: magdalena
"""

import itertools
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats

cs_sorted_fba_mc = [
    "pep",
    "icit",
    "glu__D",
    "3pg",
    # 'cyst__L',
]

# carbon sources sorted according to increasing number of strains that are
# viable in isolation (obtained with FBA)
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
    # 'cys__L'
]

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

# phylogenetic distance
df_phylo_dis = pd.read_csv("../data/phylogenetic_distance_matrix.csv")


###
# misosoup communities as list

Comms = []
for c in d_community_members:
    if c in cs_sorted_fba_mc:
        for comm in d_community_members[c]:
            Comms.append(comm)
    elif c in cs_sorted_fba_msc:
        for sp in d_community_members[c]:
            for comm in d_community_members[c][sp]:
                Comms.append(comm)


Comms.sort()  # 2424
# check numb. comms without F3M07_ and E3R11_
n_comms_F3M07__or_E3R11_ = 0
for com in Comms:
    if ("y_F3M07_" in com) or ("y_E3R11_" in com):
        n_comms_F3M07__or_E3R11_ += 1

print(
    "Total communities excluding F3M07_ and E3R11_: ",
    len(Comms) - n_comms_F3M07__or_E3R11_,
)
# 2318

Comms = list(Comms for Comms, _ in itertools.groupby(Comms))  # 777 unique communities


###
# check how many comms there are with E3R11_, F3M07_
n_F3M07_ = 0  # 38
n_F3M07 = 0  # 32
n_E3R11_ = 0  # 4
n_E3R11 = 0  # 48

for com in Comms:
    if "y_F3M07_" in com:
        n_F3M07_ += 1
    if "y_F3M07" in com:
        n_F3M07 += 1
    if "y_E3R11_" in com:
        n_E3R11_ += 1
    if "y_E3R11" in com:
        n_E3R11 += 1


###
# phylogenetic distance in misosoup communities

phyl_dis_C2 = []
phyl_dis_C3 = []
comms_not_included = 0
for comm in Comms:
    if len(comm) == 2:
        # change species names to match names used in phylo matrix
        sp1 = comm[0].split("y_")[1]
        sp2 = comm[1].split("y_")[1]
        if sp1.startswith("m_"):
            sp1 = sp1.split("m_")[1]

        if sp2.startswith("m_"):
            sp2 = sp2.split("m_")[1]

        # are species in phylo matrix?
        if (sp1 in df_phylo_dis["Unnamed: 0"].values) and (
            sp2 in df_phylo_dis["Unnamed: 0"].values
        ):
            dis = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp2][sp1].iloc[0]
            phyl_dis_C2.append(dis)
        else:
            print("Community not included in analysis: ", comm)
            comms_not_included += 1

    if len(comm) == 3:
        # change species names to match names used in phylo matrix
        sp1 = comm[0].split("y_")[1]
        sp2 = comm[1].split("y_")[1]
        sp3 = comm[2].split("y_")[1]

        if sp1.startswith("m_"):
            sp1 = sp1.split("m_")[1]
        if sp2.startswith("m_"):
            sp2 = sp2.split("m_")[1]
        if sp3.startswith("m_"):
            sp3 = sp3.split("m_")[1]

        # are species in phylo matrix?
        if (
            (sp1 in df_phylo_dis["Unnamed: 0"].values)
            and (sp2 in df_phylo_dis["Unnamed: 0"].values)
            and (sp3 in df_phylo_dis["Unnamed: 0"].values)
        ):
            dis12 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp2][sp1].iloc[0]
            dis13 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp3][sp1].iloc[0]
            dis23 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp3][sp2].iloc[0]
            phyl_dis_C3.append(np.mean([dis12, dis13, dis23]))

        else:
            print("Community not included in analysis: ", comm)
            comms_not_included += 1


###
# random communities - fundamental niche union-intersection

RandC2 = list(itertools.combinations(all_species, 2))  # 1891
RandC3 = list(itertools.combinations(all_species, 3))  # 37820

# select 500 comms at random
RandC2 = random.sample(RandC2, 500)
RandC3 = random.sample(RandC3, 500)

###
# phylogenetic distance in random communities
phyl_dis_randC2 = []
phyl_dis_randC3 = []

commsR_not_included = 0
for comm in RandC2:
    sp1 = comm[0]
    sp2 = comm[1]
    if sp1.startswith("m_"):
        sp1 = sp1.split("m_")[1]

    if sp2.startswith("m_"):
        sp2 = sp2.split("m_")[1]

    # are species in phylo matrix?
    if (sp1 in df_phylo_dis["Unnamed: 0"].values) and (
        sp2 in df_phylo_dis["Unnamed: 0"].values
    ):
        dis = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp2][sp1].iloc[0]
        phyl_dis_randC2.append(dis)
    else:
        commsR_not_included += 1


for comm in RandC3:
    sp1 = comm[0]
    sp2 = comm[1]
    sp3 = comm[2]

    if sp1.startswith("m_"):
        sp1 = sp1.split("m_")[1]
    if sp2.startswith("m_"):
        sp2 = sp2.split("m_")[1]
    if sp3.startswith("m_"):
        sp3 = sp3.split("m_")[1]

    # are species in phylo matrix?
    if (
        (sp1 in df_phylo_dis["Unnamed: 0"].values)
        and (sp2 in df_phylo_dis["Unnamed: 0"].values)
        and (sp3 in df_phylo_dis["Unnamed: 0"].values)
    ):
        dis12 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp2][sp1].iloc[0]
        dis13 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp3][sp1].iloc[0]
        dis23 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == sp3][sp2].iloc[0]
        phyl_dis_randC3.append(np.mean([dis12, dis13, dis23]))

    else:
        commsR_not_included += 1


"""
F3M07_
E3R11_
not included in phylo matrix. The analysis therefore doesnt include the phylo distance for
42 communities found with misosoup and 88 random communities
"""


###
# all possible phylo distances with the strains used

d_tree = {}
d_tree["Arcobacteraceae"] = ["A1R12"]
d_tree["Cyclobacteriaceae"] = ["E3R18"]
d_tree["Flavobacteriaceae"] = [
    "E3R01",
    "B2M17",
    "I2M19",
    "C3R15",
    "E3R09",
    "D3R19",
    "B2M06",
    "E3M18",
    "A2M03",
    "C3R17",
    "B3R18",
]
d_tree["Rhodobacteraceae"] = [
    "C2R09",
    "C3M06",
    "4E07",
    "5C01",
    "C3M10",
    "A3M17",
    "4A10",
    "A3R06",
    "4F10",
    "G2R14",
    "D2R18",
    "D2R04",
    "F2R14",
    "G3M19",
]
d_tree["Vibrionaceae"] = ["3C02", "G2R10", "1A01", "I3M07", "C3R12", "6D03", "E3R11"]
d_tree["Psychromonadaceae"] = ["6C06", "B3M02"]
d_tree["Psychrobiaceae"] = ["6D02"]
d_tree["Alteromonadaceae"] = [
    "B3R10",
    "I2R16",
    "F3M07",
    "3D05",
    "C2R02",
    "D2M02",
    "E2M01",
    "C2M11",
    "A3R12",
    "A3R16",
    "D2R05",
    "A3R04",
    "C1M14",
    "4B03",
]
d_tree["Oleiphilaceae"] = ["D2M19", "F3R11", "F3R08"]
d_tree["Halomonadaceae"] = ["B2M13"]
d_tree["Nitrincolaceae"] = ["3B05", "C1R06"]
d_tree["Saccharospirillaceae"] = ["G2M07"]
d_tree["Cellvibrionaceae"] = ["E3M09", "E3M17"]

# d_genus = {}
# d_tree['Arcobacteraceae'] = ['A1R12']
# d_tree['Cyclobacteriaceae'] = ['E3R18']
# d_tree['Flavobacteriaceae'] = ['E3R01','B2M17','I2M19','C3R15','E3R09','D3R19','B2M06','E3M18','A2M03','C3R17','B3R18']
# d_tree['Rhodobacteraceae'] = ['C2R09','C3M06','4E07','5C01','C3M10','A3M17','4A10','A3R06','4F10','G2R14','D2R18','D2R04','F2R14','G3M19']
# d_tree['Vibrionaceae'] = ['3C02','G2R10','1A01','I3M07','C3R12','6D03','E3R11']
# d_tree['Psychromonadaceae'] = ['6C06','B3M02']
# d_tree['Psychrobiaceae'] = ['6D02']
# d_tree['Alteromonadaceae'] = ['B3R10','I2R16','F3M07','3D05','C2R02','D2M02','E2M01','C2M11','A3R12','A3R16','D2R05','A3R04','C1M14','4B03']
# d_tree['Oleiphilaceae'] = ['D2M19','F3R11','F3R08']
# d_tree['Halomonadaceae'] = ['B2M13']
# d_tree['Nitrincolaceae'] = ['3B05','C1R06']
# d_tree['Saccharospirillaceae'] = ['G2M07']
# d_tree['Cellvibrionaceae']

d_distance = {}

# phylo distance between strains of same family
for fam in d_tree:
    # if there are at  least 2 strains of this family, calculate phylo distance
    if len(d_tree[fam]) > 1:
        d_distance[fam] = []
        for ix_st1 in range(len(d_tree[fam]) - 1):
            for ix_st2 in range(ix_st1 + 1, len(d_tree[fam])):
                st1 = d_tree[fam][ix_st1]
                st2 = d_tree[fam][ix_st2]

                dis12 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == st2][st1].iloc[0]
                d_distance[fam].append(dis12)


# phylo distance between members of different families
d_distance["between families"] = []

for ix_f1 in range(len(d_tree) - 1):
    fam1 = list(d_tree.keys())[ix_f1]
    for ix_f2 in range(ix_f1 + 1, len(d_tree)):
        fam2 = list(d_tree.keys())[ix_f2]
        for st1 in d_tree[fam1]:
            for st2 in d_tree[fam2]:
                dis12 = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == st2][st1].iloc[0]
                d_distance["between families"].append(dis12)


###
# plot

color_green = "#b9d7b1"
color_orange = "#eca04f"

f = plt.figure()

data = []
for k in d_distance:
    data.append(d_distance[k])


data.append(phyl_dis_randC2)
data.append(phyl_dis_C2)

ax = plt.boxplot(
    data,
    whis=[5, 95],
    patch_artist=True,
    medianprops=dict(color="black"),
    meanline=True,
    showmeans=True,
    meanprops=dict(color="black"),
)

colors = ["#D3D3D3"] * (len(d_distance.keys()) - 1) + [
    "#B0B0B0",
    "#696969",
    color_green,
]
for patch, color in zip(ax["boxes"], colors):
    patch.set_facecolor(color)

plt.xticks(
    range(1, len(data) + 1),
    labels=list(d_distance.keys()) + ["2-strain random", "2-strain misosoup"],
    rotation="vertical",
)

plt.ylabel("Phylogenetic distance", size=12, family="Arial")
plt.show()

f.savefig(
    path_results_script
    + "phylo_distance_comms_2st_rand_misosoup_tol_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


###
# plot

color_green = "#b9d7b1"
color_orange = "#eca04f"

f = plt.figure()

data = [phyl_dis_C2, phyl_dis_randC2, phyl_dis_C3, phyl_dis_randC3]
ax = plt.boxplot(
    data,
    whis=[5, 95],
    patch_artist=True,
    medianprops=dict(color="black"),
    meanline=True,
    showmeans=True,
    meanprops=dict(color="black"),
)

colors = [color_green, "grey", color_orange, "grey"]
for patch, color in zip(ax["boxes"], colors):
    patch.set_facecolor(color)

plt.xticks([1.5, 3.5], labels=["2", "3"])
plt.xlabel("Species in community", size=12, family="Arial")
plt.ylabel("Phylogenetic distance", size=12, family="Arial")
plt.show()

f.savefig(
    path_results_script
    + "phylo_distance_comms_rand_misosoup_tol_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)

# Welch's t-test
#######

# Welch’s t-test is a nonparametric univariate test that tests for a significant difference between the mean of two unrelated groups.

print("Test differences in phylogenetic distance between random and misosoup comms")

sp_comm = {}
sp_comm[2] = [phyl_dis_C2, phyl_dis_randC2]
sp_comm[3] = [phyl_dis_C3, phyl_dis_randC3]

for comm_size in [2, 3]:
    print("Welchs t-test for " + str(comm_size) + "-sp communities")

    # Perform Welch's t-test
    sample1 = sp_comm[comm_size][0]
    sample2 = sp_comm[comm_size][1]
    t_stat, p_value = stats.ttest_ind(sample1, sample2, equal_var=False)

    # One-sided test (Check if sample1 < sample2)
    if t_stat < 0:  # Check if the mean of sample1 is less than sample2
        p_value_one_sided = p_value / 2
    else:
        p_value_one_sided = 1 - (p_value / 2)

    print(f"T-statistic: {t_stat}")
    print(f"One-sided p-value: {p_value_one_sided}")

    # Interpretation
    alpha = 0.05  # Significance level
    if p_value_one_sided < alpha:
        print(
            "Reject the null hypothesis: sample1 is significantly smaller than sample2."
        )
    else:
        print(
            "Fail to reject the null hypothesis: No significant evidence that sample1 is smaller."
        )

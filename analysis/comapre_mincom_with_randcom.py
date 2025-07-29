#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:12:09 2024

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

with open(path_results_script + "d_cobra_growth.yaml", "r") as file:
    d_cobra_growth = yaml.load(file, Loader=yaml.Loader)  # Loader=yaml.SafeLoader)

d_species_fundamental_niche = {}
for s in d_cobra_growth:
    d_species_fundamental_niche[s] = []
    for c in d_cobra_growth[s]:
        if d_cobra_growth[s][c] != 0 and d_cobra_growth[s][c].objective_value > 0.00001:
            d_species_fundamental_niche[s].append(c)

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
            Comms.append(comm)
    elif c in cs_sorted_fba_msc:
        for sp in d_community_members[c]:
            for comm in d_community_members[c][sp]:
                Comms.append(comm)


Comms.sort()  # 2424
Comms = list(Comms for Comms, _ in itertools.groupby(Comms))  # 777 unique communities


###
# fundamental niche union-intersection


def fundamental_niche_dic_for_comm(comm, d_species_fundamental_niche):
    d_here = {}

    first_species = comm[0]

    if first_species.startswith("y_"):
        for y_sp in comm:
            sp = y_sp.split("y_")[1]

            if "m_" in sp:
                sp = sp.split("m_")[1]
            d_here[sp] = d_species_fundamental_niche[sp]
    else:
        for sp in comm:
            if "m_" in sp:
                sp = sp.split("m_")[1]
            d_here[sp] = d_species_fundamental_niche[sp]

    return d_here


def union_fundamental_niches(d_here):
    union = []
    for sp in d_here:
        for cs in d_here[sp]:
            if cs not in union:
                union.append(cs)

    return union


def intersection_fundamental_niches(d_here):
    intersection = []
    cs_for_every_species = True
    first_species = list(d_here.keys())[0]
    for cs in d_here[first_species]:
        for s in d_here:
            if cs not in d_here[s]:
                cs_for_every_species = False

        if cs_for_every_species == True:
            intersection.append(cs)

    return intersection


###
# MC/MCS fundamental niche union-intersection

# 496 2-species communities (tol 7)
fn_union_2sp = []
fn_intersection_2sp = []
fn_union_minus_intersection_2sp = []

# 251 3-species communities (tol 7)
fn_union_3sp = []
fn_intersection_3sp = []
fn_union_minus_intersection_3sp = []

for comm in Comms:
    # transitory dic with the fundamental niche of the species in comm
    d_here = fundamental_niche_dic_for_comm(comm, d_species_fundamental_niche)

    # union
    union = union_fundamental_niches(d_here)

    # intersection
    intersection = intersection_fundamental_niches(d_here)

    # union-intersection
    if len(comm) == 2:
        fn_union_2sp.append(len(union))
        fn_intersection_2sp.append(len(intersection))
        fn_union_minus_intersection_2sp.append(len(union) - len(intersection))
    elif len(comm) == 3:
        fn_union_3sp.append(len(union))
        fn_intersection_3sp.append(len(intersection))
        fn_union_minus_intersection_3sp.append(len(union) - len(intersection))


###
# random communities - fundamental niche union-intersection
RandC2 = list(itertools.combinations(all_species, 2))
RandC3 = list(itertools.combinations(all_species, 3))

# select 500 comms at random
RandC2 = random.sample(RandC2, 500)
RandC3 = random.sample(RandC3, 500)

fn_union_2sprand = []
fn_intersection_2sprand = []
fn_union_minus_intersection_2sprand = []

fn_union_3sprand = []
fn_intersection_3sprand = []
fn_union_minus_intersection_3sprand = []

for comm in RandC2:
    # transitory dic with the fundamental niche of the species in comm
    d_here = fundamental_niche_dic_for_comm(comm, d_species_fundamental_niche)

    # union
    union = union_fundamental_niches(d_here)

    # intersection
    intersection = intersection_fundamental_niches(d_here)

    # union-intersection
    fn_union_2sprand.append(len(union))
    fn_intersection_2sprand.append(len(intersection))
    fn_union_minus_intersection_2sprand.append(len(union) - len(intersection))

for comm in RandC3:
    # transitory dic with the fundamental niche of the species in comm
    d_here = fundamental_niche_dic_for_comm(comm, d_species_fundamental_niche)

    # union
    union = union_fundamental_niches(d_here)

    # intersection
    intersection = intersection_fundamental_niches(d_here)

    # union-intersection
    fn_union_3sprand.append(len(union))
    fn_intersection_3sprand.append(len(intersection))
    fn_union_minus_intersection_3sprand.append(len(union) - len(intersection))


###
# plot

color_green = "#b9d7b1"
color_orange = "#eca04f"


def round_to_next5(n):
    return n + (5 - n) % 5


def scatter_hist(x, y, ax, ax_histx, ax_histy, color_plot, transparency_plot):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, c=color_plot, alpha=transparency_plot)
    # ax.set_xlim((0,80))
    # ax.set_ylim((0,60))

    # now determine nice limits by hand:
    binwidth = 5
    xymax = 80
    lim = (int(xymax / binwidth) + 1) * binwidth

    bins = np.arange(0, lim, binwidth)
    histx_vals = ax_histx.hist(
        x, bins=bins, color=color_plot, ec="white", alpha=transparency_plot
    )
    ax_histx.axvline(np.mean(x), linestyle="--", linewidth=1.5, color=color_plot)

    histy_vals = ax_histy.hist(
        y,
        bins=bins,
        orientation="horizontal",
        color=color_plot,
        ec="white",
        alpha=transparency_plot,
    )
    ax_histy.axhline(np.mean(y), linestyle="--", linewidth=1.5, color=color_plot)

    # ylim = max(max(histx_vals[0]), max(histy_vals[0]))
    # ylim = round_to_next5(ylim)
    # ax_histx.set_xlim((0,80))
    # ax_histx.set_ylim((0,ylim))
    # ax_histy.set_xlim((0,ylim))
    # ax_histy.set_ylim((0,60))


###
# plot comms 2-sp
# Start with a square Figure.
fig = plt.figure(figsize=(6, 6))
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
# ax_histx.set_yticks([0,20],[0,20],size=12, family='Arial')
# ax_histy.set_xticks([0,20],[0,20],size=12, family='Arial')

# Draw the scatter plot and marginals.

# 2-sp random comm
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_2sprand]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sprand]
scatter_hist(x_data, y_data, ax, ax_histx, ax_histy, "black", 0.5)

# 2-sp MC/MSC comm
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_2sp]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sp]
scatter_hist(x_data, y_data, ax, ax_histx, ax_histy, color_green, 0.5)


ax.set_xticks(list(range(0, 100, 20)), list(range(0, 100, 20)))
ax.set_yticks(list(range(0, 100, 20)), list(range(0, 100, 20)))

ax.set_xlabel("Fundamental niche union (%)", size=12, family="Arial")
ax.set_ylabel("Fundamental niche intersection (%)", size=12, family="Arial")
plt.show()

fig.savefig(
    path_results_script
    + "fn_intersection_vs_union_comms_2sp_tol_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


###
# plot comms 3-sp
# Start with a square Figure.
fig = plt.figure(figsize=(6, 6))
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
# ax_histx.set_yticks([0,20],[0,20],size=12, family='Arial')
# ax_histy.set_xticks([0,20],[0,20],size=12, family='Arial')

# Draw the scatter plot and marginals.

# 3-sp random comm
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_3sprand]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sprand]
scatter_hist(x_data, y_data, ax, ax_histx, ax_histy, "black", 0.5)

# 3-sp MC/MSC comm
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_3sp]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sp]
scatter_hist(x_data, y_data, ax, ax_histx, ax_histy, color_orange, 0.5)


ax.set_xticks(list(range(0, 100, 20)), list(range(0, 100, 20)), size=12, family="Arial")
ax.set_yticks(list(range(0, 100, 20)), list(range(0, 100, 20)), size=12, family="Arial")

ax.set_xlabel("Fundamental niche union (%)", size=12, family="Arial")
ax.set_ylabel("Fundamental niche intersection (%)", size=12, family="Arial")
plt.show()

fig.savefig(
    path_results_script
    + "fn_intersection_vs_union_comms_3sp_tol_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


###
# plot union-intersection

data = [
    [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_2sp],
    [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_2sprand],
    [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_3sp],
    [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_3sprand],
]

f = plt.figure()
ax = plt.boxplot(data, whis=[5, 95], patch_artist=True, medianprops=dict(color="black"))

colors = [color_green, "grey", color_orange, "grey"]
for patch, color in zip(ax["boxes"], colors):
    patch.set_facecolor(color)

plt.xticks(range(1, 5), labels=["2", "2", "3", "3"])
plt.xlabel("Species in community", size=12, family="Arial")
plt.ylabel("Fundamental niche union - intersection", size=12, family="Arial")
plt.show()

f.savefig(
    path_results_script + "fn_union_minus_intersection_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# plot union and intersection separately
f = plt.figure(figsize=(8, 4), dpi=80)
gs = gridspec.GridSpec(1, 2)

##
# union
ax_0 = plt.subplot(gs[0, 0])

data = [
    [x * 100 / len(cs_sorted_fba) for x in fn_union_2sp],
    [x * 100 / len(cs_sorted_fba) for x in fn_union_2sprand],
    [x * 100 / len(cs_sorted_fba) for x in fn_union_3sp],
    [x * 100 / len(cs_sorted_fba) for x in fn_union_3sprand],
]

ax = plt.boxplot(data, whis=[5, 95], patch_artist=True, medianprops=dict(color="black"))

colors = [color_green, "grey", color_orange, "grey"]
for patch, color in zip(ax["boxes"], colors):
    patch.set_facecolor(color)

plt.xticks(range(1, 5), labels=["2", "2", "3", "3"])
plt.xlabel("Species in community", size=12, family="Arial")
plt.ylabel("Fundamental niche union", size=12, family="Arial")

##
# intersection
ax_1 = plt.subplot(gs[0, 1])

data = [
    [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sp],
    [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sprand],
    [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sp],
    [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sprand],
]

ax = plt.boxplot(data, whis=[5, 95], patch_artist=True, medianprops=dict(color="black"))

colors = [color_green, "grey", color_orange, "grey"]
for patch, color in zip(ax["boxes"], colors):
    patch.set_facecolor(color)

plt.xticks(range(1, 5), labels=["2", "2", "3", "3"])
plt.xlabel("Species in community", size=12, family="Arial")
plt.ylabel("Fundamental niche intersection", size=12, family="Arial")
plt.show()

f.savefig(
    path_results_script + "fn_union_and_intersection_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# plot symetric difference vs union/intersection
f = plt.figure(figsize=(8, 8), dpi=80)
gs = gridspec.GridSpec(2, 2)

ax_00 = plt.subplot(gs[0, 0])
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_2sp]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_2sp]
plt.plot(x_data, y_data, "o", color=color_green, alpha=0.4)
# plt.xlabel('Fundamental niche union - intersection', size=12, family='Arial')
plt.ylabel("Fundamental niche \n union - intersection (%)", size=12, family="Arial")

ax_01 = plt.subplot(gs[0, 1], sharey=ax_00)
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sp]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_2sp]
plt.plot(x_data, y_data, "o", color=color_green, alpha=0.4)
# plt.xlabel('Fundamental niche union - intersection', size=12, family='Arial')
# plt.ylabel('Fundamental niche union', size=12, family='Arial')

ax_10 = plt.subplot(gs[1, 0], sharex=ax_00)
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_3sp]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_3sp]
plt.plot(x_data, y_data, "o", color=color_orange)
# plt.xlabel('Fundamental niche union - intersection', size=12, family='Arial')
plt.xlabel("Fundamental niche \n union (%)", size=12, family="Arial")
plt.ylabel("Fundamental niche \n union - intersection (%)", size=12, family="Arial")

ax_11 = plt.subplot(gs[1, 1], sharex=ax_01, sharey=ax_10)
x_data = [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sp]
y_data = [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_3sp]
plt.plot(x_data, y_data, "o", color=color_orange, alpha=0.4)
plt.xlabel("Fundamental niche \n intersection (%)", size=12, family="Arial")
# plt.ylabel('Fundamental niche intersection', size=12, family='Arial')

plt.show()
f.savefig(
    path_results_script
    + "fn_union_minus_intersection_vs_union-intersection_tol_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


#######
# Welch's t-test
#######

# Welch’s t-test is a nonparametric univariate test that tests for a significant difference between the mean of two unrelated groups.

# union-intersection
print("Test differences in union-intersection")
print("Welchs t-test for 2-sp communities")
print(
    stats.ttest_ind(
        [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_2sp],
        [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_2sprand],
        equal_var=False,
    )
)

print("Welchs t-test for 3-sp communities")
print(
    stats.ttest_ind(
        [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_3sp],
        [x * 100 / len(cs_sorted_fba) for x in fn_union_minus_intersection_3sprand],
        equal_var=False,
    )
)

print("--------------")

# union
print("Test differences in union")
print("Welchs t-test for 2-sp communities")
print(
    stats.ttest_ind(
        [x * 100 / len(cs_sorted_fba) for x in fn_union_2sp],
        [x * 100 / len(cs_sorted_fba) for x in fn_union_2sprand],
        equal_var=False,
    )
)

print("Welchs t-test for 3-sp communities")
print(
    stats.ttest_ind(
        [x * 100 / len(cs_sorted_fba) for x in fn_union_3sp],
        [x * 100 / len(cs_sorted_fba) for x in fn_union_3sprand],
        equal_var=False,
    )
)

print("--------------")

# intersection
print("Test differences in intersection")
print("Welchs t-test for 2-sp communities")
print(
    stats.ttest_ind(
        [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sp],
        [x * 100 / len(cs_sorted_fba) for x in fn_intersection_2sprand],
        equal_var=False,
    )
)

print("Welchs t-test for 3-sp communities")
print(
    stats.ttest_ind(
        [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sp],
        [x * 100 / len(cs_sorted_fba) for x in fn_intersection_3sprand],
        equal_var=False,
    )
)

##
# tolerance 7
"""
Tolerance:  7
Welchs t-test for 2-sp communities
Ttest_indResult(statistic=11.308658424603177, pvalue=9.089043394014926e-28)
Welchs t-test for 3-sp communities
Ttest_indResult(statistic=8.053990974108574, pvalue=4.3723465296153475e-15)
"""

# tolerance 8
"""
Tolerance:  8
Welchs t-test for 2-sp communities
Ttest_indResult(statistic=9.81457896576191, pvalue=1.2923111323056841e-21)
Welchs t-test for 3-sp communities
Ttest_indResult(statistic=9.455213865909505, pvalue=2.8674783052193465e-20)
"""

# tolerance 9
"""
Tolerance:  9
Welchs t-test for 2-sp communities
Ttest_indResult(statistic=10.289010482574826, pvalue=1.646930661791367e-23)
Welchs t-test for 3-sp communities
Ttest_indResult(statistic=10.600941227740815, pvalue=1.0582251942430328e-24)
"""


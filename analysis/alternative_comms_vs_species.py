#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:13:23 2024

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

d_niche = {}

inisolation = [0] * len(all_species)
incomm_or_isolation = [0] * len(all_species)

for si, s in enumerate(all_species):
    d_niche[s] = {}
    for cs in cs_sorted_fba:
        if type(d_growth_misosoup[cs][s]) == float:
            inisolation[si] = inisolation[si] + 1
            incomm_or_isolation[si] = incomm_or_isolation[si] + 1
        elif d_growth_misosoup[cs][s] == "in community":  # find comm size
            incomm_or_isolation[si] = incomm_or_isolation[si] + 1

    d_niche[s]["Fundamental niche size (%)"] = (
        inisolation[si] * 100 / len(cs_sorted_fba)
    )
    d_niche[s]["Realized niche size (%)"] = (
        incomm_or_isolation[si] * 100 / len(cs_sorted_fba)
    )

###
# transform to pandas db and save
df_niche = pd.DataFrame.from_dict(d_niche)
df_niche.T.to_csv(
    "out/analysis/Species_fundamental_realized_niche_size.csv", index=True
)

###

niche_expansion_degree = [
    (incomm_or_isolation[s] - inisolation[s]) / len(cs_sorted_fba)
    for s in range(len(all_species))
]


# order species according to realized niche size
species_ordered_rn = [x for y, x in sorted(zip(incomm_or_isolation, all_species))]
ordered_rn = [y for y, x in sorted(zip(incomm_or_isolation, all_species))]
fn_ordered_rn = [x for y, x in sorted(zip(incomm_or_isolation, inisolation))]

# keep species with niche expansion (to make sure there are communities)
species_niche_expansion = [
    s for si, s in enumerate(all_species) if niche_expansion_degree[si] != 0
]

###
# alternative comms for species
d_alternative_comms = {}

# add alternative comms from MSC
for s in species_niche_expansion:
    d_alternative_comms[s] = {}
    for c in cs_sorted_fba_msc:
        if 2 in d_comm_size_and_number[c][s] or 3 in d_comm_size_and_number[c][s]:
            d_alternative_comms[s][c] = 0
            for com_size in d_comm_size_and_number[c][s]:
                if com_size != 1:  # sum communities of 2 or 3 species
                    d_alternative_comms[s][c] = (
                        d_alternative_comms[s][c]
                        + d_comm_size_and_number[c][s][com_size]
                    )


alternative_msc = []
for s in d_alternative_comms:
    for c in d_alternative_comms[s]:
        alternative_msc.append(d_alternative_comms[s][c])

print(
    "The number of alternative MSC for each environment and species varied between "
    + str(min(alternative_msc))
    + " and "
    + str(max(alternative_msc))
    + ", with a mean of "
    + str(np.mean(alternative_msc))
    + "+/-"
    + str(np.std(alternative_msc))
)
"The number of alternative MSC for each environment and species varied between 1 and 37, with a mean of 5.473815461346633+/-5.382611484204079"


# add alternative comms from MC
for cs in cs_sorted_fba_mc:
    if cs in d_community_members:
        # go through the comms found on that cs
        for comm in d_community_members[cs]:
            # get species in that comm
            for s in comm:
                s_name = s.split("y_")[1]
                if s_name not in d_alternative_comms:
                    d_alternative_comms[s_name] = {}
                if cs not in d_alternative_comms[s_name]:
                    d_alternative_comms[s_name][cs] = 0
                d_alternative_comms[s_name][cs] += 1


###
# transform to pandas db
df_alternative_comms = pd.DataFrame.from_dict(d_alternative_comms)
df_alternative_comms.to_csv(
    "out/analysis/Alternatice_communities_species_cs.csv", index=True
)

# order dataframe
tac_species = df_alternative_comms.sum(
    axis=0
)  # total alternative communities per species
tac_species_order = tac_species.sort_values().index
df_alternative_comms = df_alternative_comms[tac_species_order]

tac_cs = df_alternative_comms.sum(
    axis=1
)  # total alternative communities per carbon source
df_alternative_comms = df_alternative_comms.loc[tac_cs.sort_values().index]

###
# get statistics
print("Max number of alternative comms: ", np.nanmax(df_alternative_comms))
print(
    "Min number of alternative comms (ignoring cases when there are no comms): ",
    np.nanmin(df_alternative_comms),
)
print("Average number of comms: ", np.nanmean(df_alternative_comms))

###
# heatmap

# Start with a square Figure.
fig = plt.figure(figsize=(10, 6))
gs = fig.add_gridspec(
    2,
    2,
    width_ratios=((100 - 15) / 10, 15 / 10),
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
# cbar_ax = fig.add_axes([.91, .3, .03, .4])

ax = sns.heatmap(
    df_alternative_comms,
    ax=ax,
    linewidth=0.5,
    # vmin=-1, vmax=1, center=0,
    # cmap=sns.diverging_palette(20, 220, n=200),
    cmap=sns.cubehelix_palette(start=0.5, rot=-0.75, as_cmap=True),
    # cmap=sns.dark_palette("#86dad7", as_cmap=True),
    # square=True,
    # xticklabels=False,
    # cbar=False,
    cbar_ax=fig.add_subplot(gs[0, 1]),
    cbar_kws={"label": "Communities"},
)


ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histx.plot(tac_species.sort_values(), "o", color="grey")
ax_histx.set_yticks([0, 100, 200])

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax_histx.spines[s].set_visible(False)


ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
ax_histy.plot(tac_cs.sort_values().values, range(len(tac_cs)), "o", color="grey")
ax_histy.set_xticks([0, 100, 200])

# no labels
ax.tick_params(axis="x", labelbottom=False)
ax_histx.tick_params(axis="x", labelbottom=False)
ax_histy.tick_params(axis="y", labelleft=False)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax_histy.spines[s].set_visible(False)

ax.set_xlabel("Strains", size=12, family="Arial")
ax.set_ylabel("Carbon sources", size=12, family="Arial")

ax_histx.set_ylabel("Communities", size=12, family="Arial")
ax_histy.set_xlabel("Communities", size=12, family="Arial")

fig.savefig(
    path_results_script + "alternativeComms_heatmap_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# plot alternative comms vs species (each dot being a cs on which species grows in communities)
fig_all, ax_all = plt.subplots()

w = 0
x_pos = 0


# for si,s in enumerate(species_ordered_rn):
for si, s in enumerate(tac_species_order):
    if s in d_alternative_comms:
        y_pos = [d_alternative_comms[s][cs] for cs in d_alternative_comms[s]]
        y_pos_mean = np.mean(y_pos)

        # confidence interval
        ci_95 = st.norm.interval(alpha=0.95, loc=y_pos_mean, scale=st.sem(y_pos))
        ci_95min = y_pos_mean - ci_95[0]
        ci_95max = ci_95[1] - y_pos_mean
        CI_optF = np.zeros((2, 1))
        CI_optF[0, 0] = ci_95min
        CI_optF[1, 0] = ci_95max

        # ax_all.errorbar(x_pos+10, y_pos_mean, CI_optF, fmt='o', color='black',
        #                    elinewidth=3, capsize=0)
        ax_all.scatter([x_pos] * len(y_pos), y_pos, alpha=0.4, s=10, color="black")

        # total communities
        # ax_all.scatter(x_pos, sum(y_pos))

        x_pos += 25

plt.xlabel("Species", size=12, family="Arial")
plt.ylabel("Alternative communities", size=12, family="Arial")

plt.xticks([], [])

plt.show()
fig_all.savefig(
    path_results_script + "alternativeComms_vs_strains" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# how many species grow with few communities in one carbon source and with many in others?
species_with_variability_in_alt_comms = []
species_with_max10_alt_comms = []

for si, s in enumerate(tac_species_order):
    if s in d_alternative_comms:
        y_pos = [d_alternative_comms[s][cs] for cs in d_alternative_comms[s]]
        if len(y_pos) > 0:
            if min(y_pos) < 5 and max(y_pos) > 10:
                species_with_variability_in_alt_comms.append(s)
            elif max(y_pos) <= 10:
                species_with_max10_alt_comms.append(s)

print(
    "Fraction of spiecies with variability in the alternative communities: ",
    len(species_with_variability_in_alt_comms) / len(tac_species_order),
)
print(
    "Fraction of spiecies with max 10 alternative communities on every cs: ",
    len(species_with_max10_alt_comms) / len(tac_species_order),
)


###
# plot mean alternative comms vs realized niche
for si, s in enumerate(species_ordered_rn):
    if s in d_alternative_comms:
        y_pos = [d_alternative_comms[s][cs] for cs in d_alternative_comms[s]]
        y_pos_mean = np.mean(y_pos)

        # confidence interval
        ci_95 = st.norm.interval(alpha=0.95, loc=y_pos_mean, scale=st.sem(y_pos))
        ci_95min = y_pos_mean - ci_95[0]
        ci_95max = ci_95[1] - y_pos_mean
        CI_optF = np.zeros((2, 1))
        CI_optF[0, 0] = ci_95min
        CI_optF[1, 0] = ci_95max

        plt.errorbar(
            100 * ordered_rn[si] / len(cs_sorted_fba),
            y_pos_mean,
            CI_optF,
            ecolor="k",
            fmt="ko",
            markeredgecolor="white",
        )

plt.xlabel("Realized niche size", size=12, family="Arial")
plt.ylabel("Alternative communities", size=12, family="Arial")

plt.xlim((-5, 105))
# plt.xticks([], [])

plt.show()


###
# plot MEAN alternative comms vs fundamental niche
fig = plt.figure()
x_vals = []
y_vals = []
for si, s in enumerate(species_ordered_rn):
    if s in d_alternative_comms:
        y_pos = [d_alternative_comms[s][cs] for cs in d_alternative_comms[s]]
        y_pos_mean = np.mean(y_pos)

        y_vals.append(y_pos_mean)
        x_vals.append(100 * fn_ordered_rn[si] / len(cs_sorted_fba))

        # confidence interval
        ci_95 = st.norm.interval(alpha=0.95, loc=y_pos_mean, scale=st.sem(y_pos))
        ci_95min = y_pos_mean - ci_95[0]
        ci_95max = ci_95[1] - y_pos_mean
        CI_optF = np.zeros((2, 1))
        CI_optF[0, 0] = ci_95min
        CI_optF[1, 0] = ci_95max

        plt.errorbar(
            100 * fn_ordered_rn[si] / len(cs_sorted_fba),
            y_pos_mean,
            CI_optF,
            ecolor="k",
            fmt="ko",
            markeredgecolor="white",
        )

slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
if p_value < 0.05:
    if slope < 0:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#cc0b0eff"
        )
    if slope > 0:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#348bd4ff"
        )
else:
    plt.plot(
        range(100),
        [v * slope + intercept for v in range(100)],
        color="#987d82ff",
        linestyle=":",
    )

plt.xlabel("Fundamental niche size", size=12, family="Arial")
plt.ylabel("Alternative communities", size=12, family="Arial")

plt.xlim((-5, 105))
# plt.xticks([], [])

plt.show()
fig.savefig(
    path_results_script + "alternativeCommsMean_vs_fn_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# plot alternative comms vs fundamental niche

fig, ax = plt.subplots()
x_vals = []
y_vals = []
for si, s in enumerate(species_ordered_rn):
    if ordered_rn[si] != 0:  # keep only species with realized niche different from 0
        if s in d_alternative_comms:
            y_pos = [d_alternative_comms[s][cs] for cs in d_alternative_comms[s]]
            y_vals = y_vals + y_pos
            x_vals = x_vals + [100 * fn_ordered_rn[si] / len(cs_sorted_fba)] * len(
                y_pos
            )

plt.plot(x_vals, y_vals, "ko", alpha=0.4)

slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
if p_value < 0.05:
    if slope < 0:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#cc0b0eff"
        )
    if slope > 0:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#348bd4ff"
        )
else:
    plt.plot(range(100), [v * slope + intercept for v in range(100)], color="#987d82ff")

plt.xlabel("Fundamental niche size", size=12, family="Arial")
plt.ylabel("Alternative communities per carbon source", size=12, family="Arial")

plt.xlim((-5, 105))

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script + "alternativeComms_vs_fn_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# plot total alternative comms vs fundamental niche

fig, ax = plt.subplots()
x_vals = []
y_vals = []
for si, s in enumerate(species_ordered_rn):
    if ordered_rn[si] != 0:  # keep only species with realized niche different from 0
        if s in d_alternative_comms:
            alt_c = [d_alternative_comms[s][cs] for cs in d_alternative_comms[s]]
            y_vals.append(sum(alt_c))
            x_vals.append(100 * fn_ordered_rn[si] / len(cs_sorted_fba))

plt.plot(x_vals, y_vals, "ko", alpha=0.4)

slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
if p_value < 0.05:
    if slope < 0:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#cc0b0eff"
        )
    if slope > 0:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#348bd4ff"
        )
else:
    plt.plot(range(100), [v * slope + intercept for v in range(100)], color="#987d82ff")

plt.xlabel("Fundamental niche size", size=12, family="Arial")
plt.ylabel("Tolal alternative communities", size=12, family="Arial")

plt.xlim((-5, 105))

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script + "alternativeCommsTotal_vs_fn_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)

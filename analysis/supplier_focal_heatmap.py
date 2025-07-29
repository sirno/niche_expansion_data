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

path_results_script = "out/analysis/"
path_results_misosoup = "out/misosoup/"

###
# carbon sources
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

# d_community_members[cs][species], for example ['ac']['A1R12']

###
# load data on species niche size
df_niche = pd.read_csv("out/analysis/Species_fundamental_realized_niche_size.csv")

###
# order species acording to niche size
strains_ordered_fn = df_niche.sort_values("Fundamental niche size (%)")["Unnamed: 0"]
strains_ordered_rn = df_niche.sort_values("Realized niche size (%)")["Unnamed: 0"]


###
#
d_focal_supplier = {}
for cs in d_community_members:
    if cs in cs_sorted_fba_mc:
        # minimal communities (complete minimal supplying communities)
        for comm in d_community_members[cs]:
            # every strain is focal and supplier
            sp1 = comm[0].split("y_")[1]
            sp2 = comm[1].split("y_")[1]
            # sp1 is the focal species
            if sp1 not in d_focal_supplier:
                d_focal_supplier[sp1] = {}
                d_focal_supplier[sp1][sp2] = 1
            else:
                if sp2 not in d_focal_supplier[sp1]:
                    d_focal_supplier[sp1][sp2] = 1
                else:
                    d_focal_supplier[sp1][sp2] += 1

            # sp2 is the focal species
            if sp2 not in d_focal_supplier:
                d_focal_supplier[sp2] = {}
                d_focal_supplier[sp2][sp1] = 1
            else:
                if sp1 not in d_focal_supplier[sp2]:
                    d_focal_supplier[sp2][sp1] = 1
                else:
                    d_focal_supplier[sp2][sp1] += 1
    else:
        # minimal supplying communities
        for fs in d_community_members[cs]:
            if fs not in d_focal_supplier:
                d_focal_supplier[fs] = {}
            # loop through communities to get suppliers
            for comm in d_community_members[cs][fs]:
                for y_sp in comm:
                    sp = y_sp.split("y_")[1]
                    if sp != fs:
                        if sp not in d_focal_supplier[fs]:
                            d_focal_supplier[fs][sp] = 1
                        else:
                            d_focal_supplier[fs][sp] += 1


###
# transform to pandas db - focal species on the columns
df_focal_supp = pd.DataFrame.from_dict(d_focal_supplier)
df_focal_supp.to_csv("out/analysis/Focal_Supplyier.csv", index=True)


###
# heatmap
df_focal_supp_ordered_fn = df_focal_supp.reindex(columns=strains_ordered_fn)

fig, ax = plt.subplots()

ax = sns.heatmap(
    df_focal_supp_ordered_fn,
    ax=ax,
    linewidth=0.5,
    # vmin=-1, vmax=1, center=0,
    # cmap=sns.diverging_palette(20, 220, n=200),
    cmap=sns.cubehelix_palette(start=0.5, rot=-0.75, as_cmap=True),
    # cmap=sns.dark_palette("#86dad7", as_cmap=True),
    # square=True,
    xticklabels=False,
    yticklabels=False,
    # cbar=False,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    cbar_kws={"label": "Times strain acts as supplyier"},
)

ax.set_xlabel("Focal strain", size=12, family="Arial")
ax.set_ylabel("Supplier strain", size=12, family="Arial")


plt.show()
fig.savefig(
    path_results_script + "Focal_supplier_heatmap_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


fn_size = df_niche["Fundamental niche size (%)"]
rn_size = df_niche["Fundamental niche size (%)"]
niche_expansion_degree = [
    df_niche["Realized niche size (%)"].iloc[i]
    - df_niche["Fundamental niche size (%)"].iloc[i]
    for i in range(len(df_niche["Fundamental niche size (%)"]))
]

n_suppliers_all = []
for fs in df_niche["Unnamed: 0"]:
    n_suppliers = len(df_focal_supp[fs]) - sum(df_focal_supp[fs].isnull())
    n_suppliers_all.append(n_suppliers)


###
# plot num suppliers/degree of niche expansion vs fundamental niche size

x_vals = []
y_vals = []
for ix in range(len(niche_expansion_degree)):
    if niche_expansion_degree[ix] != 0:
        coef = n_suppliers_all[ix] / niche_expansion_degree[ix]
        y_vals.append(coef)
        x_vals.append(fn_size.values[ix])

        if coef > 4:
            print(n_suppliers_all[ix], niche_expansion_degree[ix])

        if n_suppliers_all[ix] < 5 and niche_expansion_degree[ix] > 10:
            print(n_suppliers_all[ix], niche_expansion_degree[ix])


fig, ax = plt.subplots()
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
plt.ylabel("Alternative suppliers/niche expansion degree", size=12, family="Arial")

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script
    + "alternativeSuppliersOverNicheExpansionDegree_vs_fn_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


###
# plot num suppliers vs fundamental niche size

x_vals = fn_size.values
y_vals = n_suppliers_all[:]


fig, ax = plt.subplots()
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
plt.ylabel("Alternative suppliers", size=12, family="Arial")

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script + "alternativeSuppliers_vs_fn_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


# remove strains with realized niche 0
x_vals = [fn_size[ix] for ix in range(len(rn_size)) if rn_size[ix] != 0]
y_vals = [n_suppliers_all[ix] for ix in range(len(rn_size)) if rn_size[ix] != 0]


fig, ax = plt.subplots()
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
plt.ylabel("Alternative suppliers", size=12, family="Arial")

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script
    + "alternativeSuppliers_vs_fn_exclude0_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


###
# plot realized niche vs. fundamental niche
fig, ax = plt.subplots()
plt.plot(fn_size.values, niche_expansion_degree, "ko", alpha=0.4)

plt.ylabel("Niche expansion degree", size=12, family="Arial")
plt.xlabel("Fundamental niche size", size=12, family="Arial")

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script + "niche_expansion_vs_fn_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# plot num suppliers vs niche expansion degree
# remove strains with realized niche 0
x_vals = niche_expansion_degree[:]
y_vals = n_suppliers_all[:]


fig, ax = plt.subplots()
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

plt.xlabel("Niche expansion degree", size=12, family="Arial")
plt.ylabel("Alternative suppliers", size=12, family="Arial")

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script
    + "alternativeSuppliers_vs_niche_expansion_degree_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)


# remove strains with realized niche 0
x_vals = [niche_expansion_degree[ix] for ix in range(len(rn_size)) if rn_size[ix] != 0]
y_vals = [n_suppliers_all[ix] for ix in range(len(rn_size)) if rn_size[ix] != 0]


fig, ax = plt.subplots()
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

plt.xlabel("Niche expansion degree", size=12, family="Arial")
plt.ylabel("Alternative suppliers", size=12, family="Arial")

# Setting the number of ticks
plt.locator_params(axis="y", nbins=4)
plt.locator_params(axis="x", nbins=4)

# Remove axes splines
for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)


plt.show()
fig.savefig(
    path_results_script
    + "alternativeSuppliers_vs_niche_expansion_degree_exclude0_"
    + str(tolerance)
    + ".pdf",
    bbox_inches="tight",
)

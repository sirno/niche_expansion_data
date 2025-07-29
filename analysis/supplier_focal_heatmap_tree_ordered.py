#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 09:59:05 2025

@author: magdalena
"""

import json

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import yaml

###
# load data

# focal-supplyier matrix
df_focal_supp = pd.read_csv("out/analysis/Focal_Supplyier.csv")
df_focal_supp.index = df_focal_supp["Unnamed: 0"]

# strains
with open("../data/Strains_used", "r") as fp:
    all_strains = json.load(fp)


# strains ordered as in tree
strains_tree_ordered = [
    "A1R12",
    "E3R18",
    "E3R01",
    "B2M17",
    "I2M19",
    "C3R15",
    "E3R09",
    "D3R19",
    "B2M06",
    "E3M18",
    "6B07",
    "C3R19",
    "A2M03",
    "C3R17",
    "B3R18",
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
    "B3M08",
    "B2R22",
    "D2R18",
    "D2R04",
    "F2R14",
    "G3M19",
    "3C02",
    "G2R10",
    "1A01",
    "I3M07",
    "C3R12",
    "6D03",
    "E3R11",
    "6C06",
    "B3M02",
    "6D02",
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
    "D2M19",
    "F3R11",
    "F3R08",
    "B2M13",
    "3B05",
    "C1R06",
    "G2M07",
    "E3M09",
    "E3M17",
]  # 64

# add m_ to strain name when needed so that it matches the names used in df_focal_supp
strains_m_tree_ordered = []
for s in strains_tree_ordered:
    if s[0].isnumeric():
        strains_m_tree_ordered.append("m_" + s)
    else:
        strains_m_tree_ordered.append(s)


strains_tree_ordered_not_in_data = [
    s for s in strains_m_tree_ordered if s not in all_strains
]
# ['m_6B07', 'C3R19', 'B3M08', 'B2R22']

###
#
df_ordered_cols = df_focal_supp.reindex(columns=strains_m_tree_ordered)
df_ordered_rows = df_ordered_cols.reindex(strains_m_tree_ordered)


###
# heatmap
fig, ax = plt.subplots(figsize=(20, 15))

ax = sns.heatmap(
    df_ordered_rows,
    ax=ax,
    linewidth=0.5,
    # vmin=-1, vmax=1, center=0,
    # cmap=sns.diverging_palette(20, 220, n=200),
    cmap=sns.cubehelix_palette(start=0.5, rot=-0.75, as_cmap=True),
    # cmap=sns.dark_palette("#86dad7", as_cmap=True),
    # square=True,
    xticklabels=1,
    yticklabels=1,
    # cbar=False,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    cbar_kws={"label": "Times strain acts as supplyier"},
)

ax.set_xlabel("Focal strain", size=12, family="Arial")
ax.set_ylabel("Supplier strain", size=12, family="Arial")


plt.show()
fig.savefig("out/analysis/Focal_supplier_heatmap_tree_order.pdf", bbox_inches="tight")


###
# count times family supplies family
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

# add m_ to strain name if needed
for fam in d_tree:
    for ist, st in enumerate(d_tree[fam]):
        if st[0].isnumeric():
            m_st = "m_" + st
            d_tree[fam][ist] = m_st

fil = open("out/analysis/Focal_supplier_heatmap_count.csv", "w")
fil.write("Supplier family,Focal family,Count,Count/strains in supp family\n")

for focal_family in d_tree:
    for sup_family in d_tree:
        sum_alt_sup = (
            df_ordered_rows[d_tree[focal_family]].loc[d_tree[sup_family]].sum().sum()
        )
        fil.write(
            sup_family
            + ","
            + focal_family
            + ","
            + str(sum_alt_sup)
            + ","
            + str(sum_alt_sup / len(d_tree[sup_family]))
            + "\n"
        )

fil.close()


###
# heatmap sum suppliers
count_times_supplier = df_ordered_rows.sum(axis=1)
count_times_supplier = np.asarray(count_times_supplier)[:, np.newaxis]

fig, ax = plt.subplots(figsize=(20, 15))

ax = sns.heatmap(
    count_times_supplier,
    ax=ax,
    linewidth=0.5,
    # vmin=-1, vmax=1, center=0,
    # cmap=sns.diverging_palette(20, 220, n=200),
    cmap=sns.cubehelix_palette(start=0.5, rot=-0.75, as_cmap=True),
    # cmap=sns.dark_palette("#86dad7", as_cmap=True),
    # square=True,
    # xticklabels=1,
    yticklabels=df_ordered_rows.index,
    # cbar=False,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    cbar_kws={"label": "Times strain acts as supplyier"},
)

# ax.set_xlabel('Focal strain', size=12, family='Arial')
ax.set_ylabel("Supplier strain", size=12, family="Arial")


plt.show()
fig.savefig("out/analysis/Supplier_heatmap_tree_order.pdf", bbox_inches="tight")


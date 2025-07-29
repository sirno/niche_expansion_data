#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:19:30 2025

@author: magdalena
"""

import json

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import yaml
from matplotlib.colors import LinearSegmentedColormap

###
# load data

# fundamental niche
df_fn = pd.read_csv("out/analysis/Fundamental_niche.csv")
df_fn.index = df_fn["Unnamed: 0"]
df_fn = df_fn.drop("Unnamed: 0", axis=1)

# realized niche with 2 strain communities
df_rn_2st_comms = pd.read_csv("out/analysis/Realized_niche_2st_comms.csv")
df_rn_2st_comms.index = df_rn_2st_comms["Unnamed: 0"]
df_rn_2st_comms = df_rn_2st_comms.drop("Unnamed: 0", axis=1)
df_rn_2st_comms = df_rn_2st_comms.replace(1, 4)

# realized niche
df_rn_3st_comms = pd.read_csv("out/analysis/Realized_niche_3st_comms.csv")
df_rn_3st_comms.index = df_rn_3st_comms["Unnamed: 0"]
df_rn_3st_comms = df_rn_3st_comms.drop("Unnamed: 0", axis=1)
df_rn_3st_comms = df_rn_3st_comms.replace(1, 6)


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


###
# fundamental niche

# order dataframe rows to match strains in tree
df_ordered_rows = df_fn.reindex(strains_tree_ordered)

###
# heatmap
fig, ax = plt.subplots(figsize=(5, 15))

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
    cbar=False,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    # cbar_kws={'label': 'Times strain acts as supplyier'}
)

ax.set_ylabel("Strains", size=12, family="Arial")
ax.set_xlabel("Carbon sources", size=12, family="Arial")

ax.set_title("Fundamental niche", size=12, family="Arial")

plt.show()
fig.savefig(
    "out/analysis/Fundamental_niche_heatmap_tree_order.pdf", bbox_inches="tight"
)


###
# realized niche with 2 strain communities

# order dataframe rows to match strains in tree
df_ordered_rows = df_rn_2st_comms.reindex(strains_tree_ordered)

###
# heatmap
fig, ax = plt.subplots(figsize=(5, 15))

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
    cbar=False,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    # cbar_kws={'label': 'Times strain acts as supplyier'}
)

ax.set_ylabel("Strains", size=12, family="Arial")
ax.set_xlabel("Carbon sources", size=12, family="Arial")

ax.set_title("Realized niche", size=12, family="Arial")

plt.show()
fig.savefig(
    "out/analysis/Realized_niche_2str_heatmap_tree_order.pdf", bbox_inches="tight"
)


###
# realized niche with 3 strain communities

# order dataframe rows to match strains in tree
df_ordered_rows = df_rn_3st_comms.reindex(strains_tree_ordered)

###
# heatmap
fig, ax = plt.subplots(figsize=(5, 15))

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
    cbar=False,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    # cbar_kws={'label': 'Times strain acts as supplyier'}
)


ax.set_ylabel("Strains", size=12, family="Arial")
ax.set_xlabel("Carbon sources", size=12, family="Arial")

ax.set_title("Realized niche", size=12, family="Arial")

plt.show()
fig.savefig(
    "out/analysis/Realized_niche_3str_heatmap_tree_order.pdf", bbox_inches="tight"
)


###
# heatmap all
df_all = df_fn + df_rn_2st_comms + df_rn_3st_comms

df_ordered_rows = df_all.reindex(strains_tree_ordered)
df_ordered_rows = df_ordered_rows.replace(1, 3)

fig, ax = plt.subplots(figsize=(5, 15))

myColors = ("#f9f5e7", "#000000B3", "#b9d7b1", "#eca04fff", "#d4c58c")
cmap = LinearSegmentedColormap.from_list("Custom", myColors, len(myColors))


ax = sns.heatmap(
    df_ordered_rows,
    ax=ax,
    linewidth=0.5,
    # vmin=-1, vmax=1, center=0,
    # cmap=sns.diverging_palette(20, 220, n=200),
    # cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True),
    # cmap=sns.dark_palette("#86dad7", as_cmap=True),
    cmap=cmap,
    # square=True,
    xticklabels=1,
    yticklabels=1,
    cbar=True,
    # cbar_ax=fig.add_subplot(gs[0, 1]),
    # cbar_kws={'label': 'Times strain acts as supplyier'}
)

vmap = {
    0: "No growth",
    3.0: "In isolation",
    4: "In 2-strain communities",
    6: "In 3-strain communities",
    10: "In 2 and 3-strain communities",
}
n = len(vmap)

# Get the colorbar object from the Seaborn heatmap
colorbar = ax.collections[0].colorbar
# The list comprehension calculates the positions to place the labels to be evenly distributed across the colorbar
r = colorbar.vmax - colorbar.vmin
colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n) + r * i / (n) for i in range(n)])
colorbar.set_ticklabels(list(vmap.values()))

ax.set_ylabel("Strains", size=12, family="Arial")
ax.set_xlabel("Carbon sources", size=12, family="Arial")

ax.set_title("Realized niche", size=12, family="Arial")

plt.show()
fig.savefig(
    "out/analysis/Fundamental_and_realized_niche_heatmap_tree_order.pdf",
    bbox_inches="tight",
)

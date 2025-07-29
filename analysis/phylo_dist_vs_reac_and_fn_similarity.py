#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:28:54 2025

@author: magdalena
"""

import json

import matplotlib.pyplot as plt
import pandas as pd
import yaml

###
# load data

# strains
with open("../data/Strains_used", "r") as fp:
    all_strains = json.load(fp)

# phylogenetic distance
df_phylo_dis = pd.read_csv("../data/phylogenetic_distance_matrix.csv")

# fundamental niche
df_fn = pd.read_csv("out/analysis/Fundamental_niche.csv")

# genes/reactions
with open("out/analysis/d_species_info_extended.yaml", "r") as f:
    d_species_info_extended = yaml.load(f, Loader=yaml.SafeLoader)

###
# fundamental niche similarity vs phylo distance
jaccard_similarity_fn = []
jaccard_similarity_reacs = []
phylo_distance = []
for is1 in range(len(all_strains) - 1):
    for is2 in range(is1 + 1, len(all_strains)):
        # pair of strains
        s1 = all_strains[is1]
        s2 = all_strains[is2]

        # change strain name if it starts with m_
        if s1.startswith("m_"):
            s1 = s1.split("m_")[1]

        if s2.startswith("m_"):
            s2 = s2.split("m_")[1]
        ###
        # fn
        fn_s1 = df_fn[df_fn["Unnamed: 0"] == s1].iloc[0][1:]
        fn_s2 = df_fn[df_fn["Unnamed: 0"] == s2].iloc[0][1:]

        # fn-intersection
        fn_intersection = len(
            [c for c in fn_s1.index if fn_s1[c] == 1 and fn_s2[c] == 1]
        )

        # fn-union
        fn_union = len([c for c in fn_s1.index if fn_s1[c] == 1 or fn_s2[c] == 1])

        # jaccard similarity
        if fn_union != 0:
            jaccard_similarity_fn.append(fn_intersection / fn_union)
        else:
            jaccard_similarity_fn.append(0)

        ###
        # reactions

        reac_s1 = d_species_info_extended[s1]["reactions"]
        reac_s2 = d_species_info_extended[s2]["reactions"]

        # reacs-intersection
        reacs_intersection = len([r for r in reac_s1 if r in reac_s2])

        # reacs-union
        reacs_union = len(reac_s1 + [r for r in reac_s2 if r not in reac_s1])

        # reacs jaccard similarity
        if reacs_union != 0:
            jc_reac = reacs_intersection / reacs_union
            jaccard_similarity_reacs.append(jc_reac)
            # if strains are identical in reacs show pair
            if jc_reac == 1:
                print(
                    "Identical pair "
                    + s1
                    + "-"
                    + s2
                    + " with "
                    + str(len(reac_s1))
                    + " reactions"
                )
        else:
            jaccard_similarity_reacs.append(0)

        ###
        # phylo distance
        phylo_d_pair = df_phylo_dis[df_phylo_dis["Unnamed: 0"] == s1][s2].iloc[0]
        phylo_distance.append(phylo_d_pair)


###
# plot reaction similarity vs phylo distance
fig, ax = plt.subplots()

plt.plot(
    phylo_distance, jaccard_similarity_reacs, "ko", alpha=0.3
)  # ,markeredgecolor ='white')

plt.ylabel("Reactions similarity (jaccard index)", size=12, family="Arial")
plt.xlabel("Phylogenetic distance", size=12, family="Arial")

# Remove axes splines
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

# Add padding between axes and labels
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=10)

# Setting the number of ticks
plt.locator_params(axis="x", nbins=4)
plt.locator_params(axis="y", nbins=4)

plt.show()
fig.savefig("out/analysis/reacs_similarity_vs_phylo_distance.pdf", bbox_inches="tight")


###
# plot fn similarity vs phylo distance
fig, ax = plt.subplots()

plt.plot(phylo_distance, jaccard_similarity_fn, "ko", alpha=0.3)

plt.ylabel("Fundamental niche similarity (jaccard index)", size=12, family="Arial")
plt.xlabel("Phylogenetic distance", size=12, family="Arial")

# Remove axes splines
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

# Add padding between axes and labels
ax.xaxis.set_tick_params(pad=5)
ax.yaxis.set_tick_params(pad=10)

# Setting the number of ticks
plt.locator_params(axis="x", nbins=4)
plt.locator_params(axis="y", nbins=4)

plt.show()
fig.savefig("out/analysis/fn_similarity_vs_phylo_distance.pdf", bbox_inches="tight")

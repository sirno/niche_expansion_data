#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:27:20 2025

@author: magdalena
"""

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

###
path_results_script = "out/analysis/"

###
# load data
# species niche size
df_niche = pd.read_csv("out/analysis/Species_fundamental_realized_niche_size.csv")

# alternative communities
df_altcomm = pd.read_csv("out/analysis/Alternatice_communities_species_cs.csv")

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
# species
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

###
# plot alternative communities vs species fundamental niche size separately for each carbon source

for cs in df_altcomm["Unnamed: 0"].values:
    data_on_cs = df_altcomm[df_altcomm["Unnamed: 0"] == cs]
    alt_comm_here = []
    fn_here = []
    for sp in [v for v in df_altcomm.columns if v != "Unnamed: 0"]:
        if not data_on_cs[sp].isnull().values[0]:
            alt_comm_here.append(data_on_cs[sp].iloc[0])
            fn_here.append(
                df_niche[df_niche["Unnamed: 0"] == sp][
                    "Fundamental niche size (%)"
                ].iloc[0]
            )

    fig, ax = plt.subplots()
    x_vals = fn_here
    y_vals = alt_comm_here
    plt.plot(x_vals, y_vals, "ko", alpha=0.4)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
    if p_value < 0.05:
        if slope < 0:
            plt.plot(
                range(100),
                [v * slope + intercept for v in range(100)],
                color="#cc0b0eff",
            )
        if slope > 0:
            plt.plot(
                range(100),
                [v * slope + intercept for v in range(100)],
                color="#348bd4ff",
            )
    else:
        plt.plot(
            range(100), [v * slope + intercept for v in range(100)], color="#987d82ff"
        )

    plt.xlabel("Fundamental niche size (%)", size=12, family="Arial")
    plt.ylabel("Alternative communities on " + cs, size=12, family="Arial")

    # Setting the number of ticks
    plt.locator_params(axis="y", nbins=4)
    plt.locator_params(axis="x", nbins=4)

    # Remove axes splines
    for s in ["top", "right"]:  # ,'top', 'bottom', 'left', 'right']:
        ax.spines[s].set_visible(False)

    plt.show()
    fig.savefig(
        path_results_script + "alternativeComs_vs_fn_on_" + cs + ".pdf",
        bbox_inches="tight",
    )

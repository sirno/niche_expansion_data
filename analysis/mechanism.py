#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 09:27:07 2024

@author: magdalena
"""

import yaml

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
    # 'cys__L',
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
# load data needed
# for tolerance in [7,8,9]:
tolerance = 7
print("Tolerance: ", tolerance)

with open(
    path_results_misosoup + "tolerance_1e-" + str(tolerance) + ".yaml", "r"
) as file:
    d_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

with open(
    path_results_script + "d_misosoup_growth_tolerance_1e-" + str(tolerance) + ".yaml",
    "r",
) as file:
    d_growth_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

with open("data/medium.yaml", "r") as file:
    d_medium = yaml.load(file, Loader=yaml.SafeLoader)


###
# functions
def focal_secretions_fun(s, comm, threshold_flux_not_zero):
    focal_secretions = []

    for reac in comm["solution"]:
        # get metabs secreted by focal
        if reac.endswith(s + "_i") and comm["solution"][reac] > threshold_flux_not_zero:
            # check that metab doesnt come in the medium
            general_reac = reac.split("_e")[0] + "_e"
            if general_reac not in d_medium["base_medium"]:
                separate = reac.split("EX_")
                separate2 = separate[1].split("_e")
                metab_secreted = separate2[0]
                focal_secretions.append(metab_secreted)

    return focal_secretions


def focal_secretions_crossfed_fun(focal_secretions, comm, threshold_flux_not_zero):
    # does anyone consume these secretions?
    focal_secretions_crossfed = []

    for sm in focal_secretions:
        reac_sm = "R_EX_" + sm + "_e"
        for reac in comm["solution"]:
            if (reac_sm in reac) and (
                comm["solution"][reac] < -1 * threshold_flux_not_zero
            ):
                focal_secretions_crossfed.append(sm)

    return focal_secretions_crossfed


#######
# get community motifs (mutualism/commensalism and who consumes the carbon source)
#######
# we classified the community flux distribution found for each minimal supplying community according to six categories:
# (i, ii) the carbon source available in the environment is consumed by the focal species,
# (iii, iv) by a supplying species
# or (v, vi) both;
# (ii, iv, vi) and where the focal species secretes something that at least one of the supplying species consumes,
# (i, iii, v) or not.
#######
d_motifs = {}
d_motifs[2] = {}  # community size 2
d_motifs[2]["i"] = []
d_motifs[2]["ii"] = []
d_motifs[2]["iii"] = []
d_motifs[2]["iv"] = []
d_motifs[2]["v"] = []
d_motifs[2]["vi"] = []
d_motifs[3] = {}  # community size 3
d_motifs[3]["i"] = []
d_motifs[3]["ii"] = []
d_motifs[3]["iii"] = []
d_motifs[3]["iv"] = []
d_motifs[3]["v"] = []
d_motifs[3]["vi"] = []

d_motifs_count = {}
d_motifs_count[2] = {}  # community size
d_motifs_count[2]["i"] = 0
d_motifs_count[2]["ii"] = 0
d_motifs_count[2]["iii"] = 0
d_motifs_count[2]["iv"] = 0
d_motifs_count[2]["v"] = 0
d_motifs_count[2]["vi"] = 0
d_motifs_count[3] = {}  # community size
d_motifs_count[3]["i"] = 0
d_motifs_count[3]["ii"] = 0
d_motifs_count[3]["iii"] = 0
d_motifs_count[3]["iv"] = 0
d_motifs_count[3]["v"] = 0
d_motifs_count[3]["vi"] = 0

threshold_flux_not_zero = 0.001

sum_communities = 0

for cs in cs_sorted_fba_msc:
    for si, s in enumerate(all_species):
        # in community?
        in_comm = 0
        if d_growth_misosoup[cs][s] == "in community":
            in_comm = 1

        if in_comm == 1:
            # analyze each comm
            for icom, comm in enumerate(d_misosoup[cs][s]):
                sum_communities += 1
                does_focal_consume_cs_value = 0
                does_suppling_consume_cs_value = 0

                # get comm size
                comm_size = len(comm["community"])

                # does the focal strain consume the carbon source?
                uptake_reac_id_focal = "R_EX_" + cs + "_e_" + s + "_i"
                if uptake_reac_id_focal in comm["solution"]:
                    if (
                        comm["solution"][uptake_reac_id_focal]
                        < -1 * threshold_flux_not_zero
                    ):
                        does_focal_consume_cs_value = 1

                # does a suppling strain consume the carbon source?
                does_suppling_consume_cs_value = 0
                suppliers = [sup.split("y_")[1] for sup in comm["community"]]
                suppliers = [sup for sup in suppliers if sup != s]

                for sup in suppliers:
                    uptake_reac_id_sup = "R_EX_" + cs + "_e_" + sup + "_i"
                    if (uptake_reac_id_sup in comm["solution"]) and (
                        comm["solution"][uptake_reac_id_sup]
                        < -1 * threshold_flux_not_zero
                    ):
                        does_suppling_consume_cs_value = 1

                # mutualistic interaction?
                focal_secretions = focal_secretions_fun(
                    s, comm, threshold_flux_not_zero
                )
                focal_secretions_crossfed = focal_secretions_crossfed_fun(
                    focal_secretions, comm, threshold_flux_not_zero
                )
                if len(focal_secretions_crossfed) != 0:
                    mutualism = 1

                # assign motif
                if (
                    does_focal_consume_cs_value == 1
                    and does_suppling_consume_cs_value == 0
                    and mutualism == 0
                ):
                    d_motifs_count[comm_size]["i"] += 1
                    d_motifs[comm_size]["i"].append([s, cs, icom])
                elif (
                    does_focal_consume_cs_value == 1
                    and does_suppling_consume_cs_value == 0
                    and mutualism == 1
                ):
                    d_motifs_count[comm_size]["ii"] += 1
                    d_motifs[comm_size]["ii"].append([s, cs, icom])
                elif (
                    does_focal_consume_cs_value == 0
                    and does_suppling_consume_cs_value == 1
                    and mutualism == 0
                ):
                    d_motifs_count[comm_size]["iii"] += 1
                    d_motifs[comm_size]["iii"].append([s, cs, icom])
                elif (
                    does_focal_consume_cs_value == 0
                    and does_suppling_consume_cs_value == 1
                    and mutualism == 1
                ):
                    d_motifs_count[comm_size]["iv"] += 1
                    d_motifs[comm_size]["iv"].append([s, cs, icom])
                elif (
                    does_focal_consume_cs_value == 1
                    and does_suppling_consume_cs_value == 1
                    and mutualism == 0
                ):
                    d_motifs_count[comm_size]["v"] += 1
                    d_motifs[comm_size]["v"].append([s, cs, icom])
                elif (
                    does_focal_consume_cs_value == 1
                    and does_suppling_consume_cs_value == 1
                    and mutualism == 1
                ):
                    d_motifs_count[comm_size]["vi"] += 1
                    d_motifs[comm_size]["vi"].append([s, cs, icom])
                else:
                    print("Didnt fit any chategory")


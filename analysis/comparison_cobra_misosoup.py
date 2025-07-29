#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:50:59 2021

@author: magdalena
"""

import yaml

# carbon sources sorted according to increasing number of strains that are
# viable in isolation (obtained with FBA)
carbon_sources = [
    "3pg",
    "cyst__L",
    "glu__D",
    "icit",
    "pep",
    "r5p",
    "akg",
    "met__L",
    "oaa",
    "f6p",
    "phe__L",
    "trp__L",
    "rmn",
    "cys__L",
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
]

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
# load cobra results
with open(path_results_script + "d_cobra_growth.yaml", "r") as file:
    d_cobra_growth = yaml.load(file, Loader=yaml.Loader)  # Loader=yaml.SafeLoader)

d_species = {}
for s in d_cobra_growth:
    d_species[s] = {}
    for c in d_cobra_growth[s]:
        d_species[s][c] = {}
        d_species[s][c]["cobra_growth"] = d_cobra_growth[s][c]

###
# load misosoup results
tolerance = 7
with open(
    path_results_script + "/d_misosoup_growth_tolerance_1e-" + str(tolerance) + ".yaml",
    "r",
) as file:
    d_growth_misosoup_7 = yaml.load(file, Loader=yaml.SafeLoader)

tolerance = 8
with open(
    path_results_script + "/d_misosoup_growth_tolerance_1e-" + str(tolerance) + ".yaml",
    "r",
) as file:
    d_growth_misosoup_8 = yaml.load(file, Loader=yaml.SafeLoader)

tolerance = 9
with open(
    path_results_script + "/d_misosoup_growth_tolerance_1e-" + str(tolerance) + ".yaml",
    "r",
) as file:
    d_growth_misosoup_9 = yaml.load(file, Loader=yaml.SafeLoader)

###
# show conditions under which there is discrepancy between misosoup (max growth) and cobra results
for c in d_growth_misosoup_7:
    for s in d_growth_misosoup_7[c]:
        if s.startswith("m_"):
            n = s[2:]
        else:
            n = s

        cobra_g = d_cobra_growth[n][c]
        if cobra_g != 0:
            cobra_g = cobra_g.objective_value
        miso_g = d_growth_misosoup_7[c][s]

        # cobra growth in isolation - misosoup doesnt
        if cobra_g > 0 and (miso_g == "no growth" or miso_g == "in community"):
            print(c, s, cobra_g, miso_g)

        # cobra doesn't grow in isolation - misosoup does
        if cobra_g == 0 and (miso_g != "no growth" and miso_g != "in community"):
            print(c, s, cobra_g, miso_g)

        # both return growth in isolation but different growth values
        if (
            cobra_g > 0
            and (miso_g != "no growth" and miso_g != "in community")
            and (type(miso_g) != float or cobra_g - miso_g > 0.001)
        ):
            print(c, s, cobra_g, miso_g)


for c in d_growth_misosoup_8:
    for s in d_growth_misosoup_8[c]:
        if s.startswith("m_"):
            n = s[2:]
        else:
            n = s

        cobra_g = d_cobra_growth[n][c]
        if cobra_g != 0:
            cobra_g = cobra_g.objective_value
        miso_g = d_growth_misosoup_8[c][s]

        # cobra growth in isolation - misosoup doesnt
        if cobra_g > 0 and (miso_g == "no growth" or miso_g == "in community"):
            print(c, s, cobra_g, miso_g)

        # cobra doesn't grow in isolation - misosoup does
        if cobra_g == 0 and (miso_g != "no growth" and miso_g != "in community"):
            print(c, s, cobra_g, miso_g)

        # both return growth in isolation but different growth values
        if (
            cobra_g > 0
            and (miso_g != "no growth" and miso_g != "in community")
            and (type(miso_g) != float or cobra_g - miso_g > 0.001)
        ):
            print(c, s, cobra_g, miso_g)


for c in d_growth_misosoup_9:
    for s in d_growth_misosoup_9[c]:
        if s.startswith("m_"):
            n = s[2:]
        else:
            n = s

        cobra_g = d_cobra_growth[n][c]
        if cobra_g != 0:
            cobra_g = cobra_g.objective_value
        miso_g = d_growth_misosoup_9[c][s]

        # cobra growth in isolation - misosoup doesnt
        if cobra_g > 0 and (miso_g == "no growth" or miso_g == "in community"):
            print(c, s, cobra_g, miso_g)

        # cobra doesn't grow in isolation - misosoup does
        if cobra_g == 0 and (miso_g != "no growth" and miso_g != "in community"):
            print(c, s, cobra_g, miso_g)

        # both return growth in isolation but different growth values
        if (
            cobra_g > 0
            and (miso_g != "no growth" and miso_g != "in community")
            and (type(miso_g) != float or cobra_g - miso_g > 0.001)
        ):
            print(c, s, cobra_g, miso_g)


###
# discrepancy between misosoup results with different tolerances
d_discrepancies = {}
for c in d_growth_misosoup_7:
    for s in d_growth_misosoup_7[c]:
        if (
            d_growth_misosoup_9[c][s] != d_growth_misosoup_7[c][s]
            or d_growth_misosoup_9[c][s] != d_growth_misosoup_8[c][s]
            or d_growth_misosoup_7[c][s] != d_growth_misosoup_8[c][s]
        ):
            # if d_growth_misosoup_7[c][s]!=d_growth_misosoup_8[c][s]:
            if c not in d_discrepancies:
                d_discrepancies[c] = []
            d_discrepancies[c].append(s)
            print(c, s, d_growth_misosoup_7[c][s], d_growth_misosoup_8[c][s])


###
# write comparison cobra-misosoup to table
fil = open(path_results_script + "/comparison_cobra_misosoup.tsv", "w")
for cs in d_growth_misosoup_7:
    fil.write("\t" + cs)

fil.write("\n")

for species in d_growth_misosoup_7[cs]:
    if species.startswith("m_"):
        n = species[2:]
    else:
        n = species
    #
    fil.write(species)

    for cs in d_growth_misosoup_7:
        if d_cobra_growth[n][cs] == 0:
            fil.write("\t" + str(d_cobra_growth[n][cs]))
        else:
            fil.write("\t" + str(d_cobra_growth[n][cs].objective_value))
        if type(d_growth_misosoup_7[cs][species]) == "str":
            fil.write(" / " + d_growth_misosoup_7[cs][species])
        else:
            fil.write(" / " + str(d_growth_misosoup_7[cs][species]))
        if type(d_growth_misosoup_8[cs][species]) == "str":
            fil.write(" / " + d_growth_misosoup_8[cs][species])
        else:
            fil.write(" / " + str(d_growth_misosoup_8[cs][species]))
        if type(d_growth_misosoup_9[cs][species]) == "str":
            fil.write(" / " + str(d_growth_misosoup_9[cs][species]))
        else:
            fil.write(" / " + str(d_growth_misosoup_9[cs][species]))

    fil.write("\n")

fil.close()

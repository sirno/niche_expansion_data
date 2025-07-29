#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 10:36:09 2022

@author: magdalena
"""

from os import listdir, makedirs
from os.path import isfile, join

import yaml
from cobra.io import read_sbml_model

# medium simulates MBM without co2 or hco3
medium_MBM_no_co2_hco3 = [
    "EX_h2o_e",
    "EX_o2_e",
    "EX_cl_e",
    "EX_na1_e",
    "EX_so4_e",
    "EX_k_e",
    "EX_mg2_e",
    "EX_ca2_e",
    "EX_nh4_e",
    "EX_pi_e",
    "EX_btn_e",
    "EX_fol_e",
    "EX_pydxn_e",
    "EX_ribflv_e",
    "EX_thm_e",
    "EX_nac_e",
    "EX_ala_B_e",
    "EX_4abz_e",
    "EX_fe3_e",
    "EX_h_e",
    "EX_cobalt2_e",
    "EX_cu2_e",
    "EX_fe2_e",
    "EX_mn2_e",
    "EX_mobd_e",
    "EX_zn2_e",
    "EX_sel_e",
    "EX_ni2_e",
    "EX_slnt_e",
    "EX_tungs_e",
    "EX_h2_e",
    "EX_cu_e",
    "EX_cbl1_e",
]

# carbon sources sorted according to increasing number of strains that are viable in isolation (obtained with FBA)
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

path_results_script = "out/"

# load the models
path_networks = "data/strains/"
onlynetworks = [
    f
    for f in listdir(path_networks)
    if (isfile(join(path_networks, f)) and f[-3:] == "xml")
]

d_species = {}
for n in onlynetworks:
    d_species[n[:-4]] = {}
    model = read_sbml_model(path_networks + n)
    d_species[n[:-4]]["model"] = model


# calculate growth in medium_MBM using cobrapy
d_cobra_growth = {}

for n in d_species:
    d_cobra_growth[n] = {}
    # set minimal medium
    for r in d_species[n]["model"].exchanges:
        if r.id in medium_MBM_no_co2_hco3:
            r.lower_bound = -1000
        else:
            r.lower_bound = 0
    # change carbon source and calculate growth
    for cs in carbon_sources:
        d_species[n][cs] = {}
        EX_cs = "EX_" + cs + "_e"
        if EX_cs in [r.id for r in d_species[n]["model"].exchanges]:
            d_species[n]["model"].reactions.get_by_id(EX_cs).lower_bound = -10

        sol_here = d_species[n]["model"].optimize()
        if sol_here.status == "optimal":
            if sol_here.objective_value >= 0.01:
                d_species[n][cs]["growth cobra"] = sol_here
                d_cobra_growth[n][cs] = sol_here
            else:
                d_species[n][cs]["growth cobra"] = 0
                d_cobra_growth[n][cs] = 0
        else:
            d_cobra_growth[n][cs] = -1

        if EX_cs in [r.id for r in d_species[n]["model"].exchanges]:
            d_species[n]["model"].reactions.get_by_id(EX_cs).lower_bound = 0

###
# save dictionary with cobra growth values
makedirs(path_results_script, exist_ok=True)
with open(join(path_results_script, "d_cobra_growth.yaml"), "w") as file:
    documents = yaml.dump(d_cobra_growth, file)

###
# fraction of strains that grow on each carbon source acording to cobra results
f_species_grow = [0] * len(carbon_sources)
for n in d_species:
    for ic, cs in enumerate(carbon_sources):
        if d_species[n][cs]["growth cobra"] != 0:
            f_species_grow[ic] += 1

f_species_grow = [f / len(d_species) for f in f_species_grow]
for ic, cs in enumerate(carbon_sources):
    print(cs, f_species_grow[ic])

"""
3pg 0.0
cyst__L 0.0
glu__D 0.0
icit 0.0
pep 0.0
r5p 0.016129032258064516
akg 0.03225806451612903
met__L 0.016129032258064516
oaa 0.04838709677419355
f6p 0.06451612903225806
phe__L 0.06451612903225806
trp__L 0.0967741935483871
rmn 0.14516129032258066
cys__L 0.16129032258064516
for 0.1935483870967742
ile__L 0.1935483870967742
rib__D 0.20967741935483872
ser__D 0.22580645161290322
acgam 0.24193548387096775
cit 0.3064516129032258
ac 0.3548387096774194
pyr 0.3870967741935484
glc__D 0.41935483870967744
pro__L 0.43548387096774194
"""

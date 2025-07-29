#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 06:50:07 2024

@author: magdalena
"""

import cobra
from cobra.io import read_sbml_model
from cobra.medium import minimal_medium
from cobra.flux_analysis import single_reaction_deletion

import random
import numpy as np
import matplotlib.pyplot as plt

#read model and print the basic
d_models = {}

for model_name in ['Rickettsia_GCF_000284195.1.xml','Rickettsia_GCF_000284195.1_gramneg.xml']:
    model = read_sbml_model(model_name)
    model.solver = 'glpk'
    d_models[model_name] = model
    
    
#super essential metabolites: exchange reactions that are essential when everything is allowed to be taken up
d_super_ess_exchanges = {}
for model_name in d_models:
    d_super_ess_exchanges[model_name] = []

for model_name in d_models:
    for ex in d_models[model_name].exchanges:
        ex.lower_bound = -1000

    super_es_metabs = single_reaction_deletion(d_models[model_name], d_models[model_name].exchanges)

    super_es_metabs.iloc[0].ids
    super_ess_exchanges = []

    for i in range(len(super_es_metabs)):
        if super_es_metabs.iloc[i].growth<0.0001:
            d_super_ess_exchanges[model_name].append(list(super_es_metabs.iloc[i].ids)[0])

    print(model_name)
    print(d_super_ess_exchanges[model_name])
    
 
    
f1, ax1 = plt.subplots()
f2, ax2 = plt.subplots()
f3, ax3 = plt.subplots()

color_green = "#b9d7b1"
color_orange = "#eca04f"
palette = [color_green, color_orange] 

d_growth_media = {}
for model_name in d_models:
    d_growth_media[model_name] = []
    
for im,model_name in enumerate(d_models):
    
    components_range = np.arange(10, len(d_models[model_name].exchanges), 10).tolist()

    number_compoents_medium = []
    growth_medium = []
    trials_with_growth = np.zeros(len(components_range))

    trials = 50
    for ir,random_components in enumerate(components_range):
        for i in range(trials):
            #create random media
            basic_medium = [r.id for r in random.sample(d_models[model_name].exchanges, random_components)]
            medium = basic_medium + [ex for ex in d_super_ess_exchanges[model_name] if not ex in basic_medium]
            #print(len(medium))
            #print(medium)
            number_compoents_medium.append(len(medium))

            #set medium
            for ex in d_models[model_name].exchanges:
                if (ex.id in medium):
                    ex.lower_bound = -1000
                else:
                    ex.lower_bound = 0

            #max growth
            sol = d_models[model_name].optimize()
            #print(sol.objective_value)
            growth_medium.append(sol.objective_value)

            if sol.objective_value > 0.0001:
                trials_with_growth[ir] += 1
                
                #save medium
                d_growth_media[model_name].append(medium)

    #plot growth vs number of components
    ax1.plot(number_compoents_medium, growth_medium, 'o', color=palette[im])
    ax1.set_xlabel('Total components in medium', size=12, family='Arial')
    ax1.set_ylabel('Growth', size=12, family='Arial')
    #plt.show()
    
    #hist growth
    ax2.hist(growth_medium, color=palette[im], alpha=0.4)
    ax2.set_xlabel('Growth')

    #plot media showing growth vs number of components
    ax3.plot(components_range, [v/trials for v in trials_with_growth], 'o', color=palette[im], label=model_name)
    ax3.set_xlabel('Components in basic medium', size=12, family='Arial')
    ax3.set_ylabel('Fraction of media showing growth', size=12, family='Arial')
    plt.legend()
    f3.savefig("FmediaGrowth_vs_Ncomponents.pdf", bbox_inches='tight')

    #plt.show()
    
 
d_metabs_ess = {}
for model_name in d_models:
    d_metabs_ess[model_name] = {}

for model_name in d_models:
    print(model_name)
    for media in d_growth_media[model_name][:]:
        #print('new medium')
        #set media
        for ex in d_models[model_name].exchanges:
            if ex.id in media:
                ex.lower_bound = -1000
            else:
                ex.lower_bound = 0
        
        #essential metabs
        ess_metabs = single_reaction_deletion(d_models[model_name], d_models[model_name].exchanges)
        for i in range(len(ess_metabs)):
            if ess_metabs.iloc[i].growth<0.0001:
                ess_metab = list(ess_metabs.iloc[i].ids)[0]
                if not ess_metab in d_metabs_ess[model_name]:
                    d_metabs_ess[model_name][ess_metab] = 1
                else:
                    d_metabs_ess[model_name][ess_metab] += 1
                    
        #print(len(d_metabs_ess[model_name]))
    



for model_name in d_metabs_ess:
    f0, ax = plt.subplots(figsize=(30, 25))
    ax.plot(range(len(d_metabs_ess[model_name])), [d_metabs_ess[model_name][m]/len(d_growth_media[model_name]) for m in d_metabs_ess[model_name]], 'o')
    ax.set_xticks(range(len(d_metabs_ess[model_name])), labels=[m for m in d_metabs_ess[model_name]], rotation='vertical', size=12, family='Arial')
    
    
###
#every essential metab
ess_metab_all = [m for m in d_metabs_ess['Rickettsia_GCF_000284195.1.xml']]+[m for m in d_metabs_ess['Rickettsia_GCF_000284195.1_gramneg.xml'] if not m in d_metabs_ess['Rickettsia_GCF_000284195.1.xml']]

f_ess_metab_Rgen = []
for m in ess_metab_all:
    if m in d_metabs_ess['Rickettsia_GCF_000284195.1.xml']:
        f_ess_metab_Rgen.append(d_metabs_ess['Rickettsia_GCF_000284195.1.xml'][m])
    else:
        f_ess_metab_Rgen.append(0)
        
f_ess_metab_Rgramneg = []
for m in ess_metab_all:
    if m in d_metabs_ess['Rickettsia_GCF_000284195.1_gramneg.xml']:
        f_ess_metab_Rgramneg.append(d_metabs_ess['Rickettsia_GCF_000284195.1_gramneg.xml'][m])
    else:
        f_ess_metab_Rgramneg.append(0)
        
#sort
f_ess_metab_Rgen_sorted = [x/len(d_growth_media['Rickettsia_GCF_000284195.1.xml']) for x,y in sorted(zip(f_ess_metab_Rgen,ess_metab_all))]
ess_metab_all_sorted = [y for x,y in sorted(zip(f_ess_metab_Rgen,ess_metab_all))]
f_ess_metab_Rgramneg_sorted = [y/len(d_growth_media['Rickettsia_GCF_000284195.1_gramneg.xml']) for x,y in sorted(zip(f_ess_metab_Rgen,f_ess_metab_Rgramneg))]

f0, ax = plt.subplots(figsize=(30, 25))
ax.plot(range(len(ess_metab_all_sorted)), f_ess_metab_Rgen_sorted, 'o', color=palette[0], label='Rickettsia_GCF_000284195.1')
ax.plot(range(len(ess_metab_all_sorted)), f_ess_metab_Rgramneg_sorted, 'x', color=palette[1], label='Rickettsia_GCF_000284195.1_gramneg')
ax.set_xticks(range(len(ess_metab_all_sorted)), labels=ess_metab_all_sorted, rotation='vertical', size=12, family='Arial')
plt.legend()
f0.savefig("freq_essential_metab.pdf", bbox_inches='tight')




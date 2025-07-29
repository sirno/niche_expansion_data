#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 11:53:05 2021

@author: magdalena
"""

#import numpy as np
import copy
import yaml
from os import listdir
#from os.path import isfile, isdir, join


"""def number_of_supplying_communities(d_output,carbon_source,species):
    number_of_communities = len(d_output[carbon_source][species])
    return(number_of_communities)"""

# def open_misosoup_file(path_results, carbon_source):
#     with open(path_results+carbon_source+'.yaml','r') as file:
#         d_output = yaml.load(file, Loader=yaml.FullLoader)
#     return(d_output)

# def communities_size_and_number(d_output,carbon_source,species):
#     number_of_communities = len(d_output[carbon_source][species])
#     size_of_communities = [0]*number_of_communities 
#     for ic,comm in enumerate(d_output[carbon_source][species]):
#         size_of_communities[ic] = len([s for s in comm if s.startswith('Growth_')])
        
#     #if community has no species, set number of communities to zero    
#     if number_of_communities==1 and size_of_communities[0]==1 and d_output[carbon_source][species][0]['Growth_'+species]==0:
#        number_of_communities=0
#        size_of_communities=[0]
        
#     return(number_of_communities,size_of_communities)

def community_members(d_misosoup,carbon_source,min_or_species):
    Comm_members = []
    for com in d_misosoup[carbon_source][min_or_species]:
        if 'community' in com:
            community_members = list(com['community'].keys())
            Comm_members.append(community_members)
            
    return(Comm_members)
            

def mc_growth_incommunity(d_misosoup,carbon_source):
    species_grow_in_MC = []
    
    Comm_members = community_members(d_misosoup,carbon_source,'min')
    for comm in Comm_members:
        for sp in comm:
            species_grow_in_MC.append(sp.split('y_')[1])
        
    return(species_grow_in_MC)
        
    
def msc_nogrowth_inisolation_incommunity(d_misosoup,carbon_source,species):
    #returns whether misosoup found that a species grows in isolation, in community or it doesnt grow in a carbon source
    growth_0 = 'Growth_'+species
    Comm_members = community_members(d_misosoup,carbon_source,species)
    
    if growth_0 in d_misosoup[carbon_source][species][0]:
        does_it_grow = 'no growth' #it does not grow
    elif len(Comm_members[0])==1:
        does_it_grow = 'in isolation'
    else:
        does_it_grow = 'in community'
    
    return(does_it_grow)
    

def growth_value_in_inisolation(d_misosoup,carbon_source,species):
    #returns misosoup growth value when species grows in isolation
    growth_0 = 'Growth_'+species
    growth_value = d_misosoup[carbon_source][species][0]['solution'][growth_0]
            
    return(growth_value)
    
    
# def growth_inisolation(d_output,carbon_source,species):
#     number_of_communities,size_of_communities,does_it_grow = nogrowth_inisolation_incommunity(d_output,carbon_source,species)
#     if does_it_grow == 'in isolation':
#         growth_inisolation_value = d_output[carbon_source][species][0]['Growth_'+species]
#     else:
#         growth_inisolation_value = 0
    
#     return(growth_inisolation_value)
    
# def get_community_members(d_output,carbon_source,species):
#     number_of_communities = len(d_output[carbon_source][species])
#     community_members = [0]*number_of_communities 
#     for ic,comm in enumerate(d_output[carbon_source][species]):
#         community_members[ic] = [s[2:] for s in comm if s.startswith('y_')]
#     return(community_members)

# #######
# #Get species and carbon sources for which misosoup finds ommunities with 2 or more members
# ####### 
# def species_carbonsource_communities(path_species, cs_sorted_fba):
#     #returns dictionary with the species and carbon sources for which misosoup found communities with 2 or more members
#     all_species = [d for d in listdir(path_species) if d[0]!='.']
#     d_communities = {}
#     for species in all_species:
#         path_results = path_species + species + '/'
#         for carbon_source in cs_sorted_fba:
#             #open misosoup output
#             d_output = open_misosoup_file(path_results, carbon_source)
#             #get community members
#             communities = get_community_members(d_output,carbon_source,species)
#             #get community reactions for each of the misosoup communities
#             for comm in communities:
#                 comm_size = len(comm)#species in community
#                 #store the species and carbon source where there are communities with 2 or more species
#                 if comm_size>1:
#                     if not species in d_communities:
#                         d_communities[species] = []
#                     if not carbon_source in d_communities[species]:
#                         d_communities[species].append(carbon_source)
#     return(d_communities)

# #######
# #Functions to get the metabolites cross-fed in a community
# ####### 
# def consumption_secretion_count(species_fluxes, threshold_flux_not_zero):
#     #count negative and positive fluxes
#     neg_count, pos_count = 0, 0
      
#     # iterating each number in list
#     for num in species_fluxes:
#         # checking condition
#         if num > threshold_flux_not_zero:
#             pos_count += 1
#         elif num < -1*threshold_flux_not_zero:
#             neg_count += 1
#     return(neg_count, pos_count)

# def find_community_crossfed_metabolites(species, d_info_community_exchanges, metabs_in_media, threshold_flux_not_zero):
#     #remove all keys that don't belong to info about exchange reactions
#     for k in d_info_community_exchanges.copy():
#         if not k.startswith('R_EX_'):
#             del d_info_community_exchanges[k]
    
#     #remove metabs that are found in the medium
#     for k in d_info_community_exchanges.copy():
#         for metab in metabs_in_media:
#             if k.startswith(metab):
#                 del d_info_community_exchanges[k]
#                 break
                
#     #create dic to store the consumption secretion of each metab by each species
#     d_consumed_secreted = {}
    
#     for k in d_info_community_exchanges:
#         if k.endswith('_i'):
#             shared_str_ends = k.find('_e_')
#             new_metab_key = k[:(shared_str_ends+2)]
#             new_species_key = k[shared_str_ends+3:-2]
#             if not new_metab_key in d_consumed_secreted:
#                 d_consumed_secreted[new_metab_key] = {}
#             d_consumed_secreted[new_metab_key][new_species_key] = d_info_community_exchanges[k]
            
#     #get cross-fed metabolites: secreted by one species (positive flux) and consumed by other (negative flux)
#     crossfed_metabolites = []
#     crossfed_metabolites_consumed_by_focal = []
#     crossfed_metabolites_secreted_by_focal = []
#     crossfed_metabolites_among_suppliers = []
    
#     for k in d_consumed_secreted:
#         species_fluxes = [d_consumed_secreted[k][s] for s in d_consumed_secreted[k]]
#         cons_count, secr_count = consumption_secretion_count(species_fluxes, threshold_flux_not_zero)
#         if cons_count>0 and secr_count>0:
#             crossfed_metabolites.append(k[5:-2])
#             if species in d_consumed_secreted[k]:
#                 if d_consumed_secreted[k][species] > threshold_flux_not_zero:
#                     crossfed_metabolites_secreted_by_focal.append(k[5:-2])
#                 elif d_consumed_secreted[k][species] < -1*threshold_flux_not_zero:
#                     crossfed_metabolites_consumed_by_focal.append(k[5:-2])                
#             else:
#                 crossfed_metabolites_among_suppliers.append(k[5:-2])

#     return(crossfed_metabolites, crossfed_metabolites_consumed_by_focal, crossfed_metabolites_secreted_by_focal, crossfed_metabolites_among_suppliers)

# #######
# #Functions for the analysis of community motifs
# ####### 
# #is the cross_feeding interaction mutualistic? (is the focal strain secreting something the suppliers take up?)
# def is_mutualistic_interaction(species, carbon_source, comm, threshold_flux_not_zero, d_media):
#     comm_here = copy.deepcopy(comm)
#     mutualism = 0
#     #remove all keys that don't belong to info about exchange reactions
#     #for k in comm.copy():
#     for k in comm_here.copy():
#         if not k.startswith('R_EX_'):
#             del comm_here[k]
                
#     #remove metabs that are found in the medium
#     #for k in comm.copy():
#     for k in comm_here.copy():    
#         for metab in d_media:
#             if k.startswith(metab):
#                 del comm_here[k]
#                 break
                
#     #find metabolites secreted by the focal strain            
#     secretions_focal = []
#     for reac in comm_here:
#         if reac.endswith('_e_'+species+'_i') and comm_here[reac]>threshold_flux_not_zero:
#             secretions_focal.append(reac)
    
#     #check if a supplier consumes anything the focal secretes
#     for metab in secretions_focal:
#         metab_general_name_ends = metab.find('_e_')
#         metab_general_name = metab[:metab_general_name_ends+3]
#         for reac in comm_here:
#             if reac.startswith(metab_general_name) and reac!=metab and comm_here[reac]<-1*threshold_flux_not_zero:
#                 mutualism = 1
#                 break
            
#     return(mutualism)

# #does the focal strain consume the carbon source from the medium?
# def does_focal_consume_cs(species, carbon_source, comm, threshold_flux_not_zero):
#     ex_cs_focal_srtain = 'R_EX_'+carbon_source+'_e_'+species+'_i'
#     does_focal_consume_cs_value = 0
#     if (ex_cs_focal_srtain in comm) and (comm[ex_cs_focal_srtain]<-1*threshold_flux_not_zero):
#         does_focal_consume_cs_value = 1
#     return(does_focal_consume_cs_value)

# #does a suppling strain consume the carbon source from the medium?
# def does_suppling_consume_cs(species, carbon_source, comm, threshold_flux_not_zero):
#     suppling_species = [s for s in comm if s.startswith('y_') and s!='y_'+species]
#     does_suppling_consume_cs_value = 0
#     for y_ss in suppling_species:
#         ss = y_ss[2:]
#         ex_cs_suppling_strain = 'R_EX_'+carbon_source+'_e_'+ss+'_i'
#         if (ex_cs_suppling_strain in comm) and (comm[ex_cs_suppling_strain]<-1*threshold_flux_not_zero):
#             does_suppling_consume_cs_value = 1
#     return(does_suppling_consume_cs_value)

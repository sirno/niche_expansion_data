#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 12:20:10 2024

@author: magdalena
"""

import matplotlib.pyplot as plt
import networkx as nx
import path
import yaml

path_results_script = "out/analysis/"
path_results_misosoup = "out/misosoup/"

###
# load data needed
# for tolerance in [7,8,9]:
tolerance = 7
print("Tolerance: ", tolerance)

with open(
    path.join(path_results_misosoup, f"tolerance_1e-{tolerance}.yaml"), "r"
) as file:
    d_misosoup = yaml.load(file, Loader=yaml.SafeLoader)

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


###
# 2 communities in glu--D, all of which have D2M02
# D2M02 takes up glu--D and secretes metabs

d_D2M02 = {}
d_D2M02["secretions_crossfed"] = {}
d_D2M02["consumptions_crossfed"] = {}

threshold_flux_not_zero = 0.001

for comm in d_misosoup["glu__D"]["min"]:
    if ("y_F3M07_" not in comm["community"]) and ("y_E3R11_" not in comm["community"]):
        D2M02_secretions = focal_secretions_fun("D2M02", comm, threshold_flux_not_zero)
        D2M02_secretions_crossfed = focal_secretions_crossfed_fun(
            D2M02_secretions, comm, threshold_flux_not_zero
        )
        for m in D2M02_secretions_crossfed:
            if m not in d_D2M02["secretions_crossfed"]:
                d_D2M02["secretions_crossfed"][m] = 1 / len(d_misosoup["glu__D"]["min"])
            else:
                d_D2M02["secretions_crossfed"][m] += 1 / len(
                    d_misosoup["glu__D"]["min"]
                )

        # D2M02 consumptions secreted by other comm members?
        comm_member = [s for s in comm["community"] if s != "y_D2M02"]
        comm_member = comm_member[0].split("y_")[1]

        comm_member_secretions = focal_secretions_fun(
            comm_member, comm, threshold_flux_not_zero
        )
        D2M02_consumptions_crossfed = focal_secretions_crossfed_fun(
            comm_member_secretions, comm, threshold_flux_not_zero
        )
        for m in D2M02_consumptions_crossfed:
            if m not in d_D2M02["consumptions_crossfed"]:
                d_D2M02["consumptions_crossfed"][m] = 1 / len(
                    d_misosoup["glu__D"]["min"]
                )
            else:
                d_D2M02["consumptions_crossfed"][m] += 1 / len(
                    d_misosoup["glu__D"]["min"]
                )


# plot

secretions_frequency = [
    d_D2M02["secretions_crossfed"][m] for m in d_D2M02["secretions_crossfed"]
]
secretions_metabs = [m for m in d_D2M02["secretions_crossfed"]]
secretions_frequency_sorted = [
    y for y, x in sorted(zip(secretions_frequency, secretions_metabs), reverse=True)
]
secretions_metabs_sorted = [
    x for y, x in sorted(zip(secretions_frequency, secretions_metabs), reverse=True)
]

consumptions_frequency = [
    d_D2M02["consumptions_crossfed"][m] for m in d_D2M02["consumptions_crossfed"]
]
consumptions_metabs = [m for m in d_D2M02["consumptions_crossfed"]]
consumptions_frequency_sorted = [
    y for y, x in sorted(zip(consumptions_frequency, consumptions_metabs), reverse=True)
]
consumptions_metabs_sorted = [
    x for y, x in sorted(zip(consumptions_frequency, consumptions_metabs), reverse=True)
]

f1 = plt.figure()
plt.plot(range(len(secretions_metabs_sorted)), secretions_frequency_sorted, "ko")
plt.xticks(
    range(len(secretions_metabs_sorted)),
    labels=secretions_metabs_sorted,
    rotation="vertical",
    size=12,
    family="Arial",
)
plt.yticks([0, 0.5, 1], [0, 0.5, 1])
plt.ylabel("Cross-feeding frequency - D2M02 secretions", size=12, family="Arial")
plt.show()
f1.savefig(
    path_results_script + "D2M02_secretions_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)

f2 = plt.figure()
plt.plot(range(len(consumptions_metabs_sorted)), consumptions_frequency_sorted, "ko")
plt.xticks(
    range(len(consumptions_metabs_sorted)),
    labels=consumptions_metabs_sorted,
    rotation="vertical",
    size=12,
    family="Arial",
)
plt.yticks([0, 0.5, 1], [0, 0.5, 1])
plt.ylabel("Cross-feeding frequency - D2M02 consumptions", size=12, family="Arial")
plt.show()
f2.savefig(
    path_results_script + "D2M02_consumptions_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)

##
# network drawing

# Create a directed graph
G = nx.DiGraph()

# Define the central node
central_node = "D2M02"

# Define outgoing edges (from central node to others)
outgoing_edges = d_D2M02["secretions_crossfed"]

# Define incoming edges (from others to central node)
incoming_edges = d_D2M02["consumptions_crossfed"]
incoming_edges["glu__D"] = 1

# write file for Alberto
fil = open(path_results_script + central_node + "_info_network.csv", "w")
for k in outgoing_edges:
    fil.write(central_node + " -> " + k + " " + str(outgoing_edges[k]) + "\n")

for k in incoming_edges:
    fil.write(central_node + " <- " + k + " " + str(incoming_edges[k]) + "\n")

fil.close()

# Add outgoing edges to the graph
for target, width in outgoing_edges.items():
    G.add_edge(central_node, target, weight=width)

# Add incoming edges to the graph
for source, width in incoming_edges.items():
    G.add_edge(source, central_node, weight=width)

# Improved layout
pos = {}
distance = 2  # Distance from the central node
side_spacing_left = -4
side_spacing_right = 3

# Identify dual-connection nodes
dual_nodes = set(incoming_edges.keys()).intersection(outgoing_edges.keys())
incoming_only = set(incoming_edges.keys()) - dual_nodes
outgoing_only = set(outgoing_edges.keys()) - dual_nodes

# Position incoming-only nodes above
for i, node in enumerate(sorted(incoming_only)):
    pos[node] = (i - len(incoming_only) / 2, distance)

# Position outgoing-only nodes below
for i, node in enumerate(sorted(outgoing_only)):
    pos[node] = (i - len(outgoing_only) / 2, -distance)

# Position dual-connection nodes on the sides
for i, node in enumerate(sorted(dual_nodes)):
    if i % 2 == 0:
        pos[node] = (side_spacing_left, (len(dual_nodes) / 2 - i) / 2)
    else:
        pos[node] = (side_spacing_right, (len(dual_nodes) / 2 - i) / 2)

# Place the central node in the center
pos[central_node] = (0, 0)

# Node sizes
node_sizes = {node: (5000 if node == central_node else 800) for node in G.nodes()}

# Improved aesthetics
fig = plt.figure(figsize=(10, 7))

# Draw nodes
nx.draw_networkx_nodes(
    G,
    pos,
    node_size=[node_sizes[node] for node in G.nodes()],  # Use dynamic node sizes
    node_color=["lightblue" if node == central_node else "grey" for node in G.nodes()],
    alpha=0.3,
    linewidths=1.5,
)

# Draw edges with smooth transitions, variable widths, and adjusted margins
edges = G.edges(data=True)
for u, v, d in edges:
    margin = 40 if v == central_node else 20
    margin_source = 20 if v == central_node else 40
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=[(u, v)],
        width=d["weight"] * 3,  # Scale widths for visibility
        connectionstyle="arc3,rad=0.2",  # Rounded edges
        edge_color="grey",
        arrowsize=15,  # Bigger arrowheads
        min_source_margin=margin_source,  # Ensure arrowheads clear source node
        min_target_margin=margin,  # Adjust dynamically for large target nodes
    )


# Draw labels
nx.draw_networkx_labels(G, pos, font_size=20, font_color="black", font_weight="bold")

# Remove the grid
plt.axis("off")

# Display the graph
plt.show()


fig.savefig(
    path_results_script + "D2M02_mechanism_network_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


###
# 13 communities in 3pg, all of which have E3R18
# E3R18 takes up 3pg and secretes metabs

d_E3R18 = {}
d_E3R18["secretions_crossfed"] = {}
d_E3R18["consumptions_crossfed"] = {}

threshold_flux_not_zero = 0.001

for comm in d_misosoup["3pg"]["min"]:
    if ("y_F3M07_" not in comm["community"]) and ("y_E3R11_" not in comm["community"]):
        E3R18_secretions = focal_secretions_fun("E3R18", comm, threshold_flux_not_zero)
        E3R18_secretions_crossfed = focal_secretions_crossfed_fun(
            E3R18_secretions, comm, threshold_flux_not_zero
        )
        for m in E3R18_secretions_crossfed:
            if m not in d_E3R18["secretions_crossfed"]:
                d_E3R18["secretions_crossfed"][m] = 1 / len(d_misosoup["3pg"]["min"])
            else:
                d_E3R18["secretions_crossfed"][m] += 1 / len(d_misosoup["3pg"]["min"])

        # E3R18 consumptions secreted by other comm members?
        comm_member = [s for s in comm["community"] if s != "y_E3R18"]
        comm_member = comm_member[0].split("y_")[1]

        comm_member_secretions = focal_secretions_fun(
            comm_member, comm, threshold_flux_not_zero
        )
        E3R18_consumptions_crossfed = focal_secretions_crossfed_fun(
            comm_member_secretions, comm, threshold_flux_not_zero
        )
        for m in E3R18_consumptions_crossfed:
            if m not in d_E3R18["consumptions_crossfed"]:
                d_E3R18["consumptions_crossfed"][m] = 1 / len(d_misosoup["3pg"]["min"])
            else:
                d_E3R18["consumptions_crossfed"][m] += 1 / len(d_misosoup["3pg"]["min"])


# plot

secretions_frequency = [
    d_E3R18["secretions_crossfed"][m] for m in d_E3R18["secretions_crossfed"]
]
secretions_metabs = [m for m in d_E3R18["secretions_crossfed"]]
secretions_frequency_sorted = [
    y for y, x in sorted(zip(secretions_frequency, secretions_metabs), reverse=True)
]
secretions_metabs_sorted = [
    x for y, x in sorted(zip(secretions_frequency, secretions_metabs), reverse=True)
]

consumptions_frequency = [
    d_E3R18["consumptions_crossfed"][m] for m in d_E3R18["consumptions_crossfed"]
]
consumptions_metabs = [m for m in d_E3R18["consumptions_crossfed"]]
consumptions_frequency_sorted = [
    y for y, x in sorted(zip(consumptions_frequency, consumptions_metabs), reverse=True)
]
consumptions_metabs_sorted = [
    x for y, x in sorted(zip(consumptions_frequency, consumptions_metabs), reverse=True)
]

f1 = plt.figure()
plt.plot(range(len(secretions_metabs_sorted)), secretions_frequency_sorted, "ko")
plt.xticks(
    range(len(secretions_metabs_sorted)),
    labels=secretions_metabs_sorted,
    rotation="vertical",
    size=12,
    family="Arial",
)
plt.yticks([0, 0.5, 1], [0, 0.5, 1])
plt.ylabel("Cross-feeding frequency- E3R18 secretions", size=12, family="Arial")
plt.show()
f1.savefig(
    path_results_script + "E3R18_secretions_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)

f2 = plt.figure()
plt.plot(range(len(consumptions_metabs_sorted)), consumptions_frequency_sorted, "ko")
plt.xticks(
    range(len(consumptions_metabs_sorted)),
    labels=consumptions_metabs_sorted,
    rotation="vertical",
    size=12,
    family="Arial",
)
plt.yticks([0, 0.5, 1], [0, 0.5, 1])
plt.ylabel("Cross-feeding frequency - E3R18 consumptions", size=12, family="Arial")
plt.show()
f2.savefig(
    path_results_script + "E3R18_consumptions_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)


##
# network drawing

# Create a directed graph
G = nx.DiGraph()

# Define the central node
central_node = "E3R18"

# Define outgoing edges (from central node to others)
outgoing_edges = d_E3R18["secretions_crossfed"]

# Define incoming edges (from others to central node)
incoming_edges = d_E3R18["consumptions_crossfed"]
incoming_edges["3pg"] = 1

# write file for Alberto
fil = open(path_results_script + central_node + "_info_network.csv", "w")
for k in outgoing_edges:
    fil.write(central_node + " -> " + k + " " + str(outgoing_edges[k]) + "\n")

for k in incoming_edges:
    fil.write(central_node + " <- " + k + " " + str(incoming_edges[k]) + "\n")

fil.close()

# Add outgoing edges to the graph
for target, width in outgoing_edges.items():
    G.add_edge(central_node, target, weight=width)

# Add incoming edges to the graph
for source, width in incoming_edges.items():
    G.add_edge(source, central_node, weight=width)

# Improved layout
pos = {}
distance = 4  # Distance from the central node
side_spacing_left = -5
side_spacing_right = 4

# Identify dual-connection nodes
dual_nodes = set(incoming_edges.keys()).intersection(outgoing_edges.keys())
incoming_only = set(incoming_edges.keys()) - dual_nodes
outgoing_only = set(outgoing_edges.keys()) - dual_nodes

# Position incoming-only nodes above
for i, node in enumerate(sorted(incoming_only)):
    pos[node] = (i - len(incoming_only) / 2, distance)

# Position outgoing-only nodes below
for i, node in enumerate(sorted(outgoing_only)):
    pos[node] = (i - len(outgoing_only) / 2, -distance)

# Position dual-connection nodes on the sides
for i, node in enumerate(sorted(dual_nodes)):
    if i % 2 == 0:
        pos[node] = (side_spacing_left, (len(dual_nodes) / 2 - i) / 2)
    else:
        pos[node] = (side_spacing_right, (len(dual_nodes) / 2 - i) / 2)

# Place the central node in the center
pos[central_node] = (0, 0)

# Node sizes
node_sizes = {node: (5000 if node == central_node else 800) for node in G.nodes()}

# Improved aesthetics
fig = plt.figure(figsize=(14, 10))

# Draw nodes
nx.draw_networkx_nodes(
    G,
    pos,
    node_size=[node_sizes[node] for node in G.nodes()],  # Use dynamic node sizes
    node_color=["lightblue" if node == central_node else "grey" for node in G.nodes()],
    alpha=0.3,
    linewidths=1.5,
)

# Draw edges with smooth transitions, variable widths, and adjusted margins
edges = G.edges(data=True)
for u, v, d in edges:
    margin = 40 if v == central_node else 20
    margin_source = 20 if v == central_node else 40
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=[(u, v)],
        width=d["weight"] * 3,  # Scale widths for visibility
        connectionstyle="arc3,rad=0.2",  # Rounded edges
        edge_color="grey",
        arrowsize=15,  # Bigger arrowheads
        min_source_margin=margin_source,  # Ensure arrowheads clear source node
        min_target_margin=margin,  # Adjust dynamically for large target nodes
    )


# Draw labels
nx.draw_networkx_labels(G, pos, font_size=15, font_color="black", font_weight="bold")

# Remove the grid
plt.axis("off")

# Display the graph
plt.show()


fig.savefig(
    path_results_script + "E3R18_mechanism_network_tol_" + str(tolerance) + ".pdf",
    bbox_inches="tight",
)

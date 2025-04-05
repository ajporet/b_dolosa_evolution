#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 06:39:01 2024

@author: Alexandra
"""

from ete3 import Tree

# Load the trees
tree1 = Tree("B_dolosa_recomb_75.tree", format=1)
tree2 = Tree("JQR_tree.tree", format=1)

# Compute the Robinson-Foulds distance and get a list of different nodes
rf, max_rf, common_edges, edges_in_tree1, edges_in_tree2, _, _ = tree1.robinson_foulds(tree2, unrooted_trees=True)
print ("Partitions in tree2 that were not found in tree1:", edges_in_tree1 - edges_in_tree2)
print ("Partitions in tree1 that were not found in tree2:", edges_in_tree2 - edges_in_tree1)

two_not_one =  list(edges_in_tree2 - edges_in_tree1)
one_not_two = edges_in_tree1 - edges_in_tree2
list_of_lists = [list(list(y) for y in s) for s in two_not_one]

# Print results
print("Robinson-Foulds Distance:", rf)
print("Max possible RF:", max_rf)
print("Common edges:", common_edges)
print("Unique edges in tree1:", edges_in_tree1)
print("Unique edges in tree2:", edges_in_tree2)

###################################
length_differences = []
topology_differences = []

# Check if the trees have the same topology (structure) using the Robinson-Foulds metric
rf_result = tree1.compare(tree2, unrooted=True)


# Traverse both trees in parallel to compare branch lengths
for node1, node2 in zip(tree1.traverse("preorder"), tree2.traverse("preorder")):
    # Check for matching leaf names
    if node1.is_leaf() and node1.name != node2.name:
        raise ValueError("Trees have mismatched leaf names; ensure they have the same labels for comparison.")
    
    # Compare branch lengths
    if abs(node1.dist - node2.dist) > 1e-6:  # Small tolerance for floating-point differences
        length_differences.append((node1.name, node1.dist, node2.dist))

# Output the results
if length_differences:
    print("\nBranches with changed lengths:")
    for name, dist1, dist2 in length_differences:
        print(f"Node: {name} | Length in Tree1: {dist1} | Length in Tree2: {dist2}")
else:
    print("\nNo significant branch length changes found.")

# Output topology differences if any were found
if topology_differences:
    print("\nTopology differences found:")
    for (node1, node2) in topology_differences:
        print(f"Tree1 node: {node1} | Tree2 node: {node2}")
else:
    print("\nNo topology differences found.")

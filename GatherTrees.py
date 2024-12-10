"""
Script by Kristine Nørtoft Kristensen
Edited by August Søndergaard
"""

import os, re

os.chdir("/home/augusthcs/Bachelorprojekt/Astral/")

for n in range(1000):
    with open(str(n)+"_random_trees.aster.tre") as f: #opening the 1000 astral trees
        tree_w_support=f.read()
    bits= re.split("(\))",tree_w_support) #removing the support values
    tree_no_support= ""
    for b in bits:
        if len(b)>1 and ":" != b[0] and bits[0] != b and bits[-1] != b:
            tree_no_support= tree_no_support+b[8:]
        else: 
            tree_no_support=tree_no_support+b
            
    with open("astral_for_mrbayes.tre", "a") as new_tree: # make a new file with alle the 1000 trees without support values
        new_tree.write(tree_no_support)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 23:58:11 2020

@author: Alexandra
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 12:59:19 2020

@author: Alexandra
"""

import scipy.io
import numpy as np


## GET LIST OF OLD LOCUS TAG GENE NAMES WITH THEIR POSITIONS 
from Bio import SeqIO
filename = "NZ_CP009795.1.gb"
locus_to_gene = dict()


for record in SeqIO.parse(filename, "genbank"):
    for f in record.features:
        if f.type == "CDS" or f.type == "tRNA" or f.type == "rRNA":
            if "old_locus_tag" in f.qualifiers:
                #genes = f.qualifiers["gene"]
                locus_tags = f.qualifiers["old_locus_tag"]
                gene_start = f.location.start.position
                gene_end = f.location.end.position
                #assert len(genes) == 1, genes
                #assert len(locus_tags) == 1, locus_tags
                locus_to_gene[locus_tags[0]] = [gene_start, gene_end]
#print("Mapped %i locus tags to genes" % len(locus_to_gene))
#locus_to_gene['PAC1_03600'] = [764732,766465] 
#locus_to_gene['PAC1_03640'] = [771842,772555] 


######### READ KEGG FILES #########
                
kegg_gene_list_file = open('KEGG_pathway_gene_table.txt', "r")

kegg_gene_list = kegg_gene_list_file.read()
kegg_gene_list = kegg_gene_list.split("\n")

path_to_gene = dict()
path_to_gene_num = dict()

for ii in range(0,len(kegg_gene_list)):
    list_entry = kegg_gene_list[ii]
    pathway_num = list_entry[5+3:5+8]
    gene_label = list_entry[18:len(list_entry)]
    if gene_label in locus_to_gene:
        genome_pos = locus_to_gene[gene_label].copy()
    
        if pathway_num not in path_to_gene:
            path_to_gene[pathway_num] = [gene_label]
            path_to_gene_num[pathway_num] = [genome_pos]
        else:
            path_to_gene[pathway_num].append(gene_label)
            path_to_gene_num[pathway_num].append(genome_pos)
        
        
########## GET KEGG DESCRIPTIONS #############
kegg_path_descriptions_file = open('KEGG_pathway_numbers.txt', "r")
kegg_path_descriptions_list = kegg_path_descriptions_file.read()
kegg_path_descriptions_list = kegg_path_descriptions_list.split("\n")

translated_path_names = dict()

for ii in range(0,len(kegg_path_descriptions_list)):
    entry =  kegg_path_descriptions_list[ii]
    if entry[0:5] in list(path_to_gene.keys()):
        translated_path_names[entry[0:5]] = entry[6:len(entry)]
    

large_grouping_index = [entry for entry in range(0,len(kegg_path_descriptions_list)) if "***" in kegg_path_descriptions_list[entry]]
large_grouping_index.append(len(kegg_path_descriptions_list))

sub_grouping_index = [entry for entry in range(0,len(kegg_path_descriptions_list)) if "%" in kegg_path_descriptions_list[entry]]
sub_grouping_index.append(len(kegg_path_descriptions_list))
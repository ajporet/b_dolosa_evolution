clear all; close all; clc; 

%% Load in KEGG data reformatted from python script and create parsable dictionary-like structures

% Load data from tree making script
load analysis_m_outputs_structures.mat

% Read in tables of KEGG annotations 
chr_1 = readtable("Functional_analysis_map_template_chr_1.xlsx");
chr_2 = readtable("Functional_analysis_map_template_chr_2.xlsx");
chr_3 = readtable("Functional_analysis_map_template_chr_3.xlsx");

% Creates a nice dict-like structure that allows on to search for the KEGG
% annotations of any gene by location or name
% These functions return:
%   chr_1_gene_location: a list of gene name, start location, end location,
%   chromosome number 
%   chr_1_KEGG_ids: all associated kegg ids with a gene
[chr_1_gene_location, chr_1_KEGG_ids] = return_gene_to_KEGG_class("Functional_analysis_map_template_chr_1.xlsx",1);
[chr_2_gene_location, chr_2_KEGG_ids] = return_gene_to_KEGG_class("Functional_analysis_map_template_chr_2.xlsx",2);
[chr_3_gene_location, chr_3_KEGG_ids] = return_gene_to_KEGG_class("Functional_analysis_map_template_chr_3.xlsx",3);

% Combine locations and kegg ids (i.e variables generated above) into one structure
all_chr_desc = [chr_1_gene_location;chr_2_gene_location;chr_3_gene_location];
all_chr_loc = str2double(all_chr_desc(:,2:3));
all_chr_KEGG_id = [chr_1_KEGG_ids;chr_2_KEGG_ids;chr_3_KEGG_ids];

% Get all unique pathways present in dolosa
unique_pathways = unique(erase([chr_1.PathwayNumber;chr_2.PathwayNumber;chr_3.PathwayNumber],char(65279)));

unfiltered_pathways = erase([chr_1.PathwayNumber;chr_2.PathwayNumber;chr_3.PathwayNumber],char(65279));
unfiltered_genes = erase([chr_1.GeneList;chr_2.GeneList;chr_3.GeneList],[char(65279),"[","]"," ","'"]);

% Make a binary matrix of (genes, kegg pathways) where 1 equals a
% matching relationship
gene_to_kegg = zeros(length(all_chr_desc),length(unique_pathways));
for ii=1:length(unfiltered_pathways)
    pathway_loc = contains(unique_pathways,unfiltered_pathways(ii));
    gene_loc = ismember(all_chr_desc(:,1),split(unfiltered_genes(ii),","));
    gene_to_kegg(gene_loc,pathway_loc) = 1;
end

%% Figure out which pathways which genes map to 

% Get mutations unique to J or Q and R
temp_SampleNames = SampleNames; temp_SampleNames(samplestoplot)=="Excluded";
J_annotation_log = sum(hasmutation(:,contains(temp_SampleNames,"P02")),2)>0;
QR_annotation_log = sum(hasmutation(:,contains(temp_SampleNames,["P06","P07"])),2)>0 & ~J_annotation_log;

J_annotation_full = annotation_full(J_annotation_log); QR_annotation_full = annotation_full(QR_annotation_log); 

% Now find the number of mutations in each KEGG pathway 

% All JQR together
hit_ = (string([annotation_full.scaffold]) == all_chr_desc(:,4)) & ([annotation_full.pos]>=all_chr_loc(:,1)) &  ([annotation_full.pos]<=all_chr_loc(:,2));
hits_per_gene = sum(hit_,2);
number_of_hits_per_pathway_JQR = sum(hits_per_gene.*gene_to_kegg);

% Only J
hit_ = (string([J_annotation_full.scaffold]) == all_chr_desc(:,4)) & ([J_annotation_full.pos]>=all_chr_loc(:,1)) &  ([J_annotation_full.pos]<=all_chr_loc(:,2));
hits_per_gene = sum(hit_,2);
number_of_hits_per_pathway_J = sum(hits_per_gene.*gene_to_kegg);

% Only QR
hit_ = (string([QR_annotation_full.scaffold]) == all_chr_desc(:,4)) & ([QR_annotation_full.pos]>=all_chr_loc(:,1)) &  ([QR_annotation_full.pos]<=all_chr_loc(:,2));
hits_per_gene = sum(hit_,2);
number_of_hits_per_pathway_QR = sum(hits_per_gene.*gene_to_kegg);
hits_in_O_ant_pathway = QR_annotation_full(sum(hit_ .* gene_to_kegg(:,54))>0);

%% Design a simulation that tracks the average number of mutations expected within each pathway assuming random distribution
% Here, we simulate droping 541 (i.e the number of mutations observed in
% total in this new dataset) across the genome. We assume every mutation
% leads to a detectable change, and do not differentiate between synonymous
% and nonsynonymous mutations. After 10,000 trials, we compare the
% distribution of mutations/KEGG pathway to the number of mutations
% empiracally observed. For subanalyses where we look at a subset of
% mutations found in patients J or Q and R, we scale this distribution
% accordingly. Ex. if we randomly found 20 mutations in a pathway after
% dropping 541 randomly, we'd scale to expect 20/(100/541) instead. 
% From there we can calculate a pseudo-p-value, 
% i.e., 1 - the probability that the number of mutations observed is
% greater than the number randomly simulated 

num_mutations_to_scatter = length(annotation_full);
rng(24680)
num_trials = 50000;
kegg_counter_simulation = nan(length(unique_pathways), num_trials);
parfor ii=1:num_trials
    random_mut_loc = randi(3406596 + 2168075 + 834424,1,num_mutations_to_scatter);
    kegg_counter_simulation(:,ii) = get_hits_per_pathway(all_chr_desc,all_chr_loc,gene_to_kegg,random_mut_loc);
    if mod(ii,num_trials/100)==0
        disp("Current iteration: " + ii)
    end

end

J_scaler = sum(J_annotation_log)/num_mutations_to_scatter; QR_scaler = sum(QR_annotation_log)/num_mutations_to_scatter;

mock_total_enrich_p_val = 1-sum(number_of_hits_per_pathway_JQR'>=kegg_counter_simulation,2)./num_trials;
mock_J_enrich_p_val = 1-sum(number_of_hits_per_pathway_J'>=(kegg_counter_simulation.*(J_scaler)),2)./num_trials;
mock_QR_enrich_p_val = 1-sum(number_of_hits_per_pathway_QR'>=(kegg_counter_simulation.*(QR_scaler)),2)./num_trials;

total_sig_hits = mock_total_enrich_p_val<(.05/(128*2));
J_sig_hits = mock_J_enrich_p_val<(.05/(128*2));
QR_sig_hits = mock_QR_enrich_p_val<(.05/(128*2));

%% Create a distance matrix of KEGG overlab by gene number
% Only in patients Q and R do we observe enrichment for specific kegg
% pathways. To determine whether this enrichment is a byproduct of
% similarity between pathways, we create an asymetric distance matrix,
% where similarity between KEGG pathways is defined as the number of common
% genes / number of total genes. The entry (ii,jj) is relative to the total
% number of genes in ii, and (jj,ii) the total number of genes in jj.

kegg_overlap = zeros(size(gene_to_kegg,2),size(gene_to_kegg,2));
for ii=1:size(gene_to_kegg,2)
    for ij=1:size(gene_to_kegg,2)
        total_alike_genes = sum(gene_to_kegg(:,ij) & gene_to_kegg(:,ii));
        kegg_overlap(ii,ij) = total_alike_genes./sum(gene_to_kegg(:,ii));        
        kegg_overlap(ij,ii) = total_alike_genes./sum(gene_to_kegg(:,ij));
    end
end

dist_loc_sig_points = contains(string(unique_pathways),unique_pathways(QR_sig_hits));

QR_distances_of_sig_pathways = kegg_overlap(dist_loc_sig_points,dist_loc_sig_points);

%% Make a csv of the pathway results 

pathway_descriptions_textfile = readlines("KEGG_pathway_numbers.txt");
pathway_descriptions = cellfun(@(x) split(x,"  "),pathway_descriptions_textfile,'UniformOutput',false);
pathway_descriptions = (pathway_descriptions(cellfun(@length, pathway_descriptions)>1));
pathway_descriptions = string([pathway_descriptions{:}]');

raw_data = [number_of_hits_per_pathway_QR', mock_QR_enrich_p_val, QR_sig_hits, number_of_hits_per_pathway_J', mock_J_enrich_p_val, J_sig_hits];
reorder_pathways  = cellfun(@(x) find(string(unique_pathways)==x), pathway_descriptions(:,1)); 

excel_file_data = [pathway_descriptions,raw_data(reorder_pathways,:)];
similarity_matrix_reordered = kegg_overlap(reorder_pathways,reorder_pathways);

similarity_matrix_reordered_sig_hits = similarity_matrix_reordered(QR_sig_hits(reorder_pathways)>0,QR_sig_hits(reorder_pathways)>0);
similarity_matrix_reordered_sig_desc = pathway_descriptions(QR_sig_hits(reorder_pathways)>0,:);

%% Functions

function hits_per_kegg_pathway = get_hits_per_pathway(all_chr_desc,all_chr_loc,gene_to_kegg,random_mut_loc)
    tailed_mut_loc = random_mut_loc;
    scaffold = 1 + (random_mut_loc>3406596) + (random_mut_loc>(3406596 + 2168075));
    tailed_mut_loc(scaffold==2) = tailed_mut_loc(scaffold==2) - 3406596;
    tailed_mut_loc(scaffold==3) = tailed_mut_loc(scaffold==3) - (3406596+2168075);
    hit_ = (string(scaffold) == all_chr_desc(:,4)) & (tailed_mut_loc>=all_chr_loc(:,1)) &  (tailed_mut_loc<=all_chr_loc(:,2));
    hits_per_gene = sum(hit_,2);
    hits_per_kegg_pathway = sum(hits_per_gene.*gene_to_kegg);
end

function [gene_locations,pathway_nums] = return_gene_to_KEGG_class(excel_file_name, chr_num)
    chr_1 = readtable(excel_file_name);
    
    gene_list_1 = cellfun(@(x) string(split(erase(x,["[","]"," ","'","﻿ ",char(65279)]), ",")), chr_1.GeneList,'UniformOutput',false);
    gene_list_1 = vertcat(gene_list_1{:});
    gene_loc_1 = cellfun(@(x) split(x,"], ["), chr_1.GeneLocationList,'UniformOutput',false);
    gene_loc_1 = cellfun(@(x) str2double(split(erase(x,["[[","]]"," ","'",char(65279)]),",")),gene_loc_1,'UniformOutput',false);
    gene_loc_1(cellfun(@(x) size(x,2), gene_loc_1)~=2) = cellfun(@(x) x', gene_loc_1(cellfun(@(x) size(x,2), gene_loc_1)~=2),'UniformOutput',false);
    
    start_and_end = cell2mat(gene_loc_1);
    gene_locations = unique([gene_list_1,start_and_end,ones(length(start_and_end),1).*chr_num],'rows');
    pathway_nums_temp = repmat({},length(gene_locations),1);
    
    for ii=1:length(gene_locations)
        kegg_nums = chr_1.PathwayNumber(contains(chr_1.GeneList,"'" + gene_locations(ii,1) + "'"));
        pathway_nums_temp(ii) = {[erase(string(kegg_nums)',char(65279))]};
    end
    
    pathway_nums = pathway_nums_temp';
end
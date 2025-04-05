clear all; close all;
load QR_breseq_output.mat
load JQR_b_dolosa_mutations.mat

sample_names_cleaned = SampleNames(samplestoplot);
hasmutation_cleaned = hasmutation(:,samplestoplot);
P067_sample_names = sample_names_cleaned(contains(sample_names_cleaned,["P06","P07"]));

% Define LPS-affecting muts from manual inspection of data 
LPS_muts= [50, 164, 188, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 265, 438,186, 187, 201];
LPS_stop_codon_mut = 189;

% Breseq indels
Breseq_names = Breseq_Vals_cleaned(:,13:end).Properties.VariableNames;
breseq_table_LPS = table2array(Breseq_Vals_cleaned(12:21,13:end));
has_an_LPS_indel_index = sum(breseq_table_LPS==1)==1;
has_an_LPS_indel = Breseq_names(has_an_LPS_indel_index);
has_an_LPS_indel = strrep(has_an_LPS_indel,'SPv','SP-v');

% Define present and mutant LPS
P02_present_LPS= ~hasmutation_cleaned(LPS_stop_codon_mut,contains(sample_names_cleaned,"P02"));
P067_mutant_LPS = hasmutation_cleaned(LPS_muts,contains(sample_names_cleaned,["P06","P07"]));
P067_mutant_LPS = (sum(P067_mutant_LPS,1)>0)' | contains(P067_sample_names, has_an_LPS_indel);

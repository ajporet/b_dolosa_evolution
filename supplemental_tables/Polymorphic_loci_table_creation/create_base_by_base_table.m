clear all; close all; clc;

load analysis_m_outputs_structures.mat
load cds_sorted.mat
load gene_table_variables.mat

temp_SampleNames = ["AU0158"; SampleNames(samplestoplot)]';

replacement_A = ["P02","P06SP-v","P07SP-v","J-LT1","J-BLB"];
replacement_B = ["J-","Q-SP-","R-SP-","J-LT-1", "J-BL-"]; 

for ij = 1:length(replacement_A)
    temp_SampleNames = replace(temp_SampleNames, replacement_A(ij), replacement_B(ij));
end



call_matrix_num = [refnti(goodpos),maNT(goodpos, samplestoplot)];
call_matrix_ATCG = strings([size(call_matrix_num)]);

genome_pos = [[annotation_full.scaffold]',[annotation_full.pos]'];

table_metadata = final_table(gene_groupings,1:8);
table_labels = ["Chromosome","Position","Gene (AK34 annotation)", "Gene (BDAG annotation)", "Annotation from NZ_CP009793.1, NZ_CP009794.1, and NZ_CP009795.1)"...
    "Accession", "Coding (CDS) or intergenic (INT)", "Number of mutations", "Number of nonsynonymous mutations"... 
    "Number of synonymous mutations"];

call_matrix_ATCG(call_matrix_num==0) = "N";
call_matrix_ATCG(call_matrix_num==1) = "A";
call_matrix_ATCG(call_matrix_num==2) = "T";
call_matrix_ATCG(call_matrix_num==3) = "C";
call_matrix_ATCG(call_matrix_num==4) = "G";

final_table_2 = [genome_pos,table_metadata,call_matrix_ATCG];
final_table_2 = [[table_labels,temp_SampleNames];final_table_2];

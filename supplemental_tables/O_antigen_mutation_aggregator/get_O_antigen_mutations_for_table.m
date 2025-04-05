clear all; close all;
load JQR_b_dolosa_mutations.mat
load cds_sorted.mat

dictionary = readtable("Burkholderia_dolosa_AU0158_2968.xlsx"); 
dct_aks = dictionary{:,6};
AK_nums_corrected = string(erase([dictionary.LocusTag],'"')); BDAG_nums = string(erase([dictionary.GeneSynonyns],'"')); 

CDS_rearranged = squeeze(struct2cell([CDS{1},CDS{2},CDS{3}]))';
CDS_AK_nums = string(CDS_rearranged(:,12)); 
CDS_AK_old_nums = string(CDS_rearranged(:,11)); 
CDS_WP_nums=(CDS_rearranged(:,6)); CDS_WP_nums(cellfun(@isempty, CDS_WP_nums))= {"NA"}; CDS_WP_nums = string(CDS_WP_nums);

CDS_bdag_nums = [];
for i =1:numel(CDS_AK_nums)
     indice_dct = find(contains(dct_aks,CDS_AK_nums(i)));
     if length(indice_dct)==0
         CDS_bdag_nums = [CDS_bdag_nums,nan];
         continue
     else
         CDS_bdag_nums  = [CDS_bdag_nums,dictionary{indice_dct,12}];
     end
 end

breseq_table = readtable("Breseq_O_ant_results.xlsx",'ReadVariableNames', true);
breseq_table_cols = breseq_table.Properties.VariableNames;
breseq_samplenames = breseq_table_cols(12:end);
replacements = ["J_","Q_","R";
                "P02", "P06", "P07"];
for i=1:3
    breseq_samplenames = strrep(breseq_samplenames, replacements(1,i), replacements(2,i));
end
breseq_samplenames = strrep(breseq_samplenames,"_","-");


sample_names_cleaned = SampleNames(samplestoplot);
hasmutation_cleaned = hasmutation(:,samplestoplot);
P067_sample_names = sample_names_cleaned(contains(sample_names_cleaned,["P06","P07"]));

% Define LPS-affecting muts from manual inspection of data 
LPS_stop_codon_mut = 189;

LPS_muts= [50, 164, 188, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 265, 438,186, 187, 201,LPS_stop_codon_mut];

% Define present and mutant LPS
columns_to_care_about = [2,3,5,6,12,34];
annotation_full_cell = squeeze(struct2cell(annotation_full))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table_of_LPS_muts = {"chr","pos","name","wp","ak","mut","affected samples","hasp6p7"};
for i=1:numel(LPS_muts)
    lps_mut = LPS_muts(i);
    affected_samples = sample_names_cleaned(hasmutation_cleaned(lps_mut,:)>0);
    has_p67 = sum(contains(affected_samples,["P06","P07"]))>0;
    mut_anno = annotation_full_cell(lps_mut,columns_to_care_about);

    table_of_LPS_muts = [table_of_LPS_muts;[mut_anno,join(affected_samples,", "),has_p67]];
end

table_of_LPS_muts = [table_of_LPS_muts, cellstr(repmat("SNV",size(table_of_LPS_muts,1),1))];
%%%%%%%%%%%%%%%%%%%%%%%%%%% Now deal with indels 
% breseq_table = readtable("Breseq_O_ant_results.xlsx",'ReadVariableNames', true);
% breseq_table_cols = breseq_table.Properties.VariableNames;

% true_LPS_indels=  Breseq_Vals_cleaned([6,10:21],:);
table_of_LPS_indels = {"chr","pos","name","wp","ak","mut","affected samples","hasp6p7"};
table_rows = {};
for i=1:size(breseq_table,1)
    ak_num = split(breseq_table{i,5}{:}," ");
    ak_num = ak_num{1};
    ind_in_cds = find(strcmp(ak_num,CDS_AK_old_nums));

    bdag_num = CDS_bdag_nums(ind_in_cds);
    akrs_num  = CDS_AK_nums(ind_in_cds);
    wp_num = CDS_WP_nums(ind_in_cds);


    table_col = table2array(breseq_table(i,12:end));
    in_samples = breseq_samplenames(table_col==1);
    has_p67 = sum(contains(in_samples,["P06","P07"]))>0;

    if length(wp_num)==0
        wp_num = "NA";
    end
    if length(akrs_num)==0
        akrs_num = "NA";
    end
    new_row = cellstr({'1',string(breseq_table{i,2}),breseq_table{i,6}{:},string(wp_num),string(akrs_num),breseq_table{i,3}{:},join(in_samples,", "),string(has_p67)});
    table_rows = [table_rows;new_row];
end
table_rows = [table_rows, cellstr(repmat("INDEL",size(table_rows,1),1))];


final_table = [table_rows;table_of_LPS_muts(2:end,:)];
bdag_col = [];
for i=1:size(final_table,1)
    tab_cds_loc = find(strcmp(CDS_AK_nums,final_table(i,5)));
    if isempty(tab_cds_loc)
        bdag_col = [bdag_col,nan];
    else
         bdag_col = [bdag_col,CDS_bdag_nums(tab_cds_loc)];
    end

end
final_table = [final_table,bdag_col'];

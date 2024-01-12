clear all; close all; clc;

%% Create a supp. table that shows the number of mutations per gene

load analysis_m_outputs_structures.mat
load cds_sorted.mat

%% Subject selection

% Turn the main mutation describing structure annotation_full into a
% workable cell array
temp_SampleNames = SampleNames(samplestoplot);
temp_hasmutation = hasmutation(:,samplestoplot);

has_QR_mut = sum(temp_hasmutation(:,contains(temp_SampleNames, ["P06", "P07"])),2)>0;
% Adds in a mutation that was mutated parallely twice, once in J and once
% unique to subject R

has_J_mut =  sum(temp_hasmutation(:,contains(temp_SampleNames, ["P02"])),2)>0;
unique_mutations_QR = has_QR_mut & ~has_J_mut;

annotation_full_cell = squeeze(struct2cell(annotation_full))';

J_or_QR_toggle = 2; % 1 for J, 0 for Q_and_$
if J_or_QR_toggle==1 % To create this spreadsheet for J, uncomment the below line
    annotation_full_cell = annotation_full_cell(has_J_mut,:);
    mut_indexer = has_J_mut;
elseif J_or_QR_toggle==0 % To create this spreadsheet for Q and R, uncomment the below line 
    annotation_full_cell = annotation_full_cell(unique_mutations_QR,:);
    mut_indexer = unique_mutations_QR;
elseif J_or_QR_toggle==2 % Make a file for all mutations
    annotation_full_cell = annotation_full_cell(:,:);
    mut_indexer = ones(size(unique_mutations_QR));
end

%

%% Create table

% In order to retrieve BDAG annotation numbers, we use this gene annotation
% spreadsheet downloaded from https://www.burkholderia.com/strain/download?
% %c1=organism&v1=Burkholderia+dolosa+AU0158+%28complete+genome%29&v2=&c2=assemblyLevel
dictionary = readtable("Burkholderia_dolosa_AU0158_2968.xlsx"); 
AK_nums_corrected = string(erase([dictionary.LocusTag],'"')); BDAG_nums = string(erase([dictionary.GeneSynonyns],'"')); 

% We also use cds_sorted.mat, a matrix of gene annotations created from the
% gb files NZ_CP009793.1, NZ_CP009794.1, and NZ_CP009795.1. We use AK
% numbers, WP numbers, and annotations from here
CDS_rearranged = squeeze(struct2cell([CDS{1},CDS{2},CDS{3}]))';
CDS_AK_nums = string(CDS_rearranged(:,12)); CDS_WP_nums=(CDS_rearranged(:,6)); CDS_WP_nums(cellfun(@isempty, CDS_WP_nums))= {"NA"}; CDS_WP_nums = string(CDS_WP_nums);
gene_annotations=(CDS_rearranged(:,3)); gene_annotations(cellfun(@isempty, gene_annotations))= {"NA"}; doubled_gene_annos = cellfun(@(x) size(x,1)>1, gene_annotations);
gene_annotations(doubled_gene_annos) = cellfun(@(x) {erase(strjoin(string(x)," "),"  ")},gene_annotations(doubled_gene_annos));
gene_annotations = string(gene_annotations);

% We group genes here 
gene_groupings = findgroups((cell2mat(annotation_full_cell(:,1)) + 1000000.*cell2mat(annotation_full_cell(:,2))));
[unique_grouping_numbers,gene_group_index] = unique(gene_groupings, 'stable');
single_gene_entries = zeros(size(gene_groupings)); single_gene_entries(gene_group_index)=1;

% Get number of syn and nonsyn mutations, as well as intergenic mutations
number_syn =  arrayfun(@(x) sum(cell2mat(annotation_full_cell(gene_groupings==x,33))=='S'), unique_grouping_numbers)';
number_nonsyn =  arrayfun(@(x) sum(cell2mat(annotation_full_cell(gene_groupings==x,33))=='N'), unique_grouping_numbers)';
number_intergen =  arrayfun(@(x) sum(isnan(cell2mat(annotation_full_cell(gene_groupings==x,19)))), unique_grouping_numbers)';
nonsyn = number_nonsyn; nonsyn(number_intergen>0) = NaN; % Remove 0's from intergenic genes
syn = number_syn; syn(number_intergen>0) = NaN; % Remove 0's from intergenic genes

% Get number of mutations total per gene
number_muts =  arrayfun(@(x) length(annotation_full_cell(gene_groupings==x,33)), unique_grouping_numbers)';

% Get nonsyn mutation identities
mutations =  arrayfun(@(x) strjoin(string([annotation_full_cell{gene_groupings==x,34}]),"  "), unique_grouping_numbers, 'UniformOutput', false)';
mutations(~cellfun(@isempty, mutations)) = cellfun(@(x) strjoin(x, "   "), mutations(~cellfun(@isempty, mutations)),'UniformOutput',false);
mutations(cellfun(@isempty, mutations)) = {""};

% Get all protein annotations for basic, non-intergenic genes (i.e AK
% number, WP number, and description) 
protein_annotations = annotation_full_cell(single_gene_entries==1,[5, 6, 12]);
protein_annotations((cellfun(@isempty,protein_annotations))) = {""};

% Translate basic, non-intergenic AK numbers into BDAG numbers
BDAG_numbers = cellfun(@(x) BDAG_nums(AK_nums_corrected==x) , protein_annotations(:,3), 'UniformOutput' , false);
BDAG_numbers((cellfun(@isempty,BDAG_numbers))) = {"NA"};
BDAG_numbers(cellfun(@(x) strcmpi(x,""), BDAG_numbers)) = {"NA"};

% Get intergenic gene description annotations  
intergenic_annotations = annotation_full_cell(single_gene_entries==1,[21,25,23,27]);
intergenic_annotations((cellfun(@isempty,intergenic_annotations))) = {""};
intergenic_desc = cellfun(@(x) string(gene_annotations(CDS_AK_nums==x)) , intergenic_annotations(:,1:2), 'UniformOutput' , false);
intergenic_desc((cellfun(@isempty,intergenic_desc))) = {""};

% Translate intergenic AK numbers into BDAG numbers and replace spaces in
% BDAG_numbers variable 
int_BDAG_numbers = cellfun(@(x) BDAG_nums(AK_nums_corrected==x) , intergenic_annotations(:,1:2), 'UniformOutput' , false);
int_BDAG_numbers(cellfun(@isempty, int_BDAG_numbers)) = {"NA"};
int_BDAG_numbers(cellfun(@(x) strcmpi(x,""), int_BDAG_numbers)) = {"NA"};
BDAG_numbers(number_intergen>0) = cellstr(int_BDAG_numbers(number_intergen>0,1) + "-" +  int_BDAG_numbers(number_intergen>0,2));

% Translate intergenic WP numbers into BDAG numbers
int_WP_numbers = cellfun(@(x) CDS_WP_nums(CDS_AK_nums==x) , intergenic_annotations(:,1:2), 'UniformOutput' , false);
int_WP_numbers(cellfun(@isempty, int_WP_numbers)) = {"NA"};
int_WP_numbers(cellfun(@(x) strcmpi(x,""), int_WP_numbers)) = {"NA"};

% Get column of CDS or INT label
CDS_or_int = repelem("CDS",length(mutations),1); CDS_or_int(number_intergen>0) = "INT";

% Replace empty regions caused by intergenic genes in the
% protein_annotations variable with correct descriptions
protein_annotations(number_intergen>0,2) = cellstr(int_WP_numbers(number_intergen>0,1) + "-" +  int_WP_numbers(number_intergen>0,2));
protein_annotations(number_intergen>0,3) = cellstr(intergenic_annotations(number_intergen>0,1) + "-" + intergenic_annotations(number_intergen>0,2));
protein_annotations(number_intergen>0,1) = cellstr(intergenic_desc(number_intergen>0,1) + "---" + intergenic_desc(number_intergen>0,2));

% Replace any blanks remaining in protein_annotations (because a gene
% translation was not present) with NAs
protein_annotations(cellfun(@(x) x=="",protein_annotations)) = {'NA'};

% Create a final table for pasting into excel. This version is a little
% messy with extra spaces, but can easily be cleaned up manually for a
% tight and well-organized excel sheet
final_table = [protein_annotations(:,3),BDAG_numbers,protein_annotations(:,1:2), cellstr(CDS_or_int), num2cell(number_muts)', num2cell(nonsyn)', num2cell(syn)', mutations'];

if J_or_QR_toggle==2
    save('gene_table_variables.mat','final_table','gene_groupings')
end

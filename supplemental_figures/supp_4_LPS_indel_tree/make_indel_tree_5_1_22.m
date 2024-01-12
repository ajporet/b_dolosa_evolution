clear all; close all;

%% Make tree with all indels 
fid = fopen('wholetree_infile.txt');
read_file = textscan(fid,'%s','delimiter','\n'); 
read_file = read_file{:}; 
header_line = read_file(1); read_file(1) = []; %removes first line
breaking_input = cellfun(@(x) split(x," "),read_file,'UniformOutput',false);
input_file_samplenames_temp = string(cellfun(@(x) x(1), breaking_input));
input_file_fasta = cellfun(@(x) x(end), breaking_input);

translating_input = string(readcell("T^toSampleName1_22_2022.xlsx"));
input_file_samplenames = cellfun(@(x) translating_input(strcmp(translating_input(:,1),x),2), input_file_samplenames_temp, 'UniformOutput', false);

indel_file = readcell("matlab_input_breseq_results");
sample_names = string(indel_file(1,14:end));

indel_matrix = string(indel_file(2:end, 14:end));
indel_matrix((indel_matrix == "?") | (indel_matrix =="Δ")) = "N"; indel_matrix(indel_matrix == "100%") = "A";
indel_matrix(indel_matrix == "1") = "A"; indel_matrix(indel_matrix == "0") = "T";

extra_phylip_indel_input = join(indel_matrix',"");

match_sample_index = cellfun(@(x) find(strcmp(sample_names,x)), input_file_samplenames, 'UniformOutput', false);

indel_addition = strings(size(1,length(input_file_samplenames)));
indel_addition(1) = join(repelem("T",87),"");

for ii=2:length(input_file_samplenames)
    if isempty(match_sample_index{ii})
        indel_addition(ii) = join(repelem("N",87),"");
    else
        indel_addition(ii) = extra_phylip_indel_input(match_sample_index{ii});
    end
end


final_phylip_input = input_file_samplenames + "	"+ input_file_fasta + indel_addition';
final_phylip_input = ["932	" + string(length(input_file_fasta{2}) + length(indel_addition{2}));final_phylip_input];

filePh = fopen('indel_tree_infile.txt','w');
fprintf(filePh,'%s\n',final_phylip_input{:});
fclose(filePh);

%% Make tree with only LPS indels 
indel_names = indel_file(2:end,7);
O_ant_indel_gene_names = ["O‑Antigen",'glycosyl transferase','GDP‑mannose'];
indel_matrix = string(indel_file(2:end, 14:end));
indel_matrix = indel_matrix(contains(indel_names,O_ant_indel_gene_names),:);

indel_matrix((indel_matrix == "?") | (indel_matrix =="Δ")) = "N"; indel_matrix(indel_matrix == "100%") = "A";
indel_matrix(indel_matrix == "1") = "A"; indel_matrix(indel_matrix == "0") = "T";

extra_phylip_indel_input = join(indel_matrix',"");

match_sample_index = cellfun(@(x) find(strcmp(sample_names,x)), input_file_samplenames, 'UniformOutput', false);

indel_addition = strings(size(1,length(input_file_samplenames)));
indel_addition(1) = join(repelem("T",87),"");

for ii=2:length(input_file_samplenames)
    if isempty(match_sample_index{ii})
        indel_addition(ii) = join(repelem("N",87),"");
    else
        indel_addition(ii) = extra_phylip_indel_input(match_sample_index{ii});
    end
end


final_phylip_input = input_file_samplenames + "	"+ input_file_fasta + indel_addition';
final_phylip_input = [["932	" + string(length(input_file_fasta{2}) + length(indel_addition{2}))];final_phylip_input];

filePh2 = fopen('indel_only_O_antigen_genes_tree_infile.txt','w');
fprintf(filePh2,'%s\n',final_phylip_input{:});
fclose(filePh2);

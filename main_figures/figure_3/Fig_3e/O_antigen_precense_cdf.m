clear; close all;

% Get isolates with LPS mutations
all_2011_mutations= readcell("Lieberman_et_al_2011_supp_table_2.xlsx");
sample_names = string(all_2011_mutations(1,[6,9:120]));
O_antigen_gene_names = ["Glycosyltransferase","Mannose-1-phosphate","O-antigen"];
O_antigen_gene_locations = find(contains(string(all_2011_mutations(2:end,4)), O_antigen_gene_names));
stop_codon_reversion_location = 232;

last_common_ancester_calls = string(all_2011_mutations(2:end,8));
sample_calls = string(all_2011_mutations(2:end,[6,9:120]));

mutation_matrix = (last_common_ancester_calls ~= sample_calls) & ~cellfun(@isempty, sample_calls);
no_mutation_matrix = (last_common_ancester_calls == sample_calls) & ~cellfun(@isempty, sample_calls);
has_reversion = mutation_matrix(stop_codon_reversion_location,:);

%has_O_antigen_indel = ["J-12-11","J-12-4","J-13-6a-$","J-13-6b","J-14-2","J-14-8"];
has_O_antigen_indel = ["J-12-11","J-13-6a-$","J-13-6b","J-14-2"];

extra_mutation = (has_reversion & sum(mutation_matrix(setdiff(O_antigen_gene_locations, stop_codon_reversion_location),:))>0) | (contains(sample_names, has_O_antigen_indel));
only_reversion = has_reversion & sum(no_mutation_matrix(setdiff(O_antigen_gene_locations, stop_codon_reversion_location),:))>=length(O_antigen_gene_locations)-1 & ~contains(sample_names, has_O_antigen_indel);
stop_codon = no_mutation_matrix(stop_codon_reversion_location,:);

%% Plot time since first sample

time_since_sampling = cellfun(@(x) split(x(3:end),"-"), sample_names, 'UniformOutput', false);
time_since_sampling = cellfun(@(x) str2num(cell2mat(regexprep(x(1),'[^0-9]',''))) + str2num(cell2mat(regexprep(x(2),'[^0-9]','')))/12, time_since_sampling(2:end));
time_since_sampling = [3.15, time_since_sampling]; % Account for AU0158 date being difficult to parse

patient_letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"];

figure('Renderer', 'painters', 'Position', [10 10 400*2 200*2])
for ii=1:length(patient_letters)
    patient_sample = contains(sample_names, patient_letters(ii));
    scatter(time_since_sampling(patient_sample & extra_mutation) - min(time_since_sampling(patient_sample)),ii.*ones(sum(patient_sample & extra_mutation),1),'r','filled','MarkerFaceAlpha',.5)
    hold on
    scatter(time_since_sampling(patient_sample & only_reversion) - min(time_since_sampling(patient_sample)),ii.*ones(sum(patient_sample & only_reversion),1),'b','filled','MarkerFaceAlpha',.5)
    hold on
    scatter(time_since_sampling(patient_sample & stop_codon) - min(time_since_sampling(patient_sample)),ii.*ones(sum(patient_sample & stop_codon),1),'k','filled','MarkerFaceAlpha',.5)
    hold on 
end

H=gca;
H.LineWidth=2;
yticks([1:length(patient_letters)])
set(gca,'FontSize',16)
set(gca,'box','off') 
set(gca,'TickDir','out');
xlabel("years from first sampled isolate in each patient")
ylabel('patient')
yticklabels(patient_letters)
xlim([-.1,8])

%% Make a cdf

patient_samp_adjust = zeros(size(time_since_sampling));
for ii=1:length(patient_letters)
    patient_sample = contains(sample_names, patient_letters(ii));
    patient_samp_adjust(patient_sample) = time_since_sampling(patient_sample) - min(time_since_sampling(patient_sample));
end

disc_data = discretize(patient_samp_adjust,-.01:1:7.99);
percent_LPS = arrayfun(@(x) mean(stop_codon(disc_data==x)), 1:8);
error_bar= zeros(2,8);

for ii=1:8
    [phat,pci] = binofit(sum(stop_codon(disc_data==ii)==1),sum(disc_data==ii));
    error_bar(2,ii) = phat - pci(1);
    error_bar(1,ii) = pci(2) - phat;
end

figure('Renderer', 'painters', 'Position', [10 10 300 300])
plot(1:8, percent_LPS,'-k', 'LineWidth', 2)
set(gca,'xtick',[])
H=gca;
H.LineWidth=3;
set(gca,'box','off') 
set(gca,'TickDir','out'); 
yticklabels([])
xlim([1,8])
yticks(0:.5:1)
xticks(1:8)
xticklabels([])


%% Load in data
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

% Seperate out array by sampling time
times = {'P06SP-v1', 'P07SP-v1', 'P06SP-v2','P07SP-v2','P06SP-v3','P07SP-v3'};
no_LPS_num_by_time = cellfun(@(x) sum(contains(P067_sample_names(P067_mutant_LPS),x)), times);
yes_LPS_num_by_time = cellfun(@(x) sum(contains(P067_sample_names(~P067_mutant_LPS),x)), times);
all_samp_by_time = cellfun(@(x) sum(contains(P067_sample_names,x)), times);

no_LPS_num_by_time = [no_LPS_num_by_time, sum(~P02_present_LPS)];
yes_LPS_num_by_time = [yes_LPS_num_by_time, sum(P02_present_LPS)];
all_samp_by_time = [all_samp_by_time, length(P02_present_LPS)];

no_LPS_num_by_time_percentage = no_LPS_num_by_time./all_samp_by_time;
yes_LPS_num_by_time_percentage = yes_LPS_num_by_time./all_samp_by_time;

errorbar_1= [];
errorbar_2= [];
for ii=1:max(size(no_LPS_num_by_time))
    [mean_val,CI] = binofit(no_LPS_num_by_time(ii),all_samp_by_time(ii),.05);
    errorbar_1 = [errorbar_1, abs(mean_val-CI(1))];
    errorbar_2 = [errorbar_2,abs(mean_val-CI(2))];
end


% Graph data
stacked_bar_data = [no_LPS_num_by_time_percentage(1:6), 0, no_LPS_num_by_time_percentage(7)]';
stacked_bar_data = [stacked_bar_data, 1-stacked_bar_data]; stacked_bar_data(7,2)=0;

figure('Renderer', 'painters', 'Position', [10 10 600 300])
yes_LPS_color = [143, 143, 143]./225;

ba = bar(stacked_bar_data,'stacked','FaceColor','flat','LineWidth',1,'EdgeColor','k');
ba(1).CData = yes_LPS_color;
ba(2).CData = [1,1,1];
hold on

errorbar([1,2,3,4,5,6,8],[no_LPS_num_by_time_percentage]',errorbar_1,errorbar_2,'.k','LineWidth',1.5,'CapSize',0)

% Graph settings
set(gca,'xtick',[])
H=gca;
H.LineWidth=3;
yticks(0:.5:1)
set(gca,'box','off') 
ylim([0,1])
set(gca,'FontSize',20)
set(gca,'TickDir','out'); 
yticklabels([])

[p_final_P06,rejectsnull_P06]  = binomial_proportion(yes_LPS_num_by_time(5),all_samp_by_time(5),yes_LPS_num_by_time(1),all_samp_by_time(1),.05);
[p_final_P07,rejectsnull_P07]  = binomial_proportion(yes_LPS_num_by_time(6),all_samp_by_time(6),yes_LPS_num_by_time(2),all_samp_by_time(2),.05);


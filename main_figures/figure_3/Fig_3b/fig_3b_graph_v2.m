clear all; close all;
load analysis_m_outputs_dnds_struct.mat

check_hasmutation_was_saved_correctly = sum(sum(hasmutation(:, samplestoplot),2)==0)==0;

SampleNames_corrected = SampleNames(samplestoplot);
hasmutation_corrected = hasmutation(:,samplestoplot);

Q_R_samples = contains(SampleNames_corrected, ["P06","P07"]);
J_samples = contains(SampleNames_corrected, ["P02"]);

all_LPS_muts= [50, 164, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 265, 438,186, 187, 201];

J_mutations = sum(hasmutation_corrected(:,J_samples),2)>0;
QR_mutations = sum(hasmutation_corrected(:,Q_R_samples),2)>0 & ~J_mutations;

fraction_J_LPS_mutations = sum(J_mutations(all_LPS_muts))/sum(J_mutations);
fraction_QR_LPS_mutations = sum(QR_mutations(all_LPS_muts))/sum(QR_mutations);

% Null model
% Locations of O-antigen genes from KEGG annotations, taken and copy and
% pasted from spreadsheets "Functional_analysis_map_template_chr_X.xlsx"

chr_2_locs=	[[1159440; 1161024], [1508907; 1509795], [1519546; 1520593], [1520579; 1521551], [1668859; 1670272], [1672140; 1673667]]';
chr_1_locs=	[[132303; 133665], [1851259; 1852144], [2287118; 2288519], [2439498; 2440521], [2456605; 2457523], [2457500; 2458541], [2463950; 2465396], [2465618; 2466509], [2466505; 2467057], [2467041; 2467935], [2467946; 2469008], [3137592; 3138753], [642228; 643260]]';		
chr_3_locs = [[390152; 391118], [476546; 477761]]';

custom_defined_O_ant_casette = [2431102,2432475;2432513,2433397;2433423,2434418;2434591,2435985;2436229,2437485;2437482,2438306;2438347,2439483;2439499,2440521;2440549,2441655;2441741,2443627;2443640,2444656;2444662,2445630;2446065,2447960;2447994,2449457;2449461,2450777;2450774,2451529;2451574,2452512;2452516,2453799;2453796,2454230;2454248,2455465;2455462,2456538;2456606,2457523;2457501,2458541;2458554,2462279;2462280,2463608;2463608,2463847;2463951,2465396;2465619,2466509;2466506,2467057;2467042,2467935;2467947,2469008;2469553,2470044;2470066,2470857;2471340,2472377;2472425,2473729;2473793,2474542;2474553,2475335];

chr_1_locs = [chr_1_locs; custom_defined_O_ant_casette];

[~,idx_1] = sort(chr_1_locs(:,1)); chr_1_locs = chr_1_locs(idx_1,:);
[~,idx_2] = sort(chr_2_locs(:,1)); chr_2_locs = chr_2_locs(idx_2,:);

merged_chr_1_locs = merge_intervals(chr_1_locs);
merged_chr_2_locs = merge_intervals(chr_2_locs);
merged_chr_3_locs = merge_intervals(chr_3_locs);

genome_size = 6409095;

chr_1_sum = sum(merged_chr_1_locs(:,2)-merged_chr_1_locs(:,1));
chr_2_sum = sum(merged_chr_2_locs(:,2)-merged_chr_2_locs(:,1));
chr_3_sum = sum(merged_chr_3_locs(:,2)-merged_chr_3_locs(:,1));

total_sum = chr_1_sum + chr_2_sum + chr_3_sum;
expected_number_of_muts = total_sum./genome_size;

%{
num_trials = 100000;
exp_fract_mut = zeros(1,num_trials);
num_muts_to_disperse = 1000;
for ii=1:num_trials
    exp_fract_mut(ii) = sum(randi(genome_size,num_muts_to_disperse,1)<=total_sum)/num_muts_to_disperse;
end
expected_number_of_muts = mean(exp_fract_mut);
%} 

% Creating error bars
errorbar_1= [];
errorbar_2= [];

[mean_val,CI] = binofit(sum(J_mutations(all_LPS_muts)),sum(J_mutations),.05);
errorbar_1 = [errorbar_1, abs(mean_val-CI(1))];
errorbar_2 = [errorbar_2,abs(mean_val-CI(2))];

[mean_val,CI] = binofit(sum(QR_mutations(all_LPS_muts)),sum(QR_mutations),.05);
errorbar_1 = [errorbar_1, abs(mean_val-CI(1))];
errorbar_2 = [errorbar_2,abs(mean_val-CI(2))];

[p_final_geno_1,~]  = binomial_proportion(sum(J_mutations(all_LPS_muts)),sum(J_mutations),sum(QR_mutations(all_LPS_muts)),sum(QR_mutations),.05);


num_J_LPS_muts = sum(J_mutations(all_LPS_muts));
num_QR_LPS_muts = sum(QR_mutations(all_LPS_muts));

binopdf(num_QR_LPS_muts,sum(QR_mutations),expected_number_of_muts) 
binopdf(num_J_LPS_muts,sum(J_mutations),expected_number_of_muts) 

p_QR_w_null_prob = sum(arrayfun(@(x) binopdf(x,sum(QR_mutations),expected_number_of_muts),num_QR_LPS_muts:sum(QR_mutations))); 
p_J_w_null_prob = sum(arrayfun(@(x) binopdf(x,sum(J_mutations),expected_number_of_muts),num_J_LPS_muts:sum(J_mutations)));

% Plotting
bar_to_plot = [fraction_J_LPS_mutations,fraction_QR_LPS_mutations,expected_number_of_muts];
figure('Renderer', 'painters', 'Position', [10 10 250 300])
ba = barh([3,2,1],bar_to_plot,'stacked','FaceColor','flat','LineWidth',1,'EdgeColor','none');
ba.CData = [126, 126, 126]./225;
hold on
errorbar(bar_to_plot(1:2),[3,2],[0 0],[0 0],errorbar_2,errorbar_1,'.k','LineWidth',1.5,'CapSize',0)
xlim([0,.16])
ylim([0,4])
set(gca,'ytick',[])
H=gca;
H.LineWidth=3;
set(gca,'FontSize',20)
set(gca,'box','off') 
xticks(linspace(0,.16,5))
set(gca,'TickDir','out'); 
xticklabels([])
xticks([0:.08:.16])

function merged_chr_2_locs = merge_intervals(chr_2_locs)

    merged_chr_2_locs = [];
    stopper = 0;
    ii=1;
    while stopper==0
        check_merge = chr_2_locs(ii,2)> chr_2_locs(:,1);
        check_merge(1:ii) = 0;
        if sum(check_merge)>0
            check_merge_indicies = find(check_merge);
            merged_chr_2_locs = [merged_chr_2_locs;[chr_2_locs(ii,1),chr_2_locs(max(check_merge_indicies),2)]];
            ii = max(check_merge_indicies);
        else
            merged_chr_2_locs = [merged_chr_2_locs;chr_2_locs(ii,:)];
            ii = ii +1;
        end
    
        if ii>length(chr_2_locs)
            stopper = 1;
        end
    end
end
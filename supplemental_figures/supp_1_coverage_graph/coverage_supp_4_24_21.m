%% Create supplemental coverage graphic
clear all; 
close all; 
load candidate_mutation_table.mat
load coveragematrix.mat

mean_cov_per_genome = mean(all_coverage_per_bp,2);

figure('Renderer', 'painters', 'Position', [10 10 300*2 150*2])
histogram(mean_cov_per_genome, 60, 'FaceColor', [126,126,126]./255)
H=gca;
H.LineWidth=2;
set(gca,'FontSize',16)
set(gca,'box','off') 
set(gca,'TickDir','out');
xlabel("coverage")
ylabel("number of isolates")

%% Get all raw data in a useable form 
clear all; close all;
file_name = "murine_model_data_machine_readable.xlsx";

sheet_names = {'Pair 1 101-317','Pair 2 107-211','Pair 3 207-311'};

dot_size = 20;

figure('Renderer', 'painters', 'Position', [10 10 600*2 150*2])
bar_locs = [[1,2,3,4];[6,7,8,9]+1;[11,12,13,14]+2];
symbols = ["o","*","o","*"];
day_1_color = [211, 213, 212]./255;
day_7_color = [126, 126, 126]./255;


for ii=1:length(sheet_names)
    counts= readtable(file_name,'Sheet',sheet_names{ii}, NumHeaderLines=1);

    %%% first trial with coutner
    inoculum = table2array(counts(2:4,1));
    counter_inoculum = table2array(counts(2:4,3));

    lung_d1_d7 = table2array([counts(9:13,1),counts(9:13,5)]);
    spleen_d1_d7 = table2array([counts(9:13,1+8),counts(9:13,5+8)]);
    
    counter_lung_d1_d7 = table2array([counts(9:13,3),counts(9:13,7)]);
    counter_spleen_d1_d7 = table2array([counts(9:13,3+8),counts(9:13,7+8)]);

    %%% second trial
    inoculum_2 = table2array(counts(2+17:4+17,1));
    
    lung_d1_d7_2 = table2array([counts(9+17:13+17,1),counts(9+17:13+17,5)]);
    spleen_d1_d7_2 = table2array([counts(9+17:13+17,1+8),counts(9+17:13+17,5+8)]);
    
    d1_means = [[lung_d1_d7(:,1);lung_d1_d7_2(:,1)],[spleen_d1_d7(:,1);spleen_d1_d7_2(:,1)]];
    d7_means = [[lung_d1_d7(:,2);lung_d1_d7_2(:,2)],[spleen_d1_d7(:,2);spleen_d1_d7_2(:,2)]];
    bar_means = [nanmean(d1_means),nanmean(d7_means)];

    bargraph = bar(bar_locs(ii,:),bar_means,'FaceColor','flat','LineWidth',1);
    bargraph.EdgeColor = 'none';


    bargraph.CData = [day_1_color;day_1_color;day_7_color;day_7_color];


    scatter_points = [d1_means,d7_means];
    hold on
    scatter(repelem(bar_locs(ii,:),5,1),scatter_points(1:5,:),'k','filled','o');
    hold on
    scatter(repelem(bar_locs(ii,:),5,1),scatter_points(6:10,:),'k','filled','s');
    hold on

    con_ints = arrayfun(@(x) confidence_int(scatter_points(:,x)), [1:4], 'UniformOutput', false);
    upper_bound = arrayfun(@(x) con_ints{x}(1), 1:4) - bar_means;
    lower_bound = bar_means - arrayfun(@(x) con_ints{x}(2), 1:4);

    er_innoc = errorbar(bar_locs(ii,:),bar_means,upper_bound,lower_bound); 
    er_innoc.Color = [0 0 0];er_innoc.LineStyle = 'none'; er_innoc.LineWidth = 1;

    disp(sheet_names{ii})
    disp("Day 1 Wilcoxen ranksum: " + ranksum(scatter_points(:,1),scatter_points(:,2)))
    disp("Day 7 Wilcoxen ranksum: " + ranksum(scatter_points(:,3),scatter_points(:,4)))

    hold on
end
    

ylim([0,1]); box off;set(gca,'TickDir','out');yticks([0:.25:1])
set(gca,'xtick',[])
H=gca;
H.LineWidth=2;
set(gca,'FontSize',15); set(gca,'Yticklabel',[]) 


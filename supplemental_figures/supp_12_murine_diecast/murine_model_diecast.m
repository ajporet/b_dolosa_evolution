%% Get all raw data in a useable form 
clear all; close all;

file_name = "machine_readable_mouse_data_9_29_21.xlsx";

lung_counts_no_triton = readtable(file_name, 'Sheet', 'lung_no_triton');
spleen_counts_no_triton = readtable(file_name, 'Sheet', 'spleen');
inoc_counts = readtable(file_name, 'Sheet', 'inoculum_counts_spread');
g2_inoc = inoc_counts.ColonyCount_Blue_NoTriton ./ (inoc_counts.ColonyCount_Blue_NoTriton + inoc_counts.ColonyCount_White_NoTriton);
g1_inoc = inoc_counts.ColonyCount_White_NoTriton ./ (inoc_counts.ColonyCount_Blue_NoTriton + inoc_counts.ColonyCount_White_NoTriton);
g2_inoc = mean(g2_inoc(1:3));
g1_inoc = mean(g1_inoc(4:6));


blue_cols = (lung_counts_no_triton.ColonyCount_Blue1 + lung_counts_no_triton.ColonyCount_Blue2)./2;
white_cols = (lung_counts_no_triton.ColonyCount_White1 + lung_counts_no_triton.ColonyCount_White2)./2;

blue_cols_sp = spleen_counts_no_triton.ColonyCount_Blue_NoTriton;
white_cols_sp = spleen_counts_no_triton.ColonyCount_White_NoTriton;

day_7_mouse_lung_g_1 = white_cols(11:15)./(white_cols(11:15)+blue_cols(11:15));
day_7_mouse_lung_g_2 = blue_cols(16:19)./(blue_cols(16:19)+white_cols(16:19));

day_7_mouse_spleen_g_1 = white_cols_sp(11:15)./(white_cols_sp(11:15)+blue_cols_sp(11:15));
day_7_mouse_spleen_g_2 = blue_cols_sp(16:19)./(white_cols_sp(16:19) + blue_cols_sp(16:19));

data_in_order_g1 = [day_7_mouse_lung_g_1,day_7_mouse_spleen_g_1];
data_in_order_g2 = [day_7_mouse_lung_g_2,day_7_mouse_spleen_g_2];

graph_means = [mean(data_in_order_g1),mean(data_in_order_g2)];

figure('Renderer', 'painters', 'Position', [10 10 200*2 200*2])
bar_locs = [1,2,4,5];
first_color = [75,75,75]./255;
second_color = [150,150,150]./255;


bargraph = bar(bar_locs,graph_means,'FaceColor','flat','LineWidth',1);
bargraph.EdgeColor = 'none';
bargraph.CData = [first_color;first_color;second_color;second_color];

hold on
scatter(repelem(bar_locs(1:2),5,1),data_in_order_g1,'k','filled','o');
hold on
scatter(repelem(bar_locs(3:4),4,1),data_in_order_g2,'k','filled','o');
hold on
%plot([0.5,2.5],[g1_inoc,g1_inoc],'k:', 'LineWidth',3);
%hold on
%plot([3.5,5.5],[g2_inoc,g2_inoc],'k:', 'LineWidth',3);
%hold on


con_ints = arrayfun(@(x) confidence_int(data_in_order_g1(:,x)), [1:2], 'UniformOutput', false);
con_ints_2 = arrayfun(@(x) confidence_int(data_in_order_g2(:,x)), [1:2], 'UniformOutput', false);
con_ints = [con_ints,con_ints_2];
upper_bound = arrayfun(@(x) con_ints{x}(1), 1:4) - graph_means;
lower_bound = graph_means - arrayfun(@(x) con_ints{x}(2), 1:4);

er_innoc = errorbar(bar_locs,graph_means,upper_bound,lower_bound); 
er_innoc.Color = [0 0 0];er_innoc.LineStyle = 'none'; er_innoc.LineWidth = 1;
ylim([0,1]); box off;set(gca,'TickDir','out');yticks([0:.25:1])
set(gca,'xtick',[])
H=gca;
H.LineWidth=2;
set(gca,'FontSize',15); set(gca,'Yticklabel',[]) 


group_1_p_val = ranksum(data_in_order_g1(:,1),data_in_order_g1(:,2));
group_2_p_val = ranksum(data_in_order_g2(:,1),data_in_order_g2(:,2));

[~,group_1_p_val_ttest] = ttest(data_in_order_g1(:,1),data_in_order_g1(:,2));
[~,group_2_p_val_ttest] = ttest(data_in_order_g2(:,1),data_in_order_g2(:,2));



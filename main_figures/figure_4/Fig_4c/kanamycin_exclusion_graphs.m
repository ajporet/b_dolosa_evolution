%% Making all raw graphs
clear all; close all; 

file_names = ["20230725_LPS dolosa THP invasion assay.xlsx", "20230801_LPS dolosa THP invasion assay.xlsx", "20230808_LPS dolosa THP invasion assay.xlsx"];

toggle_show_axes_labels = 1;

for item = file_names
    rawfile31 = readtable(item, 'Sheet',"Sheet1");
    
    percent_of_total = rawfile31.Var9;
    number_breaks = find(~isnan(percent_of_total));
    break_point = number_breaks(abs(number_breaks(1:length(number_breaks)-1) - number_breaks(2:end))>2);
    
    total_calcs = percent_of_total(number_breaks(1):break_point);
    total_calcs(5:5:end,:) = [];
    
    sample_nums = [101, 107, 207, 211, 311, 317];
    reorder = [1,6,2,4,3,5];
    % Q-101 : R-317     R-107: R-211     Q-207:R-311
    
    percent_of_total = reshape(total_calcs,[4,6]);
    percent_of_total_mean = nanmean(percent_of_total);
    locations = [1,2,4,5,7,8];
    dot_size = 20;
    
    first_color = [75,75,75]./255+.2;
    second_color = [125,125,125]./255+.15;
    third_color = [190,190,190]./255+.1;
    
    figure('Renderer', 'painters', 'Position', [10 10 150*2 150*2])
    b = bar(locations,percent_of_total_mean(reorder)); %,'facecolor',[.8,.8,.8])
    b.FaceColor = 'flat';
    b.CData([1,2],:) = repelem(first_color,2,1);
    b.CData([3,4],:) = repelem(second_color,2,1);
    b.CData([5,6],:) = repelem(third_color,2,1);
    b.EdgeColor = 'none';
    
    hold on
    scatter(locations,percent_of_total(:,reorder),dot_size,'filled','k') 
    
    box off;set(gca,'TickDir','out');
    set(gca,'xtick',[])
    H=gca;
    H.LineWidth=2;
    set(gca,'FontSize',15); 
    set(gca,'Xticklabel',[]) 
    %set(gca, 'YScale', 'log')
    rounded_y= nanmax(25*ceil(percent_of_total/25),[],'all');
    ylim([0,rounded_y])
    yticks(0:25:rounded_y);
    set(gca,'FontSize',16)
    set(gca,'box','off') 
    set(gca,'TickDir','out');
    

    if toggle_show_axes_labels==0
        yticklabels([])
        
    end

    disp(item)
    disp("Pair 1: " + ranksum(percent_of_total(:,reorder(1)),percent_of_total(:,reorder(2))));
    disp("Pair 2: " + ranksum(percent_of_total(:,reorder(3)),percent_of_total(:,reorder(4))));
    disp("Pair 3: " + ranksum(percent_of_total(:,reorder(5)),percent_of_total(:,reorder(6))));
    disp("Fold change 1: " + percent_of_total_mean(reorder(2))/percent_of_total_mean(reorder(1)));
    disp("Fold change 2: " + percent_of_total_mean(reorder(4))/percent_of_total_mean(reorder(3)));
    disp("Fold change 3: " + percent_of_total_mean(reorder(6))/percent_of_total_mean(reorder(5)));


    
end
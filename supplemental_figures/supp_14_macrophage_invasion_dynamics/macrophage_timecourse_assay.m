%% Making all raw graphs
clear all; close all; 
uptake_6 = readtable("20240806 311 207 uptake assay.xlsx",'ReadVariableNames', false);
uptake_13 = readtable("20240813 311 207 uptake assay.xlsx",'ReadVariableNames', false);

species = ["311 500","207 500","311 1000","207 1000"];
col_data = {};
for ta = {uptake_6,uptake_13}
    table_2_use = ta{:};
    first_col = table2array(table_2_use(:,1));
    cfu_col = table2array(table_2_use(:,6));
    num_row = size(first_col,1);
    time_demark = [find(contains(first_col,["min","hr"]));num_row];
    spec_loc = {};

    data_formatted = zeros(3,5,4);
    for s = species
        spec_demark = contains(first_col,s);
        spec_loc = [spec_loc,spec_demark];
    end

    
    for s =1:4
        for i =1:5
            find_points = spec_loc{s}.*linspace(1,num_row,num_row)';
            point_loc = cfu_col((find_points>=time_demark(i)) & (find_points<=time_demark(i+1)));
            data_formatted(:,i,s) = point_loc;
        end
    end
    col_data = [col_data,data_formatted];


end
times = [15,30,45,60,120];

for ii = 1:2

    data = col_data{ii};
    figure('Renderer', 'painters', 'Position', [10 10 200*2 150*2])
    lacz_mut_color = "#E48D70";lac_wt_color = "#87C7C3";
    for i=3:4
        % We only show the 1:1000 dilutions because these matched the previous
        % macrophage experiments and there were basically no difference between
        % the two (1:500 and 1:1000)
        if contains(species(i),"207")
            col = lac_wt_color;
        else
            col = lacz_mut_color;
        end
        
    
    
        x = repmat(times, 3, 1);
        y= data(:,:,i);
        scatter(reshape(x,1,15),reshape(y,1,15),'SizeData',18,'MarkerFaceColor',col,'MarkerEdgeColor',col)
        hold on
        plot(mean(x), mean(y),'Color',col,'LineWidth',2)
        hold on
        SEM = std(y)/sqrt(size(y,1));                             % Standard Error Of The Mean
        errorbar(times,mean(y),SEM,'Color',col,'LineWidth',2)
        hold on
    end
    set(gca, 'YScale', 'log');
    set(gca, 'LineWidth',2)
    box off;set(gca,'TickDir','out');
    ylim([10^3, 10^6])
    xlim([5,130])
    xticks(times)
    set(gca,'XTickLabel',times,'fontsize',14)
    disp("Table" + string(ii))
    [p,h] = ttest2(data(:,1,4),data(:,1,3));
    disp(h)
    disp("Mean dif")
    mean1 = mean(data(:,:,3));
    mean2 = mean(data(:,:,4));
    disp(log2(mean1(5)/mean1(1)))
    disp(log2(mean2(5)/mean2(1)))
    disp(mean1(1)/mean2(1))
end


    % 
    % stop
    % 
    %     
    %     figure('Renderer', 'painters', 'Position', [10 10 200*2 150*2])
    %     b = bar(1:5, mean_total_ratio(plotting_order)); %,'facecolor',[.8,.8,.8])
    %     b.FaceColor = 'flat';
    %     b.CData(:,:) = repelem([115,115,115]./255,5,1);
    %     b.EdgeColor = 'none';
    %     
    %     
    %     box off;set(gca,'TickDir','out');
    %     set(gca,'xtick',[])
    %     H=gca;
    %     H.LineWidth=2;
    %     set(gca,'FontSize',15); set(gca,'Xticklabel',[]) 
    %     set(gca,'FontSize',16)
    %     set(gca,'box','off') 
    %     set(gca,'TickDir','out');
    %     y_max = ceil(nanmax(mean_total_ratio));
    %     y_max = round(y_max*1.5);
    %     ylim([0,y_max])
    %     yticks(0:round(y_max/3):y_max);
    % 
    %     LPS_pres_samps = raw_total_ratio(:,plotting_order(1:2));
    %     LPS_abs_samps = raw_total_ratio(:,plotting_order(3:5));
    %     disp("Aggregate wilrank score: " + ranksum(LPS_pres_samps(:),LPS_abs_samps(:)));
    % 
    % 
    % end
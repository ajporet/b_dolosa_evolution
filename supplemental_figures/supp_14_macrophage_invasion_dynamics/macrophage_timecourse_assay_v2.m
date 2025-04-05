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
nums_to_use = [1,4;3,4];
innoc_0806  =[5.95E+05	1.19E+06
7.30E+05	1.46E+06];
innoc_0813 = [3.45E+05	6.90E+05
1.50E+05	3.00E+05];

innoc_0806 = [345000	300000];
innoc_0813 = [5.95E+05  7.30E+05];

innocs = {innoc_0806,innoc_0813};
figure('Renderer', 'painters', 'Position', [10 10 200*2 150*2])

trial_marker = ["diamond","o"];
actual_data = {[],[]};
for ii = 1:2
    trial_m = trial_marker(ii);
    data = col_data{ii};
    
    lacz_mut_color = "#E48D70";lac_wt_color = "#87C7C3";

    dataset_inoc = innocs{ii};
    
    count = 0;
    for i=nums_to_use(ii,:)
        count =count + 1;
        true_inoc = dataset_inoc(count);

        % We only show the 1:1000 dilutions because these matched the previous
        % macrophage experiments and there were basically no difference between
        % the two (1:500 and 1:1000)
        if contains(species(i),"207")
            col = lac_wt_color;
        else
            col = lacz_mut_color;
        end
    
        x = repmat(times, 3, 1);
        y= data(:,:,i)./true_inoc;

        actual_data(count) = {[actual_data{count};y]};

        scatter(reshape(x,1,15),reshape(y,1,15),trial_m,'SizeData',18,'MarkerFaceColor',col,'MarkerEdgeColor',col);
        hold on
    end
end    

hold on

for i=1:2
    if contains(species(i),"207")
        col = lac_wt_color;
    else
        col = lacz_mut_color;
    end
    y= actual_data{i};
    plot(times, mean(y,1),'Color',col,'LineWidth',2);
    SEM = std(y)/sqrt(size(y,1)); % Standard Error Of The Mean
    disp("YAY")
    hold on
    errorbar(times,mean(y,1),SEM,'Color',col,'LineWidth',2);
    disp(times)
    disp(SEM)
end

trial_1 = actual_data{1};
trial_2 = actual_data{2};
for i=1:5
    [p,h] = ranksum(trial_1(i,:),trial_2(i,:));
    disp(p)
end

% 
set(gca, 'YScale', 'log');
set(gca, 'LineWidth',2)
box off;set(gca,'TickDir','out');
ylim([.003  2])
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

figure('Renderer', 'painters', 'Position', [10 10  100*2 150*2])
innoc_0806 = [345000	300000];
innoc_0813 = [5.95E+05  7.30E+05];
innocs = [innoc_0806, innoc_0813];
x_val = [1 2 3 4];
trial_marker = ["diamond","diamond","o","o"];

for i=1:4
    if contains(species(i),"207")
        col = lac_wt_color;
    else
        col = lacz_mut_color;
    end
    plot([x_val(i)-.25,x_val(i)+.25], ...
        [innocs(i)-.25,innocs(i)+.25], ...
        "-" + trial_marker(i),'Color', col,'MarkerFaceColor',col,'MarkerEdgeColor',col,'LineWidth',2)
    hold on
end
set(gca, 'YScale', 'log');
set(gca, 'LineWidth',3)
box off;set(gca,'TickDir','out');
ylim([200000  800000])
xlim([0, 5])
set(gca,'XTick',[1,2,3,4])
set(gca,'xticklabel',{[]})
box off;set(gca,'TickDir','out');
set(gca,'fontsize',18)


%%%%%%%%%%%
% Define the arrays A and B
A = actual_data{1}; % Replace with your array data
B = actual_data{2}; % Replace with your array data
% Define the arrays A and B

% Calculate means
meanA = mean(A);
meanB = mean(B);

% Calculate standard error of the mean (SEM)
semA = std(A) ./ sqrt(length(A));
semB = std(B) ./ sqrt(length(B));

% Calculate the mean of the element-wise division A./B
mean_A_div_B = mean(A) ./ mean(B);

% Error propagation for f = A / B with errors in A and B:
% error_A_div_B = mean_A_div_B * sqrt((semA / meanA)^2 + (semB / meanB)^2)
error_A_div_B = mean_A_div_B .* sqrt((semA ./ meanA).^2 + (semB ./ meanB).^2);

% Output results
%ylim([1*10^5 10*10^5])

%set(gca,'YTickLabel','fontsize',14)


%300000 730000




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
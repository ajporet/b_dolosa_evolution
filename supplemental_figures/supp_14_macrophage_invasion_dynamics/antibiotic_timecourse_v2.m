clear all; close all; 
uptake_6 = readtable("antibiotic_incubation_macrophage.xlsx",'ReadVariableNames', false);
sample_order = ["P07SP-v311 (O-) 8.45E5","P07SP-v207 (O+) 1.01E6","P07SP-v311 (O-) 1.69E6","P07SP-v207 (O+) 2.02E6"];
% Use the first two columns - closest innoc match 



trial_1_v311 = table2array(uptake_6(1:4,3:5)); % The 20th
trial_1_v207 = table2array(uptake_6(1:4,6:8));

trial_2_v311 = table2array(uptake_6(6:9,3:5)); % The 27th
trial_2_v207 = table2array(uptake_6(6:9,6:8));

species = ["v311","v207"];
times = [1,2,4,6];
times = [2,4,6];
figure('Renderer', 'painters', 'Position', [10 10 200*2 150*2])

trial_marker = ["diamond","o"];
actual_data = {{trial_1_v311,trial_1_v207},{trial_2_v311,trial_2_v207}};
for ii = 1:2
    trial_m = trial_marker(ii);
    data = actual_data{ii};
    
    lacz_mut_color = "#E48D70";lac_wt_color = "#87C7C3";
    lacz_mut_color_edge = "#8A5543";lac_wt_color_edge = "#57807D";

    
    count = 0;
    for i=1:2
        count =count + 1;
        if contains(species(i),"v207")
            col = lac_wt_color;
            edge_col = lac_wt_color_edge;
            disp('y')
        else
            col = lacz_mut_color;
            edge_col = lacz_mut_color_edge;
            disp('n')
        end
    
        x = repmat(times, 3, 1)';
        y= data{i}(2:end,:)./mean(data{i}(2,:));
        disp(species(i))
        disp("trial"+ii)
        disp(y)

        scatter(x,y,trial_m,'SizeData',50,'MarkerFaceColor',col,'MarkerEdgeColor',edge_col,'LineWidth',1.5);
        hold on
    end
end    

hold on

new_data_2 = {[],[]};
for i=1:2
    trial_1 = actual_data{1}{i};
    trial_1 = trial_1(2:end,:);
    trial_2 = actual_data{2}{i};
    trial_2 = trial_2(2:end,:);


    y = [trial_1./mean(trial_1(1,:)),trial_2./mean(trial_2(1,:))]';
    new_data_2(i)= {y};

    if contains(species(i),"207")
        col = lac_wt_color;
    else
        col = lacz_mut_color;
    end

    plot(times, median(y,1),"-_",'Color',col,'LineWidth',2,'MarkerSize', 20);
    scatter(times, median(y,1),"_",'MarkerFaceColor',col,'MarkerEdgeColor',col,'SizeData',750,'LineWidth',4);

    %SEM = std(y)/sqrt(size(y,1)); % Standard Error Of The Mean
%     hold on
%     errorbar(times,mean(y,1),SEM,'Color',col,'LineWidth',2);
end
set(gca, 'LineWidth',3)
box off;set(gca,'TickDir','out');
set(gca,'XTick',times)
box off;set(gca,'TickDir','out');
set(gca,'fontsize',18)
xlim([1.5,6.5])
ylim([.5,2.5])


time_stats = zeros(1,3);
for i=1:3
    p = ranksum(new_data_2{1}(:,i),new_data_2{2}(:,i));
    time_stats(:,i)= p;
end


%%%%%% INOCULUM

innoc_311 = table2array(uptake_6(11:12,4));
innoc_217 = table2array(uptake_6(11:12,7));


figure('Renderer', 'painters', 'Position', [10 10  100*2 150*2])
innoc_0820 = [innoc_311(1)	innoc_217(1)];
innoc_0827 = [innoc_311(2)	innoc_217(2)];
innocs = [innoc_0820, innoc_0827];
x_val = [1 2 3 4];
trial_marker = ["diamond","diamond","o","o"];

species = [species,species];
for i=1:4
    if contains(species(i),"207")
        col = lac_wt_color;
        edge_col = lac_wt_color_edge;
    else
        col = lacz_mut_color;
        edge_col = lacz_mut_color_edge;
    end
    scatter([x_val(i)], ...
        [innocs(i)], ...
         200,trial_marker(i),'Color', col,'MarkerFaceColor',col,'MarkerEdgeColor',edge_col,'LineWidth',3)
    hold on
end
set(gca, 'YScale', 'log');
set(gca, 'LineWidth',3)
box off;set(gca,'TickDir','out');
set(gca,'fontsize',18)
box off;set(gca,'TickDir','out');
xlim([0, 5])
set(gca,'XTick',[1,2,3,4])
set(gca,'xticklabel',{[]})
box off;set(gca,'TickDir','out');
set(gca,'fontsize',18)
ylim([.9*min(innocs),1.1*max(innocs)])
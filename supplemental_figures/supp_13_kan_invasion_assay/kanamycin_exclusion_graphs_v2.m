%% Making all raw graphs
clear all; close all; 
rawfile31 = readtable("08_31_2020_THP_invasion.xlsx", 'Sheet',"Results");
rawfile27 = readtable("08_27_2020_THP_invasion.xlsx", 'Sheet',"Results");
rawfile20 = readtable("08_20_2020_THP_invasion.xlsx", 'Sheet',"Results");

THP_reruns = {rawfile31, rawfile27, rawfile20};
Ordering_of_samples = ["","P06SP-V101","", "P06SP-V103", "", "P07SP-V318","",  "P07SP-V317","",  "P07SP-V218","",  "P07SP-V221"];

invade_locations = [10,11,7,8,9];
total_locations = [4 5 1 2 3];

padding_value = .5;

%specific_assay = rawfile31;
%specific_assay = rawfile27;
%specific_assay = rawfile20;

for spec=THP_reruns
    specific_assay = spec{:};
    
    CFU_values = specific_assay.CFU_mL;
    total_plates= CFU_values(1:28);
    total_plates(4:4:end) = [];
    total_plates(isnan(total_plates)) = [];
    reshaped_total = (reshape(total_plates,3,6));
    mean_vals_total = mean(reshape(total_plates,3,6));
    
    invaded_plates= CFU_values(28+2:56);
    invaded_plates(4:4:end) = [];
    invaded_plates(isnan(invaded_plates)) = [];
    reshaped_invaded = (reshape(invaded_plates,3,6));
    mean_vals_invaded = mean(reshape(invaded_plates,3,6));
    
    x_vals_scatter_total = repelem(total_locations,3);
    x_means_scatter_total = repelem(total_locations,2); x_means_scatter_total(1:2:end) = x_means_scatter_total(1:2:end) - padding_value; x_means_scatter_total(2:2:end) = x_means_scatter_total(2:2:end) + padding_value;
    x_vals_scatter_invade = repelem(invade_locations,3);
    x_means_scatter_invade= repelem(invade_locations,2); x_means_scatter_invade(1:2:end) = x_means_scatter_invade(1:2:end) - padding_value; x_means_scatter_invade(2:2:end) = x_means_scatter_invade(2:2:end) + padding_value;
    
    mean_total_ratio =  mean_vals_invaded./mean_vals_total*100;
    raw_total_ratio = reshaped_invaded./reshaped_total*100;

    %% Plot without v218 because it wasn't resequenced and had its identity confirmed
    
    plotting_order = [2,1,3,4,6];
    
    dot_size = 20;
    
    figure('Renderer', 'painters', 'Position', [10 10 200*2 150*2])
    b = bar(1:5, mean_total_ratio(plotting_order)); %,'facecolor',[.8,.8,.8])
    b.FaceColor = 'flat';
    b.CData(:,:) = repelem([115,115,115]./255,5,1);
    b.EdgeColor = 'none';
    
    
    hold on
    scatter(repelem(1:5,3,1), raw_total_ratio(:,plotting_order),dot_size,'filled','k','jitter','on','jitterAmount',.075)
    
    disp(raw_total_ratio(:,plotting_order))
    box off;set(gca,'TickDir','out');
    set(gca,'xtick',[])
    H=gca;
    H.LineWidth=2;
    set(gca,'FontSize',15); set(gca,'Xticklabel',[]) 
    set(gca,'FontSize',16)
    set(gca,'box','off') 
    set(gca,'TickDir','out');
    y_max = ceil(nanmax(mean_total_ratio));
    y_max = round(y_max*1.5);
    ylim([0,y_max])
    yticks(0:round(y_max/3):y_max);

    LPS_pres_samps = raw_total_ratio(:,plotting_order(1:2));
    LPS_abs_samps = raw_total_ratio(:,plotting_order(3:5));
    disp("Aggregate wilrank score: " + ranksum(LPS_pres_samps(:),LPS_abs_samps(:)));


end
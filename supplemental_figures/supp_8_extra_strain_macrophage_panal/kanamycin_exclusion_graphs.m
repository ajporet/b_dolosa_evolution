%% Making all raw graphs
clear all; close all; 
rawfile31 = readtable("08_31_2020_THP_invasion.xlsx", 'Sheet',"Results");
rawfile27 = readtable("08_27_2020_THP_invasion.xlsx", 'Sheet',"Results");
rawfile20 = readtable("08_20_2020_THP_invasion.xlsx", 'Sheet',"Results");

THP_reruns = {rawfile31, rawfile27, rawfile20};

invade_locations = [10,11,7,8,9];
total_locations = [4 5 1 2 3];

padding_value = .5;

for file=1:length(THP_reruns)
    specific_assay = THP_reruns{file};
    
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
    
    %% Plot without v218 because it wasn't resequenced
    
    plotting_order = [2,1,3,4,6];
    
    dot_size = 20;
    
    mean_plotting_values_invade =  repelem(mean_vals_invaded(plotting_order),2);
    mean_plotting_values_total =  repelem(mean_vals_total(plotting_order),2);
    
    figure('Renderer', 'painters', 'Position', [10 10 300*2 200*2])
    b = bar([total_locations(1:5),invade_locations(1:5)], [mean_vals_total(plotting_order),mean_vals_invaded(plotting_order)]); %,'facecolor',[.8,.8,.8])
    b.FaceColor = 'flat';
    b.CData([1,2,3,6,7,8],:) = repelem([102, 194, 165]./255,6,1);
    b.CData([4,5,9,10],:) = repelem([252, 141, 98]./255,4,1);
    
    % Switched colors
    b.CData([1,2,3,6,7,8],:) = repelem([211, 213, 212]./255,6,1);
    b.CData([4,5,9,10],:) = repelem([126, 126, 126]./255,4,1);
    b.EdgeColor = 'none';
        
    hold on
    scatter(x_vals_scatter_total(1:15), reshape(reshaped_total(:,plotting_order),1,15),dot_size,'filled','k')
    hold on
    scatter(x_vals_scatter_invade(1:15), reshape(reshaped_invaded(:,plotting_order),1,15),dot_size,'filled','k')
    hold on    
    
    box off;set(gca,'TickDir','out');
    set(gca,'xtick',[])
    H=gca;
    H.LineWidth=2;
    set(gca,'FontSize',15); set(gca,'Xticklabel',[]) 
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',16)
    set(gca,'box','off') 
    set(gca,'TickDir','out');
    
end
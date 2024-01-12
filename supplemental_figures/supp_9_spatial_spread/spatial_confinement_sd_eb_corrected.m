
%% Load and clean preset values
function [number_of_unique_sites, average_site_similarity, all_site_similarity, all_sites]  = spatial_confinement_sd_eb(namelist, graph_color)
  
    % Decription: generates a graph of spatial confinement, as introduced in
    % https://www.nature.com/articles/ncomms14078 
   
    load 2020_09_19_pairwise.mat
    
    SNP_scalar = 1./0.1850;
    array_normalized = round(x2020_09_19_pairwise_array.*SNP_scalar);
    names_no_reference = x2020_09_19_pairwise_names(1:931);

    % Removing blood
    names_no_reference_NB_gen = names_no_reference(contains(names_no_reference,{'LT','SL','LG'}));
    array_norm_NB_gen = array_normalized(contains(names_no_reference,{'LT','SL','LG'}),contains(names_no_reference,{'LT','SL','LG'}));

    % Finding the specific clade/input sample indexes
    indexes_of_namelist = contains(names_no_reference_NB_gen, namelist);
    names_no_reference_NB = names_no_reference_NB_gen(indexes_of_namelist);
    array_norm_NB = array_norm_NB_gen(indexes_of_namelist,indexes_of_namelist);
    
    [number_of_unique_sites, average_site_similarity, all_site_similarity, all_sites] = get_site_diversity_stats(names_no_reference_NB,{'LT','SL','LG'});

    %% Processing Data and Graph
    only_upper_array = array_norm_NB;
    only_upper_array(tril(ones(size(array_norm_NB)))==1) = NaN;

    x_less_than = [0:34];
    prob_same_site = zeros(size(x_less_than));

    null_model_sites_to_compare = {};

    for ii=0:34
        [samp_1, samp_2] = find(only_upper_array<=ii);
        null_model_sites_to_compare = [null_model_sites_to_compare,{{samp_1,samp_2}}];
        counting_same_site = sum(all_sites(samp_1)==all_sites(samp_2));
        [phat,pci] = binofit(counting_same_site,length(samp_1),.05);
        prob_same_site(ii+1) = phat;
        confidence_lower(ii+1) = pci(1);
        confidence_upper(ii+1) = pci(2);
    end

    %% Null model
    null_site_scrambles = cell(1,1000);
    null_site_scrambles = cellfun(@(x) all_sites(randperm(length(all_sites))), null_site_scrambles,'UniformOutput', false);
    null_model_prob_same_site = zeros(size(x_less_than));
    standard_dev = zeros(1,length(null_model_sites_to_compare));

    for ik = 1:length(null_model_sites_to_compare)
        disp(ik)
        samp_1 = null_model_sites_to_compare{ik}{1};
        samp_2 = null_model_sites_to_compare{ik}{2};
        number_of_same_site = cellfun(@(x) sum(x(samp_1)==x(samp_2)), null_site_scrambles);
        sorted_prob_values = sort(number_of_same_site,'ascend');

        null_model_prob_same_site(ik) = mean(number_of_same_site./length(samp_1));
        standard_dev(ik) = std(number_of_same_site./length(samp_1));

        null_model_prob_same_site(ik) = phat;

        confidence_null_lower(ik) = sorted_prob_values(round(length(sorted_prob_values).*.025))./length(samp_1);
        confidence_null_upper(ik) = sorted_prob_values(round(length(sorted_prob_values).*.975))./length(samp_1);
    end

    proportionalized_prob_same_site = prob_same_site;

    %% Error "propogarion"

    % Originally we used a different statistic that required error
    % propogation - this was changed to another method that did not require
    % such statistics. Nevertheless, old (albeit confusing) variable names
    % remain. 

    % Experimental data error
    del_confidence_lower = prob_same_site-confidence_lower;
    del_confidence_upper = confidence_upper-prob_same_site;

    % Null model error
    del_confidence_null_lower = null_model_prob_same_site-confidence_null_lower;
    del_confidence_null_upper = confidence_null_upper-null_model_prob_same_site;


    del_confidence_lower_propogated  = del_confidence_lower;
    del_confidence_upper_propogated = del_confidence_upper;

    del_confidence_null_lower_propogated = del_confidence_null_lower;
    del_confidence_null_upper_propogated = del_confidence_null_upper;




    %% Plotting
    markersize = 10;

    % Using public functions to mimic Hattie's paper
    figure;
    ax = axes(); 
    shadedErrorBar(x_less_than, proportionalized_prob_same_site, [del_confidence_upper_propogated;del_confidence_lower_propogated],'lineProps',{"-", 'markerfacecolor',graph_color, 'color',graph_color})
    hold on
    scatter(x_less_than, proportionalized_prob_same_site,markersize,graph_color,'filled','HandleVisibility','off');
    hold on
    shadedErrorBar(x_less_than, null_model_prob_same_site, [del_confidence_null_upper_propogated;del_confidence_null_lower_propogated]);
    hold on
    scatter(x_less_than, null_model_prob_same_site,markersize, 'filled','k','HandleVisibility','off')

    xlabel('Pairwise distance ≤ d (SNVs)','FontSize',14)
    ylabel('Probability of lung site co-localization','FontSize',14)
    xlim([0 34])
    xticks([0:2:34])
    H=gca;
    NumTicks = 5;
    L = get(H,'YLim');
    set(H,'YTick',linspace(L(1),L(2),NumTicks))
    a = get(H,'XTickLabel');
    set(H,'XTickLabel',a,'fontsize',16)
    set(H,'linewidth',2)
    legend({'Site', 'Null model'})
    set(H,'TickDir','out')
    
end
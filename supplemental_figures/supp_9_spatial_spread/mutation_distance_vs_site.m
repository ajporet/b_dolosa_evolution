%% Load and clean preset values
clear;
load 2020_09_19_pairwise.mat
SNP_scalar = 1./0.1850;
array_normalized = round(x2020_09_19_pairwise_array.*SNP_scalar);
names_no_reference = x2020_09_19_pairwise_names(1:931);

% Removing blood
names_no_reference_NB = names_no_reference(contains(names_no_reference,{'LT','SL','LG'}));
array_norm_NB = array_normalized(contains(names_no_reference,{'LT','SL','LG'}),contains(names_no_reference,{'LT','SL','LG'}));

[number_of_unique_sites, average_site_similarity, all_site_similarity, all_sites] = get_site_diversity_stats(names_no_reference_NB,{'LT','SL','LG'});

%% Processing Data and Graph
only_upper_array = array_norm_NB;
only_upper_array(tril(ones(size(array_norm_NB)))==1) = NaN;

x_less_than = [0:34];
prob_same_site = zeros(size(x_less_than));
confidence_lower = zeros(size(x_less_than));
confidence_upper = zeros(size(x_less_than));

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
confidence_null_lower = zeros(size(x_less_than));
confidence_null_upper = zeros(size(x_less_than));


for ik = 1:length(null_model_sites_to_compare)
    samp_1 = null_model_sites_to_compare{ik}{1};
    samp_2 = null_model_sites_to_compare{ik}{2};
    number_of_same_site = cellfun(@(x) sum(x(samp_1)==x(samp_2)), null_site_scrambles);
    
    [phat,pci] = binofit(round(mean(number_of_same_site)),length(samp_1),.05);
    
    null_model_prob_same_site(ik) = phat;
    confidence_null_lower(ik) = pci(1);
    confidence_null_upper(ik) = pci(2);
end

proportionalized_prob_same_site = prob_same_site./null_model_prob_same_site;
%% Error propogarion

% Experimental data error
del_confidence_lower = prob_same_site-confidence_lower;
del_confidence_upper = confidence_upper-prob_same_site;

% Null model error
del_confidence_null_lower = null_model_prob_same_site-confidence_null_lower;
del_confidence_null_upper = confidence_null_upper-null_model_prob_same_site;

% Propogating error
del_confidence_lower_propogated = proportionalized_prob_same_site.*sqrt( (del_confidence_lower./prob_same_site).^2 + (del_confidence_null_lower./null_model_prob_same_site).^2);
del_confidence_upper_propogated = proportionalized_prob_same_site.*sqrt( (del_confidence_upper./prob_same_site).^2 + (del_confidence_null_upper./null_model_prob_same_site).^2);

del_confidence_null_lower_propogated = sqrt(2.*((del_confidence_null_lower./null_model_prob_same_site).^2)).*1;
del_confidence_null_upper_propogated = sqrt(2.*((del_confidence_null_upper./null_model_prob_same_site).^2)).*1;


%% Plotting
markersize = 10;

%Basic plot
figure;
plot(x_less_than, proportionalized_prob_same_site,'Color','r', 'Marker','o','MarkerFaceColor','r','MarkerSize',3);
hold on
errorbar(x_less_than, proportionalized_prob_same_site, del_confidence_lower_propogated, del_confidence_upper_propogated)

hold on 

plot(x_less_than, null_model_prob_same_site./null_model_prob_same_site,'Color','b', 'Marker','o','MarkerFaceColor','b','MarkerSize',3);
hold on
errorbar(x_less_than, null_model_prob_same_site./null_model_prob_same_site, del_confidence_null_lower_propogated, del_confidence_null_upper_propogated)


% Using public functions to mimic Hattie's paper
figure;
shadedErrorBar(x_less_than, proportionalized_prob_same_site, [del_confidence_upper_propogated;del_confidence_lower_propogated],'lineProps',{'r-'})
hold on
scatter(x_less_than, proportionalized_prob_same_site,markersize, 'filled','r','HandleVisibility','off')

hold on
shadedErrorBar(x_less_than, null_model_prob_same_site./null_model_prob_same_site, [del_confidence_null_upper_propogated;del_confidence_null_lower_propogated])
hold on
scatter(x_less_than, null_model_prob_same_site./null_model_prob_same_site,markersize, 'filled','k','HandleVisibility','off')

xlabel('Pairwise distance ≤ d (SNPs)','FontSize',14)
ylabel('Spatial confinement η(d)','FontSize',14)
xlim([0 34])
xticks([0:2:34])
yticks([0:1:5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
set(gca,'linewidth',1)
legend({'Site', 'Null'})
set(gca,'TickDir','out');
set(gca, 'XScale', 'log')

%% New graph
%Basic plot
x_less_than_nonzero = x_less_than;
x_less_than_nonzero(1) = .01; % Due to bug in patch

figure;
plot(x_less_than_nonzero, proportionalized_prob_same_site,'Color','r', 'Marker','o','MarkerFaceColor','r','MarkerSize',3);
%hold on
%plot(x_less_than, proportionalized_prob_same_site - del_confidence_lower_propogated,'r')
%hold on
%plot(x_less_than, proportionalized_prob_same_site + del_confidence_upper_propogated,'r')
patch([x_less_than_nonzero fliplr(x_less_than_nonzero)], [(proportionalized_prob_same_site - del_confidence_lower_propogated) fliplr(proportionalized_prob_same_site + del_confidence_upper_propogated)], 'r','Edgecolor','none')
alpha(0.3)
hold on 

plot(x_less_than_nonzero, null_model_prob_same_site./null_model_prob_same_site,'Color','k', 'Marker','o','MarkerFaceColor','k','MarkerSize',3);
hold on
patch([x_less_than_nonzero fliplr(x_less_than_nonzero)], [(null_model_prob_same_site./null_model_prob_same_site - del_confidence_null_lower_propogated) fliplr(null_model_prob_same_site./null_model_prob_same_site + del_confidence_null_lower_propogated)], 'k','Edgecolor','none')
alpha(0.3)

%errorbar(x_less_than, null_model_prob_same_site./null_model_prob_same_site, del_confidence_null_lower_propogated, del_confidence_null_upper_propogated)

set(gca, 'XScale', 'log')
%set(gca,'XScale',logfreq);

xlim([0 34])
xticks([0:2:34])
yticks([0:1:5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
set(gca,'linewidth',1)
legend({'Site', 'Null'})
set(gca,'TickDir','out');

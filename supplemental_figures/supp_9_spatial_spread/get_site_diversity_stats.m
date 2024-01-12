function [number_of_unique, average_difference_site, site_similarity, sites] = get_site_diversity_stats(stringlist,splitter_string)
    sites = get_site_numbers(stringlist,splitter_string);
    number_of_unique = max(size(unique(sites)));
    site_similarity = cellfun(@(x) (sum(x==sites)-1)./(length(sites)-1),sites);
    average_difference_site = mean(site_similarity);
end 
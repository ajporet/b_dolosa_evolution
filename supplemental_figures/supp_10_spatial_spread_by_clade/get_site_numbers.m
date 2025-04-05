function [list_of_site_numbers] = get_site_numbers(stringlist,splitter_string)
    sites_A = cellfun(@(x) strsplit(x,splitter_string),stringlist,'UniformOutput',false);
    sites_A = cellfun(@(x) x(2), sites_A);
    site_numbers = regexp(sites_A,'\d*','Match');
    site_numbers = [site_numbers{:}];
    site_holder = strings(size(site_numbers));
    bloodsite = cellfun(@(x) strlength(x)==3, site_numbers);
    double_delim = cellfun(@(x) strlength(x)==4, site_numbers);
    triple_delim = cellfun(@(x) strlength(x)==5, site_numbers);
    
    site_holder(bloodsite) = cellfun(@(x) "B" + x(1),site_numbers(bloodsite),'UniformOutput',false);
    site_holder(double_delim) = cellfun(@(x) x(1:2),site_numbers(double_delim),'UniformOutput',false);
    site_holder(triple_delim) = cellfun(@(x) x(1:3),site_numbers(triple_delim),'UniformOutput',false);
    list_of_site_numbers = site_holder;
end

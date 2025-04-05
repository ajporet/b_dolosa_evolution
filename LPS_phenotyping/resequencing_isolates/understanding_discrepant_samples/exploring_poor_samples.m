%% PLEASE RUN AFTER analysis.m variables have been generated

poor_samples = {["P07SP-v318", "P07SP-v318_flipped"],

["P07SP-v302_deep_well", "P07SP-v302"],

["P07SP-v207_flipped", "P07SP-v207"],

["P07SP-v221_flipped", "P07SP-v316"],

["P07SP-v205_flipped", "P07SP-v316"],

["P07SP-v316_new_seq", "P07SP-v316"],

["P06SP-v305_resequenced_2", "P06SP-v324"],

["P02BLB3e08_first_resequenced_1_5", "P02LT-0815"],

["P02LT-2407_first_resequenced_1_5", "P02LT-2403"],

["P06SP-v101_lacZ","P06SP-v101_unaltered"]};

mismatched_sample_desc = [];

for ii=1:length(poor_samples)
    pair_to_consider = poor_samples{ii};
    sample_1 = pair_to_consider(1);
    sample_2 = pair_to_consider(2);

    sample_id_1 = find(SampleNames(samplestoplot)==sample_1);
    sample_id_2 = find(SampleNames(samplestoplot)==sample_2);

    mismatched_locs = find(calls_for_treei(:,sample_id_1)~=calls_for_treei(:,sample_id_2) & calls_for_treei(:,sample_id_2)~=0 & calls_for_treei(:,sample_id_1)~=0);
    true_calls = calls_for_treei(mismatched_locs,[sample_id_1,sample_id_2]);
    coverage_bad_calls = coverage(positions_to_show_in_table(mismatched_locs),[sample_id_1,sample_id_2]);
    centile_of_coverage = [];
    for ij=1:size(coverage_bad_calls,1)
        cent_temp_1 = centile_finder(coverage(positions_to_show_in_table,sample_id_1),coverage_bad_calls(ij,1));
        cent_temp_2 = centile_finder(coverage(positions_to_show_in_table,sample_id_1),coverage_bad_calls(ij,2));
        centile_of_coverage = [centile_of_coverage;[cent_temp_1,cent_temp_2]];
    end

    mismatched_sample_desc = [mismatched_sample_desc; [repelem(sample_1,size(true_calls,1))',repelem(sample_2,size(true_calls,1))',true_calls,coverage_bad_calls,centile_of_coverage,mismatched_locs]];
end

supposed_mutations = annotation_full(double(mismatched_sample_desc(:,9)));
mismatched_sample_desc_no_repetitive_values = mismatched_sample_desc(mismatched_sample_desc(:,9)~="308",:);
supposed_mutations_no_repetitive_values = annotation_full(double(mismatched_sample_desc_no_repetitive_values(:,9)));

mismatched_sample_desc_no_repetitive_no_low_cov_values  = mismatched_sample_desc_no_repetitive_values(double(mismatched_sample_desc_no_repetitive_values(:,5))>2,:);
supposed_mutations_no_repetitive_values_no_low_cov_values  = annotation_full(double(mismatched_sample_desc_no_repetitive_no_low_cov_values(:,9)));

disp(mismatched_sample_desc_no_repetitive_no_low_cov_values)
disp(supposed_mutations_no_repetitive_values_no_low_cov_values)


% Calculatethe number of samples (if randomly grabbed) that one would
% expect to fall into a filtered out isolate

fraction_of_filtered_isolates = (987-981)/987;
all_used_samp_names = SampleNames(samplestoplot);
resequenced_isolates = contains(all_used_samp_names, "P0") & contains(all_used_samp_names, ["resequenced", "deep_well","seq","unaltered","repulled","flipped"]);
unique_resequencing = cellfun(@(x) split(x,"_",1), all_used_samp_names(resequenced_isolates), 'UniformOutput', false);
[unique_resequencing, unique_idx] = unique(cellfun(@(x) x(1),unique_resequencing));

[unique_genotypes, ~, ~] = unique(hasmutation(:,~resequenced_isolates)', 'rows');
chance_of_unique_isolate = size(unique_genotypes,1)/sum(~resequenced_isolates);
expected_number_of_unsequenced_matches = fraction_of_filtered_isolates*length(unique_resequencing)*chance_of_unique_isolate;


function centile = centile_finder(data,value)
    % from https://www.mathworks.com/matlabcentral/answers/182131-percentile-of-a-value-based-on-array-of-data
    data = data(:)';
    value = value(:);
    
    nless = sum(data < value, 2);
    nequal = sum(data == value, 2);
    centile = 100 * (nless + 0.5.*nequal) / length(data);
end


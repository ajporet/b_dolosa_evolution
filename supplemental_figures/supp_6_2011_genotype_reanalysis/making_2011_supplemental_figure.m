clear; close all;

% Get isolates with LPS mutations
all_2011_mutations= readcell("Lieberman_et_al_2011_supp_table_2.xlsx");
sample_names = string(all_2011_mutations(1,[6,9:120]));
O_antigen_gene_names = ["Glycosyltransferase","Mannose-1-phosphate","O-antigen"];
O_antigen_gene_locations = find(contains(string(all_2011_mutations(2:end,4)), O_antigen_gene_names));
stop_codon_reversion_location = 232;

last_common_ancester_calls = string(all_2011_mutations(2:end,8));
sample_calls = string(all_2011_mutations(2:end,[6,9:120]));

mutation_matrix = (last_common_ancester_calls ~= sample_calls) & ~cellfun(@isempty, sample_calls);
no_mutation_matrix = (last_common_ancester_calls == sample_calls) & ~cellfun(@isempty, sample_calls);
has_reversion = mutation_matrix(stop_codon_reversion_location,:);

%has_O_antigen_indel = ["J-12-11","J-12-4","J-13-6a-$","J-13-6b","J-14-2","J-14-8"];
has_O_antigen_indel = ["J-12-11","J-13-6a-$","J-13-6b","J-14-2"];

extra_mutation = (has_reversion & sum(mutation_matrix(setdiff(O_antigen_gene_locations, stop_codon_reversion_location),:))>0) | (contains(sample_names, has_O_antigen_indel));
only_reversion = has_reversion & sum(no_mutation_matrix(setdiff(O_antigen_gene_locations, stop_codon_reversion_location),:))>=length(O_antigen_gene_locations)-1 & ~contains(sample_names, has_O_antigen_indel);
stop_codon = no_mutation_matrix(stop_codon_reversion_location,:);

%% Get gel images
name_dictionary = string(table2cell(readtable("2011strainDict.txt")));
dates_of_sampling_temp = name_dictionary(:,2);
dates_of_sampling = NaT(size(dates_of_sampling_temp));
name_dict_samp_names = name_dictionary(:,9);

for ij =1:length(dates_of_sampling_temp)
    sample_to_parse = char(dates_of_sampling_temp(ij));
    if str2num(sample_to_parse([end-1:end]))>10
        sample_to_parse([end-3:end-2]) = "19";
    else
        sample_to_parse([end-3:end-2]) = "20";
    end
    dates_of_sampling(ij) = datetime(sample_to_parse,'InputFormat','MM/dd/uuuu');
end

gel1 = mat2gray(rgb2gray(imread('Gel Images/gel1.png')));
gel2 = mat2gray(rgb2gray(imread('Gel Images/gel2.png')));
gel3 = mat2gray(rgb2gray(imread('Gel Images/gel3.png')));

gel1_xcrop = [95.6948261924009,125.554163298302,152.095796281326,185.272837510105,223.979385610348,259.368229587712,290.333468067906,331.251818916734,361.111156022635,395.394098625707,425.253435731609,459.536378334681,498.242926434923,526.996362166532,557.961600646726,594.456345998383];
gel2_xcrop = [43.1522821576764,77.9937759336100,117.190456431535,149.854356846473,190.139834024896,219.537344398340,262.000414937759,290.309128630705,330.594605809129,370.880082987552,406.810373443983,448.184647302905,477.582157676349,514.601244813278,553.797925311203,589.728215767635];
gel3_xcrop = [94.9470691163605,131.940944881890,157.001312335958,186.835083114611,225.022309711286,265.596237970254,300.203412073491,326.457130358705,363.451006124235,401.638232720910,437.438757655293,474.432633420823,507.846456692914,536.486876640420,572.287401574803,606.894575678040];

gel1_lanes = get_lanes(gel1, gel1_xcrop, 253.7414, 587.6724);
gel2_lanes = get_lanes(gel2, gel2_xcrop, 205.0176, 560.0294);
gel3_lanes = get_lanes(gel3, gel3_xcrop, 254.7385, 605.4122);

gel1_names = translate_marker(["J1", "AU0158", "K1", "N2", "A4", "A1", "S4"],name_dictionary);
gel2_names = translate_marker(["R4", "M3", "M5", "S3", "C1", "C17", "C2"],name_dictionary);
gel3_names = translate_marker(["N1", "L3", "B12", "D2", "D10", "D11", "D6"],name_dictionary);

all_phenotyped_samps = [gel1_names; gel2_names; gel3_names];
all_samp_lanes = [gel1_lanes, gel2_lanes, gel3_lanes];

%% Catagorize gels 
one_LPS_mut_lanes = all_samp_lanes(contains(all_phenotyped_samps,sample_names(only_reversion)));
mult_LPS_muts_lanes = all_samp_lanes(contains(all_phenotyped_samps,sample_names(extra_mutation)));
no_LPS_muts_lanes = all_samp_lanes(contains(all_phenotyped_samps,sample_names(stop_codon)));

one_LPS_mut_lanes_names = all_phenotyped_samps(contains(all_phenotyped_samps,sample_names(only_reversion)));
mult_LPS_muts_lanes_names = all_phenotyped_samps(contains(all_phenotyped_samps,sample_names(extra_mutation)));
no_LPS_muts_lanes_names = all_phenotyped_samps(contains(all_phenotyped_samps,sample_names(stop_codon)));

figure;
sgtitle("Single LPS mutation: stop codon reversion")
for ii=1:length(one_LPS_mut_lanes)
    subplot(1,length(one_LPS_mut_lanes),ii)
    imshow(one_LPS_mut_lanes{ii})
    imwrite(one_LPS_mut_lanes{ii},"individual_gel_images/reversion/" + one_LPS_mut_lanes_names(ii) + ".png")

end

figure;
sgtitle("Multiple LPS mutations: stop codon + other")
for ii=1:length(mult_LPS_muts_lanes)
    subplot(1,length(mult_LPS_muts_lanes),ii)
    imshow(mult_LPS_muts_lanes{ii})
    imwrite(mult_LPS_muts_lanes{ii},"individual_gel_images/further_mutation/" + mult_LPS_muts_lanes_names(ii) + ".png")

end

figure;
sgtitle("No stop codon mutation")
for ii=1:length(no_LPS_muts_lanes)
    subplot(1,length(no_LPS_muts_lanes),ii)
    imshow(no_LPS_muts_lanes{ii})
    imwrite(no_LPS_muts_lanes{ii},"individual_gel_images/no_stop_codon/" + no_LPS_muts_lanes_names(ii) + ".png")
end

%% Plot time since first sample

time_since_sampling = cellfun(@(x) split(x(3:end),"-"), sample_names, 'UniformOutput', false);
time_since_sampling = cellfun(@(x) str2num(cell2mat(regexprep(x(1),'[^0-9]',''))) + str2num(cell2mat(regexprep(x(2),'[^0-9]','')))/12, time_since_sampling(2:end));
time_since_sampling = [3.15, time_since_sampling]; % Account for AU0158 date being difficult to parse

patient_letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"];

figure('Renderer', 'painters', 'Position', [10 10 400*2 200*2])
for ii=1:length(patient_letters)
    patient_sample = contains(sample_names, patient_letters(ii));
    scatter(time_since_sampling(patient_sample & extra_mutation) - min(time_since_sampling(patient_sample)),ii.*ones(sum(patient_sample & extra_mutation),1),'r','filled','MarkerFaceAlpha',.5)
    hold on
    scatter(time_since_sampling(patient_sample & only_reversion) - min(time_since_sampling(patient_sample)),ii.*ones(sum(patient_sample & only_reversion),1),'b','filled','MarkerFaceAlpha',.5)
    hold on
    scatter(time_since_sampling(patient_sample & stop_codon) - min(time_since_sampling(patient_sample)),ii.*ones(sum(patient_sample & stop_codon),1),'k','filled','MarkerFaceAlpha',.5)
    hold on 

    max_stop = max(time_since_sampling(patient_sample & stop_codon));
    min_edit = min([time_since_sampling(patient_sample & only_reversion),time_since_sampling(patient_sample & extra_mutation)]);

    if max_stop>min_edit
        disp(patient_letters(ii))
    end
    
end

H=gca;
H.LineWidth=2;
yticks([1:length(patient_letters)])
set(gca,'FontSize',16)
set(gca,'box','off') 
set(gca,'TickDir','out');
xlabel("years from first sampled isolate ")
ylabel('patient')
yticklabels(patient_letters)
xlim([-.1,8])

%% Plot time of sampling
time_of_sampling = cellfun(@(x) dates_of_sampling(strcmp(name_dict_samp_names,x)), sample_names(2:end));
time_of_sampling = [datetime("REDACTED",'InputFormat','MM/dd/uuuu'),time_of_sampling];
name_dict_samp_names = name_dictionary(:,9);

figure('Renderer', 'painters', 'Position', [10 10 400*2 200*2])
for ii=1:length(patient_letters)
    patient_sample = contains(sample_names, patient_letters(ii));
    scatter(time_of_sampling(patient_sample & extra_mutation),ii.*ones(sum(patient_sample & extra_mutation),1),'r','filled','MarkerFaceAlpha',.5)
    hold on
    scatter(time_of_sampling(patient_sample & only_reversion),ii.*ones(sum(patient_sample & only_reversion),1),'b','filled','MarkerFaceAlpha',.5)
    hold on
    scatter(time_of_sampling(patient_sample & stop_codon),ii.*ones(sum(patient_sample & stop_codon),1),'k','filled','MarkerFaceAlpha',.5)
    hold on 
end

H=gca;
H.LineWidth=2;
yticks([1:length(patient_letters)])
set(gca,'FontSize',16)
set(gca,'box','off') 
set(gca,'TickDir','out');
xlabel("date of sampling ")
ylabel('patient')
yticklabels(patient_letters)
time_span = max(time_of_sampling) - min(time_of_sampling);
xlim([min(time_of_sampling)-time_span*.2/16,max(time_of_sampling)+time_span*.2/16])

%% Functions
function translated_names = translate_marker(marker_names, description_2011)
    dictionary_BD_marker_temp = strrep(description_2011(:,1),'BD_','');
    dictionary_BD_marker_temp = strrep(dictionary_BD_marker_temp,'BM_','');
    
    % Removing extra 0's
    has_extra_zero_temp = arrayfun(@(x) x{:}(2) == "0", dictionary_BD_marker_temp);
    dictionary_BD_marker_temp(has_extra_zero_temp) =  erase(dictionary_BD_marker_temp(has_extra_zero_temp),'0');
    
    dictionary_BD_marker = dictionary_BD_marker_temp;

    indexes_of_marker = cellfun(@(x) find(dictionary_BD_marker == x), marker_names);
    translated_names = description_2011(:,9);
    translated_names = translated_names(indexes_of_marker);
end
function lane_image_cell_array = get_lanes(gel_image, x_crop_coor, bottom_crop, top_crop)
% All bottom crops occur at the bottom standard line. The top crop is
% present at the top standard line  
    x_crop_coor = round(x_crop_coor);
    bottom_crop = round(bottom_crop);
    top_crop = round(top_crop);
    
    lane_image_cell_array = cell(1,length(x_crop_coor)./2-1);
    for ii=[3:2:16] % This skips the standard
        gel_lane = gel_image(bottom_crop:top_crop,x_crop_coor(ii):x_crop_coor(ii+1));
        lane_image_cell_array((ii-1)/2) = {gel_lane};
    end
end
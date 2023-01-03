%% Parameters  -- IMPORTANT TO ADJUST !!!!
% Much of this will need to vary according to your reference genome,
% coverage, and particular samples

% NOTE: this analysis.m is editted for for the specific case of combining
% sequencing of the 2011 dolosa and 2020 dolosa samples. Because one is
% paired end and the other is not, both require different filtering. 

global refgenome
load('candidate_mutation_table.mat')
refgenome='Bdolosa_new';

%% Define different sample sets 
% In this code, Lieberman et al.'s data is processed seperately from the
% modern, "2020" data (date references the time of this code's writing).
% This is because the new dataset is from high coverage, modern Illumina
% sequencing, while the 2011 samples are lower covrage single end 454
% reads. Different filtering cutoffs are needed for both datasets. As such,
% we process each set of samples differently, find all relevant polymorphic
% positions, and then combine those positions into a master list for final
% tree making. 

samps_2011 = contains(SampleNames,'_');
samps_2020 = ~samps_2011;

%%  2022 study cutoffs
% Include sample in tree
min_average_coverage_to_include_sample_2020 = 12; 
include_dirty_sample_in_tree_2020 = 10; %maximum number of N's in a sample after all other filtering - removes a sample 

%Filter on a per-call basis
min_maf_for_call_2020 = .9; % maf = major allele frequency
min_cov_per_strand_for_call_2020 = 12;% cov = coverage
min_qual_for_call_2020 = 0; % qual = quality from samtools

% Cuttoffs for all positions once samples have been cleaned
max_fraction_ambigious_samples_2020 = .2; % ambiguous = a sample that failed one filter at that position, used to be .15
min_median_coverage_position_2020 = 40;

%% Lieberman et al. (2011) cutoffs
% Include sample in tree
min_average_coverage_to_include_sample_2011 = 0; 
include_dirty_sample_in_tree_2011 = 30; %maximum number of N's in a sample after all other filtering - removes a sample 

%Filter on a per-call basis
min_maf_for_call_2011 = .9; % maf = major allele frequency
min_cov_TOTAL_for_call_2011 = 20;% cov = coverage
min_qual_for_call_2011 = 0; % qual = quality from samtools

% Cuttoffs for all positions once samples have been cleaned
max_fraction_ambigious_samples_2011 = .4; 
min_median_coverage_position_2011 = 20; 

%% General sample cutoffs
promotersize=250;
min_maf_for_analysis=.80;
disp('Params set')

%% Enviornment set up 
workingdir=char(pwd);
REFGENOMEFOLDER=['~/Dropbox (MIT)/Lieberman Lab/Reference_Genomes/' refgenome];
SCRIPTSDIRECTORY = ['~/Dropbox (MIT)/Lieberman Lab/scripts'];
path(SCRIPTSDIRECTORY,path);
NTs='ATCG';

%% Remove undesired samples based on name and/or coverage
coverage=squeeze(sum(counts(1:8,:,:)));

goodsamples_2011 = mean(coverage)>min_average_coverage_to_include_sample_2011;
goodsamples_2020 = mean(coverage)>min_average_coverage_to_include_sample_2020;

goodsamples = (samps_2011.*goodsamples_2011' + samps_2020.*goodsamples_2020')==1;

SampleNames=SampleNames(goodsamples);
counts=counts(:,:,goodsamples);
Quals=Quals(:,goodsamples);
coverage=coverage(:,goodsamples);
num_samples=numel(SampleNames);

samps_2011_cleaned = contains(SampleNames,'_');
samps_2020_cleaned = ~samps_2011_cleaned;

coverage_forward_strand=squeeze(sum(counts(1:4,:,:)));
coverage_reverse_strand=squeeze(sum(counts(5:8,:,:)));
coverage_total = squeeze(sum(counts(1:8,:,:)));

Quals = -1*Quals; %use -Quals because this way higher numbers are more confident

%% Read in genome information
[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);
refnt = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p,ChrStarts));

ancnt=refnt;
[~,ancnti]=ismember(refnt,NTs);
ancnti_m=repmat(ancnti,1,num_samples);

%% Make some basic structures for finding mutations
contig_positions=p2chrpos(p, ChrStarts);
[maf, maNT, minorNT, minorAF] = div_major_allele_freq(counts);
[~,refnti]=ismember(refnt,NTs);

mutantAF=zeros(size(maNT));
mutantAF(maNT~=ancnti_m)=maf(maNT~=ancnti_m);
mutantAF(minorNT~=ancnti_m)=mutantAF(minorNT~=ancnti_m)+minorAF(minorNT~=ancnti_m); %this construction allows for positions with two different mutations

disp('Basic structures made')
%% Find positions with fixed mutations

num_samples_2011 = sum(samps_2011_cleaned);
num_samples_2020 = sum(samps_2020_cleaned);

% Find the nucleotide identity at each position.
Calls_2011=maNT(:,samps_2011_cleaned);
Calls_2020=maNT(:,samps_2020_cleaned);

% 2020 filters
Calls_2020(Quals(:,samps_2020_cleaned) < min_qual_for_call_2020 | maf(:,samps_2020_cleaned)< min_maf_for_call_2020 | coverage_forward_strand(:,samps_2020_cleaned) < min_cov_per_strand_for_call_2020 | coverage_reverse_strand(:,samps_2020_cleaned) < min_cov_per_strand_for_call_2020)=0; 
Calls_2020(sum(Calls_2020<1,2)>=(num_samples_2020*max_fraction_ambigious_samples_2020) | median(coverage(:,samps_2020_cleaned),2)<min_median_coverage_position_2020,:)=0;

% 2011 filters
Calls_2011(Quals(:,samps_2011_cleaned) < min_qual_for_call_2011 | maf(:,samps_2011_cleaned)< min_maf_for_call_2011 | coverage_total(:,samps_2011_cleaned) < min_cov_TOTAL_for_call_2011)=0; 
Calls_2011(sum(Calls_2011<1,2)>=(num_samples_2011*max_fraction_ambigious_samples_2011) | median(coverage(:,samps_2011_cleaned),2)<min_median_coverage_position_2011,:)=0;

%[MutQual, MutQualIsolates] = ana_mutation_quality(Calls,Quals) ; 
[MutQual_2011, MutQualIsolates_2011] = ana_mutation_quality(Calls_2011,Quals(:,samps_2011_cleaned)) ; 
[MutQual_2020, MutQualIsolates_2020] = ana_mutation_quality(Calls_2020,Quals(:,samps_2020_cleaned)) ; 


fixedmutation_2011=((maNT(:,samps_2011_cleaned)~=repmat(ancnti,1,num_samples_2011)) & maNT(:,samps_2011_cleaned)>0 & repmat(MutQual_2011,1,num_samples_2011)>=1 & maf(:,samps_2011_cleaned)> min_maf_for_analysis);
fixedmutation_2020=((maNT(:,samps_2020_cleaned)~=repmat(ancnti,1,num_samples_2020)) & maNT(:,samps_2020_cleaned)>0 & repmat(MutQual_2020,1,num_samples_2020)>=1 & maf(:,samps_2020_cleaned)> min_maf_for_analysis);

hasmutation= zeros(size(maNT));
hasmutation(:,samps_2011_cleaned) =fixedmutation_2011;
hasmutation(:,samps_2020_cleaned) =fixedmutation_2020;

goodpos=find(sum(hasmutation,2)>0);

disp('Finding SNPs')

%% Display table (with or without annotation)

positions_to_show_in_table=goodpos;
samples_to_show=1:numel(SampleNames);
num_contigs=max(contig_positions(:,1));
contig_lengths=[ChrStarts(2:end) GenomeLength]-ChrStarts;

% Whether or not to sort the clickable table by quality
QualSort=0; % toggle to sort the table by quality or not
SpecificSamples = 0; %manually input a sample list to observe s
SpecificSite = 0; %manualy select a site NOT WORKING 

sitetoparse = positions_to_show_in_table;

% Uncomment the following three lines if there is an annotated genome
annotations = annotate_mutations_gb(p2chrpos(p(positions_to_show_in_table),ChrStarts),REFGENOMEFOLDER) ;
annotation_full= append_annotations(annotations, ancnti(positions_to_show_in_table), Calls(positions_to_show_in_table,:), counts(:,positions_to_show_in_table,:), hasmutation(positions_to_show_in_table,:), promotersize) ; %% adds information about particular mutations observed, based on
%clickable_snp_table(annotation_full, Calls(sitetoparse,samples_to_show), counts(:,sitetoparse,samples_to_show), SampleNames(samples_to_show), ScafNames, MutQual(positions_to_show_in_table), QualSort);
hasmutation=hasmutation(goodpos,:);

chrpos=p2chrpos(p,ChrStarts);

%% Finding potential contaminants

maf_good=maf(positions_to_show_in_table,:);
maf_mean=mean(maf_good);
maf_outliers=sum(maf_good<.9);

 %% parsimony tree
 disp('Making tree')

samplestoplot=1:numel(SampleNames);
quality_positions=goodpos;
 
calls_for_treei=maNT(quality_positions,samplestoplot);
maf2 = maf(quality_positions,samplestoplot);
calls_for_treei(maf2 < min_maf_for_analysis)=0;

%checks for samples with too many N's
badsamples_temp = sum(calls_for_treei==0);
badsamples_2020 = (badsamples_temp>include_dirty_sample_in_tree_2020) & samps_2020_cleaned';
badsamples_2011 = (badsamples_temp>include_dirty_sample_in_tree_2011) & samps_2011_cleaned'; 
badsamples = badsamples_2020 | badsamples_2011;

indices_dirty_samples = find(badsamples);

samplestoplot(indices_dirty_samples) = [];
%reruns code to reform calls_for_treei

calls_for_treei=maNT(quality_positions,samplestoplot);
maf2 = maf(quality_positions,samplestoplot);
calls_for_treei(maf2 < min_maf_for_analysis)=0;

calls_for_tree=zeros(size(calls_for_treei));
calls_for_tree(calls_for_treei > 0)=NTs(calls_for_treei(calls_for_treei>0));
calls_for_tree(calls_for_treei < 1)='?';


% ADD REFERENCE AT THESE POSITIONS -- don't bother in this case because
% this is from de novo assembly
outgroup_nts = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p(quality_positions),ChrStarts));
 
TreeSampleNames= {SampleNames{samplestoplot}};
 
calls_for_tree = [outgroup_nts, calls_for_tree];
TreeSampleNames={'Reference' SampleNames{samplestoplot}};


%% Generate Parsimony tree
disp('Renaming Tree')
 [treefilename, tempnames] =generate_parsimony_tree_rename(calls_for_tree, TreeSampleNames, refgenome);
 
fprintf(1,'\nDone with tree\n');
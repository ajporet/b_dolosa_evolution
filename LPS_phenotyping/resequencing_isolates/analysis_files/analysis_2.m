%% Parameters  -- IMPORTANT TO ADJUST !!!!
% Much of this will need to vary according to your reference genome,
% coverage, and particular samples
clear all; close all;

global refgenome
load('candidate_mutation_table.mat')

refgenome='Bdolosa_new';
min_average_coverage_to_include_sample = 5; %used to be 8, min cov 9
                        % CHANGED TO 5 DUE TO NEW DATA, used to be 10

%Filter on a per-call basis
min_maf_for_call = .9; % maf = major allele frequency
min_cov_per_strand_for_call = 10;% cov = coverage
min_qual_for_call = 0;  % qual = quality from samtools

%Filter on a by-position basis (across samples)
max_fraction_ambigious_samples = .2; % ambiguous = a sample that failed one filter at that position
min_median_coverage_position = 20; % note: related to min_average_coverage_to_include_sample

promotersize=250;

%filters for using once we have a trusted position
min_maf_for_analysis=.80;
include_dirty_sample_in_tree = 1000; %maximum number of N's in a sample after all other filtering 
%above filter completely removes a sample - best to set after seeing how other filters interact with samples
disp('Params set')
%% Enviornment set up -- probably won't need to change

workingdir=char(pwd);
REFGENOMEFOLDER=['~/Dropbox (MIT)/Lieberman Lab/Reference_Genomes/' refgenome];
SCRIPTSDIRECTORY = ['~/Dropbox (MIT)/Lieberman Lab/scripts'];
path(SCRIPTSDIRECTORY,path);

NTs='ATCG';

%% Remove undesired samples based on name and/or coverage
coverage=squeeze(sum(counts(1:8,:,:)));

goodsamples = mean(coverage) > min_average_coverage_to_include_sample;

SampleNames=SampleNames(goodsamples);
counts=counts(:,:,goodsamples);
Quals=Quals(:,goodsamples);
coverage=coverage(:,goodsamples);
num_samples=numel(SampleNames);

coverage_forward_strand=squeeze(sum(counts(1:4,:,:)));
coverage_reverse_strand=squeeze(sum(counts(5:8,:,:)));

Quals = -1*Quals; %use -Quals because this way higher numbers are more confident

%% Read in genome information

[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);

refnt = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p,ChrStarts));

% Is be more complicated if refnt ~= ancnt
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

% Find the nucleotide identity at each position.
Calls=maNT;
Calls(abs(Quals) < min_qual_for_call | maf< min_maf_for_call | coverage_forward_strand < min_cov_per_strand_for_call | coverage_reverse_strand < min_cov_per_strand_for_call)=0; %Quals < min_qual_for_call |

Calls(sum(Calls<1,2)>=(num_samples*max_fraction_ambigious_samples) | median(coverage,2)<min_median_coverage_position,:)=0;

[MutQual, MutQualIsolates] = ana_mutation_quality(Calls,Quals) ; 


fixedmutation=((maNT~=repmat(ancnti,1,num_samples)) & maNT>0 & repmat(MutQual,1,num_samples)>=1 & maf> min_maf_for_analysis);

hasmutation= fixedmutation;
goodpos=find(sum(fixedmutation,2)>0);

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
clickable_snp_table(annotation_full, Calls(sitetoparse,samples_to_show), counts(:,sitetoparse,samples_to_show), SampleNames(samples_to_show), ScafNames, MutQual(positions_to_show_in_table), QualSort);
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
badsamples = sum(calls_for_treei==0);
indices_dirty_samples = find(badsamples>include_dirty_sample_in_tree);

samplestoplot(indices_dirty_samples) = [];

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
stop
% Parsimony tree
[treefilename, tempnames] =generate_parsimony_tree_rename(calls_for_tree, TreeSampleNames, refgenome);
relabeltree(treefilename,newtreefilename,tempnames, TreeSampleNames);

fprintf(1,'\nDone with tree\n');

 
%% Make a tree for each SNP location

if exist('tree_counting','dir')
    rmdir('tree_counting','s') 
end
mkdir('tree_counting')
cd('tree_counting')
 
fid=fopen('for_tree_labeling.csv','w');
fprintf(fid,'chr,pos');
for i=1:numel(TreeSampleNames)
    fprintf(fid,[',' TreeSampleNames{i}]);
end
for i=1:numel(goodpos)
    fprintf(fid,['\n' num2str(contig_positions(positions_to_show_in_table(i),1)) ',' num2str(contig_positions(positions_to_show_in_table(i),2))]);
    for j=1:size(calls_for_tree,2)
        fprintf(fid,[',' calls_for_tree(i,j)]);
    end
end
fprintf(fid,'\n');
fclose(fid);
eval(['! python2.7 '  SCRIPTSDIRECTORY '/countMutations.py ../' treefilename ' for_tree_labeling.csv'])
cd('..')


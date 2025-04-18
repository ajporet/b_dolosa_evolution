function build_mutation_table_master_fast

%% Parameters  -- IMPORTANT TO ADJUST !!!!
% Much of SCRIPTSDIRECTORYthis will need to vary according to your reference genome and
% coverage

%Input files that you need: csv sample names, strain.pileup, vcf_to_quals
%(strain.vcf file), variant.vcf (generate_positions)

folder=char(pwd);
cd('..');
masterdir=char(pwd);

SCRIPTSDIRECTORY = '/scratch/mit_lieberman/scripts';
REF_GENOME_DIRECTORY = '/scratch/mit_lieberman/reference_genomes/Bdolosa_new';
other_positions_to_consider={};

looseFQmax=-30;

loose_parameters=struct('minorfreqthreshold',.05, 'minreads_perstrand',1,...
    'maxreads_perstrand_percentile', 100,'minreads_perstrand_per_allele',2,...
    'min_bq',15,'min_mq',30, 'min_td', 10, 'max_td',90, 'max_sbp', 10,...
    'max_bqp', 255,'max_tdp',255, 'max_percent_ends', .90, 'max_percent_indels', .90, 'min_control_MAF', .01);

% Most threshold checks are strictly > or strictly <
% and many of these filters are not used in this instance

% min_td and max_td are not symmetrical relative to read length of 100
% because some reads were trimmed prior to alignment
% maxreads_perstrand_percentile is which threshold in list .01:.01:1 ...
%    e.g. 98 is 98 percentile of covered positions


%% Enviornment set up -- IMPORTANT TO DO !!!!
% Much of this will need to vary according to your computing system


RUN_ON_COMPUTING_CLUSTER = 1;  %set to zero if you aren't running on a cluster --  may take a very long time

TEMPORARYFOLDER =  [masterdir '/temp/'];


% Options for the cluster, used for HMS Orchestra LSF
% Please adjust this and send_jobs_to_cluster if on other computing cluster
jobsubmitoptions_short = '1:00:00';
jobsubmitoptions_long='12:00:00';

if ~exist(TEMPORARYFOLDER,'dir')
    mkdir(TEMPORARYFOLDER);
end

fprintf(['Usings scripts directory: ' SCRIPTSDIRECTORY  '\n']);
fprintf(1,['Genome : ' REF_GENOME_DIRECTORY]);

path(SCRIPTSDIRECTORY,path);


%% Read in csv files

cd(folder)

save('for_matlab_scripts','SCRIPTSDIRECTORY','REF_GENOME_DIRECTORY','TEMPORARYFOLDER');

SampleInfo = read_sample_names ;

SampleDirs = {} ;
SampleNames={SampleInfo(:).Sample}';
fprintf('\nNumber of lines in sampleNames: %i\n',numel(SampleNames));

for i=1:numel(SampleNames)
    SampleDirs{i} = [SampleInfo(i).ExperimentFolder '/' SampleInfo(i).AlignmentFolder ] ;
    %fprintf('\nDiversity.mat file: %i\n', [SampleDirs{i} '/diversity.mat']); %georgia's addition

end

in_outgroup=zeros(numel(SampleNames),1);
if isfield(SampleInfo,'Outgroup')
    in_outgroup=[SampleInfo(:).Outgroup];
end
    
fprintf('\nNumber of samples read from csv file: %i\n', numel(SampleNames));


fprintf('\nReading reference genome...\n');
[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REF_GENOME_DIRECTORY);



%% Summarize pileup files in an easily accesible MATLAB structure


fprintf(1,'Create diversity.mat and Create quals.mat for each sample... \n') ;

cmds = {} ;
for i=1:length(SampleDirs)
    if ~exist([SampleDirs{i} '/diversity.mat'],'file') ;
        fprintf(1,'Need to create diversity.mat for %i sample... \n',i) ; %Georgia's addition

        if RUN_ON_COMPUTING_CLUSTER == 1
            cmds{end+1} = ['matlab -r "path(' char(39) SCRIPTSDIRECTORY char(39) ',path); pileup_to_diversity_matrix_light(' char(39) SampleDirs{i} char(39) ')"'];
        else
            pileup_to_diversity_matrix_light(SampleDirs{i});
        end
    end
end

send_jobs_to_cluster(cmds,jobsubmitoptions_long,RUN_ON_COMPUTING_CLUSTER,{'.'});

%STOP

%% Summarize quality scores from vcf file in an easily accesible MATLAB structure


for i=1:length(SampleDirs)
    if ~exist([SampleDirs{i} '/quals.mat'],'file') ;
        fprintf(1,'Need to create quals.mat for %i sample... \n',i) ; %Georgia's addition

        if RUN_ON_COMPUTING_CLUSTER == 1
            cmds{end+1} = ['matlab -r "path(' char(39) SCRIPTSDIRECTORY char(39) ',path); vcf_to_quals(' char(39) SampleDirs{i} char(39) ')"'];
        else
            vcf_to_quals(SampleDirs{i});
        end
    end
    
end

%run the things
send_jobs_to_cluster(cmds,jobsubmitoptions_long,RUN_ON_COMPUTING_CLUSTER,{'.'});


%% Make a list of candiate positions

fprintf(1,'\n\nFinding positions with at least 1 fixed mutation...\n');
cp = generate_positions_aro({SampleDirs{~in_outgroup}}, {SampleNames{~in_outgroup}}', looseFQmax, RUN_ON_COMPUTING_CLUSTER, jobsubmitoptions_short, TEMPORARYFOLDER,SCRIPTSDIRECTORY);
fprintf(1,['Found ' num2str(length(cp)) ' positions where provided vcf called a fixed variant in at least one in-group sample with FQ score < ' num2str(looseFQmax) '\n']) ;

dp =[];
% fprintf(1,'\nFinding single nucleotide positions with within-sample polymorphisms...\n');
% [dp, coveragethresholds] = find_diverse_positions_no_control(loose_parameters, {SampleDirs{~in_outgroup}}, {SampleNames{~in_outgroup}}', RUN_ON_COMPUTING_CLUSTER, jobsubmitoptions_short,TEMPORARYFOLDER,SCRIPTSDIRECTORY);
% fprintf(1,'Found %i positions with within-sample polymorphism that meets loose parameters in at least 1 in-group sample \n',length(dp)) ;


op=[];
if numel(other_positions_to_consider)>0
    for i=1:numel(other_positions_to_consider)
        other=load(other_positions_to_consider{i});
        op=[op; chrpos2index(other.pos_list, ChrStarts)];
    end
end
fprintf(1,'\nConsidering %g positions previously specified \n',length(op)) ;

%% Combine positions

allp=unique([dp; op; cp;]);
p=sort(allp);
p=p(p>0); %remove any 0s
positions=p2chrpos(p,ChrStarts);


%% Get counts and mutgenvcf for snp positions

fprintf(1,'\nAcquiring detailed information at each potential position...\n');

fprintf(1,'Gathering quality scores at each candidate position ...\n');
Quals=zeros(length(p), length(SampleDirs),'int16');
for i=1:length(SampleDirs)
    fprintf(1,'Loading quals matrix for sample: %g  \n',i) ;
    load([SampleDirs{i} '/quals.mat']);
    Quals(:,i)=quals(p);
end


fprintf(1,'Gathering detailed read information from MATLAB counts matrix at each candidate position...\n');
counts=zeros(8, length(p), length(SampleDirs),'uint16');
for i=1:length(SampleDirs)
    fprintf(1,'Loading counts matrix for sample: %g  \n',i) ;
    load([SampleDirs{i} '/diversity.mat']);
    counts(:,:,i)=data(1:8,p);
end


fprintf(1,'Getting all the coverage information...\n');
[all_coverage_per_bp, ~] = get_all_coverage(SampleInfo, GenomeLength);

%% Save

save('candidate_mutation_table', 'SampleNames', 'p', 'counts', 'Quals', 'in_outgroup', '-v7.3') ;
save('coveragematrix', 'all_coverage_per_bp', '-v7.3'); 

end

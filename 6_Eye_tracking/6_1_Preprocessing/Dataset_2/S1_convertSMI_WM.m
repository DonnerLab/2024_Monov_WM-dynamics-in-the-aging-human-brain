% Convert from SMI .idf files to .asc files
clear all 
loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupilWM/1.converted/';
savepath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_preproc/yhc_preproc/2.asc_converted/';
behav_path = '/mnt/homes/home028/gmonov/SCZ/Data/working_memory/'; 

addpath '/mnt/homes/home030/aarazi/fieldtrip-20201009'  % tell Matlab where FieldTrip is 
ft_defaults

load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/yhc.mat']);
allsubj = yhc;
sess = 2; % WM performed in second session 

% Loop through individual files
for subj = 11:17 %18:length(allsubj) 
    filz = dir([[behav_path,allsubj{subj},filesep,'S2/Behaviour/'] '*.mat']); % find number of blocks 
    nblocks = length(filz); 
    for s = 1:length(sess)
        for b = 1:nblocks

            fprintf('Subj %s, session % d, block %d\n',subj,s,b)

            % Define file names
            ascFile = [loadpath,allsubj{subj},'_',num2str(sess(s)),'_',num2str(b),'_WM Samples.txt'];
            matFile = [savepath,allsubj{subj},'_',num2str(sess(s)),'_',num2str(b),'_WM.mat'];

            % Read the asc file into matlab & save
            if ~exist(matFile)  % ignores if file exists already
                tic
                asc = read_SMI_asc(ascFile);
                toc
                save(matFile,'asc')
            end
        end
    end

end 


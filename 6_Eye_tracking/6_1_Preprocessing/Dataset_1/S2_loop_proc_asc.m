clear

loadpath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/asc_files/';
savepath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/matConverted/';

addpath '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/'
addpath '/mnt/homes/home030/aarazi/fieldtrip-20201009'  % tell Matlab where FieldTrip is 
ft_defaults

%% define subjects 
% Load subject IDs
load('/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/behav_subj.mat');
hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final; 
allsubj = horzcat (hc,pat,cog_def);

% Loop through subjects
for subj = 1:length(allsubj)
    
    sdirs = dir([loadpath,allsubj{subj}]);
    sdirs = sdirs(3:end);
    
    for s = 1:1 %
                
        sesspath = [loadpath,allsubj{subj},filesep,sdirs(s+2).name,filesep];
        bnames = dir([sesspath,'*.asc']);
        subj_dir = [savepath,allsubj{subj},'/'];
        if ~exist(subj_dir,'dir'), mkdir(subj_dir); end
        
        for b = 1:length(bnames)
            
            fprintf('Subj %s, %s\n',allsubj{subj},bnames(b).name)
            
            % Define file names
            ascFile = [sesspath,filesep,bnames(b).name];
            matFile = [subj_dir,bnames(b).name(1:end-4),'.mat'];
            
            
                % Read the asc file into matlab
                asc = read_eyelink_ascNK_AU(ascFile);
                
                % Create data and events structures, and matrices of blink/saccade times
                data = asc2dat(asc); clear asc
                
                % Create refined event structure
                data = MCI_WM_trialfun(data,b);
                
                % Save output
               save(matFile,'data')
         
            clear data 
        end
    end
    clearvars -except allsubj subj loadpath savepath
end
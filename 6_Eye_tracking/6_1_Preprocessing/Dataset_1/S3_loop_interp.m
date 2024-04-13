% Interpolate pupil time-series

clear; close all; clc

user_input = false;

loadpath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/matConverted/';
savepath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/interpolated/';
 if ~exist(savepath,'dir'), mkdir(savepath); end
addpath '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/'


%% define subjects 
% Load subject IDs
load('/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/behav_subj.mat');
hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final; 
allsubj = horzcat (hc,pat,cog_def);
allsubj = {'15'}; 
% Loop through subjects
for_checking = {};
% Loop through individual files
for subj = 1:length(allsubj)
    subj_dir = [savepath,allsubj{subj},'/'];
   if ~exist(subj_dir,'dir'), mkdir(subj_dir); end
    bnames = dir([loadpath,allsubj{subj},filesep,'*.mat']);
    
    for b = 1:length(bnames)
        fprintf('Subj %s, %s\n',allsubj{subj},bnames(b).name)  % print progress
        
        % Load current asc conversion
        inFile = [loadpath,allsubj{subj},filesep,bnames(b).name];
        outFile = [subj_dir,bnames(b).name(1:end-4),'_interp.mat'];
        
       % if ~exist(outFile)
            load(inFile)
            
            % Run interpolation routine & add results to data structure
            [newpupil, newXgaze, newYgaze, newblinksmp, badsmp, h] = blink_interpolatePM(data,bnames(b).name);
            data.pupil = newpupil;
            data.Xgaze = newXgaze;
            data.Ygaze = newYgaze;
            data.newblinksmp = newblinksmp;
            data.badsmp = badsmp;
            
            pause(2)  % wait for a second to give user a better look when user input is specified to false
            
            if user_input
                % Get user input on interpolation quality if specified
                instr = input('Happy with interpolation? (y/n) ...','s');
                
                if strcmp(instr,'y')
                    save(outFile,'data')
                else
                    for_checking{end+1} = bnames(b).name(1:end-4);
                end
            else % Otherwise just save
                save(outFile,'data')
                savefig(h, [savepath,'plots',filesep,bnames(b).name(1:end-4),'.fig'],'compact')
               close all
            end
            close all
       % end
    end
    clearvars -except *path *subj for_checking user_input
end

% Reminder of files to be checked further
if ~isempty(for_checking)
    fprintf('\nThe following files need to be inspected more closely:\n')
    for f = 1:length(for_checking)
        fprintf('%s\n',for_checking{f})
    end
end
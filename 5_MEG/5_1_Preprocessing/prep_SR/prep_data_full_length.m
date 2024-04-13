%% Concatenate full length 3 and 9s trials for decoding analysis

% Script for preparing MCI WM task data for source reconstruction 
clear all
close all
addpath '/mnt/homes/home028/gmonov/meg_analysis/task_analysis/'; 

addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/comps_rejection/';
        
        comps2rej = readtable ('comps_rejection_wm','Range','A1:L64'); 
        IDs = comps2rej{1:64,1};
        
for idx_file=1:length(IDs)
    ID = IDs{idx_file,1}; 
    if strcmp(ID(end),'1') %only 1st session 
    data_new_9s = [];
    data_new_3s = [];
    loadpath =['/mnt/homes/home028/gmonov/meg_analysis/sr_separated_data/', ID] %only 1st session 
    load([loadpath, filesep, ID '_1.mat']) %load 1s data 
    data_1s = data; 
    load([loadpath, filesep, ID '_3.mat']) %load 3s data 
    data_3s = data; 
    try % try whether 9s trials are present 
    load([loadpath, filesep, ID '_9.mat']) %load 9s data 
    data_9s = data; 
    catch, data_9s = []; end 
    
    % Add first second to 3s and 9s delay 
    
    data_new_3s = data_3s;
    data_new_3s.trial(:,:,data_3s.time<=max(data_1s.time)) = [];%delete overlap with 1st second 
    data_new_3s.trial=cat(3,data_1s.trial(data_1s.trialinfo(:,4)==3 & data_1s.trialinfo(:,24)==0,:,:),data_new_3s.trial); % concatenate with first second 
    data_new_3s.time = unique([data_1s.time,data_3s.time]);
    if ~isempty(data_9s)
    data_new_9s = data_9s;
    data_new_9s.trial(:,:,data_9s.time<=max(data_1s.time)) = [];%delete overlap with 1st second 
    data_new_9s.trial=cat(3,data_1s.trial(data_1s.trialinfo(:,4)==9 & data_1s.trialinfo(:,24)==0,:,:),data_new_9s.trial); % concatenate with first second 
    data_new_9s.time = unique([data_1s.time,data_9s.time]);
    
    data_new_3s.trial = cat(1,data_new_3s.trial,data_new_9s.trial(:,:,data_new_9s.time<=max(data_new_3s.time)));
    data_new_3s.trialinfo = [data_new_3s.trialinfo;data_new_9s.trialinfo];
    data_new_3s.sampleinfo = [data_new_3s.sampleinfo;data_new_9s.sampleinfo];
    end 
    % Add first 3s of the 9s trials to the 3s delay trial data 
    
    % Save new full-length data sets 
    
       mat_name = ['/mnt/homes/home028/gmonov/meg_analysis/full_length_3_9s/' ID '/'];
    
        if 7==exist(mat_name,'dir')
            cd(mat_name)
        else
            mkdir(mat_name)
            cd(mat_name)
        end
        
         if ~isempty(data_new_3s) % Do not save when no 9s trials are there 
         data = data_new_3s;
         datastore = [ID '_3'];
         save(datastore,'data','-v7.3');
         clear data 
         end 
         if ~isempty(data_new_9s) % Do not save when no 9s trials are there 
         data = data_new_9s;
         datastore = [ID '_9'];
         save(datastore,'data','-v7.3');
         clear data 
         end 
    end   
end 
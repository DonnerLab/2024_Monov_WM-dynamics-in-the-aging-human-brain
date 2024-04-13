% 1. Script for preparing MCI WM task data for source reconstruction 

clear all
close all
addpath '/mnt/homes/home028/gmonov/meg_analysis/task_analysis/'; 

addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/comps_rejection/';
        
        comps2rej = readtable ('comps_rejection_wm','Range','A1:L64'); 
        IDs = comps2rej{1:64,1};
        
   for idx_file=1:length(IDs)
     ID = IDs{idx_file,1}; 
      % load cleaned postICA data 
     loadpath = ['/mnt/homes/home028/gmonov/meg_analysis/data_trials/', ID]
     load([loadpath, filesep, 'data_clean_postICA_' ID '.mat'])
     %Initialize structures to separate data 
     trials_1s = all_trials_cl;
     trials_1s.dimord = 'rpt_chan_time';
     trials_3s = all_trials_cl; 
     trials_3s.dimord = 'rpt_chan_time';
     trials_9s = all_trials_cl; 
     trials_9s.dimord = 'rpt_chan_time';
     
   % extract data 1st window (from -0.6s to +1.8s)
   % Save this data set 
     for l = 1:length(trials_1s.trialinfo(:,1)) %loop through every clean trials
         trials_1s.trial{1,l}(:,all_trials_cl.time{1,l}>1.8)=[]; 
         trials_1s.time{1,l}(all_trials_cl.time{1,l}>1.8)=[]; 
     end 
    
     
     % extract data: time window from +1.2s to +3.8s --> only 3s delay
     % duration trials
     trials_3s.trialinfo(all_trials_cl.trialinfo(:,4)~=3,:)=[];  %delete all trials that are not 3s delay trials 
     trials_3s.sampleinfo(all_trials_cl.trialinfo(:,4)~=3,:)=[]; 
     trials_3s.time(:,all_trials_cl.trialinfo(:,4)~=3)=[]; 
     trials_3s.trial(:,all_trials_cl.trialinfo(:,4)~=3)=[]; 
     trials_3s.eye_alltrials(all_trials_cl.trialinfo(:,4)~=3,:)=[]; 
     trials_3s.eye_trialtiming(all_trials_cl.trialinfo(:,4)~=3,:)=[]; 
     

     trials_3s.sampleinfo(trials_3s.trialinfo(:,24)>0,:)=[];  %delete all trials that had artifacts after the first second 
     trials_3s.time(:,trials_3s.trialinfo(:,24)>0)=[]; 
     trials_3s.trial(:,trials_3s.trialinfo(:,24)>0)=[]; 
     trials_3s.eye_alltrials(trials_3s.trialinfo(:,24)>0,:)=[]; 
     trials_3s.eye_trialtiming(trials_3s.trialinfo(:,24)>0,:)=[]; 
     trials_3s.trialinfo(trials_3s.trialinfo(:,24)>0,:)=[]; 
     
     for l = 1:length(trials_3s.trialinfo(:,1)) %loop through every clean 3s trial
         trials_3s.trial{1,l}(:,trials_3s.time{1,l}<1.2)=[]; %extract specific time windows of 3s trials 
         trials_3s.time{1,l}(trials_3s.time{1,l}<1.2)=[];   
         if max(trials_3s.time{1,l})>3.8 || length(trials_3s.trial{1,l})>1041
            disp error3s 
         end 
         
     end 
     
     % extract data: time window from +1.2s to +9.8s --> only 9s delay
     trials_9s.trialinfo(all_trials_cl.trialinfo(:,4)~=9,:)=[];  % delete all trials that are not 9s delay trials 
     trials_9s.sampleinfo(all_trials_cl.trialinfo(:,4)~=9,:)=[]; 
     trials_9s.time(:,all_trials_cl.trialinfo(:,4)~=9)=[]; 
     trials_9s.trial(:,all_trials_cl.trialinfo(:,4)~=9)=[]; 
     trials_9s.eye_alltrials(all_trials_cl.trialinfo(:,4)~=9,:)=[]; 
     trials_9s.eye_trialtiming(all_trials_cl.trialinfo(:,4)~=9,:)=[]; 
     

     trials_9s.sampleinfo(trials_9s.trialinfo(:,24)>0,:)=[];  % delete all trials that had artifacts after the first second 
     trials_9s.time(:,trials_9s.trialinfo(:,24)>0)=[]; 
     trials_9s.trial(:,trials_9s.trialinfo(:,24)>0)=[]; 
     trials_9s.eye_alltrials(trials_9s.trialinfo(:,24)>0,:)=[]; 
     trials_9s.eye_trialtiming(trials_9s.trialinfo(:,24)>0,:)=[]; 
     trials_9s.trialinfo(trials_9s.trialinfo(:,24)>0,:)=[]; 

     for l = 1:length(trials_9s.trialinfo(:,1)) % loop through every clean 9s trial
         trials_9s.trial{1,l}(:,trials_9s.time{1,l}<1.2)=[]; % extract specific time windows of 9s trials 
         trials_9s.time{1,l}(trials_9s.time{1,l}<1.2)=[];  
         if max(trials_9s.time{1,l})>9.8 || length(trials_9s.trial{1,l})>4161
            disp error9s
         end 
     end 
     
     % Save new data sets
     % Define path/folder to save cleaned data
  
    mat_name = ['/mnt/homes/home028/gmonov/meg_analysis/sr_separated_data_uncut/' ID '/'];
    
        if 7==exist(mat_name,'dir')
            cd(mat_name)
        else
            mkdir(mat_name)
            cd(mat_name)
        end

    datastore = ['1sDelay_' ID];
    save(datastore,'trials_1s','-v7.3');
    datastore = ['3sDelay_' ID];
    save(datastore,'trials_3s','-v7.3');
    datastore = ['9sDelay_' ID];
    save(datastore,'trials_9s','-v7.3');
    
    % Clean up
     clear trials_1s
     clear trials_3s
     clear trials_9s
     
   end 
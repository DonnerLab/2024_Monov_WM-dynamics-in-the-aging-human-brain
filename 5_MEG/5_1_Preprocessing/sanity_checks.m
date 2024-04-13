%% Checking if behavioral data and meg data match, make updated master file, mark any inconsistencies in this master file 
clear close all

% first load IDs/colors/file_info 
colorpath = '/mnt/homes/home028/gmonov/functions/colors/'
load([colorpath,'colors']);

addpath '/mnt/homes/home028/gmonov/meg_analysis/task_analysis/'; 

addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/comps_rejection/';
        
        comps2rej = readtable ('comps_rejection_wm','Range','A1:L64'); 
        IDs = comps2rej{1:64,1};
        final_checks={}; 
 for idx_file=1:length(IDs)
     ID = IDs{idx_file,1}; 
      % load cleaned postICA data 
     loadpath = ['/mnt/homes/home028/gmonov/meg_analysis/data_trials/', ID]
     load([loadpath, filesep, 'data_clean_postICA_' ID '.mat'])
     % load cleaned behavioral data 
     if strcmp(ID,'P03_1')
         loadpath =['/mnt/homes/home028/gmonov/behavioral_analysis/clean_allbehav_data/', ID(1:3),filesep,'S',ID(5)]
         load([loadpath, filesep, ID,'_clean_allbehav.mat'])
     else loadpath =['/mnt/homes/home028/gmonov/behavioral_analysis/clean_allbehav_data/', ID(1:2),filesep,'S',ID(4)]
          load([loadpath, filesep, ID,'_clean_allbehav.mat'])
     end 
     % delete first block 
     if strcmp(ID,'43_2')
         block_count = allbehav(:,11)
         allbehav(block_count==1,:)=[]; 
     end 
   
     % Check whether early block terminations other inconsistencies of trialcounts between meg and bhav data set were handled correctly 
     % If yes the same number of trials that were either marked with 0 or 1
     % have been initially extracted from the meg data sets 
     if numel(all_trials_cl.file_info.trialex(all_trials_cl.file_info.trialex==1))+numel(all_trials_cl.file_info.trialex(all_trials_cl.file_info.trialex==0))+numel(all_trials_cl.file_info.trialex(all_trials_cl.file_info.trialex==6)) == all_trials_cl.alltrials_count
         final_checks.termination = 0; 
     else final_checks.termination = 1; 
     end 
     
     % Check whether RT/wrong button artifacts are removed in the same way 
     for c=1:length(all_trials_cl.file_info.trialex(:,1)) %loop through blocks
        if strcmp(ID,'43_2')
            block_trials = allbehav(allbehav(:,11)==c+1,:); 
            else block_trials = allbehav(allbehav(:,11)==c,:); 
        end 
         clean_trials = 0; 
         for t=1:length(all_trials_cl.file_info.trialex(c,:))%loop through trials
            
             if all_trials_cl.file_info.trialex(c,t) == 1 || all_trials_cl.file_info.trialex(c,t) == 6 || all_trials_cl.file_info.trialex(c,t) == 3  
                 clean_trials(end+1,1) = t; 
             end 
         end 
         clean_trials(clean_trials==0)=[]; 
         if clean_trials == block_trials(:,10) %same trialno in cleaned behavioral data as marked in trialex?
            final_checks.trialex_check{c} = 0;
         else final_checks.trialex_check{c} = 1;
         end 
             
     end 
     
     for p=1:length(all_trials_cl.trialinfo(:,1)) %loop thorugh trials in cleaned meg data set 
         %Check if delay duration and length of the trials indicated via
         %sample times match!
         trial_dur = all_trials_cl.sampleinfo(p,2)-all_trials_cl.sampleinfo(p,1)
         if all_trials_cl.trialinfo(p,4) == 1 && ~isempty(intersect(trial_dur,950:970)) % allow for +/-10 samples tolerance 
          final_checks.delayduration_check{p} = 0;  
         elseif all_trials_cl.trialinfo(p,4) == 3 && ~isempty(intersect(trial_dur,1750:1770))
          final_checks.delayduration_check{p} = 0;  
         elseif all_trials_cl.trialinfo(p,4) == 9 && ~isempty(intersect(trial_dur,4150:4170))
          final_checks.delayduration_check{p} = 0; 
         else final_checks.delayduration_check{p} = 1;  
         end 
         
         %Check if saved trigger times match as well 
         if all_trials_cl.trialinfo(p,16)-all_trials_cl.trialinfo(p,14)+(0.6*400+0.3*400) == trial_dur %Probe onset - mem onset +/-0.6s
          final_checks.delayduration_trigger{p} = 0; 
         else final_checks.delayduration_trigger{p} = 1;  
         end 
         %Check whether saved trigger times exceeds actual trial timings
         %find last trigger 
         last_trigger = all_trials_cl.trialinfo(p,16);
         if last_trigger+0.3*400 == all_trials_cl.sampleinfo(p,2)
            final_checks.last_trigger{p} = 0;
         else final_checks.last_trigger{p} = 1;
         end 
         
        %Check whether trigger timing makes sense
        trigger_timings = all_trials_cl.trialinfo(p,14:end-1);
        trigger_timings(isnan(trigger_timings)) = []; 
        
        if length(trigger_timings)==6 %There should be 6 triggers related to each trial 
           final_checks.trigger_count{p} = 0;
        else final_checks.trigger_count{p} = 1;
        end 
        
        if trigger_timings == sort(trigger_timings) 
            final_checks.trigger_timings{p} = 0;
        elseif ismember(1,trigger_timings ~= sort(trigger_timings)) 
            trigger_timings_new = sort(trigger_timings); 
             swap = 4;
             trigger_timings_new([swap;swap+1]) = trigger_timings_new([swap+1 swap]);
             if trigger_timings_new == trigger_timings
             final_checks.trigger_timings{p} = 0;
             end 
        else final_checks.trigger_timings{p} = 1;
        end 
            
         %Check whether hand responses match (left/-1 = different, right/+1 =
         %same) 
         if all_trials_cl.trialinfo(p,6)==-1 && ~isnan(all_trials_cl.trialinfo(p,19)) && isnan(all_trials_cl.trialinfo(p,18)) && isnan(all_trials_cl.trialinfo(p,20))
            final_checks.response{p} = 0;
         elseif all_trials_cl.trialinfo(p,6)==1 && ~isnan(all_trials_cl.trialinfo(p,18)) && isnan(all_trials_cl.trialinfo(p,19)) && isnan(all_trials_cl.trialinfo(p,20))
            final_checks.response{p} = 0;
         else final_checks.response{p} = 1;
         end  
         %Check whether feedback cues match 
        if all_trials_cl.trialinfo(p,7)==1 && ~isnan(all_trials_cl.trialinfo(p,21)) && isnan(all_trials_cl.trialinfo(p,22)) && isnan(all_trials_cl.trialinfo(p,23))
            final_checks.feedback{p} = 0;
         elseif all_trials_cl.trialinfo(p,7)==0 && ~isnan(all_trials_cl.trialinfo(p,22)) && isnan(all_trials_cl.trialinfo(p,21)) && isnan(all_trials_cl.trialinfo(p,23))
            final_checks.feedback{p} = 0;
         else final_checks.feedback{p} = 1;
         end  
            
     end 
         %Check whether all wrong button presses/bad feedback are gone 
         if ~ismember(99,all_trials_cl.trialinfo(:,6))%wrong button in behavioral file 
             final_checks.wrongbutton{1} = 0;
         else final_checks.wrongbutton{1} = 1;
         end 
         if isnan(all_trials_cl.trialinfo(1:end,20))%Response bad
             final_checks.wrongbutton{2} = 0;
         else final_checks.wrongbutton{2} = 1;
         end 
         if isnan(all_trials_cl.trialinfo(1:end,23))%Feedback bad
             final_checks.wrongbutton{3} = 0;
         else final_checks.wrongbutton{3} = 1;
         end 
                  %Check whether all RT/button press failures have been
                  %removed properly 
         if ~ismember(0,all_trials_cl.trialinfo(:,12))%marked with 0 due to RT /wrong button in trialex should be gone  
             final_checks.RT{1} = 0;
         else final_checks.RT{1} = 1;
         end 
     
        %CREATE NEW MASTERFILE  
        file_info=all_trials_cl.file_info;
        if strcmp(ID,'43_2') %load cleaned behavior again with 1st block 
          loadpath =['/mnt/homes/home028/gmonov/behavioral_analysis/clean_allbehav_data/', ID(1:2),filesep,'S',ID(4)]
          load([loadpath, filesep, ID,'_clean_allbehav.mat'])
        end 
     
        %Add cleaned behavior 
        file_info.clean_behavior=allbehav;
        % Add info about subjects with poor performance --> potential
        % exclusion 
        if mean(allbehav(:,6))<0.6 %accuracy lower than 60% 
           file_info.low_accuracy=1; 
        else file_info.low_accuracy=0;
        end
        %How many trials in meg data set? 
        file_info.meg_all_trials = all_trials_cl.alltrials_count; 
        %How many completely clean trials in meg data set? 
        artifactual_trials = all_trials_cl.trialinfo(:,24); 
        file_info.clean_trials=length(artifactual_trials(artifactual_trials==0));
        %How many trials that are clean at least up to 1s delay duration 
        file_info.clean_1s=length(artifactual_trials);
        
        if ismember(6,file_info.trialex)
           file_info.cog_art = 1
        else file_info.cog_art = 0
        end 
    cd '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/metal_artifacts/'
    exc_metal = isfile(strcat(ID,'_metal.mat')); 

    if strcmp(ID,'46_1') % a lot of artifacts as well but not reaching the threshold of 60% of all samples, however no clean trials left when rigorously correcting for metal artifacts  
        exc_metal = 1; 
    end
    file_info.metal_artifacts=exc_metal;
    
    %Check whether saved ocular artifacts are in correct alignement 
    
    
    %Indicate whether final check shows any irregularites 
    if  sum(final_checks.termination) == 0 &&...
        sum(cell2mat(final_checks.trialex_check)) == 0 &&...
        sum(cell2mat(final_checks.delayduration_check)) == 0 &&...
        sum(cell2mat(final_checks.delayduration_trigger)) == 0 &&...
        sum(cell2mat(final_checks.trigger_timings)) == 0 &&...
        sum(cell2mat(final_checks.last_trigger)) == 0 &&...
        sum(cell2mat(final_checks.trigger_count)) == 0 &&...
        sum(cell2mat(final_checks.response)) == 0 &&...
        sum(cell2mat(final_checks.wrongbutton)) == 0 &&...
        sum(cell2mat(final_checks.feedback)) == 0 &&...
        sum(cell2mat(final_checks.RT)) == 0
         file_info.final_checks = 0; 
    else file_info.final_checks = 1;
    end 
    file_info.final_check_specifics = final_checks; 
    file_info.ID = ID; 
    file_info.artifacts.metal = all_trials_cl.metal_count; 
    file_info.artifacts.muscle = all_trials_cl.muscle_count; 
    file_info.artifacts.RTbuttonpress = all_trials_cl.cog_count; 
    file_info.artifacts.jump = all_trials_cl.jump_count; 
    file_info.artifacts.headmovement = all_trials_cl.HM_count; 
    file_info.artifacts.ocular = all_trials_cl.eye_trialtiming; 
    file_info.meg_sampleinfo = all_trials_cl.sampleinfo; 
    file_info.meg_trialinfo = all_trials_cl.trialinfo; 
    file_info.blockBounds = all_trials_cl.blockBounds; 
     % load cleaned postICA data for resting state
     loadpath =['/mnt/homes/home028/gmonov/meg_analysis/data_restingstate/', ID]
     load([loadpath, filesep, 'rs_data_clean_postICA_' ID '.mat'])
     rs_duration  = length(data_cl.trial{1,1}(1,:))/400; %bring back to seconds 
     edge_count  = length(data_cl.edges(data_cl.edges==1)); 
    file_info.restingstate.duration = rs_duration; 
    file_info.restingstate.edgecount = edge_count; 
     % Add info how things change with z-threshold=10 for muscle artifacts 
     loadpath =['/mnt/homes/home028/gmonov/meg_analysis/data_trials_10/', ID]
     load([loadpath, filesep, 'data_clean_postICA_' ID '.mat'])
     artifactual_trials = all_trials_cl.trialinfo(:,24); 
    file_info.thres10.clean1s = length(artifactual_trials);
    file_info.thres10.clean_all = length(artifactual_trials(artifactual_trials==0));
    file_info.thres10.muscletrials = all_trials_cl.muscle_count; 
    % Merge all file_infos together in one overview file 
    datasets_overview(idx_file)=file_info;   
   clear file_info
   clear final_checks
   clear rs_duration
   clear edge_count

 end 
     
save(['/mnt/homes/home028/gmonov/meg_analysis/datasets_overview1.mat'], 'datasets_overview')

     

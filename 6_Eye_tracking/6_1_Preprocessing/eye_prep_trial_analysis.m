% Prepare eye-tracking data for decoding analysis (performed in Python) 
% Extraxt X- and Y-gaze positions that will serve as decoding features 
% Run checks to verify alignment with behavioural data 
% Sampling rate can be determined (new_sr). 
% Native sampling rate of Dataset 1: 1000 Hz; Dataset 2: 250 Hz
% Gina Monov, UKE, January 2024 

%% Define subject IDs/paths etc. 

clear all 
close all

load('/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat');
load('/Users/ginamonov/Servers/mountpoint1/meg_analysis/datasets_overview.mat') % trial exclusion
hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final;
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']); %YHC IDs
mci_subj = horzcat (hc,pat,cog_def);
allsubj = horzcat (mci_subj,yhc); 

few_trials_subj = {}; %save info which subjects have fewer than 60 trials for analysis 
%Exclude subjects due to extremely bad data quality that made preprocessing
%impossible / missing data files 
% allsubj = mci_subj; 
allsubj = allsubj(~ismember(allsubj,{'15','37'})); 

loadpath = '/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/interpolated/';
loadpath_yhc = '/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/eye_preproc/yhc_preproc/4.interpolated/';
behav_path_mci = '/Users/ginamonov/Servers/mountpoint1/behav_eye_data/';
behav_path_yhc = '/Users/ginamonov/Servers/mountpoint1/SCZ/Data/working_memory/';

savepath = '/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/data4decoding/';
%% Loop over subjects 
error_subj = {}; 
for subj = 1:length(allsubj)
   % find out trial count for each subject 
   if ismember(allsubj{subj},mci_subj)
    for iii = 1:length(datasets_overview) 
        if datasets_overview(iii).sess == 1 && strcmp(datasets_overview(iii).subj,allsubj{subj})
            trialex = datasets_overview(iii).trialex;  
        end 
    end 
    elseif ismember(allsubj{subj},yhc)
       %Create trialex for yhc as saved in datasets_overview for subjects
       %of MCI study 
       allbehav = [];
       trialex = []; 
      files = dir([[behav_path_yhc,allsubj{subj},filesep,'S2',filesep,'Behaviour',filesep] '*mat']);

      for f = 1:length(files) 
          load([files(f).folder,filesep,files(f).name])
          Behav(1:end,10) = f; 
          allbehav = [allbehav;Behav]; 
          
      end 
         allbehav(:,11) = ones(length(allbehav),1); 
         R = allbehav(:,7);
         allbehav(R>=(mean(R)+4*std(R)),11) = 0;
         R = allbehav(:,7);
         allbehav(R<=0.2,11) =0;
         choice = allbehav(:,5);
         allbehav(choice==99,11)=0;
         
         trialex = ones(length(files),63); 
         
         for f = 1:length(files) 
             if length(allbehav(allbehav(:,10)==f)) ~= 63
                trialex(f,1:length(allbehav(allbehav(:,10)==f))) = allbehav(allbehav(:,10) == f,11); 
                trialex(f,length(allbehav(allbehav(:,10)==f))+1:63) = nan; 
             else 
             trialex(f,:) = allbehav(allbehav(:,10) == f,11); 
             end 
         end 

   end  
         clear Behav allbehav R choice files
    
    
   data_all_blocks.stim_loc = [];  % initialize data of all blocks 
   data_all_blocks.events = [];  
   data_all_blocks.Ygaze = []; 
   data_all_blocks.Xgaze = []; 
   data_all_blocks.eccentricity = []; 
   data_all_blocks.polar_angle = []; 
   data_all_blocks.times = []; 
   data_all_blocks.trials = []; 
   stim_loc = []; % initialize stimulus location across all blocks 
    % check how many blocks exist in the folder 
    if ismember(allsubj{subj},mci_subj)
         blocks = dir([[loadpath,allsubj{subj},filesep] '*.mat']);
    elseif ismember(allsubj{subj},yhc)  
        blocks = dir([[loadpath_yhc,allsubj{subj}] '*.mat']);
    end 
    
    if strcmp(allsubj{subj},'01') || strcmp(allsubj{subj},'02') 
        blocks(1) = []; %Exclude first block because it was interrupted after a few trials / performance issues 
    end 
    
    subj_path = [savepath,allsubj{subj},'/']; 
    if ~exist(subj_path,'dir'), mkdir(subj_path); end 
    
    for b = 1:length(blocks)   % loop over the blocks 
        block_trialex = trialex(b,:)';
        if ismember(allsubj{subj},mci_subj)
            load([loadpath,allsubj{subj},filesep,blocks(b).name]); % load interpolated eye tracking data 
        elseif ismember(allsubj{subj},yhc)
            load([loadpath_yhc,filesep,blocks(b).name]); % load interpolated eye tracking data 
        end 
         % load corresponding behavioral file 
         if ismember(allsubj{subj},mci_subj)
              %Exceptions
              if strcmp(allsubj{subj},'01') || strcmp(allsubj{subj},'02') || strcmp(allsubj{subj},'27') %27: only second eyetracking block with usable data 
                load([behav_path_mci,allsubj{subj},'/S1/Behaviour/',allsubj{subj},'_1_',num2str(b+1)]); % exlude the first block 
              elseif strcmp(allsubj{subj},'25')
                  load([behav_path_mci,allsubj{subj},'/S1/Behaviour/02C_1_',num2str(b+1)]); % exlude the first block, correct ID
              elseif strcmp(allsubj{subj},'24')
                  load([behav_path_mci,allsubj{subj},'/S1/Behaviour/01C_1_',num2str(b)]); % correct ID
              else  load([behav_path_mci,allsubj{subj},'/S1/Behaviour/',allsubj{subj},'_1_',num2str(b)]); % load behavioral data 
              end 
         elseif ismember(allsubj{subj},yhc) %load younger subjects' data from other study 
               load([behav_path_yhc,allsubj{subj},filesep,'S2',filesep,'Behaviour/',allsubj{subj},'_2_',num2str(b)])
         end 
    %% run a few checks, mark and remove excluded trials and save the sample stimulus location 
    if strcmp(allsubj{subj},'17') && b == 3 || strcmp(allsubj{subj},'22') && b == 3 || strcmp(allsubj{subj},'25') && b == 3 || strcmp(allsubj{subj},'30') && b == 3  || strcmp(allsubj{subj},'26') && b == 3
        data.event(end,:) = []; % exclude only the last trial because it was not saved in behavioral file when the block was terminated earlier 
    end 
       data.event(:,end+1) = Behav(:,1); % add sample stimulus location
       acc = Behav(:,6); 
       acc(acc==0) = 2; 
       acc(isnan(acc)) = 3; 
       
       if sum(acc == data.event(:,9)) ~= length(acc) % Check whether behavioral and eye tracking data are aligned
          error_subj = horzcat(error_subj,allsubj{subj}); % save subjects for which this is not the case to fix this 
       end 
       % exclude trials as for the behavioral analysis
       data.event(block_trialex == 0,:) = []; 
       clear Behav acc block_trialex
       
   
    %% Rearrange data (y-gaze and x-gaze) for trial-based analysis (sample onset to 1s delay, i.e. 1.5 s after ~1500 samples (rate 1000 Hz) 
    % & Mark bad samples for quality-based trial exclusion
    
       old_sr = data.fsample; % original sampling rate
       new_sr = 250; % new sampling rate; choose either 100, 250 or 1000 Hz 
 
       bad_samples_idx = {}; 
       bdsmp = unique(data.badsmp); 
       exc_prop = 0.6; % proportion threshold of bad samples for trial exclusion 
       exc_cons = 0.3*data.fsample; % threshold for acceptable consecutive artifacts 
       edges = [0.1,0.1]; % seconds to increase the trial interval to extract in the first step 
       edge_samples = edges * old_sr; % compute how many samples the edges will be denepnding on the sampling rate of the original data 
       ds_window = 0.012; % averaging window +/- time point of interest for downsampling (12 ms to ensure this exact window can be captured in 1000 Hz and 250 Hz recordings)
       
       times2extract = 0-edges(1):1/old_sr:1.5+edges(2); 
       
       trial_exc = zeros(length(data.event(:,1)),1);
       reu = []; 
       for t = 1:length(data.event(:,1))
           mem_time = data.event(t,2); 
           tr_samples = mem_time:mem_time+(1.5*data.fsample)-1;
           pad = [tr_samples(1)-edge_samples(1);tr_samples(end)+edge_samples(2)]; % add some data for baseline period and after end of 1s delay
           
           bad(t,:) = ismember(tr_samples,bdsmp); %find all samples within a trial that had to be interpolated 
           bad_samples_idx{1,t} = tr_samples(ismember(tr_samples,bdsmp)); % save affected samples 
           bad_samples_idx{2,t} = diff(bad_samples_idx{1,t}); % mark consecutive bad samples (diff == 1) 
           n_bad(t) = sum(bad(t,:)); 
           
           if n_bad(t) > 1.5*data.fsample*exc_prop % if total amount of interpolated samples exceeds predefined proportion, mark trial for exclusion 
               trial_exc(t) = 1; 
           end 
           new_c = 0; % counter for consecutuve artifactual samples 
           cons_art = []; 
           % Mark trials with prolonged artifacts for exclusion 
           if ~isempty(bad_samples_idx{1,t}) 
              for eu = 1:length(bad_samples_idx{2,t})  
                  if bad_samples_idx{2,t}(eu) == 1
                      new_c = new_c+1; 
                  elseif bad_samples_idx{2,t}(eu) ~= 1
                      cons_art(end+1,1) = new_c; % save number of samples in a consecutive artifact window 
                      new_c = 0; % reset 
                  end 
                  if eu == length(bad_samples_idx{2,t})  && bad_samples_idx{2,t}(eu) == 1 % save artifact ending 
                     cons_art(end+1,1) = new_c;  
                  end 
              end 
              if isempty(cons_art) && new_c ~= 0
                  cons_art = new_c; 
              end 
              
              if sum(cons_art>=exc_cons) >= 1 %find out whether there is a prolonged artifact
                  trial_exc(t) = 1;
        
              end 
           end 
               
           n_bad(t) = sum(bad(t,:)); 
           X_gaze(t,:) = data.Xgaze(pad(1,1):pad(2,1)); 
           Y_gaze(t,:) = data.Ygaze(pad(1,1):pad(2,1));
           if std(X_gaze(t,:)) < 10^-6 % find trials where corneal reflection was lost entirely 
                trial_exc(t) = 1; 
                reu(end+1,1) = t;
           end

           clear mem_time tr_samples
       end 
       
      
       times = data.times(1:1.5/(1/old_sr)); 
       clear t
       if sum(trial_exc) ~= 0
         X_gaze(trial_exc==1,:) = []; 
         Y_gaze(trial_exc==1,:) = []; 
         data.event(trial_exc==1,:) = []; 
       end 
    %% Downsampling 
   if new_sr == 100
    times_100 = -0.05:1/new_sr:1.5+0.05; 
    dps = length(times_100); 
    
    X_gaze_100 = zeros(length(data.event(:,1)),dps); 
    Y_gaze_100 = zeros(length(data.event(:,1)),dps); 
    
     times2extract = round(times2extract,3); 
     times_100 = round(times_100,2); 
     
    for t = 1:length(data.event(:,1)) %loop over trials 
       for rs = 1:dps % loop over (downsampled) time points
            [minValue,idx] = min(abs(times2extract-times_100(rs))); 
            test(rs) = minValue; 
            idx = idx-ds_window*old_sr:1:idx+ds_window*old_sr; % find indices 
            X_gaze_100(t,rs) = mean(X_gaze(t,idx)); 
            Y_gaze_100(t,rs) = mean(Y_gaze(t,idx)); 
            clear idx fr
        end 
    end 
    
%     Overwrite with downsampled versions 
      clear Y_gaze X_gaze
      Y_gaze = Y_gaze_100; 
      X_gaze = X_gaze_100; 
      times = times_100; 
   else 
      X_gaze_new = []; 
      Y_gaze_new = [];  
      ds = old_sr/new_sr; 
          for t = 1:length(data.event(:,1)) %loop over trials 
                 X_gaze_new(t,:) = downsample(X_gaze(t,:),ds); 
                 Y_gaze_new(t,:) = downsample(Y_gaze(t,:),ds); 
          end  
          times = downsample(times2extract(1:end-1),ds); 
          clear Y_gaze X_gaze
          Y_gaze = Y_gaze_new; 
          X_gaze = X_gaze_new; 
         
   end 
       
    %% Compute polar angle and eccentricty from y- and x-gaze position 
        if ismember(allsubj{subj},mci_subj) 
            res = [1920 1080];
            w =  51.9347; 
            dist = 60; 
        elseif ismember(allsubj{subj},yhc) 
            res = [1680 1050];
            w = 47.5; 
            dist = 52; 
        end 
        pa = []; 
        ecc = []; 
       for tr = 1:length(data.event(:,1)) %loop over trials 
          for ti = 1:length(times) %loop over times 
               [pa_x,ecc_x]=gaze_pixel2polar(X_gaze(tr,ti),Y_gaze(tr,ti),res,w,dist); 
               pa(tr,ti) = pa_x; 
               ecc(tr,ti) = ecc_x; 
               clear pa_x ecc_x
          end 
       end 
       
 
    %% Fusion of data across all blocks
      data_all_blocks.events = [data_all_blocks.events;data.event]; 
      data_all_blocks.Ygaze  = [data_all_blocks.Ygaze;Y_gaze]; 
      data_all_blocks.Xgaze  = [data_all_blocks.Xgaze;X_gaze];
      data_all_blocks.eccentricity  = [data_all_blocks.eccentricity;ecc];
      data_all_blocks.polar_angle  = [data_all_blocks.polar_angle;pa];
      data_all_blocks.times  = times; 
      data_all_blocks.stim_loc  = data_all_blocks.events(:,end); 
      
      % Clean up
      clearvars -except blocks b allsubj *subj *path yhc datasets_overview trialex data_all_blocks error_subj *path_yhc *path_mci times new_sr
    end 
    %% create .trials (trials x channels for decoding x time)
    data_all_blocks.trials(:,:,1) = data_all_blocks.eccentricity;
    data_all_blocks.trials(:,:,2) = data_all_blocks.polar_angle;
    data_all_blocks.trials = permute(data_all_blocks.trials,[1,3,2]); 
    trial_data = data_all_blocks.trials; 
    stim_loc = data_all_blocks.events(:,end); 
    data_all_blocks.stim_loc = stim_loc; 
    times = times'; 
    
    data_all_blocks.trials_x_y(:,:,1) = data_all_blocks.Xgaze;
    data_all_blocks.trials_x_y(:,:,2) = data_all_blocks.Ygaze;
    data_all_blocks.trials_x_y = permute(data_all_blocks.trials_x_y,[1,3,2]); 
    trial_data_x_y = data_all_blocks.trials_x_y; 
    
 %% SAVE 

     save([savepath,'stim_loc/',allsubj{subj},'_stim_loc_', num2str(new_sr), '.mat'],'stim_loc')  % save to-be-decoded variable (stimulus location) separately 
     save([savepath,'trials4decoding/',allsubj{subj},'_trial_data_', num2str(new_sr), '.mat'],'trial_data','times','new_sr','stim_loc')  % save trials x channels for decoding x time separately 
     save([savepath,'trials4decoding_x_y/',allsubj{subj},'_trial_data_x_y_', num2str(new_sr), '.mat'],'trial_data_x_y','times','new_sr','stim_loc')  % save trials x channels for decoding x time separately 
     save([subj_path,filesep,allsubj{subj},'_eye_trials_', num2str(new_sr), '.mat'],'data_all_blocks')  % save data_all_blocks 
     if length(stim_loc) < 60 
        few_trials_subj = horzcat(few_trials_subj,allsubj{subj}); 
     end 
% Clean up 
     clear data_all_blocks trialex trial_data times blocks b

end 

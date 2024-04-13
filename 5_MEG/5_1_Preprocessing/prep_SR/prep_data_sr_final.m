% 2. Script for preparing MCI WM task data for source reconstruction 
clear all
close all
addpath '/mnt/homes/home028/gmonov/meg_analysis/task_analysis/'; 

addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/comps_rejection/';
        
        comps2rej = readtable ('comps_rejection_wm','Range','A1:L64'); 
        IDs = comps2rej{1:64,1};
        %find out what is the minimal/maximal sample count for each delay duration
        min1s = nan; 
        min3s = nan; 
        min9s = nan; 
        maxtime1s = nan; 
        maxtime3s = nan; 
        maxtime9s = nan; 
        
  for idx_file=1:length(IDs)
     ID = IDs{idx_file,1}; 
      % load cleaned postICA data 
     loadpath =['/mnt/homes/home028/gmonov/meg_analysis/sr_separated_data_uncut/', ID]
     load([loadpath, filesep, '1sDelay_' ID '.mat'])
     load([loadpath, filesep, '3sDelay_' ID '.mat'])
     load([loadpath, filesep, '9sDelay_' ID '.mat'])
    for  l = 1:length(trials_1s.trialinfo(:,1)) 
        if trials_1s.trialinfo(l,4)==1 % only 1s delay trials 
          min1s(end+1,1) = length(trials_1s.trial{1,l}(1,:)); 
          maxtime1s(end+1,1)=max(trials_1s.time{1,l}(1,:));
        end 
    end 
    
    for  l = 1:length(trials_3s.trialinfo(:,1)) 
          min3s(end+1,1) = length(trials_3s.trial{1,l}(1,:)); 
          maxtime3s(end+1,1)=max(trials_3s.time{1,l}(1,:));
    end 
    
    for  l = 1:length(trials_9s.trialinfo(:,1)) 
          min9s(end+1,1) = length(trials_9s.trial{1,l}(1,:)); 
          maxtime9s(end+1,1)=max(trials_9s.time{1,l}(1,:));
    end 
    
  end 
  
  %delete nan at the beginning 
   min1s(isnan(min1s)) = []; 
   min3s(isnan(min3s)) = []; 
   min9s(isnan(min9s)) = []; 
   
   maxtime1s(isnan(maxtime1s)) = []; 
   maxtime3s(isnan(maxtime3s)) = []; 
   maxtime9s(isnan(maxtime9s)) = []; 

   
   % Check if all lengths make sense first 
   
   check1s = min1s-(400*(0.6+1.8)); 
   check3s = min3s-(400*(3.8-1.2)); 
   check9s = min9s-(400*(9.8-1.2)); 
   
   checks = [check1s;check3s;check9s];
   % checks == 1 is the exact sample count 
   for i = 1:length(checks) 
       if checks(i)>11 || checks(i)<-9
           disp error
       end
   end 
  
   
   mintime1s = min(maxtime1s);
   mintime3s = min(maxtime3s);
   mintime9s = min(maxtime9s);
      
   % find out how many samples must be cut from actual end 
   cut1s = (1.8-mintime1s)*400;
   cut3s = (3.8-mintime3s)*400;
   cut9s = (9.8-mintime9s)*400;
    
   for idx_file=1:length(IDs)
     ID = IDs{idx_file,1}; 
      % load cleaned postICA data 
     loadpath =['/mnt/homes/home028/gmonov/meg_analysis/sr_separated_data_uncut/', ID]
     load([loadpath, filesep, '1sDelay_' ID '.mat'])
     load([loadpath, filesep, '3sDelay_' ID '.mat'])
     load([loadpath, filesep, '9sDelay_' ID '.mat'])

   % cut 1s delay duration at to the sample number that exists for every
   % trial/every participant 
 
     for l = 1:length(trials_1s.trialinfo(:,1)) %loop through every clean trials
         
         trials_1s.trial{1,l}(:,trials_1s.time{1,l}>mintime1s)=[]; 
         trials_1s.time{1,l}(trials_1s.time{1,l}>mintime1s)=[];
         trials_1s.sampleinfo(l,2) = trials_1s.sampleinfo(l,1)+(960-cut1s); % adapt sampleinfo to actually saved times and trial data 

     end 

   % cut 3s delay duration at to the sample number that exists for every
   % trial/every participant 
 
     for l = 1:length(trials_3s.trialinfo(:,1)) %loop through 3s clean trials
         trials_3s.trial{1,l}(:,trials_3s.time{1,l}>mintime3s)=[]; 
         trials_3s.time{1,l}(trials_3s.time{1,l}>mintime3s)=[]; 
     end 
         trials_3s.sampleinfo(:,2) = trials_3s.sampleinfo(:,1)+(1760-cut3s); %adjust end of trials to be equal 
         trials_3s.sampleinfo(:,1) = trials_3s.sampleinfo(:,1)+((0.6+1.2)*400); 
         
    % cut 9s delay duration at to the sample number that exists for every
    % trial/every participant 
 
     for l = 1:length(trials_9s.trialinfo(:,1)) %loop through 9s clean trials
         trials_9s.trial{1,l}(:,trials_9s.time{1,l}>mintime9s)=[]; 
         trials_9s.time{1,l}(trials_9s.time{1,l}>mintime9s)=[]; 
     end 
         trials_9s.sampleinfo(:,2) = trials_9s.sampleinfo(:,1)+(4160-cut9s); % adjust end of trials to be equal 
         trials_9s.sampleinfo(:,1) = trials_9s.sampleinfo(:,1)+((0.6+1.2)*400); 
         
         %Adjust data sets 1s
         tlen=[];
          for t = 1:length(trials_1s.time)
             tlen(t) = length(trials_1s.time{t}); % get lengths of all trials (in samples)
          end

          min_tr = find(tlen==min(tlen),1,'first');

          data_mne = struct; % build data structure for importing to mne
          data_mne.time = trials_1s.time{min_tr};
          data_mne.label = trials_1s.label;
          data_mne.dimord = 'rpt_chan_time';
          data_mne.trialinfo = trials_1s.trialinfo;
          data_mne.trialInfoLabel = trials_1s.trialInfoLabel;
          data_mne.sampleinfo = trials_1s.sampleinfo;
          data_mne.blockBounds = trials_1s.blockBounds;
          data_mne.grad = trials_1s.grad;
          data_mne.eye_trialtiming = trials_1s.eye_trialtiming;
          data_mne.eye_alltrials = trials_1s.eye_alltrials;
          data_mne.file_info = trials_1s.file_info;
          data_mne.elec = trials_1s.elec;

           trial_mne = [];
           
        for t = 1:length(trials_1s.trial)
           trial_mne(t,:,:) = trials_1s.trial{t}(:,1:length(data_mne.time));
        end
       clear alltrials

        data_mne.trial = trial_mne;

        trials_1s = data_mne;
      clear data_mne
      
               %Adjust data sets 3s
     if ~isempty(trials_3s.time) % Do not save when no 3s trials are there 
         tlen=[];
          for t = 1:length(trials_3s.time)
             tlen(t) = length(trials_3s.time{t}); % get lengths of all trials (in samples)
          end

          min_tr = find(tlen==min(tlen),1,'first');

          data_mne = struct; % build data structure for importing to mne
          data_mne.time = trials_3s.time{min_tr};
          data_mne.label = trials_3s.label;
          data_mne.dimord = 'rpt_chan_time';
          data_mne.trialinfo = trials_3s.trialinfo;
          data_mne.trialInfoLabel = trials_3s.trialInfoLabel;
          data_mne.sampleinfo = trials_3s.sampleinfo;
          data_mne.blockBounds = trials_3s.blockBounds;
          data_mne.grad = trials_3s.grad;
          data_mne.eye_trialtiming = trials_3s.eye_trialtiming;
          data_mne.eye_alltrials = trials_3s.eye_alltrials;
          data_mne.file_info = trials_3s.file_info;
          data_mne.elec = trials_3s.elec;

           trial_mne = [];
           
        for t = 1:length(trials_3s.trial)
           trial_mne(t,:,:) = trials_3s.trial{t}(:,1:length(data_mne.time));
        end
       clear alltrials

        data_mne.trial = trial_mne;

        trials_3s = data_mne;
      clear data_mne
     end 
      
               %Adjust data sets 9s
       if ~isempty(trials_9s.time) % Do not execute when no 9s trials are there 
         tlen=[];
          for t = 1:length(trials_9s.time)
             tlen(t) = length(trials_9s.time{t}); % get lengths of all trials (in samples)
          end

          min_tr = find(tlen==min(tlen),1,'first');

          data_mne = struct; % build data structure for importing to mne
          data_mne.time = trials_9s.time{min_tr};
          data_mne.label = trials_9s.label;
          data_mne.dimord = 'rpt_chan_time';
          data_mne.trialinfo = trials_9s.trialinfo;
          data_mne.trialInfoLabel = trials_9s.trialInfoLabel;
          data_mne.sampleinfo = trials_9s.sampleinfo;
          data_mne.blockBounds = trials_9s.blockBounds;
          data_mne.grad = trials_9s.grad;
          data_mne.eye_trialtiming = trials_9s.eye_trialtiming;
          data_mne.eye_alltrials = trials_9s.eye_alltrials;
          data_mne.file_info = trials_9s.file_info;
          data_mne.elec = trials_9s.elec;

           trial_mne = [];
           
        for t = 1:length(trials_9s.trial)
           trial_mne(t,:,:) = trials_9s.trial{t}(:,1:length(data_mne.time));
        end
       clear alltrials

        data_mne.trial = trial_mne;

        trials_9s = data_mne;
      clear data_mne
     end 
      % Ajust sampleinfo to timings in raw data instead of relative to
      % block onset
      if strcmp(ID,'43_2') % here no triggers were recoreded for the first block, change blocknums back to make it consistent with the other datasets 
          trials_1s.trialinfo(:,11) = trials_1s.trialinfo(:,11)-1; 
          trials_3s.trialinfo(:,11) = trials_3s.trialinfo(:,11)-1; 
          trials_9s.trialinfo(:,11) = trials_9s.trialinfo(:,11)-1; 
      end 
 
      for w=1:max(trials_1s.trialinfo(:,11)) %last block
          for l = 1:length(trials_1s.sampleinfo(:,1))
              if trials_1s.trialinfo(l,11) == w
               trials_1s.sampleinfo(l,:) = trials_1s.sampleinfo(l,:)+trials_1s.blockBounds(w,1); 
              end 
          end
      end 
       for w=1:max(trials_3s.trialinfo(:,11))
          for l = 1:length(trials_3s.sampleinfo(:,1))
              if trials_3s.trialinfo(l,11) == w
               trials_3s.sampleinfo(l,:) = trials_3s.sampleinfo(l,:)+trials_3s.blockBounds(w,1); 
              end 
          end
       end
      
       for w=1:max(trials_9s.trialinfo(:,11))
          for l = 1:length(trials_9s.sampleinfo(:,1))
              if trials_9s.trialinfo(l,11) == w
               trials_9s.sampleinfo(l,:) = trials_9s.sampleinfo(l,:)+trials_9s.blockBounds(w,1); 
              end 
          end
       end 
       
       % Add unique trial ID 
       
       for z = 1:length(trials_1s.trialinfo(:,1))
           uniqueID = []; 
           uniqueID = [num2str(trials_1s.trialinfo(z,1)),num2str(trials_1s.trialinfo(z,11)),num2str(trials_1s.trialinfo(z,13))];
           trials_1s.trialinfo(z,25) = str2num(uniqueID); 
       end 
       for z = 1:length(trials_3s.trialinfo(:,1))
           uniqueID = []; 
           uniqueID = [num2str(trials_3s.trialinfo(z,1)),num2str(trials_3s.trialinfo(z,11)),num2str(trials_3s.trialinfo(z,13))];
           trials_3s.trialinfo(z,25) = str2num(uniqueID); 
       end 
       for z = 1:length(trials_9s.trialinfo(:,1))
           uniqueID = []; 
           uniqueID = [num2str(trials_9s.trialinfo(z,1)),num2str(trials_9s.trialinfo(z,11)),num2str(trials_9s.trialinfo(z,13))];
           trials_9s.trialinfo(z,25) = str2num(uniqueID); 
       end 
       
            % Save new data sets
            % Define path/folder to save cleaned data
  
       mat_name = ['/mnt/homes/home028/gmonov/meg_analysis/sr_separated_data/' ID '/'];
    
        if 7==exist(mat_name,'dir')
            cd(mat_name)
        else
            mkdir(mat_name)
            cd(mat_name)
        end
         data = trials_1s;
         datastore = [ID '_1'];
         save(datastore,'data','-v7.3');
         clear data
         if ~isempty(trials_3s.time) % Do not save when no 3s trials are there 
         data = trials_3s;
         datastore = [ID '_3'];
         save(datastore,'data','-v7.3');
         clear data
         end
         if ~isempty(trials_9s.time) % Do not save when no 9s trials are there 
         data = trials_9s;
         datastore = [ID '_9'];
         save(datastore,'data','-v7.3');
         clear data 
         end 
    
   end 
         
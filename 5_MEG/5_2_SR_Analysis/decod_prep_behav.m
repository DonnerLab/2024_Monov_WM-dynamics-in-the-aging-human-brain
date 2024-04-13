% Save the sample stimulus location for decoding analyses 
% Gina Monov, UKE, 2022

clear all
close all
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
delays = [1,3,9]; 
for d = 1:length(delays) 
for f = 1:length(datasets_overview) 
    mem_loc = [];  
    tIDs = []; 
    ID = []; 
    trial_meta_del = []; 
    if datasets_overview(f).sess == 1
       datasets_overview(f).meg_trialinfo(:,end+1)=datasets_overview(f).meg_trial_ids_1s(:,1); %add trial ids to the trialinfo
       
       if delays(d) == 1
       tIDs = datasets_overview(f).meg_trial_ids_1s; 
       elseif delays(d) == 3
       tIDs = [datasets_overview(f).meg_trial_ids_3s;datasets_overview(f).meg_trial_ids_9s]; 
       elseif delays(d) == 9 
       tIDs = datasets_overview(f).meg_trial_ids_9s; 
       end 
       
       trial_meta_del = datasets_overview(f).meg_trialinfo(ismember(datasets_overview(f).meg_trialinfo(:,25),tIDs),:); % trialinfo with trials in current delay duration 
          % reorder trial_meta_del
           Indices = zeros(length(tIDs),1);
         for k = 1:numel(tIDs)
           Indices(k) = find(tIDs(k)==trial_meta_del(:,25));
         end
          trial_meta_del = trial_meta_del(Indices,:);
          
       mem_loc = trial_meta_del(:,2);   
           
       ID = cellstr(datasets_overview(f).ID); 
    save(['/mnt/homes/home028/gmonov/meg_analysis/Decoding/behav4decoding/', ID{1,1}, '_', num2str(delays(d)), '.mat'], 'mem_loc', 'tIDs') 
    end

end 



end 

   % Check certain things from datasets overview file and add a few
   % important infos 
   
    clear close all 

    
    addpath('/Users/ginamonov/Servers/mountpoint1/functions/');
    load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat'])
    loadpath = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'
    loadpath2 = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_analysis/clean_WM_data'
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/all_subj_mci.mat']);
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/all_patients.mat']);
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/cog_deficits.mat']);
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/hc1.mat']);
    
     
loadpath3 = '/Users/ginamonov/Servers/mountpoint1//meg_analysis/';
load([loadpath3,'datasets_overview1.mat']);

for l=1:length(datasets_overview)
    %add info about groups 
    if ismember (datasets_overview(l).subj,hc1)
        datasets_overview(l).group = 'hc'
    elseif ismember (datasets_overview(l).subj,all_pat)
        datasets_overview(l).group = 'mci'; 
    elseif ismember (datasets_overview(l).subj,all_cog_def)
        datasets_overview(l).group = 'cog_def'; 
    else datasets_overview(l).group = 'exc';
    end
    % add unique trial ids to check if they match the meg datasets
    log_ind = zeros(length(datasets_overview(l).meg_trialinfo(:,1)),1)
    for z = 1:length(datasets_overview(l).meg_trialinfo(:,1))
        if strcmp(datasets_overview(l).ID,'43_2')
            uniqueID = [num2str(datasets_overview(l).meg_trialinfo(z,1)),num2str(datasets_overview(l).meg_trialinfo(z,11)-1),num2str(datasets_overview(l).meg_trialinfo(z,13))];
        else
            uniqueID = [num2str(datasets_overview(l).meg_trialinfo(z,1)),num2str(datasets_overview(l).meg_trialinfo(z,11)),num2str(datasets_overview(l).meg_trialinfo(z,13))];
        end
        datasets_overview(l).meg_trial_ids_1s(z,1) = str2num(uniqueID);
     if datasets_overview(l).meg_trialinfo(z,4)==3 && datasets_overview(l).meg_trialinfo(z,24)==0 % delay period and no artifacts after the first second 
         log_ind(z) = 3;
     elseif datasets_overview(l).meg_trialinfo(z,4)==9 && datasets_overview(l).meg_trialinfo(z,24)==0
          log_ind(z) = 9;
     end
     clear uniqueID
    end
    % Add trial IDs fpr 3s and 9s only
    datasets_overview(l).meg_trial_ids_3s = datasets_overview(l).meg_trial_ids_1s(log_ind==3);
    datasets_overview(l).meg_trial_ids_9s = datasets_overview(l).meg_trial_ids_1s(log_ind==9);
end

save(['/Users/ginamonov/Servers/mountpoint1//meg_analysis/datasets_overview.mat'], 'datasets_overview') % Final master file which will be used for analyses 
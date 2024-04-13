% Pulls source reconstructed data for each subject and delay duration
% Saves the TFRs of power modulation
% Computes neural markers across-trial variability, within-trial
% variability, and power modulation during delay 
% Gina Monov, UKE, 2022 

function AVE_basic_TFR(nsubj)
addpath('/mnt/homes/home028/gmonov/functions/')
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
load('/mnt/homes/home028/gmonov/meg_analysis/all_glasser_parcels.mat')

new_ROIs={};
for g = 1:length(ROIs) 
    new_ROIs{g} = replace(ROIs{g},'_','-'); 
end 
ROIs = new_ROIs; 

all_glasser = {}; 
for z = 1:length(ROIs) 
    all_glasser(z) = ROIs(1,z); 
end 

% Set complete glasser to 'yes' to save the TFR based neural markers % Set
% to 'no' if you want to save the data averaged across the parecels within
% one parcel (includes JW M1 Hand area)
complete_glasser = 'yes'; % just put yes if data of all 180 parcels should be pulled 
for d = 1:3 % loop over delay durations 
    delays = [1,3,9];
    delay = delays(d); 
    
% % run a check on behavior and trialex
% for o = 1:length(datasets_overview) 
%     if length(datasets_overview(o).trialex(datasets_overview(o).trialex==1))+...
%           length(datasets_overview(o).trialex(datasets_overview(o).trialex==6))+...
%           length(datasets_overview(o).trialex(datasets_overview(o).trialex==3))  == length(datasets_overview(o).clean_behavior)
%         x = 1;
%     else disp(datasets_overview(o).ID);
%     end
% end 


%Pull all subjects with (usable) meg data 
allsubj={};
for f = 1:length(datasets_overview) 
    if datasets_overview(f).sess == 1
    allsubj = horzcat(allsubj,datasets_overview(f).subj);
    end
end 



%Pull number of sessions
nsess = [];
for f = 1:length(allsubj) 
    subj_sess = []; 
    for l = 1:length(datasets_overview) 
       if strcmp(allsubj{f},datasets_overview(l).subj)
          subj_sess = [subj_sess;datasets_overview(l).sess];
       end 
    end 
    nsess = [nsess,max(subj_sess)];
end 

subj = allsubj{nsubj};

proc = 'F';
type = 'Averaged';  % possibilities: Averaged, Lateralized, Pair

% Add as a last cluster all glasser parcels to save them for the whole
% cortex plots, but don't save in the end 
parcels={{'V1'},...
    {'V2', 'V3', 'V4'},...
    {'V6', 'V3A', 'V7', 'IPS1', 'V3B', 'V6A'},...
    {'MST', 'LO1', 'LO2', 'MT','V4t', 'FST', 'LO3', 'V3CD', 'PH'},...
    {'V8', 'FFC', 'PIT', 'VMV1', 'VMV3', 'VMV2', 'VVC'},...
    {'7Pm', '7AL', '7Am', '7PL', '7PC', 'LIPv', 'VIP', 'MIP', 'LIPd', 'AIP'},...
    {'6a', '6d'},...
    {'FEF', 'PEF'},...
    {'6v', '6r'},...
    {'EC', 'PreS', 'H', 'PeEc', 'PHA1', 'PHA3', 'TF', 'PHA2'},...
    {'PCV', '7m', 'POS1', 'v23ab', 'd23ab', '31pv',  '31pd', '31a', 'RSC', 'POS2', 'ProS', 'DVT', '23d', '23c'},...
    {'SFL', '8Av', '8Ad', '8BL', '9p', '8C', 'p9-46v', '46', 'a9-46v', '9-46d', '9a', 'i6-8', 's6-8'},...
    {'JWG_lr_M1'},...
    all_glasser}; 

names={'V1', 'Early_visual', 'Dorsal_visual', 'MT+_complex', 'Ventral_visual', 'Superior_parieta', 'PMd', 'Eye_fields', 'PMv',...
    'Medial_temporal', 'Posterior_cingulate', 'Dorsolateral_prefrontal', 'M1','all_glasser'};

for i=1:length(names)
    clusters(i).name=names{i};
    clusters(i).roi=parcels{i};
    clusters(i).parcel_data={};
end


if strcmp(complete_glasser,'yes')
megpath = ['/mnt/homes/home028/gmonov/meg_analysis/source_reconstruction/agg_complete_glasser/'];
clusters(1:end-1) = []; 
elseif strcmp(complete_glasser,'no') % refers to an old version of the analysis; only used to etract JWs M1 hand area for the paper plots 
megpath = ['/mnt/homes/home028/gmonov/meg_analysis/source_reconstruction/agg/'];
clusters(end) = []; 
end 
savepath = ['/mnt/homes/home028/gmonov/meg_analysis/task_analysis/TFR/output_ave_TFR/'];


% ==================================================================
% LOAD DATA
% ==================================================================
% loop through sessions and load data

% pull dataset index in the overview structure 
   for m=1:length(datasets_overview)
         if strcmp(datasets_overview(m).subj,subj) && datasets_overview(m).sess ==1
            overview_idx = m; 
        end
   end
   
for nclust = 1:length(clusters)
for parcels = 1:length(clusters(nclust).roi) %loop over the parcels in the cluster 
  
    %only considering 1st session
    %Loading data for delay duration of interest
    fname = [megpath,'S',subj,'_SESS1_' num2str(delay) '_',proc,'_agg.hdf'];
    fprintf('Loading %s...\n',fname)
    
    % load info about source rec dataset structure
    info = h5info(fname);  % pull information about file structure
    error_parcel = []; % initialize variable in case there is no data for a parcel 
    
    % pull times & frequencies vectors, trial IDs
    
        freqs = sort(cellfun(@str2num,[{info.Groups(1).Groups(1).Datasets(:).Name}]));
        try
         times = h5readatt(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(1)),'/'],'cols_time')';
        catch 
          clusters(nclust).parcel_data{1,parcels}=nan; %set to nan when there is no data for this parcel 
          error_parcel = 1; 
        end
        
    if isempty(error_parcel) % only when there is data...
       times = round(times,2);   % can be very minor error in timings which can throw off indexing - rounding to remove this
  
    
    % construct data matrix for this cluster (trials*freqs*times)
    parcel_data = [];
    for f = 1:length(freqs)
        data_in = h5read(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(f)),'/'])';
        
        mne_ids = double(h5readatt(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(f)),'/'],'rows_trial'))'; % pull trial IDs (transformed to deal with earlier bug)
       
        % get the correct order of trial ids
        if delay == 1 
       sorted_ids = datasets_overview(overview_idx).meg_trial_ids_1s';
       elseif delay == 3
           sorted_ids = datasets_overview(overview_idx).meg_trial_ids_3s';
       elseif delay == 9
           sorted_ids = datasets_overview(overview_idx).meg_trial_ids_9s';
        end 

       
        if mne_ids(:,1)~=sorted_ids(:,1)  % if trials aren't in the same order, sort
            disp order_is_incorrect
           Indices = zeros(length(sorted_ids),1);
         for k = 1:numel(sorted_ids)
           Indices(k) = find(sorted_ids(k)==mne_ids);
         end
            data_in = data_in(Indices,:);
        end
        parcel_data(:,f,:) = data_in;
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Load the 9s data for the 3s trials as well%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if delay == 3 && ~strcmp(subj,'01') && ~strcmp(subj,'02')% no clean data for 9s delay for this subject 
        
    fname = [megpath,'S',subj,'_SESS1_9_',proc,'_agg.hdf'];
    fprintf('Loading %s...\n',fname)
    
    % load info about source rec dataset structure
    info = h5info(fname);  % pull information about file structure
    %error_parcel = []; % initialize variable in case there is no data for a parcel 
    
    % pull times & frequencies vectors, trial IDs
    
        freqs = sort(cellfun(@str2num,[{info.Groups(1).Groups(1).Datasets(:).Name}]));
        
        times9 = h5readatt(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(1)),'/'],'cols_time')';
    
        times9 = round(times9,2);   % can be very minor error in timings which can throw off indexing - rounding to remove this
  
   
    % construct data matrix for this cluster (trials*freqs*times)
    
    parcel_data_9s = [];
    for f = 1:length(freqs)
        data_in = h5read(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(f)),'/'])';
        
        mne_ids = double(h5readatt(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(f)),'/'],'rows_trial'))'; % pull trial IDs (transformed to deal with earlier bug)
       
        % get the correct order of trial ids

           sorted_ids = datasets_overview(overview_idx).meg_trial_ids_9s';
       
        if mne_ids(:,1)~=sorted_ids(:,1)  % if trials aren't in the same order, sort
            disp order_is_incorrect
         Indices = zeros(length(sorted_ids),1);
         for k = 1:numel(sorted_ids)
           Indices(k) = find(sorted_ids(k)==mne_ids);
         end
            data_in = data_in(Indices,:);
        end
        parcel_data_9s(:,f,:) = data_in;
        
    end
    %trim to 3s delay duration 
    parcel_data_9s = parcel_data_9s(:,:,times9<=max(times));
    
    %concatenate 3s and 9s trials 
    parcel_data = cat(1,parcel_data,parcel_data_9s);
    
    end 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%load the first seconds from the 3s and 9s trials & concatenate%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if delay > 1 
        
        
    %Loading first second
    fname = [megpath,'S',subj,'_SESS1_1_',proc,'_agg.hdf'];
    fprintf('Loading %s...\n',fname)
    
    % load info about source rec dataset structure
    info = h5info(fname);  % pull information about file structure 
    
    % pull times & frequencies vectors, trial IDs
    
        freqs = sort(cellfun(@str2num,[{info.Groups(1).Groups(1).Datasets(:).Name}]));
     
        times1 = h5readatt(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(1)),'/'],'cols_time')';
        
        times1 = round(times1,2);   % can be very minor error in timings which can throw off indexing - rounding to remove this
  
    
    % construct data matrix for this cluster (trials*freqs*times)
    
    parcel_data_1s = [];
    for f = 1:length(freqs)
        data_in = h5read(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(f)),'/'])';
        
        mne_ids = double(h5readatt(fname,['/',type,'/',char(clusters(nclust).roi(parcels)),'/',num2str(freqs(f)),'/'],'rows_trial'))'; % pull trial IDs (transformed to deal with earlier bug)
       
        % get the correct order of trial ids

           sorted_ids = datasets_overview(overview_idx).meg_trial_ids_1s';
       
        if mne_ids(:,1)~=sorted_ids(:,1)  % if trials aren't in the same order, sort
            disp order_is_incorrect
           Indices = zeros(length(sorted_ids),1);
         for k = 1:numel(sorted_ids)
           Indices(k) = find(sorted_ids(k)==mne_ids);
         end
            data_in = data_in(Indices,:);
        end
        parcel_data_1s(:,f,:) = data_in;
        
    end
    
    %merge times 
    trl_timing = [times1, times];
    trl_timing = unique(trl_timing); 
    %Delete trials that don't match with the trials of interest
    % extract index of these trial ids 
    if delay == 3
        if strcmp(subj,'02')
       all_trials = datasets_overview(overview_idx).meg_trial_ids_3s; %one 9s trial - not source reconstructed!
        else
        all_trials = [datasets_overview(overview_idx).meg_trial_ids_3s;datasets_overview(overview_idx).meg_trial_ids_9s];
        end
        s_idx = ismember(datasets_overview(overview_idx).meg_trial_ids_1s,all_trials);
    elseif delay == 9 
       s_idx = ismember(datasets_overview(overview_idx).meg_trial_ids_1s,datasets_overview(overview_idx).meg_trial_ids_9s); 
       all_trials = datasets_overview(overview_idx).meg_trial_ids_9s; 
    end 
    
    %Remove trials
    parcel_data_1s = parcel_data_1s(s_idx,:,:);
    %Bring them in the correct order
    ex_trials =  datasets_overview(overview_idx).meg_trial_ids_1s(s_idx); 
    Indices = zeros(length(all_trials),1);
         for k = 1:numel(all_trials)
           Indices(k) = find(all_trials(k)==ex_trials);
         end
    parcel_data_1s = parcel_data_1s(Indices,:,:);
    %Remove overlapping data 
    parcel_data_1s = parcel_data_1s(:,:,times1<=1.4);
    parcel_data = parcel_data(:,:,times>1.4);
    
    %concatenate data
    parcel_data = cat(3,parcel_data_1s,parcel_data);
    
    %rename trl_timing for saving times 
    times = trl_timing; 
    else all_trials = datasets_overview(overview_idx).meg_trial_ids_1s; 
        
    end 
    
    % Extract data separately for correct and incorrect trials
    
        accuracy = zeros(length(all_trials(:,1)),1);
        sorted_ids = all_trials(:,1)'; 
        
        
         datasets_overview(overview_idx).meg_trialinfo(:,end+1)=datasets_overview(overview_idx).meg_trial_ids_1s(:,1); %add trial ids to the trialinfo
        % extract trialinfo only for selected trials in the correct delay
        % period 
        trial_meta_del = datasets_overview(overview_idx).meg_trialinfo(ismember(datasets_overview(overview_idx).meg_trialinfo(:,25),sorted_ids),:); % trialinfo with trials in current delay duration 
          % reorder trial_meta_del to match the meg data
           Indices = zeros(length(sorted_ids),1);
         for k = 1:numel(sorted_ids)
           Indices(k) = find(sorted_ids(k)==trial_meta_del(:,25));
         end
          trial_meta_del = trial_meta_del(Indices,:);
          
       accuracy = trial_meta_del(:,7); % Extract if response was correct
       data_correct = parcel_data; %initialize parcel data only correct trials 
       data_error = parcel_data; %initialize parcel data only error trials 
       data_correct(accuracy==0,:,:) = []; % remove incorrect trials 
       data_error(accuracy==1,:,:) = []; % remove correct trials 
       if delay == 1
       % Extract data for within and across trial-variability    
        across_trial_var = parcel_data(:,freqs>=5 & freqs <=35, times>=0.5 & times <1.5); 
        across_trial_var = mean(across_trial_var,2); 
        across_trial_var = squeeze(mean(across_trial_var,3));
        clusters(nclust).across_trial_var{1,parcels} = var(across_trial_var); 
        across_trial = var(across_trial_var); 
        
        within_trial_var = parcel_data(:,freqs>=5 & freqs <=35, times>=0.5 & times <1.5); 
        within_trial_var = squeeze(mean(within_trial_var,2)); 
        within_trial_var = var(within_trial_var,0,2); 
        clusters(nclust).within_trial_var{1,parcels} = mean(within_trial_var);
        within_trial = mean(within_trial_var); 
        
        % Extract power in delay period 
        
        del_power = parcel_data(:,freqs>=5 & freqs <=35, times>=0.5 & times <1.5); 
        del_power = mean(del_power,2); 
        del_power = squeeze(mean(del_power,3)); 
        clusters(nclust).del_power{1,parcels} = mean(del_power); 
        mean_TFR = mean(del_power); 
        % save this data separately for all the glasser parcels in order to pull it easier
        if strcmp(clusters(nclust).name,'all_glasser')
        save(['/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/mean_TFR/',subj, '_', clusters(nclust).roi{1,parcels}, '.mat'] , 'mean_TFR');
        save(['/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/within_trial/',subj, '_', clusters(nclust).roi{1,parcels}, '.mat'], 'within_trial');
        save(['/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/across_trial/',subj, '_', clusters(nclust).roi{1,parcels}, '.mat'], 'across_trial'); 
        end 
       end 
    clusters(nclust).parcel_data{1,parcels}=mean(parcel_data,1); %save mean of the data within this parcel
    clusters(nclust).data_correct{1,parcels}=mean(data_correct,1);
    clusters(nclust).data_error{1,parcels}=mean(data_error,1);
    clusters(nclust).trial_count=size(parcel_data,1);%save how many trials contributed
    clear parcel_data error_parcel
    end
end
end

if strcmp(complete_glasser,'yes')
   save(['/mnt/homes/home028/gmonov/meg_analysis/task_analysis/TFR/output_ave_TFR_all_glasser/',subj,'_' num2str(delay),'_',proc,'_AVE_TFR_all_glasser.mat'],'times','freqs','clusters')
end 

% Save mean within entire cluster
%First delete the last row that only contains all the remaining Glasser
%parcels 

if strcmp(complete_glasser,'no')
 
for u = 1:length(clusters) 
    all_parcels = []; 
    for g=1:length(clusters(u).parcel_data)
        if ~isnan(clusters(u).parcel_data{1,g})
          all_parcels = [all_parcels;clusters(u).parcel_data{1,g}]; 
        end
    end 
    clusters(u).averaged = mean(all_parcels,1);
end 
   

for u = 1:length(clusters) 
    all_parcels = []; 
    for g=1:length(clusters(u).data_correct)
        if ~isnan(clusters(u).data_correct{1,g})
          all_parcels = [all_parcels;clusters(u).data_correct{1,g}]; 
        end
    end 
    clusters(u).correct_averaged = mean(all_parcels,1);
end 

for u = 1:length(clusters) 
    all_parcels = []; 
    for g=1:length(clusters(u).data_error)
        if ~isnan(clusters(u).data_error{1,g})
          all_parcels = [all_parcels;clusters(u).data_error{1,g}]; 
        end
    end 
    clusters(u).error_averaged = mean(all_parcels,1);
end 

% Save the across and within trial vaiability 
if delay == 1
   for u = 1:length(clusters) 
    all_parcels = []; 
    for g=1:length(clusters(u).across_trial_var)
        if ~isnan(clusters(u).across_trial_var{1,g})
          all_parcels = [all_parcels;clusters(u).across_trial_var{1,g}]; 
        end
    end 
    clusters(u).all_across_tr_var = all_parcels; % save data for every single parcel 
    clusters(u).across_trial_var = mean(all_parcels,1);
   end  
   
   for u = 1:length(clusters) 
    all_parcels = []; 
    for g=1:length(clusters(u).within_trial_var)
        if ~isnan(clusters(u).within_trial_var{1,g})
          all_parcels = [all_parcels;clusters(u).within_trial_var{1,g}]; 
        end
    end 
    clusters(u).all_within_tr_var = all_parcels; % save data for every single parcel 
    clusters(u).within_trial_var = mean(all_parcels,1);
   end  
   
   for u = 1:length(clusters) 
    all_parcels = []; 
    for g=1:length(clusters(u).del_power)
        if ~isnan(clusters(u).del_power{1,g})
          all_parcels = [all_parcels;clusters(u).del_power{1,g}]; 
        end
    end 
    clusters(u).all_del_power = all_parcels; % save data for every single parcel 
    clusters(u).del_power = mean(all_parcels,1);
   end     
    
    
end 



savename = [savepath,subj,'_' num2str(delay),'_',proc,'_AVE_TFR.mat'];


save(savename,'times','freqs','clusters')
clearvars -except d nsubj datasets_overview
end

end 
end


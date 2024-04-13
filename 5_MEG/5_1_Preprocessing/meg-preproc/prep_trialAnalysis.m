function [] = prep_trialAnalysis(idx_file)
%PREP_TRIALANALYSIS prepares data for a trial based analysis
%   Data is split into trials, artifacts get removed
%   Output is cleaned data split into trials of one file 


    % Go to file_info folder and pull information from there 
    cd /mnt/homes/home028/gmonov/meg_analysis/Dataset_info/final_masterfile
   
    % Load important information about files
    % files contains the complete names of all files that must be processed
    % info_EL_blocks: matrix with 5 columns
    % col 1 - first block to process with regard to behavioral data (normally 1)
    % col 2 - last block to process with regard to behavioral data (normally 2-4)
    % col 3 - number of blocks in file (normally 3)
    % col 4 - session (1 or 2)
    % col 5 - eyelink data first block (1: usable, 0: use veog instead)
    % col 6 - eyelink data second block
    % col 7 - eyelink data third block (if nan 
    loadpath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/final_masterfile/';
    load([loadpath,'Master_file.mat']);
    info_EL_blocks = zeros(length(file_info),7);
    files = {zeros(length(file_info),1)};
    for f = 1:length(file_info)
        info_EL_blocks(f,1)=file_info(f).blocknums(1);
        info_EL_blocks(f,2)=file_info(f).blocknums(end);
        info_EL_blocks(f,3)=length(file_info(f).blocknums);
        info_EL_blocks(f,4)=file_info(f).sess;
        info_EL_blocks(f,5)=file_info(f).ELex(1);
        info_EL_blocks(f,6)=file_info(f).ELex(2);
        if length(file_info(f).blocknums)==4
            info_EL_blocks(f,7:8)=file_info(f).ELex(3:4);
        elseif length(file_info(f).blocknums)==3
            info_EL_blocks(f,7)=file_info(f).ELex(3);
            info_EL_blocks(f,8)=nan;
        else info_EL_blocks(f,7:8)=nan;
        end 
       
        files{f,1}=file_info(f).name; 
    end 
        
    filein = files{idx_file};
    ID = [file_info(idx_file).subj,'_',num2str(file_info(idx_file).sess)];
    current_file_info = file_info(idx_file);
    
    %Check whether dataset contains a lot of metal artifacts 
    cd '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/metal_artifacts/'
    exc_metal = isfile(strcat(ID,'_metal.mat')); 
    if strcmp(ID,'46_1')
        exc_metal = 1; 
    end 
    addpath '/mnt/homes/home032/aarazi/fieldtrip-20201009'
    addpath '/mnt/homes/home028/gmonov/meg_data/'
    addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/prep_trialAnalysis/'
    addpath '/mnt/homes/home028/gmonov/functions/'
    ft_defaults


    fprintf('\n\n ---------------- \n PROCESSING FILE %s...\n ---------------- \n\n', ID)
    
    %% Construct trial matrix
    
    % =====================================================================
    %   1   Construct trial matrix
    % =====================================================================
    % go to data folder
    cd /mnt/homes/home028/gmonov/meg_data/

    fprintf('\n ---------- 1 Constructing trial matrix ----------\n')
    
    % Specify cfg
    cfgin = [];
    cfgin.idx = idx_file; 
    cfgin.ID = ID;
    cfgin.subj = file_info(idx_file).subj;
    cfgin.session = file_info(idx_file).sess;
    cfgin.fileNum = file_info(idx_file).rec;
    cfgin.dataset = filein;
    cfgin.trialdef.prestim = 0.6; % 0.6 s second offset for TF baseline
    cfgin.trialdef.poststim = 0.3; % add 0.3 s after Probe-onset
    cfgin.trialfun = 'trialbasedfun_MCI';
    cfgin = ft_definetrial(cfgin);
    
    % Needed cfgin stored in cfgin.trl
    cfgtrial = cfgin.trl;
    clear cfgin
    cfgtrial.alltrl = cfgtrial.trl;
    cfgtrial = rmfield(cfgtrial,'trl');
    
    % Downsample originally created trial matrix
    cfgtrial.alltrl(:,1:3) = round(cfgtrial.alltrl(:,1:3)/1200*400);
    cfgtrial.alltrl(:,17:26) = round(cfgtrial.alltrl(:,17:26)/1200*400);
    % Downsample the block bounds
    cfgtrial.blockBound_trl = round(cfgtrial.blockBound_trl/1200*400);


    %% Loop through blocks

    firstBl = 1;
    
   
    loop_vec = firstBl:info_EL_blocks(idx_file,3);

    for block = loop_vec
    
        % Consider the actual block number for behavioral data (1-8)
        actualblock = (block-1) + info_EL_blocks(idx_file,1);


        fprintf('\n\n ---------------- \n Loop through BLOCK #%d...\n ---------------- \n\n', block)

        % =====================================================================
        %   2   DEFINE BLOCK AND LOAD RELEVANT CHANNELS 
        % =====================================================================

        fprintf('\n ---------- 2 Defining block and relevant channels ----------\n')

        % Specify cfg
        cfgbl = [];
        cfgbl.dataset = filein;
        cfgbl.block = block; 
        cfgbl.ID = ID;
        cfgbl.trialdef.prestim = 10; % Add 10 seconds before and after the block
        cfgbl.trialdef.poststim = 10;
        cfgbl.trialfun = 'trialfun_MCI_continuous';
        cfgbl = ft_definetrial(cfgbl);
        cfgbl.continuous = 'yes'; % read in data as continuous
        cfgbl.channel = {'meg','EEG001','EEG002','EEG003','EEG057','EEG058','EEG059','HLC*','UADC*'};
        data = ft_preprocessing(cfgbl);
        fsample_old = data.fsample;
        sampleinfo_old = data.sampleinfo; 
        
        % =====================================================================
        %   3    REMOVE LINE NOISE
        % =====================================================================

        fprintf('\n ---------- 3 Filtering out line noise (50 Hz) and its harmonics -----------\n')

        cfg             = [];
        cfg.bsfilter    = 'yes';
        cfg.bsfreq      = [49 51; 99 101; 149 151];
        data            = ft_preprocessing(cfg, data);
        
        % =====================================================================
        %   4   HIGH PASS FILTER (cutoff 0.1 Hz)
        % =====================================================================

        fprintf('\n ---------- 4 High pass filtering (cutoff 0.1 Hz) -----------\n')

        cfg          = [];
        cfg.hpfilter = 'yes';
        cfg.hpfreq   = 0.1;  
        cfg.hpfiltord = 3;
        cfg.hpfilttyoe = 'firws';
        data = ft_preprocessing(cfg,data);
        
        % =====================================================================
        %   5    DOWN SAMPLE DATA TO 400 HZ
        % =====================================================================

        fprintf('\n ---------- 5 Resampling block -----------\n')

        cfgres.resample = 'yes';
%         cfgres.fsample = 1200;
        cfgres.resamplefs = 400;
        cfg.detrend = 'no';
        data = ft_resampledata(cfgres, data);
        sampleinfo_new = [round(sampleinfo_old(1))./fsample_old.*cfgres.resamplefs...
                          round(sampleinfo_old(1))./fsample_old.*cfgres.resamplefs+length(data.time{1}-1)]; 

        % =====================================================================
        %   6    LOAD ICA weights and IDs of artifactual components
        % =====================================================================

        fprintf('\n ---------- 6 Load ICA weights and IDs of artifactual components -----------\n')
        name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/comp_ICA/' ID];
        load([name '/comp_' ID '.mat' ]);
        addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/comps_rejection/';
        
        comps2rej = readtable ('comps_rejection_wm','Range','A1:L64'); 
        for rejloop=1:height(comps2rej)
            if strcmp(table2array(comps2rej(rejloop,1)),{ID})
               rejComps = table2array(unique(comps2rej(rejloop,2:end)));
               rejComps(isnan(rejComps))=[]; 
            end 
        end 
       
        
        % Remove artifactual components
        cfg = [];
        cfg.component = rejComps;
        data_cl = ft_rejectcomponent(cfg, comp, data);
        
        % Separate current block ! ACTUALBLOCK CORRECT?
        cfgtrial.trl = cfgtrial.alltrl(cfgtrial.alltrl(:,14)==actualblock,:);
        
        % Substract block onset from all sample numbers (except for offset) in the trial matrix
        %if block start hasnt been recoreded the block bound might be
        %negative 
        if cfgtrial.blockBound_trl(block,1)<0
            cfgtrial.trl(:,1:2) = cfgtrial.trl(:,1:2) + (cfgtrial.blockBound_trl(block,1));
            cfgtrial.trl(:,17:26) = cfgtrial.trl(:,17:26) + (cfgtrial.blockBound_trl(block,1));
        else
        cfgtrial.trl(:,1:2) = cfgtrial.trl(:,1:2) - (cfgtrial.blockBound_trl(block,1));
        cfgtrial.trl(:,17:26) = cfgtrial.trl(:,17:26) - (cfgtrial.blockBound_trl(block,1));
        end
        % Segment block's data into trials
        trials = ft_redefinetrial(cfgtrial, data_cl);
        
        % Load previously created artifact matrices (aligned to block onset!)
        name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/artifact_identification/preICA_artifactMatrices/' ID];
 
        artifacts = ['Artifacts_' ID '_Block_' num2str(block)];

        load([name '/' artifacts]);
        
        % Mark trials that overlap with an artifact
        % Keep trials with delay duration >1s, when artifacts happened
        % later in time (a trial with 1 s delay duration lasts for 2.7s
        % (1080 samples at sampling rate == 400) 
        % from memorandum onset to probe onset including 0.6s pre/ 0.3s
        % post memorandum  onset 
        % stimulus 
        trls2remove = [];
        eye_art = []; % marking trials with saccades/blinks instead of removing them 
        eye_art_blink=[]; 
        eye_art_sacc=[]; 
        %Initialize variables for counting the amount of trials that contain a
        %certain artifact 
        HM_count = []; 
        jump_count = []; 
        muscle_count = []; 
        metal_count = []; 
        cog_count = [];
        all_trials_count = length(cfgtrial.alltrl(:,1)); 
        artifactual_trials = zeros(length(trials.sampleinfo),1); %only trials marked with 0 are clean from any artifacts
        % how many trials per block? 
            trbl1=cfgtrial.alltrl((cfgtrial.alltrl(:,14)==1)); 
            trbl2=cfgtrial.alltrl((cfgtrial.alltrl(:,14)==2)); 
            trbl3=cfgtrial.alltrl((cfgtrial.alltrl(:,14)==3));

        for j = 1:length(trials.sampleinfo)
            
            % Mark trials that have negative sample times for removal
            if trials.sampleinfo(j,1)<0
               trls2remove = [trls2remove j];
            end 
            % Mark trials with head movements
            for i = 1:size(artifact_headM,1)
                if ~isempty(intersect([artifact_headM(i,1):artifact_headM(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))%overlap between trial and artifact time window?
                    if min(intersect([artifact_headM(i,1):artifact_headM(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))-trials.sampleinfo(j,1)>960 %keep trials when everything up to 1s delay duration is clean 
                       artifactual_trials(j) = artifactual_trials(j)+1; %add 1 if for every artifactual sequence that is detected in this trial  
                       HM_count = [HM_count j]; 
                    else trls2remove = [trls2remove j];
                         artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                         HM_count = [HM_count j]; 
                    end 
                end
            end
            % Mark trials with jumps
            for i = 1:size(artifact_Jump,1)
                if ~isempty(intersect([artifact_Jump(i,1):artifact_Jump(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))%overlap between trial and artifact time window?
                    if min(intersect([artifact_Jump(i,1):artifact_Jump(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))-trials.sampleinfo(j,1)>960 %keep trials when everything up to 1s delay duration is clean 
                       artifactual_trials(j) = artifactual_trials(j)+1; %add 1 if for every artifactual sequence that is detected in this trial  
                       jump_count = [jump_count j];
                    else trls2remove = [trls2remove j];
                         artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                         jump_count = [jump_count j]; 
                    end 
                end
            end
            % Mark trials with muscle artifacts
            for i = 1:size(artifact_Muscle,1)
                if ~isempty(intersect([artifact_Muscle(i,1):artifact_Muscle(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))%overlap between trial and artifact time window?
                    if min(intersect([artifact_Muscle(i,1):artifact_Muscle(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))-trials.sampleinfo(j,1)>960 %keep trials when everything up to 1s delay duration is clean 
                       artifactual_trials(j) = artifactual_trials(j)+1; %add 1 if for every artifactual sequence that is detected in this trial   
                       muscle_count = [muscle_count j]; 
                    else trls2remove = [trls2remove j];
                         artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                         muscle_count = [muscle_count j]; 
                    end 
                end
            end

            %Mark trials with metal artefacts (only if dataset isn't full
            %of them) 
           if exc_metal == 0 
            for i = 1:size(artifact_metal,1)
                if ~isempty(intersect([artifact_metal(i,1):artifact_metal(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))%overlap between trial and artifact time window?
                    if min(intersect([artifact_metal(i,1):artifact_metal(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))-trials.sampleinfo(j,1)>960 %keep trials when everything up to 1s delay duration is clean 
                       artifactual_trials(j) = artifactual_trials(j)+1; %add 1 if for every artifactual sequence that is detected in this trial  
                       metal_count = [metal_count j]; 
                    else trls2remove = [trls2remove j];
                         artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                         metal_count = [metal_count j]; 
                    end 
                end
            end
           end 
            
%             Mark trials with cognitive artefacts according to trialex
%             matrix in file_info 

            if actualblock==1
               if cfgtrial.alltrl(j,15)==0
                  trls2remove = [trls2remove j];
                  artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                  cog_count = [cog_count j]; 
               end
            elseif  actualblock==2
                if cfgtrial.alltrl(j+length(trbl1),15)==0
                  trls2remove = [trls2remove j];
                  artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                  cog_count = [cog_count j]; 
               end
            elseif actualblock==3
                if cfgtrial.alltrl(j+length(trbl1)+length(trbl2),15)==0
                  trls2remove = [trls2remove j];
                  artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                  cog_count = [cog_count j]; 
               end
            elseif actualblock ==4
                if cfgtrial.alltrl(j+length(trbl1)+length(trbl2)+length(trbl3),15)==0
                  trls2remove = [trls2remove j];
                  artifactual_trials(j) = nan; % set to nan if trial will later be removed (artifacts already for the time of shortest delay duration)
                  cog_count = [cog_count j]; 
               end
            end 
           
        end
        
        HM_count_all(block) = length(unique(HM_count)); 
        jump_count_all(block) = length(unique(jump_count)); 
        muscle_count_all(block) = length(unique(muscle_count)); 
        metal_count_all(block) = length(unique(metal_count)); 
        cog_count_all(block) = length(unique(cog_count)); 
        %add info about trials that contain artifacts only after 1 s delay
        %duration 
        trials.trialinfo(:,end+1)=artifactual_trials; 
        % Some trials needed to be removed because of multiple artifacts
        trls2remove = unique(trls2remove); % Remove douplicates
        numtrial = 1:length(trials.sampleinfo);
        trls2keep = numtrial(~ismember(numtrial,trls2remove));
        
        cfg = [];
        cfg.trials = trls2keep;
          % Remove these trials
        trials = ft_selectdata(cfg, trials);
        eye_artifacts_vec = zeros(size(cfg.trials)); 
        eye_artifacts = {}; % initialize the final cell array where the first column gives information whether an eye related artifact has occured in this trial
        % 0 = no eye-related artifact 
        % 1 = blink and saccade artifacts
        % 2 = only saccade(s)
        % 3 = only blink(s)
        % 2nd column: saccade times 3rd column: blink times 
                    %Mark trials with saccades
                    for j=1:length(trials.sampleinfo)
                        for i = 1:size(artifact_saccade,1)
                            if ~isempty(intersect([artifact_saccade(i,1):artifact_saccade(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))
                                eye_art = [eye_art j];
                                eye_art_sacc = [eye_art_sacc,[j; artifact_saccade(i,1);artifact_saccade(i,2)]];
                            end
                        end
                        % Mark trials with blinks
                        for i = 1:size(artifact_Ygaze,1)
                            if ~isempty(intersect([artifact_Ygaze(i,1):artifact_Ygaze(i,2)],[trials.sampleinfo(j,1):trials.sampleinfo(j,2)]))
                                eye_art = [eye_art j];
                                eye_art_blink = [eye_art_blink, [j; artifact_Ygaze(i,1);artifact_Ygaze(i,2)]];
                            end
                        end
                    end
                    
        for l=1:length(eye_artifacts_vec)
            if ~isempty (eye_art_blink) && ismember(l,eye_art_blink(1,:)) && ~isempty (eye_art_sacc) && ismember(l,eye_art_sacc(1,:)) 
                eye_artifacts_vec(l) = 1;
            elseif ~isempty (eye_art_sacc) && ismember(l,eye_art_sacc(1,:))
                eye_artifacts_vec(l) = 2;
            elseif ~isempty (eye_art_blink) &&ismember(l,eye_art_blink(1,:))
                eye_artifacts_vec(l) = 3;
            end 
            eye_artifacts{l,1} = {eye_artifacts_vec(l)};
            eye_artifacts{l,2} = {}; %specific saccade times
            eye_artifacts{l,3} = {}; %specific blink times 
        end 
          if isempty(eye_art_sacc) disp 'No Saccades'
          else 
             for l=1:length([eye_artifacts{:,1}])
               for k=1:length(eye_art_sacc(1,:))
                  if eye_art_sacc(1,k) == l 
                 eye_artifacts{l,2} = [eye_artifacts{l,2}, [eye_art_sacc(2,k); eye_art_sacc(3,k)]];
                  end 
               end 
             end
          end  
      if isempty(eye_art_blink) disp 'No Blinks'
      else
        for l=1:length([eye_artifacts{:,1}])
        for k=1:length(eye_art_blink(1,:))
              if eye_art_blink(1,k) == l 
                 eye_artifacts{l,3} = [eye_artifacts{l,3}, [eye_art_blink(2,k); eye_art_blink(3,k)]];
              end 
        end 
        end 
      end 
        
         if block ==1
            eye_artifacts_1 = eye_artifacts;
         elseif block == 2
          eye_artifacts_2 = eye_artifacts;     
         elseif block == 3
             eye_artifacts_3 = eye_artifacts;
         elseif block ==4
             eye_artifacts_4 = eye_artifacts;
         end 
         
        
        % Concatenate blocks
        if block == 1 
            old_trials = trials;
        else
            cfg = [];
            old_trials = ft_appenddata(cfg, old_trials, trials);
        end
    end

    % Define path/folder to save cleaned data
  
    mat_name = ['/mnt/homes/home028/gmonov/meg_analysis/data_trials/' ID '/'];
    
        if 7==exist(mat_name,'dir')
            cd(mat_name)
        else
            mkdir(mat_name)
            cd(mat_name)
        end
    
    % Rename
    all_trials_cl = old_trials;
    if block == 1
       all_trials_cl.eye_artifacts = eye_artifacts_1;
    elseif block == 2
      all_trials_cl.eye_artifacts = {eye_artifacts_1,eye_artifacts_2};
    elseif block == 3
      all_trials_cl.eye_artifacts = {eye_artifacts_1,eye_artifacts_2,eye_artifacts_3};
    elseif block == 4
      all_trials_cl.eye_artifacts = {eye_artifacts_1,eye_artifacts_2,eye_artifacts_3, eye_artifacts_4};
    end 
    all_trials_cl.trialInfoLabel = cfgtrial.trialInfoLabel(:,4:end);
    all_trials_cl.eye_alltrials = vertcat(all_trials_cl.eye_artifacts{1,1:end});
% including blink/saccade times relative to trial onset 
all_trials_cl.eye_trialtiming = all_trials_cl.eye_alltrials;
for o = 1:length(all_trials_cl.eye_trialtiming(:,1))
    if all_trials_cl.eye_trialtiming{o,1}{1,1}  == 1 % both saccades and blinks in this trial 
        % realigning times of saccades to trial onset
        for sacc=1:length(all_trials_cl.eye_trialtiming{o,2})
            for start_end = 1:2
                all_trials_cl.eye_trialtiming{o,2}{1,sacc}(start_end,:) = all_trials_cl.eye_alltrials{o,2}{1,sacc}(start_end,:)-all_trials_cl.sampleinfo(o,1);
            end
        end
        for blink=1:length(all_trials_cl.eye_trialtiming{o,3})
            for start_end = 1:2
                all_trials_cl.eye_trialtiming{o,3}{1,blink}(start_end,:) = all_trials_cl.eye_alltrials{o,3}{1,blink}(start_end,:)-all_trials_cl.sampleinfo(o,1);
            end
        end
        
    elseif all_trials_cl.eye_trialtiming{o,1}{1,1}  == 3
        for blink=1:length(all_trials_cl.eye_trialtiming{o,3})
            for start_end = 1:2
                all_trials_cl.eye_trialtiming{o,3}{1,blink}(start_end,:) = all_trials_cl.eye_alltrials{o,3}{1,blink}(start_end,:)-all_trials_cl.sampleinfo(o,1);
            end
        end
    elseif all_trials_cl.eye_trialtiming{o,1}{1,1}  == 2
        for sacc=1:length(all_trials_cl.eye_trialtiming{o,2})
            for start_end = 1:2
                all_trials_cl.eye_trialtiming{o,2}{1,sacc}(start_end,:) = all_trials_cl.eye_alltrials{o,2}{1,sacc}(start_end,:)-all_trials_cl.sampleinfo(o,1);
            end
        end
    end
end 
    


    all_trials_cl.blockBounds = cfgtrial.blockBound_trl;
    
    %add info how many trials contained a sprecific artifact
    
    all_trials_cl.HM_count = sum(HM_count_all); 
    all_trials_cl.jump_count = sum(jump_count_all); 
    all_trials_cl.muscle_count = sum(muscle_count_all); 
    all_trials_cl.metal_count = sum(metal_count_all); 
    all_trials_cl.cog_count = sum(cog_count_all); 
    all_trials_cl.alltrials_count = all_trials_count; 
    % add label for additional column on whether there are artifacts after
    % 960 samples (== trials with a 1s delay duration)
    all_trials_cl.trialInfoLabel(1,end+1)={'artifacts'}; 
    % add file_info row of current data set for sanity checks 
    all_trials_cl.file_info = current_file_info; 
    %correct that actual 1st block didn't contain triggers 
    if strcmp(ID,'43_2')
       all_trials_cl.trialinfo(:,11)=all_trials_cl.trialinfo(:,11)+1;
       all_trials_cl.file_info.blocknums=[2,3];
    end 
    
    clear old_trials;
    
    % Save concatenated data file
    datastore = ['data_clean_postICA_' ID];
    save(datastore,'all_trials_cl','-v7.3');

end


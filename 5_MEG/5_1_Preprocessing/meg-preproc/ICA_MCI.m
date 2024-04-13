function [] = ICA_MCI(idx_file)
% PERFORM ICA and save ICA weightings as well as data/artifact matrices 
    % - Realtive to block onset!
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
    % col 7 - eyelink data third block (if nan's the block doesn't exist
    % for this participant) 
    % col 8 - eyelink data fourth block 
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
    
    sample_counter = 0;
    trial_data = [];
    
    addpath '/mnt/homes/home032/aarazi/fieldtrip-20201009'
    addpath '/mnt/homes/home028/gmonov/meg_data/'
    addpath '/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/'
    addpath '/mnt/homes/home028/gmonov/functions/'
    ft_defaults
        % Go to data folder
   cd /mnt/homes/home028/gmonov/meg_data/


    %% Loop thorugh blocks
    %   Loop trough each block separately
    %   Load data of relevant channels, downsampling to 400Hz, filter out 50Hz
    %   and harmonics (line noise), load artifact matrices and remove head
    %   movements, muscle and jump artifacts



    firstBl = 1;

    for block = firstBl:info_EL_blocks(idx_file,3)

        cd /mnt/homes/home028/gmonov/meg_data/  

        fprintf('\n\n ---------------- \n PREPARING AND PERFORMING ICA FILE %s, BLOCK #%d...\n ---------------- \n\n', ID, block)

        % =====================================================================
        %   1   DEFINE BLOCK AND LOAD RELEVANT CHANNELS 
        % =====================================================================

        fprintf('\n ---------- 1 Defining block and relevant channels ----------\n')

        % Specify cfg
        cfgin = [];
        cfgin.dataset = filein;
        cfgin.ID = ID;
        cfgin.block = block;
        cfgin.trialdef.prestim = 10; % Add 10 seconds before and after the block
        cfgin.trialdef.poststim = 10;
        cfgin.trialfun = 'trialfun_MCI_continuous';
        cfgin = ft_definetrial(cfgin);
        cfgin.continuous = 'yes'; % read in data as continuous
        cfgin.channel = {'meg','EEG001','EEG002','EEG003','EEG057','EEG058','EEG059','HLC*','UADC*'}; %added EEG059 --> in MCI study ECG channel
        data = ft_preprocessing(cfgin);
        fsample_old = data.fsample;
        sampleinfo_old = data.sampleinfo; 

        % =====================================================================
        %   2    DOWN SAMPLE DATA TO 400 HZ
        % =====================================================================

        fprintf('\n ---------- 2 Resampling block -----------\n')

        cfgres.resample = 'yes';
%         cfgres.fsample = 1200;
        cfgres.resamplefs = 400;
        cfg.detrend = 'no';
        data = ft_resampledata(cfgres, data);
        sampleinfo_new = [round(sampleinfo_old(1))./fsample_old.*cfgres.resamplefs...
                          round(sampleinfo_old(1))./fsample_old.*cfgres.resamplefs+length(data.time{1}-1)]; 

        % =====================================================================
        %   3    REMOVE LINE NOISE
        % =====================================================================

        fprintf('\n ---------- 3 Filtering out line noise (50 Hz) and its harmonics -----------\n')

        cfg             = [];
        cfg.bsfilter    = 'yes';
        cfg.bsfreq      = [49 51; 99 101; 149 151];
        data            = ft_preprocessing(cfg, data);

        % =====================================================================
        %   4   HIGH PASS FILTER (cutoff 1 Hz)
        % =====================================================================

        fprintf('\n ---------- 4 High pass filtering (cutoff 1 Hz) -----------\n')

        cfg          = [];
        cfg.hpfilter = 'yes';
        cfg.hpfreq   = 1;
        data = ft_preprocessing(cfg,data);

        % =====================================================================
        %   5   REMOVE ARTIFACTS (HEAD MOVEMENT, SQUID JUMPS, MUSCLE ARTIFACTS)
        % =====================================================================

        fprintf('\n ---------- 5 Remove artifacts -----------\n')

        % Load previously created artifact matrices (aligned to block onset!)
        name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/artifact_identification/preICA_artifactMatrices/' ID];
        artifacts = ['Artifacts_' ID '_Block_' num2str(block) '.mat'];
        load([name '/' artifacts]);
        %______________check metal artifacts ___ GM____________
        
        % Percentage of samples containing metal artifacts of all samples
        % 1st step: Finding out number of samples in each block 
        blocksamples = file_info(idx_file).blocktimes(block,:);
        blocksamples = blocksamples*400; %from seconds to samples
        number_samples = blocksamples(1,2)-blocksamples(1,1)+8000; %compute total number of samples within a block 
       % account for datasets with no metal artifacts (to avoid dividing by 0)
        if isempty(artifact_metal)
            percentage_metal_all(block,1) = 0; 
        else
        %create matrix in which zeros mean there is no artifact detected
        %and  1 means there is an artifact detected in this sample 

        metal = zeros(int32(number_samples),1); 
        
        for l = 1:length(artifact_metal(:,1))
            metal(artifact_metal(l,1):artifact_metal(l,2)) = 1; 
        end 
        % How many samples contain metal artifacts? 
        samples_w_metal(block)=length(metal(metal==1)); 
        percentage_metal_all (block,1) = (samples_w_metal(block)/number_samples)*100; 
        end 
        
        %________________________________________________________
        % Concatenate artifact matrices
        if percentage_metal_all (block,1) >= 60 
            artifact_HMJ = [artifact_headM; artifact_Muscle; artifact_Jump]; % discard metal artifacts when a large part of samples contain them 
            if block==length(file_info(idx_file).blocktimes(:,1)) && mean(percentage_metal_all) >= 60 
            save(['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/metal_artifacts/', ID, '_metal.mat'],'percentage_metal_all'); % save info which datasets contain a lot of metal artifacts
            end
        else
        artifact_HMJ = [artifact_headM; artifact_Muscle; artifact_Jump; artifact_metal]; 
        end 
        % Sort according to time stamp
        artifact_HMJ = sort(artifact_HMJ,'ascend');

        % Merge consecutive artifacts into 1 if they're Xms together
        coalesce = 1;
        if size(artifact_HMJ,1)>1
            cmsmp = artifact_HMJ(1,:);
            for b = 1:size(artifact_HMJ,1)-1
                if artifact_HMJ(b+1,1) - cmsmp(end,2) < coalesce*data.fsample
                    cmsmp(end,2) = artifact_HMJ(b+1,2);
                else
                    cmsmp(end+1,:) = artifact_HMJ(b+1,:);
                end
            end
            artifact_HMJ = cmsmp; clear cmsmp
        end
        

         % Initialize vectors of artifactual samples from each modality
        artsmp_Ygaze = zeros(1,max(data.sampleinfo)); % initial for vEOG
        artsmp_Sacc = artsmp_Ygaze; % Initial for sacc
        artsmp_HMJ = artsmp_Ygaze; % Initial for HMJ
       

        % Mark bad samples 
        for b = 1:size(artifact_HMJ,1)
            artsmp_HMJ(artifact_HMJ(b,1):artifact_HMJ(b,2)) = 1;
        end
        for b = 1:size(artifact_Ygaze,1)
            artsmp_Ygaze(artifact_Ygaze(b,1):artifact_Ygaze(b,2)) = 1;
        end
        for b = 1:size(artifact_saccade,1)
            artsmp_Sacc(artifact_saccade(b,1):artifact_saccade(b,2)) = 1;
        end
       

        % For blinks and saccades, keep only samples where HMJ will be removed
        artsmp_V = artsmp_Ygaze(artsmp_HMJ~=1);
        artsmp_S = artsmp_Sacc(artsmp_HMJ~=1);

        % Redefine start and end points for each artifactual epoch
        % New Ygaze artifact vector
        artifact_blink = []; artflag = 0; i=find(artsmp_V==1,1,'first');
        while i<length(artsmp_V)
            if artflag == 0 && artsmp_V(i)==1
                artifact_blink(end+1,1) = i;
                artflag = 1;
            elseif artflag == 1 && artsmp_V(i)==0
                artifact_blink(end,2) = i-1;
                artflag = 0;
            end
            i = i+1;
        end
        if artflag==1, artifact_blink(end,2) = i; end  

        % New saccade artifact vectors
        artifact_sacc = []; artflag = 0; i=find(artsmp_S==1,1,'first');
        while i<length(artsmp_S)
            if artflag == 0 && artsmp_S(i)==1
                artifact_sacc(end+1,1) = i;
                artflag = 1;
            elseif artflag == 1 && artsmp_S(i)==0
                artifact_sacc(end,2) = i-1;
                artflag = 0;
            end
            i = i+1;
        end
        if artflag==1, artifact_sacc(end,2) = i; end  

        % Delete head movements, muscle and jump artifacts
        [data,indNanTrl] = delete_artifact_Timescale(artifact_HMJ,data);

        % Adjust the data.sampleinfo to correspond to new trial size
        sampleinfo_oldHMJ = data.sampleinfo; % Remember what sampleinfo was previous to HMJ artifact deletion
        data.sampleinfo = [1 length(data.trial{1})];

        % =====================================================================
        %   6  CONCATENATE CLEANED DATA AND REALIGNED BLINK/SACCADE MATRICES
        %      ACROSS BLOCKS
        % =====================================================================

        fprintf('\n ---------- 6 Concatenating cleaned data -----------\n')

        artifact_HMJ_all(block).HMJ = artifact_HMJ;
        artifact_blink_all(block).blink = artifact_blink;
        artifact_sacc_all(block).sacc = artifact_sacc;
        sinfo_old_all(block).sold = sampleinfo_old;
        fs_old_all(block).fs = fsample_old;
        sinfo_new_all(block).snew = sampleinfo_new;
        sinfo_HMJ_all(block).sHMJ = sampleinfo_oldHMJ;

        fprintf('\n ---------- Adjusting samples -----------\n')
        % Offset the artifact vectors by the number of samples in the previous
        % block (so that they're all aligned when concatenated)

        if block == 1 
            data_struct = data;
            artifact_HMJ_fin = artifact_HMJ;
            artifact_blink_fin = artifact_blink;
            artifact_sacc_fin = artifact_sacc;
            sinfo_adjusted_all(block).sinfo = data.sampleinfo;
        elseif block > 1 % Adjust/Offest all blocks > 1
            artifact_HMJ_fin = [artifact_HMJ_fin; artifact_HMJ + sample_counter];
            artifact_blink_fin = [artifact_blink_fin; artifact_blink + sample_counter];
            artifact_sacc_fin = [artifact_sacc_fin; artifact_sacc + sample_counter];
            sinfo_adjusted_all(block).sinfo = data.sampleinfo + sample_counter;
        end

        sample_counter = [sample_counter + max(data.sampleinfo)];
        trial_data = [trial_data data.trial{:}];

        fprintf('\n ---------- Cleaning up -----------\n')
        clear data artifact_Ygaze artifact_saccade artifact_Muscle 
        clear artifact_Jump artifact_HMJ artifact_blink artifact_sacc
        clear artsmp_blink artsmp_Sacc artsmp_HMJ artsmp_V artsmp_S artsmp_Ygaze
        clear blinks metal blocksamples number_samples 
    end

        % =====================================================================
        %   7  SAVE CLEANED AND CONCATENATED DATA & REALIGNED ARTIFACT MATRICES
        % =====================================================================

        fprintf('\n ---------- 7 Saving data -----------\n')

        % Create data structure with cleaned & concatenated data across blocks
        data = data_struct;
        data.sampleinfo = [1 sample_counter];
        data.time = {linspace(0, sample_counter/data.fsample, sample_counter)};
        data.trial = {trial_data};

        % Re-name artifact matrices to simpler names
        artifact_HMJ = artifact_HMJ_fin;
        artifact_saccade = artifact_sacc_fin;
        artifact_blink = artifact_blink_fin;
        % correct if trials and sampleinfo aren't properly aligned
        if length(data.trial{1,1}(1,:))~= data.sampleinfo(2)
            data.sampleinfo(2) = length(data.trial{1,1}(1,:));
            data.time{1,1}(:,data.sampleinfo(2)+1:end)=[];
        end 
       
        % Define path/folder to save cleaned data
        cd('/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/data_clean/') 
        name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/data_clean/' ID];

        if 7==exist(name,'dir')
            cd(name)
        else
            mkdir(name)
            cd(name)
        end
        
        % Save concatenated data file
        datastore = ['data_clean_' ID];
        save(datastore,'data','-v7.3');

        % Save the artifacts (aligned)
        artstore = ['all_artifacts_' ID];
        save(artstore,'artifact_HMJ','artifact_blink','artifact_saccade')

        % Save other file information (just in case!)
        sstore = ['all_sampleinfo_' ID];
        save(sstore,'sinfo_old_all', 'fs_old_all', 'sinfo_new_all',...
        'sinfo_HMJ_all','sinfo_adjusted_all')

        % Cleaning up before running the ICA
        fprintf('\n ---------- Cleaning up -----------\n')
        clearvars -except data ID filein name


    %%
    % =====================================================================
    %   8  ICA
    % =====================================================================

    fprintf('\n ---------- 8 Performing ICA -----------\n')

    try

        cfg = [];
        cfg.channel = 'MEG';
        dataMEG = ft_selectdata(cfg,data);

        % Run the ICA
        %cfg.numcomponent = 25;
        comp = ft_componentanalysis(cfg, dataMEG);

        % Define path/folder to save comp
        cd('/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/comp_ICA/') 
        name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/comp_ICA/' ID];

        if 7==exist(name,'dir')
            cd(name)
        else
            mkdir(name)
            cd(name)
        end

        % Save comp
        savecomp = ['comp_' ID];
        save(savecomp,'comp','-v7.3')

    catch err

        fprintf(['\n ----------Error performing ICA Subject #' ID ' -----------\n'])
        cd('/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/comp_ICA/')
        fid=fopen('logfile_ICA','a+');
        c = clock;
        fprintf(fid,sprintf('\n\n\n\nNew entry for %s at %i/%i/%i %i:%i\n\n\n\n',...
            ID,fix(c(1)),fix(c(2)),fix(c(3)),fix(c(4)),fix(c(5))))

        fprintf(fid,'%s',err.getReport('extended','hyperlinks','off'))

        fclose(fid)    
    end
end
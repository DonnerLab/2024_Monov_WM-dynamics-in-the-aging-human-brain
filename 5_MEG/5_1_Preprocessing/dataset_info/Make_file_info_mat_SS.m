% FILE INFO SCRIPT 2: ADD BLOCK/EVENT INFORMATION AND SAVE EOG/EL PLOTS

function Make_file_info_mat_SS(s)

% FieldTrip 
addpath '/mnt/homes/home032/aarazi/fieldtrip-20201009'  % tell Matlab where FieldTrip is
ft_defaults

% Path/subject information
subj = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '35' '36' 'P03' '37' '38' '39' '40' '41' '42' '43' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54'}; % 01 and 02 are 'C' subjects where MEG dataset has different naming scheme - should be 24 and 25, respectively
megpath = '/mnt/homes/home028/gmonov/meg_data/';
behavpath = '/mnt/homes/home028/gmonov/behav_eye_data/';
savepath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/';
 
ntrials = 63;

% Load initial master file
load([savepath,'Master_file.mat'])

% Fill in some more detailed informtion (trial exceptions, etc.)

for f = 1:length(file_info)
    
    if strcmp({file_info(f).subj},subj{s})  
        
        fprintf('Processing file %s...\n',file_info(f).name)
        
        if strcmp({file_info(f).subj},subj{s})==1 && file_info(f).sess==1 % 1st session
            file_info(f).trialex = ones(length(file_info(f).blocknums),ntrials);
            
            % Session 1 -Behavior
            % loading behavioral data
            fullpath = [behavpath,subj{s},filesep,'S1',filesep,'Behaviour',filesep];
            files = dir([fullpath,filesep,'*.mat'])
            behav1_1=[]; % behavior session 1, block 1
            load([fullpath,filesep,files(1).name]);
            behav1_1 = [behav1_1; Behav.']; % transpose behavioral data matrix
            behav1_2=[]; % behavior session 1, block 2
            load([fullpath,files(2).name]);
            behav1_2 = [behav1_2; Behav.'];
            if length(file_info(f).blocknums)>2;
                behav1_3=[]; % behavior session 1, block 3
                load([fullpath,files(3).name]);
                behav1_3 = [behav1_3; Behav.'];
            end
            if length(file_info(f).blocknums)>3;
                behav1_4=[]; % behavior session 1, block 4 --> only exists for '25' % '10'
                load([fullpath,files(4).name]);
                behav1_4 = [behav1_4; Behav.'];
            end
            
            % Initial trial exceptions: 
            % excluding trials within a block (too short or too long
            % RT, wrong button press, or no trial (block ended earlier)
            % calculating RT for all blocks within one session (R)

            if  length(file_info(f).blocknums)== 4
                R = [behav1_1(7,:), behav1_2(7,:), behav1_3(7,:), behav1_4(7,:)];
            elseif length(file_info(f).blocknums)== 1
                R = behav1_1(7,:);
            elseif length(file_info(f).blocknums)== 2
                R = [behav1_1(7,:), behav1_2(7,:)];
            else R = [behav1_1(7,:), behav1_2(7,:), behav1_3(7,:)];
            end
            
            % Session 1, Block 1
            RT = behav1_1(7,:);
            behav1_1(:,RT>=(mean(R)+4*std(R))) = 999; % RT too long
            behav1_1(:,RT<=0.2)=999; % RT too fast
            choice = behav1_1(5,:); % wrong button press
            behav1_1(:,choice==99)=999;
            marks1_1 = ones(1,ntrials);
            differ = ntrials-length(behav1_1); %block ended earlier
            marks1_1(end-(0:(differ-1)))=2; %set last trials in an earlier ended block to nan's later
            marks1_1(1,behav1_1(1,:)==999) = 0;% set all 999 marked trials to zero trials in this block
            
            % Session 1 Block 2
            RT = behav1_2(7,:);
            behav1_2(:,RT>=(mean(R)+4*std(R))) = 999;
            behav1_2(:,RT<=0.2)=999;
            choice = behav1_2(5,:);
            behav1_2(:,choice==99)=999;
            marks1_2 = ones(1,ntrials);
            differ = ntrials-length(behav1_2);
            marks1_2(end-(0:(differ-1)))=2;
            marks1_2(behav1_2(1,:)==999)=0;
            
            % Session 1, Block 3
            if length(file_info(f).blocknums)>2
                RT = behav1_3(7,:);
                behav1_3(:,RT>=(mean(R)+4*std(R))) = 999;
                behav1_3(:,RT<=0.2)=999;
                choice = behav1_3(5,:);
                behav1_3(:,choice==99)=999;
                marks1_3 = ones(1,ntrials);
                differ = ntrials-length(behav1_3);
                marks1_3(end-(0:(differ-1)))=2;
                marks1_3(behav1_3(1,:)==999)=0;
            end
            % Session 1, Block 4 
            if length(file_info(f).blocknums)>3
                RT = behav1_4(7,:);
                behav1_4(:,RT>=(mean(R)+4*std(R))) = 999;
                behav1_4(:,RT<=0.2)=999;
                choice = behav1_4(5,:);
                behav1_4(:,choice==99)=999;
                
                marks1_4 = ones(1,ntrials);
                differ = ntrials-length(behav1_4);
                marks1_4(end-(0:(differ-1)))=2;
                marks1_4(behav1_4(1,:)==999)=0;
            end
            
            %Concatenating marks
            
            if  length(file_info(f).blocknums)== 4
                marks1 = [marks1_1; marks1_2; marks1_3; marks1_4];
            elseif length(file_info(f).blocknums)== 1
                marks1 = marks1_1;
            elseif length(file_info(f).blocknums)== 2
                marks1 = [marks1_1; marks1_2];
            else  marks1 = [marks1_1; marks1_2; marks1_3];
            end
            
            % Transforming marks into trial exclusions
            
                file_info(f).trialex(marks1==0)= 0;
                file_info(f).trialex(marks1==2)= nan;
            
        end
        
        % Session 2 - Behavior
        %___________________________________
        
        if strcmp(file_info(f).subj,subj{s}) && file_info(f).sess==2 % 2nd session
            file_info(f).trialex = ones(length(file_info(f).blocknums),ntrials);  % initializing trial exception matrix (block*trial, 1=present, 0=missing)
            
            fullpath = [behavpath,subj{s},filesep,'S2',filesep,'Behaviour',filesep];
            files = dir([fullpath,filesep,'*.mat']);
            behav2_1=[]; % behavior session 2, block 1
            load([fullpath,files(1).name]);
            behav2_1 = [behav2_1; Behav.']; %transposing behavioral matrix
            behav2_2=[]; % behavior session 2, block 2
            load([fullpath,files(2).name]);
            behav2_2 = [behav2_2; Behav.'];
            behav2_3=[]; % behavior session 2, block 3
            if length(file_info(f).blocknums) > 2
                load([fullpath,files(3).name]);
                behav2_3 = [behav2_3; Behav.'];
            end
            % R is reaction time of all blocks concatenated, RT stands for
            % each block 
            clear R
            if  length(file_info(f).blocknums)== 4
                R = [behav2_1(7,:), behav2_2(7,:), behav2_3(7,:), behav2_4(7,:)];
            elseif length(file_info(f).blocknums)== 1
                R = behav2_1(7,:);
            elseif length(file_info(f).blocknums)== 2
                R = [behav2_1(7,:), behav2_2(7,:)];
            else R = [behav2_1(7,:), behav2_2(7,:), behav2_3(7,:)];
            end
            
            % Session 2, Block 1
            RT = behav2_1(7,:);
            behav2_1(:,RT>=(mean(R)+4*std(R))) = 999; % RT too long
            behav2_1(:,RT<=0.2)=999; % RT too fast
            choice = behav2_1(5,:); % wrong button press
            behav2_1(:,choice==99)=999;
            
            marks2_1 = ones(1,ntrials);
            differ = ntrials-length(behav2_1); %block ended earlier
            marks2_1(end-(0:(differ-1)))=2; %set trials in an earlier ended block to nan's
            marks2_1(1,behav2_1(1,:)==999) = 0;% set all 999 marked trials to zero trials in this block
            
            % Session 2, Block 2
            RT = behav2_2(7,:);
            behav2_2(:,RT>=(mean(R)+4*std(R))) = 999;
            behav2_2(:,RT<=0.2)=999;
            choice = behav2_2(5,:);
            behav2_2(:,choice==99)=999;
            
            marks2_2 = ones(1,ntrials);
            differ = ntrials-length(behav2_2);
            marks2_2(end-(0:(differ-1)))=2;
            marks2_2(behav2_2(1,:)==999)=0;
            
            % Session 2, Block 3
            if length(file_info(f).blocknums)>2 
                RT = behav2_3(7,:);
                behav2_3(:,RT>=(mean(R)+4*std(R))) = 999;
                behav2_3(:,RT<=0.2)=999;
                choice = behav2_3(5,:);
                behav2_3(:,choice==99)=999;
                
                marks2_3 = ones(1,ntrials);
                differ = ntrials-length(behav2_3);
                marks2_3(end-(0:(differ-1)))=2;
                marks2_3(behav2_3(1,:)==999)=0;
            end
            
            if length(file_info(f).blocknums)==1
                marks2 = marks2_1; 
            elseif length(file_info(f).blocknums)== 2
                marks2 = [marks2_1; marks2_2];
            elseif length(file_info(f).blocknums)== 3
                marks2 = [marks2_1; marks2_2; marks2_3];
            end
            
                file_info(f).trialex(marks2==0)= 0;
                file_info(f).trialex(marks2==2)=nan;
           
            
        end

        % Mark blocks that were flagged by the experimenter during data collection for any reason ('6')
        % so that these can be easily identified later 
        %____________________________________________________________________
        % Excluding 1st block -->'01': 1st block cancelled after 2 trials; '02': participant doesn't
        % understand button boxes in 1st block; '25': 1st block cancelled because
        % there were many mistakes, participant didn't understand task
        % rules 
        
        if strcmp(file_info(f).subj,'01') && file_info(f).sess==1 || strcmp(file_info(f).subj,'25') && file_info(f).sess==1
            for t=1:length(file_info(f).trialex(1,:))
                if ~isnan(file_info(f).trialex(1,t))
                    file_info(f).trialex(1,t) = 0;
                end
            end
        end
        
        
        if strcmp(file_info(f).subj,'02') && file_info(f).sess==1 || strcmp(file_info(f).subj,'41') && file_info(f).sess==1
            for t=1:length(file_info(f).trialex(1,:))
                if ~isnan(file_info(f).trialex(1,t)) && file_info(f).trialex(1,t) > 0
                    file_info(f).trialex(1,t) = 6;
                end
            end
        end
        
        
        if strcmp(file_info(f).subj,'33') && file_info(f).sess==1 % '33': wrong button many times in the third block
            for t=1:length(file_info(f).trialex(3,:))
                if ~isnan(file_info(f).trialex(3,t)) && file_info(f).trialex(3,t) > 0
                    file_info(f).trialex(3,t) = 6;
                end
            end
        end
        
        if strcmp(file_info(f).subj,'09') && file_info(f).sess==1 % '09': participant doesn't understand task
            for t=1:length(file_info(f).trialex(1:2,:))
                for loop=1:2
                 if ~isnan(file_info(f).trialex(loop,t)) && file_info(f).trialex(loop,t) > 0
                    file_info(f).trialex(loop,t) = 6;
                 end
                end 
            end
        end
        
        if strcmp(file_info(f).subj,'09') && file_info(f).sess==2  % '09': participant forgets task rules in the 3rd block
            for t=1:length(file_info(f).trialex(3,:))
                if ~isnan(file_info(f).trialex(3,t)) && file_info(f).trialex(3,t) > 0
                    file_info(f).trialex(3,t) = 6;
                end
            end
        end
        if strcmp(file_info(f).subj,'10')  % '10': participant doesn't understand task
            for t=1:length(file_info(f).trialex(1:4,:))
                for loop = 1:4
                if ~isnan(file_info(f).trialex(loop,t)) && file_info(f).trialex(loop,t) > 0
                    file_info(f).trialex(loop,t) = 6;
                end
                end 
            end
        end
        if strcmp(file_info(f).subj,'52') % '52': participant doesn't understand task
            for t=1:length(file_info(f).trialex(1:2,:))
                for loop = 1:2
                if ~isnan(file_info(f).trialex(loop,t)) && file_info(f).trialex(loop,t) > 0
                    file_info(f).trialex(loop,t) = 6;
                end
                end
            end 
        end
         if strcmp(file_info(f).subj,'41')&& file_info(f).sess==2  % '41': participant lapses often, falls asleep sometimes, forgets task rules
            for t=1:length(file_info(f).trialex(1:2,:))
                for loop = 1:2
                if ~isnan(file_info(f).trialex(loop,t)) && file_info(f).trialex(loop,t) > 0
                    file_info(f).trialex(loop,t) = 6;
                end
                end 
            end
        end
 
        % Read meg events for this file 
        event = ft_read_event([megpath,file_info(f).name]);
 
        real_events = [];
        for e = 1:length(event)
            if strcmp(event(e).type,'UPPT001')
                real_events(end+1,1) = event(e).value;      % event type
                real_events(end,2) = event(e).sample/1200;  % event time (in secs relative to rec onset)
            end
        end
        
        if isempty(real_events)==1
            fprintf('No events recorded in %s...\n',file_info(f).name)
            file_info(f).eventsM = 'no events recorded';
            file_info(f).blocktimes = nan;
            file_info(f).trialex = 'no events recorded';
            file_info(f).rest = 0;
            file_info(f).resttime = nan;
            
            if strcmp(file_info(f).subj,subj{s})&& file_info(f).sess == 1
                save([savepath,filesep,'ss_mat_files',filesep,subj{s},'_1.mat'],'file_info')
            else
                save([savepath,filesep,'ss_mat_files',filesep,subj{s},'_2.mat'],'file_info')
            end
        else
            % testing which triggers are in the dataset 
        which_triggers=unique(real_events(:,1));
        
            % remove repeating adjacent events (can happen in particular for '1's (block starts))
            coalesce = 0.06;  % distance, in seconds, between events of same type to merge
            
            reps = find(diff(real_events(:,1))==0); % find event repetitions
            reptimes = diff(real_events(:,2)); reptimes = reptimes(reps); % compute time differences between repetitions
            for_c=[];
            for t = 1:length(reps)
                if abs(reptimes(t)) < coalesce
                    for_c(end+1) = reps(t)+1;  % marking repetition to be removed
                end
            end
            real_events(for_c,:)=[];  % remove repetitions
            
            % dealing with exceptions
            % 24: no associated behavioral file for very first block, just
            % deleting events here --> no need to change more stuff in the
            % next step 
            % 01: first block cancelled after a few trials, but behav file
            % exists (taken care of in make_final_file_info) 
            
            if strcmp(file_info(f).subj,'24')||strcmp(file_info(f).subj,'01')&&file_info(f).sess==1  % 1st block is a mess/consists only of few trials, and there might be no associated behavioural file - just removing events here
                starts = find(real_events(:,1)==1);
                real_events(1:(starts(2)-1),:) = [];
                clear starts
            end
            %'25': for first block there seems to be no associated
            %behavioral file, for the second block, which only consists of
            %a few trials events should also be deleted
            if  strcmp(file_info(f).subj,'25')
                starts = find(real_events(:,1)==1);
                real_events(1:(starts(3)-1),:) = [];
                clear starts
            end 
            
            % add events to master file
            file_info(f).eventsM = real_events;                                    % MEG events (for later debugging)
            
            % Getting times of resting block
            if ~isempty(real_events(:,1)==3) % 2 is trigger for start of resting block but also end of task block (terrible design), so need to be careful; 3 is end of resting block
                rest_end = find(real_events(:,1)==3);
                rest_start = rest_end; while real_events(rest_start,1)~=2, rest_start=rest_start-1; end
                file_info(f).resttime = [real_events(rest_start,2) real_events(rest_end,2)];
            else
                twos = find(real_events(:,1)==2);  % if there's no '3', then taking rest start to be '2' which occurs soonest after a preceding '2' (i.e. block end)
                rest_start = twos(find(diff(twos)==min(diff(twos)))+1);
                file_info(f).resttime = [rest_start nan];
            end
            file_info(f).rest = 1; %resting block present(1=yes;0=no), assuming that rest bloc is there                                            % resting block present (1=yes, 0=no) - here I assume that a resting block is present in all datasets
            
            % Segment meg events into blocks
            bstarts = find(real_events(:,1)==1)'+1;  % get first logged event in each block
            if real_events(1,1)~=1 || bstarts(1)>20, bstarts = [1 bstarts]; end  % adding early block start in cases where this was missed in meg
            bstarts(end+1) = size(real_events,1)+1;
            
            bends = find(real_events(:,1)==2)';  % get last logged event in each block
            bends(bends==rest_start) = [];
            bends = bends-1;
            
            % sometimes bstarts and bends are repeated
            repeated_bstarts=[];
            repeated_bends=[];
            repeated_bstarts=find(diff(bstarts(1,:))<=2)
            repeated_bends=find(diff(bends(1,:))<=2)
            for k=1:length(repeated_bends)
                bends(repeated_bends(1,k)+1)=0.99;
            end
            bends(bends==0.99)=[];
            for p=1:length(repeated_bstarts)
                bstarts(repeated_bstarts(1,p))=0.99;
            end
            bstarts(bstarts==0.99)=[];
            % for participant 29 there is a bstart recorded
            % which seems to be in the middle of the block--> just removing
            % this event
            if strcmp(file_info(f).subj,'29')&&file_info(f).sess==1
                
                if bstarts(3)<bends(2)
                    bstarts(3)=[];
                end
            end
            
            % add block start/end times in secs relative to recording start
            file_info(f).blocktimes = [];
            %if start of block event is not in there event 1 (first logged
            %event is the start of block) --> set to 2
            if bstarts(1)==1
                bstarts_c=bstarts
                bstarts_c(1)=2;
                file_info(f).blocktimes = [real_events((bstarts_c(1:end-1)-1)',2) real_events((bends+1)',2)];
            else  file_info(f).blocktimes = [real_events((bstarts(1:end-1)-1)',2) real_events((bends+1)',2)];
            end
            % Loop through blocks, organise meg events, load behaviour, and compare
            file_info(f).check=[];
            file_info(f).check(1:length(bends)) = 0;
            if length(file_info(f).blocknums) ~= length(bstarts)-1 || length(bstarts)~=length(bends)+1  % if number of blocks in MEG and behaviour don't match
                file_info(f).check(1:length(bends)) = 1;
            else
                for b = 1:length(bstarts)-1  % loop through blocks
                    cevents = real_events(bstarts(b):bends(b),:);  % pulling only events for this block
                    soaM = diff(cevents(cevents(:,1)==11,2));  % calculate stimulus onset asynchronies from MEG events
                    if strcmp(file_info(f).subj,'24')
                        load([behavpath,file_info(f).subj,'/S1/Behaviour/','01C_',num2str(file_info(f).sess),'_',num2str(file_info(f).blocknums(b)),'.mat'])
                    elseif strcmp(file_info(f).subj,'25')
                        load([behavpath,file_info(f).subj,'/S1/Behaviour/','02C_',num2str(file_info(f).sess),'_',num2str(file_info(f).blocknums(b)),'.mat'])
                    elseif file_info(f).sess == 1
                        load([behavpath,file_info(f).subj,'/S1/Behaviour/',file_info(f).subj,'_',num2str(file_info(f).sess),'_',num2str(file_info(f).blocknums(b)),'.mat'])
                    else load([behavpath,file_info(f).subj,'/S2/Behaviour/',file_info(f).subj,'_',num2str(file_info(f).sess),'_',num2str(file_info(f).blocknums(b)),'.mat'])   % load behaviour
                    end
                    %selecting trialex trials that are not yet marked as
                    %nans due to early block termination
                    select_trialex=~isnan(file_info(f).trialex(b,:));
                    soaB = diff(Behav(:,8));  % calculate stimulus onset asynchronies from behavioural files
                    if length(soaM)~=length(soaB)  % if there's a different number of trials, test if this is due to missing trials @ start/end
                        lag = length(soaB)-length(soaM);  % lag<0 means more trials in meg, lag>0 means more trials in behaviour (should never really be former...)
                        if lag > 0
                            if isempty(find(abs(soaM-soaB((1+lag):end))>0.1))  % trials missing @ start of meg
                               file_info(f).trialex(b,1:lag) = 3;
                                  
                            elseif isempty(find(abs(soaB(1:end-lag)-soaM)>0.1))   % trials missing @ end of meg
                                file_info(f).trialex(b,(end-lag+1):end) = 4; % GM: added since otherwise one additional trial to the lag would have been removed 
                            else
                                file_info(f).check(b) = 2;
                            end
                        elseif lag == -1  % when block terminated by user, there's 1 trial more in MEG
                            if isempty(find(abs(soaM(1:length(soaB))-soaB)>0.1))
                                file_info(f).trialex(b,(length(soaB)+2):end) = 4;  % if all other trials in both are matched for SOA, just mark extra one in MEG for exclusion; GM: added +2 since soa is numerically one less than actual trials and same as above 
                            else file_info(f).check(b) = 3;  % otherwise mark for checking
                            end
                        else
                            file_info(f).check(b) = 3;
                        end
                    elseif ~isempty(find(abs(soaM-soaB)>0.1))   % if there's an SOA discrepancy greater than 0.1s, mark for checking
                        file_info(f).check(b) = 4;
                    end
                    
                end
            end
           
            % Load meg data for visualization of EOG, ECG & Eyelink data
            plot_time = 60;  % amount of time over which to plot data (in seconds)
            
            cfg=[];
            cfg.dataset = [megpath,file_info(f).name];
            cfg.event = event;
            cfg.continuous = 'yes';
            cfg.demean = 'yes';
            cfg.channel = {'EEG057' 'EEG058' 'EEG059' 'UADC001' 'UADC002' 'UADC003' 'UADC004'};
            data = ft_preprocessing(cfg);
            
            good_s = zeros(1,size(data.trial{1},2));
            for b = 1:length(bstarts)-1
                good_s(round(real_events(bstarts(b),2)*1200):(round(real_events(bstarts(b),2)*1200)+5*60*1200)) = 1;  % only using from block start to plus 5 minutes per block
            end
            if length(good_s)>size(data.trial{1},2), good_s = good_s(1:size(data.trial{1},2)); end
            r = corrcoef(data.trial{1}(:,logical(good_s))');  % calculating correlation coefficients between channels
            
            if length(bstarts)>4, pb = 4; else pb = 1; end  % pick block to plot
            % see if ylim works in second block of participant 02
            if strcmp(file_info(f).subj,'02')&&file_info(f).sess==1
                pb=2;
            end 
                
            cols = {'k','k','k','b','b','b','b'};
            figure('Position',[80, 100, 1500, 1000]),
            for c = 1:length(cfg.channel)
                subplot(length(cfg.channel),1,c), hold on
                plot((1:length(data.trial{1}(c,round(real_events(bstarts(pb),2)*1200):round((real_events(bstarts(pb),2)+plot_time)*1200))))./1200,...
                    data.trial{1}(c,round(real_events(bstarts(pb),2)*1200):round((real_events(bstarts(pb),2)+plot_time)*1200)),cols{c})
                title(sprintf('%s, R = %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',cfg.channel{c},r(c,:))), xlim([0 plot_time])
                if min(data.trial{1}(c,round(real_events(bstarts(pb),2)*1200):round((real_events(bstarts(pb),2)+plot_time)*1200))) == max(data.trial{1}(c,round(real_events(bstarts(pb),2)*1200):round((real_events(bstarts(pb),2)+plot_time)*1200)))
                    ylim([-1 1]);
                else
                ylim([min(data.trial{1}(c,round(real_events(bstarts(pb),2)*1200):round((real_events(bstarts(pb),2)+plot_time)*1200))) max(data.trial{1}(c,round(real_events(bstarts(pb),2)*1200):round((real_events(bstarts(pb),2)+plot_time)*1200)))])
                end 
                if c==length(cfg.channel), xlabel('Time relative to 1st block onset (s)'), end
            end
            
            print([savepath,filesep,'60s_plots',filesep,file_info(f).name,'.png'],'-dpng');
            close all
            
           % Correct if the early terminated trials were overwritten,
           % important for trial analysis 
           if file_info(f).sess==1
            file_info(f).trialex(marks1==2)=nan;
           elseif file_info(f).sess==2
            file_info(f).trialex(marks2==2)=nan;
           end 
            
            % Save subject-specific version

            if strcmp(file_info(f).subj,subj{s})&& file_info(f).sess == 1
                save([savepath,filesep,'ss_mat_files',filesep,subj{s},'_1.mat'],'file_info')
                
            elseif strcmp(file_info(f).subj,subj{s})&& file_info(f).sess == 2
                save([savepath,filesep,'ss_mat_files',filesep,subj{s},'_2.mat'],'file_info')
            end

        end
          
    end
end
end
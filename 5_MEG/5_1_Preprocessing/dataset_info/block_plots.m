% Plotting of entire blocks of eye data to make sure that
% the quality is good, on which to base the decision whether one should use the EOG or 
% eye tracker data for the preprocessing 

function block_plots(s)
 
% FieldTrip stuff
addpath '/mnt/homes/home032/aarazi/fieldtrip-20201009'  % FieldTrip 
addpath '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/file_info_scripts/'
ft_defaults
 
% Path/subject stuff
subj = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '35' '36' 'P03' '37' '38' '39' '40' '41' '42' '43' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54'};   % 01 and 02 are 'C' subjects where MEG dataset has different naming scheme - should be 24 and 25, respectively
megpath = '/mnt/homes/home028/gmonov/meg_data/';
behavpath = '/mnt/homes/home028/gmonov/behav_eye_data/';
savepath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/';

ntrials = 63;

% Load master file
load([savepath,'Master_file.mat'])

for f = 1:length(file_info);
    if strcmp(file_info(f).subj,subj{s})
    
    fprintf('Processing file %s...\n',file_info(f).name)
    
    %Read meg events for this file 
    event = ft_read_event([megpath,file_info(f).name]);
    real_events = [];
    
    for e = 1:length(event)
        if strcmp(event(e).type,'UPPT001')
            real_events(end+1,1) = event(e).value;      % event type
            real_events(end,2) = event(e).sample/1200;  % event time (in secs relative to rec onset)
        end
    end
    
if isempty(real_events)==1
    disp 'no events recorded'
else 
    
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
   if  strcmp(file_info(f).subj,'24')||strcmp(file_info(f).subj,'25')||strcmp(file_info(f).subj,'01')&&file_info(f).sess==1 % 1st block is a mess/consists only of few trials, and there might be no associated behavioural file - just removing events here
       starts = find(real_events(:,1)==1);
       real_events(1:(starts(2)-1),:) = [];
       clear starts
    end
    
    % add events to master file
    file_info(f).eventsM = real_events;                                    % MEG events (for later debugging)

    % Getting times of resting block
    if ~isempty(real_events(:,1)==3) % 2 is trigger for start of resting block but also end of task block, so need to be careful; 3 is end of resting block
        rest_end = find(real_events(:,1)==3);
        rest_start = rest_end; while real_events(rest_start,1)~=2, rest_start=rest_start-1; end
        file_info(f).resttime = [real_events(rest_start,2) real_events(rest_end,2)];
    else
        twos = find(real_events(:,1)==2);  % if there's no '3', then taking rest start to be '2' which occurs soonest after a preceding '2' (i.e. block end)
        rest_start = twos(find(diff(twos)==min(diff(twos)))+1);
        file_info(f).resttime = [rest_start nan];
    end
    
    file_info(f).rest = 1;                                             % resting block present (1=yes, 0=no) - here I assume that a resting block is present in all datasets
    
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
            % for participant 29 there is a block start recorded
            % which seems to be in the middle of the block --> just removing
            % this event
            if strcmp(file_info(f).subj,'29')&&file_info(f).sess==1
                if bstarts(3)<bends(2)
                    bstarts(3)=[];
                end 
            end 
    
     % add block start/end times in secs relative to recording start
    file_info(f).blocktimes = [];
    
     % if start of block event is not in there event 1 (first logged
     % event is the start of block) --> set to 2
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
            soaB = diff(Behav(:,8));  % calculate stimulus onset asynchronies from behavioural files
            
            if length(soaM)~=length(soaB)  % if there's a different number of trials, test if this is due to missing trials @ start/end
                lag = length(soaB)-length(soaM);  % lag<0 means more trials in meg, lag>0 means more trials in behaviour (should never really be former...)
                if lag > 0
                    if isempty(find(abs(soaM-soaB((1+lag):end))>0.1))  % trials missing @ start of meg
                        file_info(f).trialex(b,1:lag) = 0;
                    elseif isempty(find(abs(soaB(1:end-lag)-soaM)>0.1))   % trials missing @ end of meg
                        file_info(f).trialex(b,(end-lag):end) = 0;
                    else
                        file_info(f).check(b) = 2;
                    end
                elseif lag == -1  % when block terminated by user, there's 1 trial more in MEG
                    if isempty(find(abs(soaM(1:length(soaB))-soaB)>0.1))
                        file_info(f).trialex(b,(length(soaB)+1):end) = 0;  % if all other trials in both are matched for SOA, just mark extra one in MEG for exclusion
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
    cfg=[];
    cfg.dataset = [megpath,file_info(f).name];
    cfg.event = event;
    cfg.continuous = 'yes';
    cfg.demean = 'yes';
    cfg.channel = {'EEG057';'EEG058';'EEG059';'UADC001';'UADC002';'UADC003';'UADC004'};
    data = ft_preprocessing(cfg);
    
      good_s = zeros(1,size(data.trial{1},2));
            for b = 1:length(bstarts)-1
                good_s(round(real_events(bstarts(b),2)*1200):(round(real_events(bstarts(b),2)*1200)+5*60*1200)) = 1;  
            end  % only using from block start to plus 5 minutes per block
 
    if length(good_s)>size(data.trial{1},2), good_s = good_s(1:size(data.trial{1},2)); end
    r = corrcoef(data.trial{1}(:,logical(good_s))');  % calculating correlation coefficients between channels
    for pb = 1:length(bstarts(1:end-1))   % loop through participants blocks
        plot_time = real_events(bends(pb),2)-real_events(bstarts(pb),2)
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
    print([savepath,filesep,'block_plot',filesep,num2str(pb),filesep,num2str(pb),'_',file_info(f).name,'.png'],'-dpng');
    end
    close all
    end 

end
    end 
end 
end
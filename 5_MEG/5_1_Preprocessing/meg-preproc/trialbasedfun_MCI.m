function [cfgin] = trialbasedfun_MCI(cfgin)
% TRIALBASEDFUN_MCI trial definition function
%   cfgin.trialfun = 'trialfun_MCI'  
%   for trial based analysis of working memory task 


% Path behavioral data 

    dir_behav = '/mnt/homes/home028/gmonov/behav_eye_data/';


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
  
    
    % Pull out important information about dataset to load
    subject = cfgin.subj;
    session = cfgin.session;
    fileNum = cfgin.fileNum;
    ID = cfgin.ID;
    idx_file = find(strcmp(files, cfgin.dataset));
    filein = files{idx_file};
    resttime = file_info(idx_file).resttime*1200; % resample timing back from seconds 
    trialex = file_info(idx_file).trialex';
    

%%
% =====================================================================
%   INFORMATION FROM BEHAVIORAL DATA
% =====================================================================

fprintf('\n\n ---------------- \n Processing behavioral information...\n ---------------- \n\n');

% Check the run 

firstblock = info_EL_blocks(idx_file,1);
lastblock = info_EL_blocks(idx_file,2);
nblock = info_EL_blocks(idx_file,1)-1; 

numFile = fileNum;     
 % List files in these directories
      dir_behaviour = [dir_behav, subject, filesep,  'S', num2str(session), filesep, 'Behaviour',filesep]; 
      files_behav = dir([dir_behaviour,'*.mat']);
block=[];
trials_info = [];

% Exceptions Session 1 
if session == 1
%Excluding 1st block because it has been interrupted 
    if strcmp(subject,'01')||strcmp(subject,'25')
        files_behav(1)=[];
    end
    
end 

% Exceptions Session 2
if session == 2
        if strcmp(subject,'43') % no triggers for the first block in second session 
        files_behav(1)=[];
        end
end 
    allbehav=[];
    
    
for block = firstblock:lastblock % Loop through files/blocks 
    missing_trials = 0; 
    % Define specific path
     load([dir_behaviour, files_behav(block).name]);% Load behavioural data
     Behav(1:length(Behav),end+1)= block; % including information which block in the Behav Matrix
         %add trialcounts --> Behavioral file always starts at first trial
     for trialno=1:length(Behav(:,1))
         Behav(trialno,12)=trialno;
     end
     trialex_block=trialex(:,block);

      for trialexloop=1:length(trialex_block)
         if trialex_block(trialexloop)==3 || isnan(trialex_block(trialexloop)) % trials missing in meg file at start of a block 
          missing_trials(end+1)= trialexloop; %set to nan when these trials only exist in behavioral dataset!
         end 
      end
       missing_trials(missing_trials==0)=[]; 
       select_trials=~ismember(Behav(1:end,12),missing_trials(1:end));  
       Behav = Behav(select_trials,1:end);
       trialex_block(trialex_block==3)=[]; 
       trialex_block(isnan(trialex_block))=[]; 
       Behav(:,11)=trialex_block;
       allbehav = [allbehav;Behav];
end 

 
   trials_info=allbehav;
    

%%
% =====================================================================
%   INFORMATION FROM MEG DATA
% =====================================================================

    % _____________________________________________________________________
    %   TRIGGERS OF INTEREST

trg.block_start = 1;  % start of block
trg.block_end = 2;    % end of block
trg.rest_start = 2; %start resting state
trg.rest_end = 3; %end resting state 

trg.mem_on = 11;   % onset of memoranda
trg.mem_off = 12;  % offset of memoranda

trg.probe_on = 31;    % onset of probe
trg.probe_off = 32;   % offset of probe

trg.resp_left = 41;    % 'left' response
trg.resp_right = 42;   % 'right' response
trg.resp_bad = 43;     % bad response (either double press, or task-irrelevant button)

trg.fb_correct = 51;   % feedback for hit or correct rejection
trg.fb_error = 52;     % feedback for false alarm, miss, mislocalization, or premature response
trg.fb_bad = 53;       % feedback for bad responses

trg.rest_on = 61;    % onset of rest break period
trg.rest_off = 62;   % offset of rest break period


    % Variable for trigger indices
    trg_idx.mem_on = [];               
    trg_idx.mem_off = [];               
    trg_idx.probe_on = [];          
    trg_idx.probe_off = [];            
    trg_idx.resp_left = [];                   
    trg_idx.resp_right = []; 
    trg_idx.resp_bad = []; 
    trg_idx.fb_correct = []; 
    trg_idx.fb_error = []; 
    trg_idx.fb_bad = []; 

   
    trl = []; % function output, matrix to store trial information    
    % _____________________________________________________________________


cd '/mnt/homes/home028/gmonov/meg_data/';
% Load events within specified block times from file_info and 


% Read header and event information

hdr = ft_read_header(cfgin.dataset);
% making clearer what is resting state--> setting the start of resting
% state to 4 instead of 2 
events = file_info(idx_file).eventsM; % events are already conveniently saved in file_info 
events(:,2)=events(:,2).*1200; %resampling from seconds 
events(events(:,2)==resttime(1),1)=4;
trg.rest_start = 4; %correcting for this 
trgVal = events(:,1); 
trgIdx = events(:,2);


% Find events with triggers of interest and save into trg_idx structure
fields = fieldnames(trg_idx);

for trg_type = 1:numel(fields)
    for trg_event = 1:length(trgVal)
        if ismember(trgVal(trg_event),trg.(fields{trg_type}))
            trg_idx.(fields{trg_type})(end+1,1) = trgIdx(trg_event);
        end
    end
end


% block start and block end are repetive sometimes (already taken account
% for in file_info, therefore replacing the information about block timings
% in trg_idx structure)
blocktimes = file_info(idx_file).blocktimes;
blocktimes = blocktimes.*1200;
trg_idx.block_start=blocktimes(:,1);
trg_idx.block_end=blocktimes(:,2);
trg_idx.resp = [trg_idx.resp_right;trg_idx.resp_left;trg_idx.resp_bad];
trg_idx.fb = [trg_idx.fb_correct;trg_idx.fb_error;trg_idx.fb_bad];
trg_idx.resp = sort(trg_idx.resp); 
trg_idx.fb = sort(trg_idx.fb); 

% _________________________________________________________________________
% !!!   CHECK TRIGGERS & DEAL WITH EXCEPTIONS   !!!
fprintf('\n\n ---------------- \n Check triggers...\n ---------------- \n\n');


trg_idx.resp =  trg_idx.resp(trg_idx.resp >= trg_idx.mem_on(1));
trg_idx.resp_left =  trg_idx.resp_left(trg_idx.resp_left >= trg_idx.mem_on(1));
trg_idx.resp_right =  trg_idx.resp_right(trg_idx.resp_right >= trg_idx.mem_on(1));
trg_idx.resp_bad =  trg_idx.resp_bad(trg_idx.resp_bad >= trg_idx.mem_on(1));

trg_idx.resp =  trg_idx.resp(trg_idx.resp <= trg_idx.fb(end));
trg_idx.resp_left =  trg_idx.resp_left(trg_idx.resp_left <= trg_idx.fb(end));
trg_idx.resp_right =  trg_idx.resp_right(trg_idx.resp_right <= trg_idx.fb(end));
trg_idx.resp_bad =  trg_idx.resp_bad(trg_idx.resp_bad <= trg_idx.fb(end));

startBlockSamp = [trg_idx.block_start];
endBlockSamp = [trg_idx.block_end];


% Check whether there were more than one responses/feedbacks on a trial, if
% that's the case only keep the first response/feedback

respnum =[]; 

if length(trg_idx.mem_on)<length(trg_idx.resp)
    %1. identify trial (considering mem_on) where this has happened 
    

        for o=1:length(trg_idx.mem_on)-1
                    
                    
                    respcue = trg_idx.resp(trg_idx.resp>trg_idx.mem_on(o)); %&& trg_idx.resp(trg_idx.resp<trg_idx.mem_on(o+1)))
                    respcue(respcue>=trg_idx.mem_on(o+1))=[]; 
                    respnum = length (respcue); 
                    if respnum >1
       
                   % respcue(2:respnum+1) = trg_idx.resp(trg_idx.resp>trg_idx.mem_on(o))% && trg_idx.resp<trg_idx.mem_on(o+1)); %list all responses
                    
                    respcue = sort(respcue); 
                    respcue2keep = respcue(respcue>trg_idx.probe_on(o)); %only keep responses that were submitted after probe onset 
                    respcue2keep = respcue2keep(1); %only keep first response after probe onset 
                    respcue2remove=respcue;
                    respcue2remove(respcue2remove==respcue2keep)=[]; % remove the first response from removal matrix 
                    trg_idx.resp(ismember(trg_idx.resp,respcue2remove)) = []; 
                    trg_idx.resp_left(ismember(trg_idx.resp_left,respcue2remove)) = [];
                    trg_idx.resp_right(ismember(trg_idx.resp_right,respcue2remove)) = [];
                    trg_idx.resp_bad(ismember(trg_idx.resp_bad,respcue2remove)) = [];
       
   
                   end
        end 
end
 %Still too many responses? Consider what happened at the end 
if length(trg_idx.mem_on)<length(trg_idx.resp)
   respnum = numel(trg_idx.resp(trg_idx.resp>trg_idx.mem_on(end)))
   respcue = [-1]; 
   if respnum >1
       
       respcue(2:respnum+1) = trg_idx.resp(trg_idx.resp>trg_idx.mem_on(end)); 
       respcue(respcue==-1)=[]; 
       respcue = sort(respcue); 
       respcue2keep = respcue(respcue>trg_idx.probe_on(end)); 
       respcue2keep = respcue2keep(1); %only keep first response after probe onset 
       respcue2remove=respcue;
       respcue2remove(respcue2remove==respcue2keep)=[]; % remove the first response from removal matrix 
       trg_idx.resp(ismember(trg_idx.resp,respcue2remove)) = []; 
       trg_idx.resp_left(ismember(trg_idx.resp_left,respcue2remove)) = [];
       trg_idx.resp_right(ismember(trg_idx.resp_right,respcue2remove)) = [];
       trg_idx.resp_bad(ismember(trg_idx.resp_bad,respcue2remove)) = [];
   end 
end 
 
%Do the same for repetitive feedback cues

fbnum =[]; 

% %testing the removal of feeback
% trg_idx.fb=[trg_idx.fb(1)-1;trg_idx.fb];
% trg_idx.fb_correct=[trg_idx.fb(1);trg_idx.fb_correct];

if length(trg_idx.mem_on)<length(trg_idx.fb)
    %1. identify trial (considering mem_on) where this has happened 
    

        for o=1:length(trg_idx.mem_on)-1
                    
                   
                    fbcue = trg_idx.fb(trg_idx.fb>trg_idx.mem_on(o)); %&& trg_idx.fb(trg_idx.fb<trg_idx.mem_on(o+1)))
                    fbcue(fbcue>=trg_idx.mem_on(o+1))=[]; 
                    fbnum = length (fbcue); 
                    if fbnum >1
       
                   % fbcue(2:fbnum+1) = trg_idx.fb(trg_idx.fb>trg_idx.mem_on(o))% && trg_idx.fb<trg_idx.mem_on(o+1)); %list all fbonses
                    
                    fbcue = sort(fbcue); 
                    fbcue2keep = fbcue(fbcue>trg_idx.resp(o)); %only keep fbonses that were submitted after probe onset 
                    fbcue2keep = fbcue2keep(1); %only keep first fbonse after probe onset 
                    fbcue2remove=fbcue;
                    fbcue2remove(fbcue2remove==fbcue2keep)=[]; % remove the first fbonse from removal matrix 
                    trg_idx.fb(ismember(trg_idx.fb,fbcue2remove)) = []; 
                    trg_idx.fb_correct(ismember(trg_idx.fb_correct,fbcue2remove)) = [];
                    trg_idx.fb_error(ismember(trg_idx.fb_error,fbcue2remove)) = [];
                    trg_idx.fb_bad(ismember(trg_idx.fb_bad,fbcue2remove)) = [];
       
   
                    end
        end 
end
 %Still too many feedback cues? Consider what happened at the end 
if length(trg_idx.mem_on)<length(trg_idx.fb)
   fbnum = numel(trg_idx.fb(trg_idx.fb>trg_idx.mem_on(end)))
   fbcue = [-1]; 
   if fbnum >1
       
       fbcue(2:fbnum+1) = trg_idx.fb(trg_idx.fb>trg_idx.mem_on(end)); 
       fbcue(fbcue==-1)=[]; 
       fbcue = sort(fbcue); 
       fbcue2keep = fbcue(fbcue>trg_idx.probe_on(end)); 
       fbcue2keep = fbcue2keep(1); %only keep first feedback after probe onset 
       fbcue2remove=fbcue;
       fbcue2remove(fbcue2remove==fbcue2keep)=[]; % remove the first feedback from removal matrix 
       trg_idx.fb(ismember(trg_idx.fb,fbcue2remove)) = []; 
       trg_idx.fb_correct(ismember(trg_idx.fb_correct,fbcue2remove)) = [];
       trg_idx.fb_error(ismember(trg_idx.fb_error,fbcue2remove)) = [];
       trg_idx.fb_bad(ismember(trg_idx.fb_bad,fbcue2remove)) = [];
   end 
end 

%__________________________________________________________________________

% Initialize counters
ntrial = 0;         % Counter for trials within dataset
% currBlock = 0;      % End sample of current block

fprintf('\n\n ---------------- \n Loop through blocks...\n ---------------- \n\n');
% Loop through the blocks

% for start_block = trg_idx.block_start
for block_loop=1:length(trg_idx.block_start)
    err = 0;
    block_trl = []; 
    ntrial_meg = 0;
    nblock = nblock+1;
   
    
    fprintf('\n\n ---------------- \n Loop through blocks #%d\n ---------------- \n\n', nblock);
 
     
    % Extract blocks from behavioral data corresponding to the current file  
    trials_behavBlock = trials_info(trials_info(:,10)==nblock,:);   % Selecting trials for this block
    
    % Create vectors containing only samples from current meg block   
    MemOn_block = trg_idx.mem_on(find(trg_idx.mem_on <= trg_idx.block_end(nblock) & trg_idx.mem_on >= trg_idx.block_start(nblock)));   
    MemOff_block = trg_idx.mem_off(find(trg_idx.mem_off <= trg_idx.block_end(nblock) & trg_idx.mem_off >= trg_idx.block_start(nblock)));  
    ProbeOn_block =  trg_idx.probe_on(find(trg_idx.probe_on <= trg_idx.block_end(nblock) & trg_idx.probe_on >= trg_idx.block_start(nblock)));  
    ProbeOff_block = trg_idx.probe_off(find(trg_idx.probe_off <= trg_idx.block_end(nblock) & trg_idx.probe_off >= trg_idx.block_start(nblock)));  
    RespRight_block = trg_idx.resp_right(find(trg_idx.resp_right <= trg_idx.block_end(nblock) & trg_idx.resp_right >= trg_idx.block_start(nblock)));  
    RespLeft_block = trg_idx.resp_left(find(trg_idx.resp_left <= trg_idx.block_end(nblock) & trg_idx.resp_left >= trg_idx.block_start(nblock)));
    RespBad_block = trg_idx.resp_bad(find(trg_idx.resp_bad <= trg_idx.block_end(nblock) & trg_idx.resp_bad >= trg_idx.block_start(nblock)));
    FbCorrect_block = trg_idx.fb_correct(find(trg_idx.fb_correct <= trg_idx.block_end(nblock) & trg_idx.fb_correct >= trg_idx.block_start(nblock)));
    FbError_block = trg_idx.fb_error(find(trg_idx.fb_error <= trg_idx.block_end(nblock) & trg_idx.fb_error >= trg_idx.block_start(nblock)));
    FbBad_block = trg_idx.fb_bad(find(trg_idx.fb_bad <= trg_idx.block_end(nblock) & trg_idx.fb_bad >= trg_idx.block_start(nblock)));
    resp_block = [RespRight_block;RespLeft_block;RespBad_block]; % creating vector for any response 
    resp_block = sort(resp_block);
    Fb_block = [FbCorrect_block;FbError_block;FbBad_block]; % creating vector for any response 
    Fb_block = sort(Fb_block);
    %remove triggers when block end was forced and trial has not ended
    %properly 
    if MemOn_block(end) > Fb_block(end) %end of entire trial ist the feedback here!
    resp_block(resp_block>MemOn_block(end))=[]; 
    RespRight_block(RespRight_block>MemOn_block(end))=[];
    RespLeft_block(RespLeft_block>MemOn_block(end))=[];
    RespBad_block(RespBad_block>MemOn_block(end))=[];
    Fb_block(Fb_block>MemOn_block(end))=[]; 
    FbCorrect_block(FbCorrect_block>MemOn_block(end))=[];
    FbError_block(FbError_block>MemOn_block(end))=[];
    FbBad_block(FbBad_block>MemOn_block(end))=[]; 
    ProbeOn_block(ProbeOn_block>MemOn_block(end))=[];
    ProbeOff_block(ProbeOff_block>MemOn_block(end))=[];
    MemOff_block(MemOff_block>MemOn_block(end))=[];
    MemOn_block(end) = [];% Remove last element wherever necessary, i.e. triggers that are part of unfinished trial 
    end

    
    % Setting the sample from MEG data into "trial_info"
    trl_behavIdx = length(trials_behavBlock);
    

        loop_vec = length(MemOn_block):-1:1;
    
   
    fprintf('\n\n ---------------- \n Loop through trials...\n ---------------- \n\n');
    % is there a memorandum onset, but the trial has been terminated, so no
    % corresponding behavioral data is avalable? 
    if length(MemOn_block)>length(allbehav(allbehav(:,10)==nblock))
       MemOn_block(end-1:end)=[]; 
    end 
    % Create trials based on the trialOn trigger
    for i = 1:length(MemOn_block)
        ntrial_meg = ntrial_meg + 1;
        fprintf('\n... Loop through trial %d', ntrial_meg);
        
        
        % Determine where the trial starts with respect to the event
        if isfield(cfgin.trialdef, 'prestim')
            trlOff = round(-cfgin.trialdef.prestim*hdr.Fs);
            trlStart = trgIdx(trgIdx==MemOn_block(i)) + trlOff; % Shift trial start 
        else
            trlOff = -hdr.Fs;
            trlStart = trgIdx(trgIdx==MemOn_block(i));
        end
        
        % Determine where the trial ends (default response cue changend to Probe Onset!)
        k = ProbeOn_block(ntrial_meg);
        if isfield(cfgin.trialdef, 'poststim')
           trlEnd = trgIdx(trgIdx==k) + round(cfgin.trialdef.poststim*hdr.Fs);  % Shift trial end
        else
            trlEnd = trgIdx(trgIdx==k);
        end
                % Determine where the feedback was given
        k = Fb_block(ntrial_meg);
        
           Feedback_time = trgIdx(trgIdx==k); 

      
    
   % GM: Extend info about event times for later 
        event_times = zeros (1,length(fields)); %initialize matrix for specific times 
        event_times(event_times==0)=nan; %set all values to nans 
        
        
for z=1:length(MemOn_block(:,1))
    if MemOn_block(z)>=trlStart && MemOn_block(z)<=Feedback_time
       event_times(1) =  MemOn_block(z);
    end 
end 
for z=1:length(MemOff_block(:,1))
    if MemOff_block(z)>=trlStart && MemOff_block(z)<=Feedback_time
       event_times(2) =  MemOff_block(z);
    end 
end 
for z=1:length(ProbeOn_block(:,1))
    if ProbeOn_block(z)>=trlStart && ProbeOn_block(z)<=Feedback_time
       event_times(3) =  ProbeOn_block(z);
    end 
end 
for z=1:length(ProbeOff_block(:,1))
    if ProbeOff_block(z)>=trlStart && ProbeOff_block(z)<=Feedback_time
       event_times(4) =  ProbeOff_block(z);
    end 
end 
for z=1:length(RespRight_block(:,1))
    if RespRight_block(z)>=trlStart && RespRight_block(z)<=Feedback_time
       event_times(5) =  RespRight_block(z);
    end 
end 
for z=1:length(RespLeft_block(:,1))
    if RespLeft_block(z)>=trlStart && RespLeft_block(z)<=Feedback_time
       event_times(6) =  RespLeft_block(z);
    end 
end 
if ~isempty(RespBad_block)
for z=1:length(RespBad_block(:,1))
    if RespBad_block(z)>=trlStart && RespBad_block(z)<=Feedback_time
       event_times(7) =  RespBad_block(z);
    end 
end 
end
if ~isempty(FbCorrect_block)
for z=1:length(FbCorrect_block(:,1))
    if FbCorrect_block(z)>=trlStart && FbCorrect_block(z)<=Feedback_time
       event_times(8) =  FbCorrect_block(z);
    end 
end 
end
if ~isempty(FbError_block)
for z=1:length(FbError_block(:,1))
    if FbError_block(z)>=trlStart && FbError_block(z)<=Feedback_time
       event_times(9) =  FbError_block(z);
    end 
end 
end 
if ~isempty(FbBad_block)
for z=1:length(FbBad_block(:,1))
    if FbBad_block(z)>=trlStart && FbBad_block(z)<=Feedback_time
       event_times(10) =  FbBad_block(z);
    end 
end 
end 
 
        % Concatenate information 
        new_trial = [trlStart trlEnd trlOff session trials_behavBlock(ntrial_meg,:) event_times]; 
        block_trl = [block_trl; new_trial];    
        
    end % End of meg trial loop
    
    trl = [trl; block_trl];  
       
end % End of meg block loop


% Construct matrix containing start and end points of each block (+-10 seconds)
offset_mat = -round(10*hdr.Fs); 

% Shift trial start
startSamples_mat = startBlockSamp - round(10*hdr.Fs); % Shift trial start
endSamples_mat = endBlockSamp + round(10*hdr.Fs); % Shift trial end


cfgin.trl = trl;

cfgin.trialInfoLabel = {'startSample','endSample','offset','session','location_memorandum', 'location_probe',...
    'delay','ttype','response','correct_false','RT','idk','eyetracking','block','trialex','trialnumber'...
    'MemOn', 'MemOff', 'ProbeOn', 'ProbeOff', 'RespRight', 'RespLeft', 'RespBad', 'FbCorrect', 'FbError', 'FbBad'};
cfgin.blockBound_trl = [startSamples_mat endSamples_mat offset_mat*ones(length(startSamples_mat),1)];

cfgin.event = events;

end % End of entire function


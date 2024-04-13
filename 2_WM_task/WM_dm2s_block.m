% TASK: remember spatial location of a checkerboard presented somewhere
% along isoeccentric line in lower half of visual field, and say whether it
% matches location of a probe stimulus after 1 of 3 different delays.


sca, clear
commandwindow;  % stops task response keys being written into task script

testing = 0;  % set to 1 for debugging version: brief stimuli, few trials
language = 'german';  % switch between 'german' and 'english'


%% Path details
cd /home/userm/MCI/Task_SpatialOnly_MEG/
datadir = '/home/userm/MCI/Data/';


%% Basic PTB setup
options = setup;
if strcmp(options.do_trigger,'yes'), addpath(['matlabtrigger',filesep]), else addpath(['faketrigger',filesep]), end
trigger_enc = setup_trigger;
Screen('Preference', 'SkipSyncTests', 0 );


%% Subject details
if ~testing
    subj = input('Subject initials? ', 's');
    sess = input('Session number? ', 's');
else
    subj = 'test';
    sess = '1';
    options.do_trigger = 'no';
end

behavpath = [datadir,subj,filesep,'S',sess,filesep,'Behaviour',filesep];
ELpath = [datadir,subj,filesep,'S',sess,filesep,'Eyetracking',filesep];

if ~isdir([datadir,subj,filesep,'S',sess])  % making new directories for this session if they don't exist, and setting block number to 1
      mkdir(behavpath)
      mkdir(ELpath)
      b = '1';  fprintf('\nNew subject/session combination...\n')
      WaitSecs(2);
else
    b = count_blocks([behavpath,'*.mat']);  % otherwise, counting number of existing files to determine current block number
    if isempty(b), b='1'; end   % catch just in case there are no saved behavioural files (if script errored before b1 end and now rerunning)
end
fprintf('\nCurrent block number is %s...\n',b)
WaitSecs(2);

% Suppress keyboard input to command line
ListenChar(2)


%% Fixed task settings
% Trial stuff
delays     = [1 3 9];    % delay priod durations
tpd        =      21;    % number of trials per delay
p_match    =    0.33;    % proportion of trials on which sample matches probe
p_near     =    0.33;    % proportion of trials on which probe is presented at nearest location  probe

% Timing (all in seconds)
timing.mem       =     0.5;    % memoranda presentation time
timing.probe     =     0.5;    % probe presentation time
timing.prefb     =     0.1;    % time between response and associated feedback
timing.fb        =    0.75;    % duration of feedback (for each of two consecutive tones)
timing.ITI       =     4.0;    % minimum time between feedback and next-trial onset
timing.rest      =     3.0;    % amount of fixed ITI time where blinks/rest is allowed (inclusive of feedback time)
timing.ITIdist   =     1.5;    % maximum increment to time between feedback and next-trial onset (determines upper bound on uniform distribution)

timing.s_flicker =       0;    % flicker the checkerboard? [1=yes, 0=no]
timing.s_refresh =     0.1;    % time between switch in checkerboard polarity within each sample (i.e. checker flicker rate, if desired)

% Fixation/cue appearance
stim.s_r             =    round(6.0*options.ppd);  % radial offset of stimuli (in d.v.a converted to pixels)
stim.fix_t           =    round(0.2*options.ppd);  % thickness of fixation cross arms
stim.fix_l           =    round(0.8*options.ppd);  % length of fixation cross arms
stim.fix_c_active    =             [255 255 255];  % colour of fixation cross for active trial period (WHITE)
stim.fix_c_go        =             [128 153 199];  % colour of fixation cross for response cueing (BLUEISH)
stim.txt_yoffset     =    round(1.0*options.ppd);  % y-offset of feedback/rest period text

% Checkerboard appearance - locations are pre-set to be @ 4 equidistant points along isoeccentric line form fixation, left & right, +/-45deg from midline
nlocs = 12;
angles = [0:(180/(nlocs+1)):180];  % possible stimulus polar angles (adding 2 so that most extreme memoranda can still be flanked either side by a near non-match)
anglesU = angles(2:end-1);

% Checkerboard patch settings
cb.size     =     2.8;   % size of one side of square within which checkerboard will be drawn (d.v.a.)
cb.freq     =     1.0;   % spatial frequency (cycles per d.v.a.)

% If only testing, overwrite some of the above
if testing == 1
    tpd = 3;
    timing.mem = 0.75;
    timing.probe = 0.75;
    timing.rest = 1.5;
    timing.ITI = 2.0;
    timing.ITIdist = 0.3;
end


%% Initialize Psychtoolbox and create textures
setup_ptb;
timing.ifi = options.frameDur;     % inter-frame interval (1/sampling rate)

[gabor.tex,gabor.rect] = createCircularChecker(window, cb.size, cb.freq, options.ppd, 0);
gabor.rect = gabor.rect-(gabor.rect(4)/2);    % making sure subsequent coordinates are always wrt screen center


try
    %% Present intro screen
    block_instructions(b,window,windowRect,language);
    
    
    %% Generate gabor locations & orientations for all trials
    % Seed random number generator
    rng('default')
    rng('shuffle')
    
    % Construct matrix of MATCH sample stimuli (ensuring at least one of each sample stimulus per delay if possible)
    matchS = nan(round(tpd*p_match),length(delays)); matchD = [];
    for d = 1:length(delays)
        if round(tpd*p_match) == length(anglesU)
            matchS(:,d) = anglesU';
        elseif round(tpd*p_match) > length(anglesU)
            tempmatchS = repmat(anglesU',floor(round(tpd*p_match)/length(anglesU)),1);
            tempmatchS = [tempmatchS; randsample(anglesU,round(mod(round(tpd*p_match),length(anglesU))))'];
            matchS(:,d) = tempmatchS;
        else matchS(:,d) = randsample(anglesU,round(tpd*p_match))';
        end
        matchD = [matchD; ones(size(matchS,1),1).*delays(d)];
    end
    matchS = reshape(matchS,length(matchD),1);
    
    % Construct matrix of NEAR-NO-MATCH sample stimuli (ensuring at least one of each sample stimulus per delay if possible)
    nearS = nan(round(tpd*p_near),length(delays)); nearD = [];
    for d = 1:length(delays)
        if round(tpd*p_near) == length(anglesU)
            nearS(:,d) = anglesU';
        elseif round(tpd*p_near) > length(anglesU)
            tempnearS = repmat(anglesU',floor(round(tpd*p_near)/length(anglesU)),1);
            tempnearS = [tempnearS; randsample(anglesU,round(mod(round(tpd*p_near),length(anglesU))))'];
            nearS(:,d) = tempnearS;
        else nearS(:,d) = randsample(anglesU,round(tpd*p_near))';
        end
        nearD = [nearD; ones(size(nearS,1),1).*delays(d)];
    end
    nearS = reshape(nearS,length(nearD),1);
    
    % Construct matrix of FAR-NO-MATCH sample stimuli (ensuring at least one of each sample stimulus per delay if possible)
    p_far = 1-p_match-p_near;
    farS = nan(round(tpd*p_far),length(delays)); farD = [];
    for d = 1:length(delays)
        if round(tpd*p_far) == length(anglesU)
            farS(:,d) = anglesU';
        elseif round(tpd*p_far) > length(anglesU)
            tempfarS = repmat(anglesU',floor(round(tpd*p_far)/length(anglesU)),1);
            tempfarS = [tempfarS; randsample(anglesU,round(mod(round(tpd*p_far),length(anglesU))))'];
            farS(:,d) = tempfarS;
        else farS(:,d) = randsample(anglesU,round(tpd*p_far))';
        end
        farD = [farD; ones(size(farS,1),1).*delays(d)];
    end
    farS = reshape(farS,length(farD),1);
    
    % Construct matrices of PROBE stimuli
    matchP = matchS;
    
    nearP = nearS;
    for i = 1:size(nearP,1)
        c_angles = [find(angles==nearS(i))-1 find(angles==nearS(i))+1];
        nearP(i) = angles(randsample(c_angles,1));
    end
    
    farP = farS;
    for i = 1:size(farP,1)
        no_angle = find(angles==farS(i))-1:find(angles==farS(i))+1;
        c_angles = find(~ismember(1:length(angles),no_angle));
        farP(i) = angles(randsample(c_angles,1));
    end
    
    % Construct full matrix of stimulus types [SAMPLE ANGLE; PROBE ANGLE; DELAY; TRIAL-TYPE] (trial-types: 1= match, 2=near no-match, 3=far no-match)
    stimIn = [[matchS matchP matchD ones(size(matchP,1),1)]; [nearS nearP nearD ones(size(nearP,1),1).*2]; [farS farP nearD ones(size(farP,1),1).*3]];
    
    % Randomize order
    shuforder = randperm(size(stimIn,1));
    stimIn = stimIn(shuforder,:);
    
    
    %% Transform dot positions from polar to cartesian (x,y) coordinates
    [X,Y] = RectCenter(windowRect);
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn(:,1)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    [XstimP,YstimP] = pol2cart(deg2rad(stimIn(:,2)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimP = XstimP+ X; YstimP = YstimP+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS XstimP YstimP];
        
    
    %% Start Eyelink recording
    if strcmp(options.et,'yes');
        Eyelink('StartRecording');
        WaitSecs(0.1);
        Eyelink('message', 'Start recording Eyelink');
    end
    
    
    %% Countdown and first fixation point
    if strcmp(language,'german')
        str1 = 'Bereit halten!';
    elseif strcmp(language,'english')
        str1 = 'Ready to go!';
    end
    cd_sec = 3;
    while cd_sec>0
        str2 = num2str(cd_sec);
        DrawFormattedText(window,str1,'center',Y-25,[255 255 255]);
        DrawFormattedText(window,str2,'center',Y+25,[255 255 255]);
        Screen('Flip', window);
        cd_sec=cd_sec-1;
        WaitSecs(0.96);
    end
    
    FixCross = [X-stim.fix_t,Y-stim.fix_l,X+stim.fix_t,Y+stim.fix_l; X-stim.fix_l,Y-stim.fix_t,X+stim.fix_l,Y+stim.fix_t];
    
    f_ps = [X-stim.fix_t/2, Y-stim.fix_l/2, X+stim.fix_t/2, Y+stim.fix_l/2;...  % middle vertical bar of fixation cross
        X-stim.fix_l/2, Y-stim.fix_t/2, X,              Y+stim.fix_t/2;...  % left horizontal bar of fixation cross
        X             , Y-stim.fix_t/2, X+stim.fix_l/2, Y+stim.fix_t/2]';   % right horizontal bar of fixation cross
    f_cs = [stim.fix_c_active' stim.fix_c_active' stim.fix_c_active'];  % initializing colours of all fixation cross parts (changes often)
    f_ps_in = [X-stim.fix_t/2; Y-stim.fix_t/2; X+stim.fix_t/2; Y+stim.fix_t/2];  % position of inner fixation circle
    f_cs_in = [127; 127; 127];  % colour of inner/outer fixation circles (doesn't change)
    
    Screen('FillRect', window, f_cs, f_ps);   % drawing fixation cross
    Screen('FillOval', window, f_cs_in, f_ps_in);   % drawing overlaid inner fixation point
    vbl = Screen('Flip', window);   % present fixation
    
    
    %% Loop through trials
    trigger(trigger_enc.block_start);  % trigger to mark start of block
    if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_start)); end
    
    on_time = [];  % onset time for first trial is specified within trial script
    Behav = zeros(size(stimIn,1),size(stimIn,2)+5);  % Behav = [stimIn(:,1:4) resp ACC RT start_time sacc]
    tRefresh = zeros(size(stimIn,1),8);  % tRefresh = [start_time t_mem_off t_probe_on t_probe_off t_resp t_fb_on t_rest_on t_rest_off]
    for t = 1:size(stimIn,1);
        
        % Presenting trial number at the bottom of the eyetracker display
        if strcmp(options.et,'yes');
            Eyelink('command', 'record_status_message "TRIAL %d/%d"', t, size(stimIn,1));
            Eyelink('message', 'TRIALID %d', t);
        end
        
        % Variable task parameters
        varopts.on_time =  on_time;                          % controls onset time of impending trial - fed back from trial function
        varopts.stimM   =  [stimIn(t,1) stimInCoord(t,1:2)]; % angle and X/Y coordinates of memorandum
        varopts.stimP   =  [stimIn(t,2) stimInCoord(t,3:4)]; % angle and X/Y coordinates of memorandum
        varopts.delay   =  stimIn(t,3);                      % delay duration
        varopts.kbqdev  =  options.kbqdev;                   % keyboard info
        varopts.ppd     =  options.ppd;                      % pixels per degree
        
        % Run trial
        [Behav(t,[1:3 5:end]),tRefresh(t,:),on_time] = WM_dm2s_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, t, strcmp(options.et,'yes'), language);
        Behav(t,4) = stimIn(t,4);
    end
    
    WaitSecs(3);
    trigger(trigger_enc.block_end);  % trigger to mark end of block
    if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_end)); end
    
    
    %% Save behavioral data
    fprintf('\nSaving behavioural data to %s\n', behavpath)
    save([behavpath,subj,'_',sess,'_',num2str(b),'.mat'],'Behav','stimInCoord','tRefresh')
    
    
    %% Calculate accuracy from current + previous blocks and display to subject
    end_block_screen(str2double(b),[behavpath,subj,'_',sess,'_'],window,windowRect,language);
    
    
    %% Save Eyelink data
    if strcmp(options.et,'yes');
        fprintf('Saving EyeLink data to %s\n', ELpath)
        eyefilename = fullfile(ELpath,options.edfFile);
        Eyelink('CloseFile');
        Eyelink('WaitForModeReady', 500);
        try
            status = Eyelink('ReceiveFile', options.edfFile, eyefilename);
            disp(['File ' eyefilename ' saved to disk']);
        catch
            warning(['File ' eyefilename ' not saved to disk']);
        end
        Eyelink('StopRecording');
    end
    
    
    %% Exit
    ListenChar(0)
    sca
    
    
catch ME  % if any errors are encountered or user requests quit, clean up and exit
    
    
    WaitSecs(2);
    trigger(trigger_enc.block_end);  % trigger to mark end of block
    if strcmp(options.et,'yes'); Eyelink('message', num2str(trigger_enc.block_end)); end
    
    if exist('Behav','var'),  % saving behaviour if some trials have been run
        if ~isempty(find(Behav(:,1)>0))
            Behav=Behav(Behav(:,1)>0,:); tRefresh=tRefresh(1:size(Behav,1),:); stimInCoord=stimInCoord(1:size(Behav,1),:);
            save([behavpath,subj,'_',sess,'_',num2str(b),'.mat'],'Behav','stimInCoord','tRefresh')
        end
    end
    if strcmp(options.et,'yes');
        Eyelink('CloseFile');
        if exist('Behav','var'),  % saving EL file if some trials have been run
            if ~isempty(find(Behav(:,1)>0))
                fprintf('Saving EyeLink data to %s\n', ELpath)
                eyefilename = fullfile(ELpath,options.edfFile);
                Eyelink('WaitForModeReady', 500);
                try
                    status = Eyelink('ReceiveFile', options.edfFile, eyefilename);
                    disp(['File ' eyefilename ' saved to disk']);
                catch
                    warning(['File ' eyefilename ' not saved to disk']);
                end
            end
        end
        Eyelink('StopRecording');
    end
    ListenChar(0)
    sca
    rethrow(ME)    
    
end



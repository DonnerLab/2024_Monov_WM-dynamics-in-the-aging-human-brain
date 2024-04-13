% TASK: remember spatial location of a checkerboard presented somewhere
% along isoeccentric line in lower half of visual field, and say whether it
% matches location of a probe stimulus after 1 of 3 different delays.

4
sca, clear
commandwindow;  % stops task response keys being written into task script

language = 'german';  % switch between 'german' and 'english'
gaze_fb  = 1;         % 1 if negative feedback should be given for fixation breaks, 0 if not


%% Path details
cd /home/userm/MCI/Task_SpatialOnly_MEG/


%% Basic PTB setup
options = setup;
if strcmp(options.do_trigger,'yes'), addpath(['matlabtrigger',filesep]), else addpath(['faketrigger',filesep]), end
trigger_enc = setup_trigger;
Screen('Preference', 'SkipSyncTests', 0 )

% Force eye-tracking on ('yes') or off ('no')
options.et = 'yes';

% Suppress keyboard input to command line
ListenChar(2)

% Useless details to pass to EL
if strcmp(options.et,'yes')
    subj = 'trn';  % giving stock subj/sess/block ID for training - EL file won't be saved anyway
    sess = '1';
    b = '1';
end

%% Fixed task settings
% Timing (all in seconds)
timing.mem       =     1.2;    % memoranda presentation time
timing.probe     =     1.2;    % probe presentation time
timing.prefb     =     0.1;    % time between response and associated feedback
timing.fb        =    0.75;    % duration of auditory feedback (for each of two consecutive tones)
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

% Checkerboard appearance
nlocs = 12;
angles = [0:(180/(nlocs+1)):180];  % possible stimulus polar angles (adding 2 so that most extreme memoranda can still be flanked either side by a near non-match)
anglesU = angles(2:end-1);

% Checkerboard patch settings
cb.size     =     2.8;   % size of one side of square within which checkerboard will be drawn (d.v.a.)
cb.freq     =     1.0;   % spatial frequency (cycles per d.v.a.)


%% Initialize Psychtoolbox and create textures
setup_ptb;
timing.ifi = options.frameDur;     % inter-frame interval (1/sampling rate)

[gabor.tex,gabor.rect] = createCircularChecker(window, cb.size, cb.freq, options.ppd, 0);
gabor.rect = gabor.rect-(gabor.rect(4)/2);    % making sure subsequent coordinates are always wrt screen center


try
    %% Start Eyelink recording
    if strcmp(options.et,'yes');
        Eyelink('StartRecording');
        WaitSecs(0.1);
        Eyelink('message', 'Start recording Eyelink');
    end
    
    %% Run short FIXATION TRAINING- fixate while checkers are presented
    [X,Y] = RectCenter(windowRect);
    
    if strcmp(language,'german')
        str1 = sprintf('Willkommen! Lassen Sie uns zuerst ein wenig Übung mit der Fixierung');
        str2 = sprintf('des Kreuzes in der Mitte des Bildschirmes kriegen.');

        str3 = sprintf('Wir werden Ihnen einige kurze Objekte zeigen, so wie die,');
        str4 = sprintf('die Sie auch während der Aufgabe zu sehen kriegen.');
        
        str5 = sprintf('Alles was Sie jetzt machen müssen ist es,');
        str6 = sprintf('dass Sie das Kreuz in der Mitte des Bildschirms direkt anschauen,');
        str7 = sprintf('und nur blinzeln wenn kein Objekt angezeigt wird.');
        
        str8 = sprintf('Schauen Sie NICHT auf die Objekte, wenn diese erscheinen.');
        str9 = 'Drücken Sie einen beliebigen Knopf zum starten...';
    elseif strcmp(language,'english')
        str1 = sprintf('Welcome! Let''s first get some practice with keeping fixation on the central cross.');
        str2 = sprintf('We will show several brief objects, like the ones you will see during the real task.');
        str3 = sprintf('All you need to do here, however, is keep looking at the cross in the middle of the screen,');
        str4 = sprintf('and only blink in the interval between objects.');
        str5 = sprintf('Do NOT look at the objects when they are presented.');
        str6 = 'Press any button to start...';
        str7 = ''; str8 = ''; str9 = '';
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-300,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y-260,[255 255 255]);
    
    DrawFormattedText(window,str3,'center',Y-160,[255 255 255]);
    DrawFormattedText(window,str4,'center',Y-120,[255 255 255]);
    
    DrawFormattedText(window,str5,'center',Y-20,[255 255 255]);
    DrawFormattedText(window,str6,'center',Y+20,[255 255 255]);
    DrawFormattedText(window,str7,'center',Y+60,[255 255 255]);
    
    DrawFormattedText(window,str8,'center',Y+140,[255 255 255]);
    DrawFormattedText(window,str9,'center',Y+220,[255 255 255]);
    Screen('Flip', window);
    
    % Wait for response
    WaitSecs(0.5);
    keyIsDown = 0;
    while ~keyIsDown
        [keyIsDown, ~, ~] = KbCheck;  % checking for response
    end
    
    % Draw random selection of stimulus locations
    stimIn = angles(randperm(length(angles),6));
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn'),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS];
    
    for t = 1:length(stimIn)
        % Variable task parameters
        varopts.stim   =  [stimIn(t) stimInCoord(t,1:2)]; % angle and X/Y coordinates of memorandum
        varopts.kbqdev  =  options.kbqdev;                 % keyboard info
        varopts.ppd     =  options.ppd;                      % pixels per degree
        
        % Run trial
        [tBehav] = WM_fixation_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, strcmp(options.et,'yes'));
        
        % Present message if fixation was broken
        if tBehav(end) > 0 && gaze_fb
            % Specify text
            if strcmp(language,'german')
                str1 = 'Sie haben vom Fixationskreuz weggeschaut!';
                str2 = 'Denken Sie daran, es ist sehr wichtig direkt das Kreuz anzuschauen,';
                str3 = 'und NICHT zu den Objekten zu schauen, wenn diese erscheinen.';
                str4 = 'Drücken Sie einen beliebigen Knopf zum fortfahren...';
            elseif strcmp(language,'english')
                str1 = 'You looked away from the fixation cross!';
                str2 = 'Remember, it''s very important to look directly at the cross,';
                str3 = 'and DO NOT look at the stimuli when they appear on screen.';
                str4 = 'Press any button to continue...';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-200,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-60,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y,[255 255 255]);
            DrawFormattedText(window,str4,'center',Y+200,[255 255 255]);
            Screen('Flip', window);
            
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
        end
    end
        
    %% Present FIRST example trial -- FAR NON-MATCH
    % Specify text
    if strcmp(language,'german')
        str1 = sprintf('OK. Lassen Sie uns nun ein paar Übungen der Aufgabe machen.');
        str2 = sprintf('Alle Objekte werden länger auf dem Bildschirm sein,');
        str3 = sprintf('sodass Sie sich an die Aufgabe gewöhnen können.');
        str4 = 'Drücken Sie einen beliebigen Knopf um mit dem ersten Beispiel zu beginnen...';
    elseif strcmp(language,'english')
        str1 = sprintf('OK. Let''s do a few practice trials of the task.');
        str2 = sprintf('Here, the objects will be presented on screen for a prolonged');
        str3 = sprintf('time to allow you to learn the sequence of events.');
        str4 = 'Press any button to see the first example trial...';
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-200,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y-60,[255 255 255]);
    DrawFormattedText(window,str3,'center',Y,[255 255 255]);
    DrawFormattedText(window,str4,'center',Y+200,[255 255 255]);
    Screen('Flip', window);
    
    % Wait for response
    WaitSecs(0.5);
    keyIsDown = 0;
    while ~keyIsDown
        [keyIsDown, ~, ~] = KbCheck;  % checking for response
    end
        
    % Specify stimulus details for this trial
    stimIn = [angles(end-3) angles(2) 1 3];
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn(:,1)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    [XstimP,YstimP] = pol2cart(deg2rad(stimIn(:,2)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimP = XstimP+ X; YstimP = YstimP+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS XstimP YstimP];
    
    % Variable task parameters
    varopts.stimM   =  [stimIn(1) stimInCoord(1:2)]; % angle and X/Y coordinates of memorandum
    varopts.stimP   =  [stimIn(2) stimInCoord(3:4)]; % angle and X/Y coordinates of memorandum
    varopts.delay   =  stimIn(3);                      % delay duration
    varopts.kbqdev  =  options.kbqdev;                 % keyboard info
    varopts.ppd     =  options.ppd;                      % pixels per degree
    
    % Run trial
    cresp = 0;
    while ~cresp
        [tBehav] = WM_dm2s_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, strcmp(options.et,'yes'));
        if tBehav(5) ~= 1
            % Specify text
            if strcmp(language,'german')
                str1 = sprintf('Falsch!');
                str2 = sprintf('Das erste und das zweite Objekt waren an VERSCHIEDENEN Orten.');
                str3 = 'Drücken Sie einen beliebigen Knopf um das Beispiel zu wiederholen...';
            elseif strcmp(language,'english')
                str1 = sprintf('Incorrect!');
                str2 = sprintf('The first and second objects were in DIFFERENT locations.');
                str3 = 'Press any button to try that trial again...';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+200,[255 255 255]);
            Screen('Flip', window);
            
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
            
        elseif tBehav(end) > 0 && gaze_fb
            if strcmp(language,'german')
                str1 = 'Sie haben vom Fixationskreuz weggeschaut!';
                str2 = 'Denken Sie daran, es ist sehr wichtig,';
                str3 = 'dass Sie während den Versuchen direkt das Kreuz anschauen,';
                str4 = 'und NICHT auf die Objekte schauen, wenn diese erscheinen.';
                str5 = 'Drücken Sie einen beliebigen Knopf um diesen Versuch zu wiederholen,';
                str6 = 'schauen Sie diesmal die ganze Zeit auf das Kreuz.';
            elseif strcmp(language,'english')
                str1 = 'You looked away from the fixation cross!';
                str2 = 'Remember, it''s very important to look directly at the cross during each trial,';
                str3 = 'and DO NOT look at the stimuli when they appear on screen.';
                str4 = 'Press any button to do that trial again, and this time look at the cross throughout...';
                str5 = '';
                str6 = '';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+0,[255 255 255]);
            DrawFormattedText(window,str4,'center',Y+40,[255 255 255]);
            DrawFormattedText(window,str5,'center',Y+200,[255 255 255]);
            DrawFormattedText(window,str6,'center',Y+240,[255 255 255]);
            Screen('Flip', window);
            
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
            
        else cresp = 1;
            break
        end
    end
    
    
    %% Present SECOND example trial -- MATCH
    % Specify text
    if strcmp(language,'german')
        str1 = sprintf('Gut gemacht! Die Objekte waren an VERSCHIEDENEN Orten.');
        str2 = sprintf('Lassen Sie uns ein weiteres Beispiel versuchen.');
        str3 = 'Drücken Sie einen beliebigen Knopf zum starten...';
    elseif strcmp(language,'english')
        str1 = sprintf('Well done! The objects were in DIFFERENT locations.');
        str2 = sprintf('Let''s try another example.');
        str3 = 'Press any button to see it...';
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
    DrawFormattedText(window,str3,'center',Y+200,[255 255 255]);
    Screen('Flip', window);
    
    % Wait for response
    WaitSecs(0.5);
    keyIsDown = 0;
    while ~keyIsDown
        [keyIsDown, ~, ~] = KbCheck;  % checking for response
    end
    
    % Specify stimulus details for this trial
    stimIn = [angles(4) angles(4) 1.5 1];
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn(:,1)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    [XstimP,YstimP] = pol2cart(deg2rad(stimIn(:,2)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimP = XstimP+ X; YstimP = YstimP+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS XstimP YstimP];
    
    % Variable task parameters
    varopts.stimM   =  [stimIn(1) stimInCoord(1:2)]; % angle and X/Y coordinates of memorandum
    varopts.stimP   =  [stimIn(2) stimInCoord(3:4)]; % angle and X/Y coordinates of memorandum
    varopts.delay   =  stimIn(3);                      % delay duration
    varopts.kbqdev  =  options.kbqdev;                 % keyboard info
    varopts.ppd     =  options.ppd;                      % pixels per degree
    
    % Run trial
    cresp = 0;
    while ~cresp
        [tBehav] = WM_dm2s_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, strcmp(options.et,'yes'));
        if tBehav(5) ~= 1
            % Specify text
            if strcmp(language,'german')
                str1 = sprintf('Falsch!');
                str2 = sprintf('Das erste und das zweite Objekt waren am SELBEN Ort.');
                str3 = 'Drücken Sie einen beliebigen Knopf um das Beispiel zu wiederholen...';
            elseif strcmp(language,'english')
                str1 = sprintf('Incorrect!');
                str2 = sprintf('The first and second objects were in the SAME location.');
                str3 = 'Press any button to try that trial again...';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+200,[255 255 255]);
            Screen('Flip', window);
            
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
            
        elseif tBehav(end) > 0 && gaze_fb
            if strcmp(language,'german')
                str1 = 'Sie haben vom Fixationskreuz weggeschaut!';
                str2 = 'Denken Sie daran, es ist sehr wichtig,';
                str3 = 'dass Sie während den Versuchen direkt das Kreuz anschauen,';
                str4 = 'und NICHT auf die Objekte schauen, wenn diese erscheinen.';
                str5 = 'Drücken Sie einen beliebigen Knopf um diesen Versuch zu wiederholen,';
                str6 = 'schauen Sie diesmal die ganze Zeit auf das Kreuz.';
            elseif strcmp(language,'english')
                str1 = 'You looked away from the fixation cross!';
                str2 = 'Remember, it''s very important to look directly at the cross during each trial,';
                str3 = 'and DO NOT look at the stimuli when they appear on screen.';
                str4 = 'Press any button to do that trial again, and this time look at the cross throughout...';
                str5 = '';
                str6 = '';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+0,[255 255 255]);
            DrawFormattedText(window,str4,'center',Y+40,[255 255 255]);
            DrawFormattedText(window,str5,'center',Y+200,[255 255 255]);
            DrawFormattedText(window,str6,'center',Y+240,[255 255 255]);
            Screen('Flip', window);
                        
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
            
        else cresp = 1;
            break
        end
    end
    
    
    %% Present THIRD example trial -- NEAR NON-MATCH
    % Specify text
    if strcmp(language,'german')
        str1 = sprintf('Gut gemacht! Die Objekte waren am SELBEN Ort.');
        str2 = sprintf('Lassen Sie uns ein drittes Beispiel versuchen.');
        str3 = 'Drücken Sie einen beliebigen Knopf zum starten...';
    elseif strcmp(language,'english')
        str1 = sprintf('Well done! The objects were in the SAME location.');
        str2 = sprintf('Let''s try a third example.');
        str3 = 'Press any button to see it...';
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
    DrawFormattedText(window,str3,'center',Y+200,[255 255 255]);
    Screen('Flip', window);
    
    % Wait for response
    WaitSecs(0.5);
    keyIsDown = 0;
    while ~keyIsDown
        [keyIsDown, ~, ~] = KbCheck;  % checking for response
    end
    
    % Specify stimulus details for this trial
    stimIn = [angles(end-3) angles(end-4) 1.5 2];
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn(:,1)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    [XstimP,YstimP] = pol2cart(deg2rad(stimIn(:,2)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimP = XstimP+ X; YstimP = YstimP+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS XstimP YstimP];
    
    % Variable task parameters
    varopts.stimM   =  [stimIn(1) stimInCoord(1:2)]; % angle and X/Y coordinates of memorandum
    varopts.stimP   =  [stimIn(2) stimInCoord(3:4)]; % angle and X/Y coordinates of memorandum
    varopts.delay   =  stimIn(3);                      % delay duration
    varopts.kbqdev  =  options.kbqdev;                 % keyboard info
    varopts.ppd     =  options.ppd;                      % pixels per degree
    
    % Run trial
    cresp = 0;
    while ~cresp
        [tBehav] = WM_dm2s_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, strcmp(options.et,'yes'));
        if tBehav(5) ~= 1
            % Specify text
            if strcmp(language,'german')
                str1 = sprintf('Falsch!');
                str2 = sprintf('Obwohl das erste und das zweite Objekt nahe beieinander lagen,');
                str3 = sprintf('waren sie an VERSCHIEDENEN Orten.');
                str4 = 'Drücken Sie einen beliebigen Knopf um das Beispiel zu wiederholen...';
            elseif strcmp(language,'english')
                str1 = sprintf('Incorrect!');
                str2 = sprintf('Although the first and second objects were very close to each other,');
                str3 = sprintf('they were actually in DIFFERENT locations.');
                str4 = 'Press any button to try that trial again...';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+20,[255 255 255]);
            DrawFormattedText(window,str4,'center',Y+200,[255 255 255]);
            Screen('Flip', window);
            
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
            
        elseif tBehav(end) > 0 && gaze_fb
            if strcmp(language,'german')
                str1 = 'Sie haben vom Fixationskreuz weggeschaut!';
                str2 = 'Denken Sie daran, es ist sehr wichtig,';
                str3 = 'dass Sie während den Versuchen direkt das Kreuz anschauen,';
                str4 = 'und NICHT auf die Objekte schauen, wenn diese erscheinen.';
                str5 = 'Drücken Sie einen beliebigen Knopf um diesen Versuch zu wiederholen,';
                str6 = 'schauen Sie diesmal die ganze Zeit auf das Kreuz.';
            elseif strcmp(language,'english')
                str1 = 'You looked away from the fixation cross!';
                str2 = 'Remember, it''s very important to look directly at the cross during each trial,';
                str3 = 'and DO NOT look at the stimuli when they appear on screen.';
                str4 = 'Press any button to do that trial again, and this time look at the cross throughout...';
                str5 = '';
                str6 = '';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+0,[255 255 255]);
            DrawFormattedText(window,str4,'center',Y+40,[255 255 255]);
            DrawFormattedText(window,str5,'center',Y+200,[255 255 255]);
            DrawFormattedText(window,str6,'center',Y+240,[255 255 255]);
            Screen('Flip', window);
                        
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
            
        else cresp = 1;
            break
        end
    end
    
    
    %% Present FOURTH example trial -- MATCH w/ LONG DELAY
    % Specify text
    if strcmp(language,'german')
        str1 = sprintf('Sehr gut! Die Objekte waren an VERSCHIEDENEN Orten.');
        str2 = sprintf('Lassen Sie uns ein letztes Beispiel machen. Diesmal ist die Zeit zwischen');
        str3 = 'dem ersten und dem zweiten Objekt länger - das kann auch in der Aufgabe passieren.';
        str4 = 'Es ist besonders schwer die ganze Zeit das Kreuz anzuschauen wenn die Pause';
        str5 = 'so lang ist, aber versuchen Sie bitte Ihr bestes.';
        str6 = 'Drücken Sie einen beliebigen Knopf zum starten...';
    elseif strcmp(language,'english')
        str1 = sprintf('Very good! The objects were in DIFFERENT locations.');
        str2 = sprintf('Let''s do a final example. For this one, the delay between');
        str3 = 'the first and second objects will be longer, which can happen in the task.';
        str4 = 'It''s particularly difficult to keep looking at the cross for the entirety of';
        str5 = 'these long delays, but please try your best to do so.';
        str6 = 'Press any button to see it...';
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-200,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y-100,[255 255 255]);
    DrawFormattedText(window,str3,'center',Y-60,[255 255 255]);
    DrawFormattedText(window,str4,'center',Y+20,[255 255 255]);
    DrawFormattedText(window,str5,'center',Y+60,[255 255 255]);
    DrawFormattedText(window,str6,'center',Y+200,[255 255 255]);
    Screen('Flip', window);
    
    % Wait for response
    WaitSecs(0.5);
    keyIsDown = 0;
    while ~keyIsDown
        [keyIsDown, ~, ~] = KbCheck;  % checking for response
    end
    
    % Specify stimulus details for this trial
    stimIn = [angles(2) angles(2) 7.5 1];
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn(:,1)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    [XstimP,YstimP] = pol2cart(deg2rad(stimIn(:,2)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimP = XstimP+ X; YstimP = YstimP+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS XstimP YstimP];
    
    % Variable task parameters
    varopts.stimM   =  [stimIn(1) stimInCoord(1:2)]; % angle and X/Y coordinates of memorandum
    varopts.stimP   =  [stimIn(2) stimInCoord(3:4)]; % angle and X/Y coordinates of memorandum
    varopts.delay   =  stimIn(3);                      % delay duration
    varopts.kbqdev  =  options.kbqdev;                 % keyboard info
    varopts.ppd     =  options.ppd;                      % pixels per degree
    
    % Run trial
    cresp = 0;
    while ~cresp
        [tBehav] = WM_dm2s_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, strcmp(options.et,'yes'));
        if tBehav(5) ~= 1
            % Specify text
            if strcmp(language,'german')
                str1 = sprintf('Falsch!');
                str2 = sprintf('Das erste und das zweite Objekt waren am SELBEN Ort.');
                str3 = 'Drücken Sie einen beliebigen Knopf um das Beispiel zu wiederholen...';
            elseif strcmp(language,'english')
                str1 = sprintf('Incorrect!');
                str2 = sprintf('The first and second locations were in fact the SAME.');
                str3 = 'Press any button to try that trial again...';
            end
            
            % Display text
            DrawFormattedText(window,str1,'center',Y-120,[255 255 255]);
            DrawFormattedText(window,str2,'center',Y-40,[255 255 255]);
            DrawFormattedText(window,str3,'center',Y+200,[255 255 255]);
            Screen('Flip', window);
            
            % Wait for response
            WaitSecs(0.5);
            keyIsDown = 0;
            while ~keyIsDown
                [keyIsDown, ~, ~] = KbCheck;  % checking for response
            end
        else cresp = 1;
            break
        end
    end
    
    
    %% Move onto sequence of practice trials with real timing
    timing.mem       =     0.5;    % memoranda presentation time
    timing.probe     =     0.5;    % probe presentation time
    
    % Specify text
    if strcmp(language,'german')
        str1 = sprintf('Hervorragend! Die Objekte waren am SELBEN Ort.');
        str2 = sprintf('Als nächstes werden wir ein paar mehr Beispiele machen, in denen die');
        str3 = 'Objekte kürzer auf dem Bildschirm sein werden. So wird auch die Aufgabe sein.';
        str4 = 'Sie sollten ab jetzt versuchen erst nach Ihrer Antwort, während der';
        str5 = 'Text auf dem Bildschirm "Jetzt Blinzeln" zeigt, zu blinzeln.';
        str6 = 'Drücken Sie einen beliebigen Knopf zum Starten...';
    elseif strcmp(language,'english')
        str1 = sprintf('Excellent! The objects were in the same location.');
        str2 = sprintf('Now we will do a few more practice trials, this time with the objects');
        str3 = 'on screen for shorter. This is how they will appear for the main task.';
        str4 = 'You should also now try to blink only after responding,';
        str5 = 'when the text on screen notifies you to do so.';
        str6 = 'Press any button to begin...';
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-200,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y-80,[255 255 255]);
    DrawFormattedText(window,str3,'center',Y-30,[255 255 255]);
    DrawFormattedText(window,str4,'center',Y+80,[255 255 255]);
    DrawFormattedText(window,str5,'center',Y+130,[255 255 255]);
    DrawFormattedText(window,str6,'center',Y+240,[255 255 255]);
    Screen('Flip', window);
    
    % Wait for response
    WaitSecs(0.5);
    keyIsDown = 0;
    while ~keyIsDown
        [keyIsDown, ~, ~] = KbCheck;  % checking for response
    end
    
    % Specify stimulus details for this trial
    stimIn = [angles(7)     angles(7)     3.0 1;
              angles(3)     angles(4)     1.0 2;
              angles(9)     angles(3)     3.0 3;
              angles(end-1) angles(end-1) 9.0 1;
              angles(5)     angles(5)     1.0 1;
              angles(8)     angles(7)     9.0 2];
    [XstimS,YstimS] = pol2cart(deg2rad(stimIn(:,1)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimS = XstimS+ X; YstimS = YstimS+Y;  % reference to screen center
    [XstimP,YstimP] = pol2cart(deg2rad(stimIn(:,2)),stim.s_r);  % get XY coords relative to origin of zero [0deg=right-most point on circle, 90deg=lower most point, moving CCW]
    XstimP = XstimP+ X; YstimP = YstimP+Y;  % reference to screen center
    stimInCoord = [XstimS YstimS XstimP YstimP];
  
    
    %% Countdown and first fixation point
    if strcmp(language,'german')
        str1 = 'Bereit machen!';
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
    
    on_time = [];  % onset time for first trial is specified within trial script
    Behav = zeros(size(stimIn,1),size(stimIn,2)+5);  % Behav = [stimIn(:,1:4) resp ACC RT start_time sacc]
    tRefresh = zeros(size(stimIn,1),8);  % tRefresh = [start_time t_mem_off t_probe_on t_probe_off t_resp t_fb_on t_rest_on t_rest_off]
    for t = 1:size(stimIn,1);
        
        % Variable task parameters
        varopts.on_time =  on_time;                          % controls onset time of impending trial - fed back from trial function
        varopts.stimM   =  [stimIn(t,1) stimInCoord(t,1:2)]; % angle and X/Y coordinates of memorandum
        varopts.stimP   =  [stimIn(t,2) stimInCoord(t,3:4)]; % angle and X/Y coordinates of memorandum
        varopts.delay   =  stimIn(t,3);                      % delay duration
        varopts.kbqdev  =  options.kbqdev;                 % keyboard info
        varopts.ppd     =  options.ppd;                      % pixels per degree
        
        % Run trial
        [Behav(t,[1:3 5:end]),tRefresh(t,:),on_time] = WM_dm2s_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, t, strcmp(options.et,'yes'), language);
        Behav(t,4) = stimIn(t,4);
    end
    
    WaitSecs(3);
    trigger(trigger_enc.block_end);  % trigger to mark end of block
    
    
    
    %% Exit
    % Specify text
    if strcmp(language,'german')
        str1 = sprintf('Sie haben das Training beendet.');
        str2 = sprintf('Gut gemacht!');
    elseif strcmp(language,'english')
        str1 = sprintf('You have now finished the training.');
        str2 = sprintf('Well done!');
    end
    
    % Display text
    DrawFormattedText(window,str1,'center',Y-40,[255 255 255]);
    DrawFormattedText(window,str2,'center',Y+80,[255 255 255]);
    Screen('Flip', window);
    WaitSecs(3)
    
    % Stop Eyelink
    if strcmp(options.et,'yes');
        Eyelink('CloseFile');
        Eyelink('WaitForModeReady', 500);
        Eyelink('StopRecording');
    end
    
    % Return to cmd line
    ListenChar(0)
    sca
    
    
catch ME  % if any errors are encountered or user requests quit, clean up and exit
    ListenChar(0)
    sca
    rethrow(ME)
end



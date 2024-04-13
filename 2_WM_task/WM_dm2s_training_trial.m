% Presents one trial of luminance patch signal detection task.
%
% Two circular patches are presented in the lower hemifield with luminances
% that fluctuate independently over time. At some point during the trial,
% the luminance of one patch may increase (on average) for some
% pre-specified duration. Trial ends with written feedback on accuracy.
%
% Outputs:
%    Behav = [signalPos signalDur signalOn signalStr RTabs RTrel ACCabs ACCdet start_time]
%    tRefresh = timing of all flips during trial
%    onset_nplus1 = specifies time for next stimulus presentation

function [Behav] = WM_dm2s_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, EL)

% Response parameters   %%%%%%%%%%%%% TO DO: COUNTER-BALANCE RESPONDING HAND ACROSS SUBJECTS
respL = '1';
respR = '4';
quitkey = 'ESCAPE';

flush_kbqueues(varopts.kbqdev);
% FlushEvents('keyDown');    %%%%%%%%%%% PM: this may mess things up......

% Get pixel coordinate of screen center
[X,Y] = RectCenter(windowRect);

% Getting eye that's being measured (matters for retrieving online gaze data)
if EL, eyeused = Eyelink('EyeAvailable'); check_saccade(eyeused, X, Y, varopts.ppd); end  % and flush EL buffer
sacc=0; % and initializing saccade counter

% Fixation cross positions and colours
f_ps = [X-stim.fix_t/2, Y-stim.fix_l/2, X+stim.fix_t/2, Y+stim.fix_l/2;...  % middle vertical bar of fixation cross
        X-stim.fix_l/2, Y-stim.fix_t/2, X,              Y+stim.fix_t/2;...  % left horizontal bar of fixation cross
        X             , Y-stim.fix_t/2, X+stim.fix_l/2, Y+stim.fix_t/2]';   % right horizontal bar of fixation cross

f_cs = [stim.fix_c_active' stim.fix_c_active' stim.fix_c_active'];  % initializing colours of all fixation cross parts (changes often)

f_ps_in = [X-stim.fix_t/2; Y-stim.fix_t/2; X+stim.fix_t/2; Y+stim.fix_t/2];  % position of inner fixation circle
f_cs_in = [127; 127; 127];  % colour of inner/outer fixation circles (doesn't change)

% Present active fixation
Screen('FillRect', window, f_cs, f_ps);
Screen('FillOval', window, f_cs_in, f_ps_in);
vbl = Screen('Flip', window);
trigger(trigger_enc.rest_off);  % trigger to mark onset of baseline period
if EL, Eyelink('message', num2str(trigger_enc.rest_off));  check_saccade(eyeused, X, Y, varopts.ppd); end

varopts.on_time = vbl + 1.75;  % onset time will be fixed @ 1.75s for training trials

WaitSecs(0.75);  if EL, check_saccade(eyeused, X, Y, varopts.ppd); end  % giving extra grace period of 0.5s where fixation break is allowed

% Present memorandum
g_ps = OffsetRect(gabor.rect, varopts.stimM(2), varopts.stimM(3));   % set new position of checkerboard
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(1), [], g_ps);  % drawing checkerboard
Screen('FillRect', window, f_cs, f_ps);  % draw fixation cross
Screen('FillOval', window, f_cs_in, f_ps_in);  % draw overlaid inner fixation point
start_time = Screen('Flip', window, varopts.on_time);  % set stimulus to be flipped at specified trial onset time
trigger(trigger_enc.mem_on);
if EL, Eyelink('message', num2str(trigger_enc.mem_on)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

if timing.s_flicker % Flicker the checker if desired
    cmask = 2;
    for f = 1:ceil((timing.mem/timing.s_refresh)-1)  % flicker checker
        Screen('DrawTexture', window, gabor.tex(cmask), [], g_ps);
        Screen('FillRect', window, f_cs, f_ps);  % draw fixation cross
        Screen('FillOval', window, f_cs_in, f_ps_in);  % draw overlaid inner fixation point
        Screen('Flip', window, start_time+(timing.s_refresh*f)-timing.ifi*0.5);
        cmask = find([1 2]~=cmask);
    end
end

Screen('FillRect', window, f_cs, f_ps);  % draw fixation only
Screen('FillOval', window, f_cs_in, f_ps_in);
t_mem_off = Screen('Flip', window, start_time + timing.mem - timing.ifi*0.5);
trigger(trigger_enc.mem_off);
if EL, Eyelink('message', num2str(trigger_enc.mem_off)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

if EL, WaitSecs(varopts.delay/2); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end  % add in extra saccade check halfway through delay period

% Present probe
g_ps = OffsetRect(gabor.rect, varopts.stimP(2), varopts.stimP(3));   % set new position of checkerboard
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('DrawTexture', window, gabor.tex(1), [], g_ps);  % drawing checkerboard
Screen('FillRect', window, f_cs, f_ps);  % draw fixation cross
Screen('FillOval', window, f_cs_in, f_ps_in);  % draw overlaid inner fixation point

flush_kbqueues(varopts.kbqdev);
% FlushEvents('keyDown'); % prepare for response checking
keyIsDown = false;

t_probe_on = Screen('Flip', window, t_mem_off + varopts.delay - timing.ifi*0.5);
trigger(trigger_enc.probe_on);
if EL, Eyelink('message', num2str(trigger_enc.probe_on)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

% Begin response checking
probe_removed = 0;  % flag to indicate whether probe has been removed yet

while ~keyIsDown
    [keyIsDown, firstPress] = check_kbqueues(varopts.kbqdev);
    % [keyIsDown, tdown, firstPress] = KbCheck;
    
    if keyIsDown  % logging response type
        tdown = GetSecs;
        RT = tdown-t_probe_on;  % logging RT relative to response cue onset
        keys = KbName(firstPress);  % retrieving string variable containing currently pressed key(s)
        if iscell(keys)
            resp = 99;  % in case of a double-press...having this as first 'if' test means it takes absolute precedence
            trigger(trigger_enc.resp_bad);  % trigger to mark a bad response
            if EL, Eyelink('message', num2str(trigger_enc.resp_bad)); end
        else
            switch keys
                case quitkey  % user requests quit
                    throw(MException('EXP:Quit', 'User request quit'));
                case {respL, 'LeftArrow', '1!', 'LeftControl'}
                    resp = -1;
                    trigger(trigger_enc.resp_left);  % trigger to mark a left (NO-MATCH) response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_left)); end
                case {respR, 'RightArrow', '4$', 'RightControl'}
                    resp = 1;
                    trigger(trigger_enc.resp_right);  % trigger to mark a right (MATCH) response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_right)); end
                otherwise
                    resp = 99;  % in case any button other than task relevant ones is pressed
                    trigger(trigger_enc.resp_bad);  % trigger to mark a bad response
                    if EL, Eyelink('message', num2str(trigger_enc.resp_bad)); end
            end
        end
    elseif ~probe_removed && GetSecs-t_probe_on+timing.ifi*0.5 > timing.probe   % remove probe if sufficient time has elasped
        f_cs = [stim.fix_c_go' stim.fix_c_go' stim.fix_c_go'];  % Go cue via fixation change
        Screen('FillRect', window, f_cs, f_ps);     % load just fixation into buffer for quick loading during response checking
        Screen('FillOval', window, f_cs_in, f_ps_in);
        t_probe_off = Screen('Flip', window);
        trigger(trigger_enc.probe_off);
        if EL, Eyelink('message', num2str(trigger_enc.probe_off)); end
        probe_removed = 1;
    end
end

if ~probe_removed  % remove probe if it hasn't been removed (i.e. if response was executed before scheduled probe offset)
    f_cs = [stim.fix_c_active' stim.fix_c_active' stim.fix_c_active'];
    Screen('FillRect', window, f_cs, f_ps);     % load just fixation into buffer for quick loading during response checking
    Screen('FillOval', window, f_cs_in, f_ps_in);
    t_probe_off = Screen('Flip', window, t_probe_on + timing.probe - timing.ifi*0.5);
    trigger(trigger_enc.probe_off);
    if EL, Eyelink('message', num2str(trigger_enc.probe_off)); end
end


% Present brief fixation-only before feedback
f_cs = [stim.fix_c_active' stim.fix_c_active' stim.fix_c_active'];
Screen('FillRect', window, f_cs, f_ps);     % load just fixation into buffer for quick loading during response checking
Screen('FillOval', window, f_cs_in, f_ps_in);
t_resp = Screen('Flip', window);


% Classify response accuracy
if (varopts.stimM(1) == varopts.stimP(1))
    cresp = 1;
else cresp = -1;
end

if resp == 99                            % bad press
    ACC = nan;
elseif resp == cresp                     % correct response
    ACC = 1;
else                                     % incorrect response
    ACC = 0;
end

% Concatenate final output variable
Behav = [varopts.stimM(1) varopts.stimP(1) varopts.delay resp ACC RT start_time sum(sacc)];

WaitSecs(0.5);

end

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

function Behav = WM_fixation_training_trial(window, windowRect, timing, stim, gabor, varopts, trigger_enc, EL)

timing.mem = 0.5;  % hardcoding stimulus presentation time to be 0.5s

% Response parameters   %%%%%%%%%%%%% TO DO: COUNTER-BALANCE RESPONDING HAND ACROSS SUBJECTS
quitkey = 'ESCAPE';

flush_kbqueues(varopts.kbqdev);
% FlushEvents('keyDown');    %%%%%%%%%%% PM: this may mess things up......

% Get pixel coordinate of screen center
[X,Y] = RectCenter(windowRect);

% Getting eye that's being measured (matters for retrieving online gaze data)
if EL, eyeused = Eyelink('EyeAvailable'); check_saccade(eyeused, X, Y, varopts.ppd); end  % also flush EL buffer
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
if EL, Eyelink('message', num2str(trigger_enc.rest_off));  [~,~,~] = Eyelink('GetQueuedData'); end

varopts.on_time = vbl + 1.2 + rand;  % onset time will be fixed uniformly distributed between 1.5 and 2.5s (taking into account 0.3s delay after stim offset)

if EL, nulls = check_saccade(eyeused, X, Y, varopts.ppd); end  % flush EL queue

% Present stim
g_ps = OffsetRect(gabor.rect, varopts.stim(2), varopts.stim(3));   % set new position of checkerboard
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
Screen('Flip', window, start_time + timing.mem - timing.ifi*0.5);
trigger(trigger_enc.mem_off);
if EL, Eyelink('message', num2str(trigger_enc.mem_off)); sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end

WaitSecs(0.3);

if EL, sacc(end+1) = check_saccade(eyeused, X, Y, varopts.ppd); end  % check for saccade

% Concatenate final output variable
Behav = [sum(sacc)];

end

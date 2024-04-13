function trig = setup_trigger()

trig.zero = 0;
trig.width = 0.005; 

trig.block_start = 1;  % start of block
trig.block_end = 2;    % end of block

trig.mem_on = 11;   % onset of memoranda
trig.mem_off = 12;  % offset of memoranda

trig.probe_on = 31;    % onset of probe
trig.probe_off = 32;   % offset of probe

trig.resp_left = 41;    % 'left' response
trig.resp_right = 42;   % 'right' response
trig.resp_bad = 43;     % bad response (either double press, or task-irrelevant button)

trig.fb_correct = 51;   % feedback for hit or correct rejection
trig.fb_error = 52;     % feedback for false alarm, miss, mislocalization, or premature response
trig.fb_bad = 53;       % feedback for bad responses

trig.rest_on = 61;    % onset of rest break period
trig.rest_off = 62;   % offset of rest break period

end

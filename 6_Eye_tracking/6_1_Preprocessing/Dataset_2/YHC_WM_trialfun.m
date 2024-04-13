 % Creates trial matrix 

function data = YHC_WM_trialfun(data,b)

event = data.event;
respdirs = [1 2 nan];  % left, right, bad resp
accs = [1 0 nan];  % correct, error, bad resp

% Pull indices of stimulus onset events
trialinds = find(cellfun(@(x) strcmp(x,'TRIALID'),{event.type}));  % TRIALID marks start of each trial

% Loop through each trial
trl = zeros(length(trialinds),12);  % initializing trial content matrix
smp = nan(length(trialinds),10);  % initializing sample times matrix


for i = 1:length(trialinds)
    
    j = trialinds(i);
    
    deets = tokenize(event(j).value);   % parsing event details
    trial_num = str2double(deets{find(strcmp(deets,'TRIALID'))+1});  % TRIAL NUMBER
    
    k = j+1;  % event counter
    s = 1;    % sample counter
    while sum(strcmp(event(k).type,{'TRIALID','2'}))==0
       if strcmp(event(k).type,'11')
          mem_on = event(k).sample;
          s = s+1; 
       elseif strcmp(event(k).type,'12')
          mem_off = event(k).sample;
       elseif strcmp(event(k).type,'31')
          probe_on = event(k).sample; 
       elseif strcmp(event(k).type,'32')
          probe_off = event(k).sample; 

       elseif strcmp(event(k).type,'41') || strcmp(event(k).type,'42') || strcmp(event(k).type,'43') %any response
          resp = event(k).sample; 
   

            if strcmp(event(k).type,'41')
              left_right_bad = 1; %left
             elseif strcmp(event(k).type,'42')
               left_right_bad = 2; %right
            elseif strcmp(event(k).type,'43')
             left_right_bad = 3; %bad
            end 

       elseif strcmp(event(k).type,'51') || strcmp(event(k).type,'52') || strcmp(event(k).type,'53') %any feedback
          fb = event(k).sample; 

             if strcmp(event(k).type,'51')
               correct_error_bad = 1; %correct
              elseif strcmp(event(k).type,'52')
                correct_error_bad = 2; %error
              elseif strcmp(event(k).type,'53')
              correct_error_bad = 3; %bad 
             end 

       elseif strcmp(event(k).type,'61')
          rest_start = event(k).sample; 
       elseif strcmp(event(k).type,'62')
          rest_end = event(k).sample; 
       end 

        k = k+1;
        
    end
    
    % Add to trial matrix
    trl(i,:) = [trial_num mem_on mem_off probe_on probe_off resp left_right_bad fb correct_error_bad rest_start rest_end b];
end

data.event = trl;  % add new event matrix into data structure (replacing old one)
data.eventsmp = smp;  % add sample times matrix into data structure

% Discard all data after block end marker
endind = find(cellfun(@(x) strcmp(x,'2'),{event.type}));  % '2' marks the end of each block
endind = event(endind).sample;

data.Xgaze          = data.Xgaze(1:endind);
data.Ygaze          = data.Ygaze(1:endind);
data.pupilL          = data.pupilL(1:endind);
data.pupilR          = data.pupilR(1:endind);
data.times          = data.times(1:endind);
data.sampleinfo     = [1 length(data.times)];
if ~isempty(data.blinksmp), data.blinksmp = data.blinksmp(data.blinksmp(:,2)<=length(data.times),:); end
if ~isempty(data.saccsmp), data.saccsmp   = data.saccsmp(data.saccsmp(:,2)<=length(data.times),:); end
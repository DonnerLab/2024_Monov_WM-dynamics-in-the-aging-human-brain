function [trl, event] = trialfun_MCI_continuous(cfgin)
%TRIALFUN_MCI_CONTINUOUS trial definition function
%  for continuous data
 
loadpath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/final_masterfile/';
load([loadpath,'Master_file.mat']);
 
    files = {zeros(length(file_info),1)};
    for f = 1:length(file_info)
 
        files{f,1}=file_info(f).name; 
    end 
idx_file = find(strcmp(files, cfgin.dataset));
ID = cfgin.ID;

% TrgVal = 1 - Start of block
% TrgVal = 2 - End of block

% read header and event information

hdr = ft_read_header(cfgin.dataset);
event = ft_read_event(cfgin.dataset);
 
% Search for trigger events

trgSample = [event(find(strcmp('UPPT001',{event.type}))).sample]';
trgTime = ([event(find(strcmp('UPPT001',{event.type}))).sample]')/1200; %event times in seconds relative to recording onset 

    
% Find specific triggers/samples for start and end of block
startSample = trgSample(find(trgTime==file_info(idx_file).blocktimes(cfgin.block,1)));
endSample = trgSample(find(trgTime==file_info(idx_file).blocktimes(cfgin.block,2)));


% Check whether start triggers are missing; add first start sample
if length(startSample) < length(endSample)
    startSample = [1+round(cfgin.trialdef.prestim*hdr.Fs); startSample];
end

% Check whether stop triggers are missing (interrupted block); add end
% sample
if length(startSample) > length(endSample)
    endSample = [endSample; trgSample(end)-round(cfgin.trialdef.poststim * hdr.Fs)];
end

% Check duration of block and discard if < 1 minutes
for i = 1:length(startSample)
    if (endSample(i)-startSample(i))/hdr.Fs < 60 % 60 seconds 
        startSample(i) = 0; % Set the start/end trigger samples that are not useable to zero to remove them afterwards
        endSample(i) = 0;
    end
end

startSample(startSample == 0) = [];
endSample(endSample == 0) = [];

if isfield(cfgin.trialdef, 'prestim')  
    trlOff = round(-cfgin.trialdef.prestim*hdr.Fs);
    startSample = max([ones(length(startSample),1),startSample + trlOff],[],2);
end


if isfield(cfgin.trialdef, 'poststim')
    endSample = min([event(end).sample*ones(length(endSample),1), endSample + round(cfgin.trialdef.poststim * hdr.Fs)],[],2);
end


% define trial matrix   trial x M (M = start sample; endSample; offset)
trl = [startSample, endSample, trlOff*ones(length(startSample),1)];

end
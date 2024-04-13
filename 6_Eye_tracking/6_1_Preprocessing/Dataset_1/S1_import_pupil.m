%importing .edf files as separate .asc files for events and samples

%% clear contents
clear 
close all
clc

%% define subjects 
% Load subject IDs
load('/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/behav_subj.mat');
hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final; 
subj = horzcat (hc,pat,cog_def);

%% add paths and set folder structure
for s = 1:length(subj)
%root directory for this project

%data folders
rawdir     = ['/mnt/homes/home028/gmonov/behav_eye_data/',subj{s},'/S1/Eyetracking/']; %folder where the .edf files are stored
wrtdir     = ['/mnt/homes/home028/gmonov/eye_tracking_data_analysis/asc_files/',subj{s},'/']; %sub-folder where the .asc files are saved to
edf2ascdir = '/mnt/homes/home030/aarazi/Toolboxes/Pupil_code-main/functions/edf2asc/'; %folder that contains the conversion program



%if the folders where the data are written out to don't yet exist, this creates them
if ~exist(wrtdir,'dir'); mkdir([wrtdir 'events']); mkdir([wrtdir 'samples']); mkdir([wrtdir 'samples_events']); end



%% get files to process

filz = dir([rawdir '*.edf']);
if strcmp(subj{s},'25')
    filz(1) = []; 
end 
%% loop over files

for fi = 1:length(filz)
    
    %% define output files, and check if they are already there
    
    edffile = filz(fi).name;

    
    
    %% convert EDF to ASCII file and put in correct folder
    
    cd(rawdir)
    
    %run the conversion to .asc
    disp(['working on ' edffile ' events']);
    system(['/mnt/homes/home030/aarazi/Toolboxes/Pupil_code-main/functions/edf2asc/edf2asc-linux -ns ' edffile ' ' edffile(1:end-4) '_e']); %-ns for no samples, so these are just the events

    
    disp(['working on ' edffile ' samples']);
    system(['/mnt/homes/home030/aarazi/Toolboxes/Pupil_code-main/functions/edf2asc/edf2asc-linux -ne ' edffile ' ' edffile(1:end-4) '_s']); %-ne for no events, so these are just the samples
    
    disp(['working on ' edffile ' samples_events']);
    system(['/mnt/homes/home030/aarazi/Toolboxes/Pupil_code-main/functions/edf2asc/edf2asc-linux ' edffile ' ' edffile(1:end-4) '_se']); % samples and events           

    
    
    movefile([edffile(1:end-4) '_e.asc'], [wrtdir 'events/'] )
    movefile([edffile(1:end-4) '_s.asc'], [wrtdir 'samples/'])
    movefile([edffile(1:end-4) '_se.asc'], [wrtdir 'samples_events/'])
   
end
 clearvars -except subj s
end



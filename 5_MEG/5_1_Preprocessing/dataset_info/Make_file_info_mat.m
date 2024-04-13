% FILE INFO SCRIPT 1: INITIALIZE MASTER FILE FOR MCI STUDY

clear, close all
clc

% Path/subject stuff
subj = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '01C' '02C' '26' '27' '28' '29' '30' '31' '32' '33' '35' '36' 'P03' '37' '38' '39' '40' '41' '42' '43' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54'};   
% 01 and 02 are 'C' subjects where MEG dataset has different naming scheme
% - should be '24' and '25', respectively; '15': MRI is named Pilot04
megpath = '/mnt/homes/home028/gmonov/meg_data/';
behavpath = '/mnt/homes/home028/gmonov/behav_eye_data/';
savepath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/';

% Get dataset names
filenames = dir([megpath,'*.ds']);

% Create final structure to contain all relevant file information
file_info = struct('name',{},'subj',{},'sess',{},'rec',{},'blocknums',{},'blocktimes',[],'trialex',{},...
    'rest',{},'resttime',{},'cVEOG',{},'cHEOG',{},'cECG',{},'cXgaze',{},'cYgaze',{},'cPupil',{},'ELex',{},'eventsM',{},'check',{});

% Fill in some basic stuff
for f = 1:length(filenames)
     if sum(strcmp(subj,filenames(f).name(1:3)))>0
        cdir = dir([megpath,filenames(f).name]);
        if sum([cdir.bytes])/1e09 > 0.1                                    % if this dataset has actual data recorded (and isn't a result of an MEG acq crash)
            file_info(end+1).name = filenames(f).name;                     % dataset name
            if strcmp('01C',filenames(f).name(1:3))
                file_info(end).subj = '24';
            elseif strcmp('02C',filenames(f).name(1:3))
                file_info(end).subj = '25';
            else file_info(end).subj = filenames(f).name(1:3);             % subject ID
            end
            
            if ~isnan(str2double(filenames(f).name(5)))
                file_info(end).sess = str2double(filenames(f).name(5));
                                                                           % session number in case of patients
            else file_info(end).sess = 1;                                  % if healthy control, set session number to 1
            end
                                                                      
            file_info(end).rec = str2double(filenames(f).name(end-3));     % recording run
        end
    
    elseif sum(strcmp(subj,filenames(f).name(1:2)))>0                      % if this is a subject in the main cohort
        cdir = dir([megpath,filenames(f).name]);
        
        if sum([cdir.bytes])/1e09 > 0.1                                    % if this dataset has actual data recorded (and isn't a result of an MEG acq crash)
            file_info(end+1).name = filenames(f).name;                     % dataset name
            file_info(end).subj = filenames(f).name(1:2);                  % subject ID
            if ~isnan(str2double(filenames(f).name(4)))
                file_info(end).sess = str2double(filenames(f).name(4));    % session number in case of patients
            else file_info(end).sess = 1;                                  % if healthy control, set session number to 1
            end
            if strcmp('11_MCI_20190923_01.ds',filenames(f).name)           % second session not named properly for participant 11
                file_info(end).sess = 2;
            end 
            file_info(end).rec = str2double(filenames(f).name(end-3));     % recording run
        end
     end
end

% Get task blocks per participant
% Session 1
 
for s = 1:length(subj)
          bls=[];
   
       filenames = dir([behavpath,subj{s},filesep,'S1/Behaviour/*.mat']);
       for b = 1:length(filenames)
              bls(end+1) = str2double(filenames(b).name(end-4));
       end
        for i=1:length(file_info); 
             if strcmp(file_info(i).subj,subj{s})==1 && file_info(i).sess==1
              file_info(i).blocknums = bls;
             end
        end 
end 
            
% Session 2

 for s = 1:length(subj)
    bls=[]; 
    try
    filenames = dir([behavpath,subj{s},filesep,'S2/Behaviour/*.mat']);
    catch 
    end
    for b = 1:length(filenames)
           bls(end+1) = str2double(filenames(b).name(end-4));
    end
    for i=1:length(file_info); 
             if strcmp(file_info(i).subj,subj{s})==1 && file_info(i).sess==2
              file_info(i).blocknums = bls;
             end
    end 
          
end 

           
% Marking Eyelink data as useable/unuseable

for f = 1:length(file_info)
    file_info(f).ELex = ones(size(file_info(f).blocknums));
    
    % '53' eyetracking data in block 3 couldn't be saved due to lacking
    % disk space
    if strcmp(file_info(f).subj,'53') 
        file_info(f).ELex = [1,1,0];
    end
    
    % Mark participants with protocolled quality problems in the eye tracking data in the first session
    if strcmp(file_info(f).subj,'04') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end  
    
    if strcmp(file_info(f).subj,'07') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 

    if strcmp(file_info(f).subj,'15') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'16') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'18') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'31') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
     if strcmp(file_info(f).subj,'39') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'42') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end  
     
    if strcmp(file_info(f).subj,'27') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
     if strcmp(file_info(f).subj,'39') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
   if strcmp(file_info(f).subj,'43') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end  
    
   if strcmp(file_info(f).subj,'45') && file_info(f).sess==1
        file_info(f).ELex = [0,1,1];
    end  
    
     if strcmp(file_info(f).subj,'47') && file_info(f).sess==1
        file_info(f).ELex = [0 1];
    end 
    
    if strcmp(file_info(f).subj,'53') && file_info(f).sess == 1
        file_info(f).ELex = [0,1,1];
    end  
    
    if strcmp(file_info(f).subj,'52') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'54') && file_info(f).sess==1
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    
    
   % Mark participants with protocolled quality problems in the eye tracking data in the second session
    if strcmp(file_info(f).subj,'01') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums)); 
    end 
    
    if strcmp(file_info(f).subj,'02') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'04') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'05') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end  
    
    if strcmp(file_info(f).subj,'07') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'11') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'16') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'29') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
    
    if strcmp(file_info(f).subj,'42') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end  
       
    if strcmp(file_info(f).subj,'43') && file_info(f).sess==2
        file_info(f).ELex = zeros(size(file_info(f).blocknums));
    end 
end

save([savepath,'Master_file.mat'],'file_info')
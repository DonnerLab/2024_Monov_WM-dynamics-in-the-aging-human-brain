% FILE INFO SCRIPT 3: ADD EOG/EL CHANNEL INFO

clear close all

% Path/subject stuff
subj = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '35' '36' 'P03' '37' '38' '39' '40' '41' '42' '43' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54'};   % 01 and 02 are 'C' subjects where MEG dataset has different naming scheme - should be 24 and 25, respectively
loadpath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/';
filepath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/ss_mat_files/';
savepath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/final_masterfile/';

% Load master file
load([loadpath,'Master_file.mat'])

% Generate cell array containing dataset-specific EOG/EL channels and add
% these to master file
chanmat = {'EEG057' 'EEG058' 'EEG059' 'UADC002' 'UADC003' 'UADC004'};

file_infoM = file_info;

% Load subject-specific file info structures and add to master file
for s = 1:length(subj)
          
   
    for f = 1:length(file_infoM);
        if strcmp(file_infoM(f).subj,subj{s})&&file_infoM(f).sess==1
            load([filepath,subj{s},'_1.mat'])
            file_infoM(f).blocktimes = file_info(f).blocktimes;
            file_infoM(f).trialex = file_info(f).trialex;
            file_infoM(f).rest = file_info(f).rest;
            file_infoM(f).resttime = file_info(f).resttime;
            file_infoM(f).eventsM = file_info(f).eventsM;
            file_infoM(f).check = file_info(f).check;
        end
        
        if strcmp(file_info(f).subj,subj{s})&&file_info(f).sess==2
            load([filepath,subj{s},'_2.mat'])
            file_infoM(f).blocktimes = file_info(f).blocktimes;
            file_infoM(f).trialex = file_info(f).trialex;
            file_infoM(f).rest = file_info(f).rest;
            file_infoM(f).resttime = file_info(f).resttime;
            file_infoM(f).eventsM = file_info(f).eventsM;
            file_infoM(f).check = file_info(f).check;
        end
         
        end 
    
    end
 
   
    
for f= 1:length(file_infoM)
   
        file_infoM(f).cVEOG = chanmat(1);
        file_infoM(f).cHEOG = chanmat(2);
        file_infoM(f).cECG = chanmat(3);
        file_infoM(f).cXgaze = chanmat(4);
        file_infoM(f).cYgaze = chanmat(5);
        file_infoM(f).cPupil = chanmat(6);
   %exceptions for eye data
   if strcmp(file_infoM(f).subj,'01')||strcmp(file_infoM(f).subj,'02')||strcmp(file_infoM(f).subj,'03')||strcmp(file_infoM(f).subj,'05')&&file_infoM(f).sess==1||strcmp(file_infoM(f).subj,'06')&&file_infoM(f).sess==1
       file_infoM(f).cXgaze = chanmat(4);
       file_infoM(f).cYgaze = chanmat(6);
       file_infoM(f).cPupil = chanmat(5);
   end 
  
   %exceptions for EOG channels --> HEOG and VEOG swapped
   if strcmp(file_infoM(f).subj,'01')&&file_infoM(f).sess==1||strcmp(file_infoM(f).subj,'11')&&file_infoM(f).sess==1||strcmp(file_infoM(f).subj,'41')&&file_infoM(f).sess==1||strcmp(file_infoM(f).subj,'45')||strcmp(file_infoM(f).subj,'49')
       file_infoM(f).cVEOG = chanmat(2);
       file_infoM(f).cHEOG = chanmat(1);
   end 
   % marking data where eyetracking shouldn't be used for preprocessing
   % because of quality issues / data not recorded due to technical
   % problems 
   if strcmp(file_infoM(f).subj,'11')&&file_infoM(f).sess==2||strcmp(file_infoM(f).subj,'15')||strcmp(file_infoM(f).subj,'36')||strcmp(file_infoM(f).subj,'37')||strcmp(file_infoM(f).subj,'38')
       file_infoM(f).cXgaze = 'no eyetracking';
       file_infoM(f).cYgaze = 'no eyetracking';
       file_infoM(f).cPupil = 'no eyetracking';
   end 
   % Final exception or inclusion based on looking at data
   % marking Eyelink data for subjects as bad(/good) after looking at block plots
   if strcmp(file_infoM(f).subj,'36')||strcmp(file_infoM(f).subj,'05')&&file_infoM(f).sess==2|| strcmp(file_infoM(f).subj,'08')||strcmp(file_infoM(f).subj,'21')||strcmp(file_infoM(f).subj,'45')&&file_infoM(f).sess==1||strcmp(file_infoM(f).subj,'47') %eye data isn't good enough for artifact rejection 
   file_infoM(f).ELex = zeros(size(file_infoM(f).ELex));
   end 
   if strcmp(file_infoM(f).subj,'01')&&file_infoM(f).sess==2||strcmp(file_infoM(f).subj,'05')&&file_infoM(f).sess==2||strcmp(file_infoM(f).subj,'18')||strcmp(file_infoM(f).subj,'29')&&file_infoM(f).sess==2 % data doesn't actually look that bad
   file_infoM(f).ELex = ones(size(file_infoM(f).ELex));
   end 
    if strcmp(file_infoM(f).subj,'10')&&file_infoM(f).sess==1% data in VEOg looks better to reflect blinks 
   file_infoM(f).ELex = zeros(size(file_infoM(f).ELex));
   end 

   if strcmp(file_infoM(f).subj,'16')&&file_infoM(f).sess==1 %only in the first block veog seems better than eyetracking data  
       file_infoM(f).ELex = [0 1 1]; 
   end
   
   if strcmp(file_infoM(f).subj,'03')&&file_infoM(f).sess==1 % calibration in the first block bad according to protocol
       file_infoM(f).ELex(1) = 0;
   end 
   if strcmp(file_infoM(f).subj,'36')||strcmp(file_infoM(f).subj,'37')||strcmp(file_infoM(f).subj,'38')%no proper eyetracking data recorded
       file_infoM(f).ELex = zeros(size(file_infoM(f).ELex));
   end
    if strcmp(file_infoM(f).subj,'39')%no proper eyetracking data recorded only in the first block 
       file_infoM(f).ELex = [0 1 1];
    end
   if strcmp(file_infoM(f).subj,'02')&& file_infoM(f).sess == 1%second block with interruptions in eyelink data 
       file_infoM(f).ELex = [1 0];
   end
    if strcmp(file_infoM(f).subj,'31')&& file_infoM(f).sess == 1%second block with interruptions in eyelink data, other blocks ok (many blinks, but consistent with EOG) 
       file_infoM(f).ELex = [1 0 1];
    end
       if strcmp(file_infoM(f).subj,'54')%eyelink data looks good 
       file_infoM(f).ELex = [1 1 1];
   end
   
   % Correcting stuff for exceptional datasets ~
   % 01-1: only 2 blocks because first dataset was cancelled early
   % 25: very first dataset doesn't have an associated behavioral file,
   % second one was cancelled after a few trials, events have been deleted 
   if strcmp(file_infoM(f).subj,'01')&&file_infoM(f).sess==1||strcmp(file_infoM(f).subj,'25')
       file_infoM(f).blocknums(end)=[]; % not counting the first block in blocknums 
       file_infoM(f).trialex(1,:)=[]; %removing the unused blocks from the trail exclusion matrices
       file_infoM(f).check=zeros(1,length((file_infoM(f).blocknums))); %setting the checks here to zero
       file_infoM(f).ELex(1)=[]; %deleting the eyelink exclusion, not using the first block anyways
   end  


   % no events recorded in for '42' and '43' in 2nd session (42: entire session, 43: only first block)
   if strcmp(file_infoM(f).subj,'43')&&file_infoM(f).sess==2
       if isnan(file_infoM(f).blocktimes)
       file_infoM(f).blocknums=[1]; % only first block in this recording without triggers 
       elseif ~isnan(file_infoM(f).blocktimes)
       file_infoM(f).blocknums=[1,2]; %only block 2 and 3 in this recording 
       file_infoM(f).trialex(1,:)=[]; %removing block without triggers from this recording 
       file_infoM(f).check=zeros(1,length((file_infoM(f).blocknums))); %setting the checks here to zero
       file_infoM(f).ELex=[0,0]; %deleting the eyelink exclusion, not using the first block anyways
       end  
   end
end 

% Save complete master file
file_info = file_infoM;
save([savepath,'Master_file.mat'],'file_info')
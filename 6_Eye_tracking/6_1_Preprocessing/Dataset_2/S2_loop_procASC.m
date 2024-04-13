% Convert asc-files to .mat 

clear all
close all

loadpath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_preproc/yhc_preproc/2.asc_converted/';
savepath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_preproc/yhc_preproc/3.mat_converted/';

addpath '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/'
addpath '/mnt/homes/home030/aarazi/fieldtrip-20201009'  % tell Matlab where FieldTrip is 
ft_defaults

load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/yhc.mat']);
allsubj = yhc;
for subj = 1:length(allsubj)
      % get files
        files = dir([[loadpath,allsubj{subj}] '*.mat']);

        % Loop through individual files
        errfiles={}; trial_counts=[];
        for f = 1:length(files)
            fprintf('File %s...\n',files(f).name)
            inFile = [loadpath,files(f).name];
            outFile = [savepath,files(f).name];

                load(inFile)
                    tic
                    % Create data and events structures, and matrices of blink/saccade times
                    data = ASC2dat(asc);

                    % Create refined event structure
                   data = YHC_WM_trialfun(data,f)
                    % Log [ID sess block ntrials]
                    trial_counts(end+1,1:4) = [str2double(files(f).name(2:4)) ...
                        str2double(files(f).name(6)) ...
                        str2double(files(f).name(8:end-4)) ...
                        size(data.event,1)];

                    % Save output
                    save(outFile,'data')

                    toc

        end

end 

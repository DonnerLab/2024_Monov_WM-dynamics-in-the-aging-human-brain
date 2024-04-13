% Interpolate pupil and gaze time-series

loadpath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_preproc/yhc_preproc/3.mat_converted/';
savepath = '/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_preproc/yhc_preproc/4.interpolated/';

% get list of all subject IDs
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/yhc.mat']);
allsubj = yhc;

for s = 1:length(allsubj)

      % pull files for only current subject
        subj = allsubj{s};

        files = dir([loadpath,subj,'*.mat']);

        errfiles={};
        for f = 1:length(files)
            fprintf('File %s...\n',files(f).name)
            inFile = [loadpath,files(f).name];
            outFile = [savepath,files(f).name];
            
            load(inFile)
            
        
                % Run interpolation routine & add results to data structure
                [newpupil, newXgaze, newYgaze, newblinksmp, badsmp, params, h, bestchan] = blink_interpolate_1chan(data,files(f).name);
                data.pupil = newpupil;
                data.Xgaze = newXgaze;
                data.Ygaze = newYgaze;
                data.newblinksmp = newblinksmp;
                data.badsmp = badsmp;
                data.interp_params = params;
                data.bestchan = bestchan;
                
                data = rmfield(data,{'pupilL','pupilR'});
                
                save(outFile,'data')
                savefig(h, [savepath,'plots',filesep,files(f).name(1:end-4),'.fig'],'compact')
                close all

        end
        
        if ~isempty(errfiles)
            save([savepath,'bad_files',filesep,subj,'_bad_files.mat'],'errfiles')
        end
        clearvars -except allsubj subj *path 
end
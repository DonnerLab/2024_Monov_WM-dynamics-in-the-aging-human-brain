% Save sample stimulus location as bilateral (<90 d.v.a.: 0 and >90
% d.v.a.: 1) 
% Gina Monov, UKE, 2023

clear all
close all

% Load subject IDs
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
% Loop over subjects
subj = horzcat(behav_mci_final,behav_hc_final,behav_cog_def_final); 
delays = [1,3,9]; 

for s=1:length(subj) 
    
    for d = 1:3 %loop over delay durations
    load(['/Users/ginamonov/Servers/mountpoint1/meg_analysis/Decoding/behav4decoding/',subj{s},'_1_',num2str(delays(d))]) % load exact sample stimulus (mem_loc) 
    bi_mem_loc = mem_loc; 
    bi_mem_loc(mem_loc<90) = 0; 
    bi_mem_loc(mem_loc>90) = 1; 
    clear mem_loc
    mem_loc = bi_mem_loc; 
    clear bi_mem_loc
    
   save(['/Users/ginamonov/Desktop/bi_mem_loc4decoding/',subj{s},'_1_',num2str(delays(d)),'.mat'],'mem_loc','tIDs') %save trial IDs and bilateral indication of sample stimulus in a new .mat file for the bilateral decoding 
    
   clear mem_loc tIDs
    end 
end 
% Extract for each subject specific behavioral data like miss rate and fa rate for later analysis (like the computation of SDT measures) 
% Gina Monov, UKE, 2022
   
clear all
close all
   
delays=[1,3,9];


load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);

subj = horzcat (behav_all_final,yhc);
loadpath1 = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'
loadpath2 = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_analysis/clean_WM_data/'
addpath('/Users/ginamonov/Servers/mountpoint1/functions/')
datamat = zeros(length(subj),length(delays));
fa=zeros(length(subj),length(delays));
CP_all=zeros(length(subj),3);

% loop through Subjects 
for s = 1:length(subj);
    if s <=length(behav_all_final)
    fullpath = [loadpath1, subj{s},filesep,'S1',filesep,subj{s}];
    
    load([fullpath,'_1_clean_allbehav.mat']);
    else    fullpath = [loadpath2, subj{s}];
    
    load([fullpath,'clean_allbehav.mat']);
    end 
    
        % exclude exceptional block due to performance
    if strcmp(subj{s},'02') 
       blockcount = allbehav(:,11);
       allbehav(blockcount == 1,:)=[]; 
    end 
    

    
    % pull behaviour
    delays = unique(allbehav(:,3));
    D = allbehav(:,3); T = allbehav(:,4); A = allbehav(:,6); 
    

   
    miss_rate = []; 
    nfa_rate = [];
    fa_rate = [];  
    H = []; 
    
    
    miss_rate = 1-(nansum(A(T==1)))./length(A(T==1));
    nfa_rate = 1-(nansum(A(T==2)))./length(A(T==2));
    ffa_rate = 1-(nansum(A(T==3)))./length(A(T==3));
    fa_rate = mean([nfa_rate,ffa_rate]); 
    
    % Compute SDT metrics
    H = -miss_rate +1; 
    
    
    save(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/miss-nfa/',subj{s},'.mat'] , 'fa_rate','nfa_rate','ffa_rate','H','miss_rate');
    
end 
% 1. Run Anova on group differences in fitted model parameters (YHC, OHC, MCI)
% 2. Compute Bayes Factor (evidence for null-hypothesis, i.e. group has no
%    effect on threshold)
%    Uses code from bayesFactor toolbox (Bart Krekelberg (2023). bayesFactor (https://github.com/klabhub/bayesFactor), GitHub. Retrieved December 2, 2023.)
% 3. Compare differences of group differences between threshold and
%    stochasticity parameters using two-sided permutation test
% Gina Monov, UKE, 2024 

clear all
close all
modelpath_mci = '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/sM+bound+theta';
modelpath_scz = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/sM+bound+theta';
addpath '/Users/ginamonov/Servers/mountpoint1/functions/';
% Load subject IDs
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
hc=behav_hc_final;
pat=behav_mci_final;
cog_def=behav_cog_def_final;
subj = horzcat(yhc,hc,pat,cog_def); 
g = [zeros(length(yhc),1);ones(length(hc),1);ones(length(pat),1)*2;ones(length(cog_def),1)*3]; 


for s = 1:length(subj) 
    
    if sum(ismember(subj{s},yhc))==1
       fullpath = [modelpath_scz,filesep, subj{s}];
    else fullpath = [modelpath_mci,filesep, subj{s}];
    end
         
    load([fullpath,'.mat']);
     
    % Pull fitted paramaters 
    noise(s,1) = pm(1); 
    criterion(s,1) = pm(2); 
    lapse(s,1) = pm(3); 
    clear pm 
    
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/lapse+noise/',subj{s},'_lapse_noise.mat']); 
    lapse_noise (s,1) = lapse_noise_single_subj; 
    clear lapse_noise_single_subj
end 

% Run one-way ANOVA for group effects (YHC, OHC, MCI) 
[p_criterion,anova_criterion]=anova1(criterion(g<3),g(g<3));
[p_noise,anova_noise]=anova1(noise(g<3),g(g<3));
[p_lapse,anova_lapse]=anova1(lapse(g<3),g(g<3));
[p_lapse_noise,anova_lapse_noise]=anova1(lapse_noise(g<3),g(g<3));
close all

%% Inititalize groups and group-specific vectors 
gg = g; 
gg(gg>0)=1; % young vs older subjects vector 

l_n_old = lapse_noise(gg==1); 
l_n_young = lapse_noise(gg==0); 

thres_old = criterion(gg==1); 
thres_ohc = criterion(g==1); 
thres_mci = criterion(g==2); 
thres_young = criterion(gg==0); 

%% Bayes factor 

[bf10,pValue] = bf_ttest2(thres_young,thres_old,'tail','both'); 
bf10_thres_old_young = bf10; 
clear bf10
bf10_p_thres_old_young = pValue; 
clear pValue 

[bf10,pValue] = bf_ttest2(thres_young,thres_mci,'tail','both'); 
bf10_thres_mci_young = bf10; 
clear bf10
bf10_p_thres_mci_young = pValue; 
clear pValue 

[bf10,pValue] = bf_ttest2(thres_young,thres_ohc,'tail','both'); 
bf10_thres_ohc_young = bf10; 
clear bf10
bf10_p_thres_ohc_young = pValue; 
clear pValue 

[bf10,pValue] = bf_ttest2(l_n_young,l_n_old,'tail','both'); 
bf10_l_n = bf10; 
clear bf10
bf10_p_l_n = pValue; 
clear pValue 

%% Compare the difference of the threshold/stochasticity parameters difference of older vs. younger subjects 

nperms = 10000; % Number of permutations 

comp_diffs_old_vs_young = permcompare_diffs(lapse_noise(gg==0),criterion(gg==0),lapse_noise(gg==1),criterion(gg==1),nperms); 

comp_diffs_mci_vs_young = permcompare_diffs(lapse_noise(g==0),criterion(g==0),lapse_noise(g==2),criterion(g==2),nperms); 

comp_diffs_ohc_vs_young = permcompare_diffs(lapse_noise(g==0),criterion(g==0),lapse_noise(g==1),criterion(g==1),nperms); 

comp_diffs_ohc_vs_mci = permcompare_diffs(lapse_noise(g==1),criterion(g==1),lapse_noise(g==2),criterion(g==2),nperms); 
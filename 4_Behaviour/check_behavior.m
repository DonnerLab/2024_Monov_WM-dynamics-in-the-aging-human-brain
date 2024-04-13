% Script for checking individuals subjects' behavior (WM accuracy) and the
% performance exclusion criterion (accuracy < 60% correct) 
% and to create a plot that shows this for each group of subjects 
% --> Extended Data Figure 1-1C
% Plotting: Violin plot toolbox: Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
% https://github.com/bastibe/Violinplot-Matlab, DOI:
% 10.5281/zenodo.4559847) 
% Compares group differences in overall task accuracy for subjects that pass the
% criterion (Permutation test: https://github.com/rudyvdbrink/Tools/tree/master/statistics) 
% Gina Monov, UKE, 2023

clear all 
close all

addpath('/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/history_bias/')
addpath('/Users/ginamonov/Servers/mountpoint1/functions/');
addpath('/Users/ginamonov/Servers/mountpoint1/functions/Violinplot-Matlab-master/');
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat'])
loadpath = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'
loadpath2 = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_analysis/clean_WM_data'
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);

hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final; 

% Add the excluded subjects with low performance to IDs
cog_def{end+1} = '52'; 
pat{end+1} = '09'; 
pat{end+1} = '10'; 
mci_subj = horzcat(hc,pat,cog_def);
subj = horzcat(mci_subj,yhc); 
    
n=[length(hc) length(pat) length(cog_def)];
nn = [1,length(hc); length(hc)+1, length(hc)+length(pat); length(hc)+length(pat)+1, length(mci_subj); length(mci_subj)+1,length(subj)]; 
 

    for s = 1:length(subj) 
        if s<=length(mci_subj)
        fullpath = [loadpath,filesep, subj{s},filesep,'S1',filesep];
        load([fullpath,subj{s},'_1_clean_allbehav.mat']);
               if strcmp(subj{s},'02') % exclude 1st block due to performance <50%
                 blockcount = allbehav(:,11);
                 mean(allbehav(blockcount==1,6))
                 allbehav(blockcount == 1,:)=[]; 
               end 
        accuracy(s,1) = mean(allbehav(:,6)); 
        n_trials(s,1) = length(allbehav(:,1)); 

        else
        fullpath = [loadpath2,filesep];
        load([fullpath,subj{s},'clean_allbehav.mat']);
        accuracy(s,1) = mean(allbehav(:,6)); 
        n_trials(s,1) = length(allbehav(:,1)); 
        end 
    end 
    
    mycolors(1).type = colors.sky; 
    myedgecolors(1).type = colors.teal; 
    mycolors(2).type = colors.rosered; 
    myedgecolors(2).type = colors.lipstick; 
    mycolors(3).type = colors.tortilla; 
    myedgecolors(3).type = colors.peanut; 
    mycolors(4).type = colors.yellow; 
    myedgecolors(4).type = colors.gold; 
    
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 8; % figure width
fig_h = 3; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];

% Violin plots 
acc.YHC = accuracy(nn(4,1):nn(4,2)); 
acc.OHC = accuracy(nn(1,1):nn(1,2)); 
acc.MCI = accuracy(nn(2,1):nn(2,2)); 
acc.UNC = accuracy(nn(3,1):nn(3,2)); 
violinplot(acc,{'YHC','OHC','MCI','UNC'},'GroupOrder',{'YHC','OHC','MCI','UNC'},'ViolinColor',[colors.yellow;colors.sky;colors.rosered;colors.peanut])
yline(0.6,'--','Color',colors.black,'LineWidth',1); 
set(gca,'box','off','TickDir','out')
ylabel('WM accuracy') 

% Group-wise comparison of error rates  
[h, p_YHC_OHC, diff, diff_null]=permtest2(acc.YHC,acc.OHC,10000);
[h, p_OHC_MCI, diff, diff_null]=permtest2(acc.OHC,acc.MCI(acc.MCI>0.6),10000);
[h, p_YHC_MCI, diff, diff_null]=permtest2(acc.YHC,acc.MCI(acc.MCI>0.6),10000);

savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/accuracy_supp1.fig']) % Extended Data Figure 1-1C
close all

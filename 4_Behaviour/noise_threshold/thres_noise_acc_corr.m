% Correlate WM accuracy with memory noise and threshold parameters
% Plot the results for all subjects, YHC, OHC, MCI separately 
% --> Extended Data Figure 3-2B
% Compare these correlations between the groups 
% Gina Monov, UKE, 2023

clear all
close all
corr_thres_noise = 'no'; %select if you want to correlate threshold and noise 

%load subj IDs, colors, model fits etc. 
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat']);

addpath('/Users/ginamonov/Servers/mountpoint1/functions/')

load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
modelpath_mci = '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/sM+bound+theta';
modelpath_scz = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/sM+bound+theta';
loadpath1 = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'; 
loadpath2 = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_analysis/clean_WM_data/'; 

% Load IDs 
hc=behav_hc_final;
mci=behav_mci_final;
cog_def=behav_cog_def_final;
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
subj = horzcat (yhc,hc,mci,cog_def);
n=[length(yhc) length(hc) length(mci) length(cog_def)];
nn = [1,length(subj); 1,length(yhc); length(yhc)+1,length(yhc)+length(hc); length(yhc)+length(hc)+1, length(yhc)+length(hc)+length(mci);length(yhc)+length(hc)+length(mci)+1, length(yhc)+length(hc)+length(mci)+length(cog_def)]; 
noise = []; 
lapse = []; 
lapse_noise = []; 
criterion = []; 
WMacc = []; 
PC1 = []; 

groups={'All','YHC','OHC','MCI'}; 

for s = 1:length(subj) 
    
    if ~ismember(subj{s},yhc)
        load([modelpath_mci,filesep,subj{s},'.mat'])
        noise(s,1)= pm(1);
        criterion(s,1)=pm(2); 
        lapse(s,1) = pm(3);
        clear pm 
        clear all_behav
        load([loadpath1,filesep,subj{s},filesep,'S1',filesep,subj{s},'_1_clean_allbehav.mat']);
           if strcmp(subj{s},'02')
              blockcount = allbehav(:,11);
              allbehav(blockcount == 1,:)=[]; 
           end 
        WMacc(s,1)= mean(allbehav(:,6)); 

        
           
    else
        load([modelpath_scz,filesep,subj{s},'.mat'])
        noise(s,1)= pm(1);
        criterion(s,1)=pm(2); 
        lapse(s,1) = pm(3);
        clear pm 
        clear all_behav
        load([loadpath2,filesep,subj{s},'clean_allbehav.mat']);
        WMacc(s,1)= mean(allbehav(:,6)); 
        
    end 
    
    load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/lapse+noise/' subj{s}, '_lapse_noise.mat']);
    lapse_noise(s,1) = lapse_noise_single_subj;
    clear lapse_noise_single_subj
end 


% Specify parameters to correlate with WM accuracy 
corr_param1_all = {'noise','criterion'}; % lapse_noise noise lapse WMacc PC1 criterion 
corr_param1_titles = {'Noise \sigma_{mem}','Threshold \delta'};

 
corr_param2 = 'WMacc'; % Define which paramter to correlate with
corr_param2_titles = {'WM accuracy'}; % give this a title 
corr_type = 'Pearson'; % Choose 'Pearson' or 'Spearman'


if strcmp(corr_thres_noise,'yes')% Select for looking at group-wise noise-threshold correlations
corr_param1_all = {'criterion'}; 
corr_param1_titles = {'Threshold \delta'};

 
corr_param2 = 'noise'; % Define which paramter to correlate with
corr_param2_titles = {'Noise \sigma_{mem}'}; % give this a title 
corr_type = 'Pearson'; % Choose 'Pearson' or 'Spearman'


end 
rows = 'Rows'; 
pw = 'Pairwise'; 
tt = 'type'; 
mec = 'MarkerEdgeColor'; 
mfc = 'MarkerFaceColor'; 
color = 'color';
lw = 'linewidth'; 
fs = 'Fontsize'; 

 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 15.5; % figure width
fig_h = 8; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


mycolors(1).type = colors.jetblack; 
mycolors(2).type = colors.canary; 
mycolors(3).type = colors.sky; 
mycolors(4).type = colors.rosered; 
myedgecolors(1).type = colors.black;
myedgecolors(2).type = colors.honey; 
myedgecolors(3).type = colors.teal; 
myedgecolors(4).type = colors.lipstick; 

scattersize = 15; 

% Make plot where each row contains the correlations for a single parameter
% to correlate and each column contains the groups (all, yhc, ohc, mci) 

counts = 0; 

for u = 1:length(corr_param1_all) 
    corr_param1 = corr_param1_all{u}; 
    
    for uu = 1:length(groups) 
       counts = counts+1;  
       eval(['s',num2str(counts),'=subplot(length(corr_param1_all),length(groups),counts), hold on']) 
       axis square
       set(gca,'fontname','Helvetica','LineWidth',0.5, 'Fontsize',7,'TickDir','out','box','off')
       eval(['[r,pval] = corr(',corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2,'(nn(uu,1):nn(uu,2)),rows,pw,tt,corr_type)']); 
       eval(['p = polyfit(', corr_param1,'(nn(uu,1):nn(uu,2)),', corr_param2,'(nn(uu,1):nn(uu,2)),1)']); 
       eval(['f = polyval(p,',corr_param1,'(nn(uu,1):nn(uu,2)))']); 
       eval(['l', num2str(counts), '=scatter(', corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2, '(nn(uu,1):nn(uu,2)),scattersize,mfc,mycolors(uu).type,mec,myedgecolors(uu).type),hold on'])
       if strcmp(corr_type,'Pearson') && pval <0.05
           eval(['sl', num2str(uu), '=plot(', corr_param1,'(nn(uu,1):nn(uu,2)),f,lw,1,color,mycolors(uu).type),hold on'])
       end 
      eval(['xlim([min(',corr_param1,'),','max(',corr_param1,')]);'])
      eval(['ylim([min(',corr_param2,'),','max(',corr_param2,')]);'])
      if strcmp(corr_thres_noise,'yes')
          ylabel(corr_param2_titles,fs,7)
      else
           ylabel('WM accuracy',fs,7)
      end 
      
      if uu == 2
      xlabel(corr_param1_titles{u},fs,7)
      end 
    if pval<10^-4
      title(['r=',num2str(round(r,2)),', p<10^{-4}'],fs,7); 
     elseif pval<0.05
        title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],fs,7); 
    else  title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],fs,7,'fontweight','normal'); 
    end 
    
    end 
end 

% Statictical test for correlation with threshold 

nperms = 10000;

corr_diffs_young_vs_old = permcompare_corrs(criterion(1:length(yhc)),WMacc(1:length(yhc)),criterion(length(yhc)+1:length(subj)),WMacc(length(yhc)+1:length(subj)),nperms); 

corr_diffs_young_vs_mci = permcompare_corrs(criterion(1:length(yhc)),WMacc(1:length(yhc)),criterion(nn(4,1):nn(4,2)),WMacc(nn(4,1):nn(4,2)),nperms); 

corr_diffs_young_vs_ohc = permcompare_corrs(criterion(1:length(yhc)),WMacc(1:length(yhc)),criterion(nn(3,1):nn(3,2)),WMacc(nn(3,1):nn(3,2)),nperms); 

corr_diffs_mci_vs_ohc = permcompare_corrs(criterion(nn(4,1):nn(4,2)),WMacc(nn(4,1):nn(4,2)),criterion(nn(3,1):nn(3,2)),WMacc(nn(3,1):nn(3,2)),nperms); 

% Perform correlational analyses of behavioural measures derived from the WM 
% task with PC1 score derived from the CERAD test battery
% Statistically compare the correlations betweeen OHC and MCI subgroups 
% Produces Figure 4 / Extended Data Figure 4-1A-B (depending on specified
% input in corr_param1_all)
% Gina Monov, UKE, 2023

clear all 
close all
clc

%load colors/cerad data/IDs
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat']);

load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')
addpath('/Users/ginamonov/Servers/mountpoint1/functions/')

load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);


% Specify parameters to correlate with PC1
corr_param1_all = {'WMacc','lapse_noise','criterion','noise','lapse'}; % lapse_noise noise lapse WMacc PC1 criterion 
corr_param1_titles = {'WM accuracy','Noise \sigma_{mem} + Lapse \theta','Threshold \delta','Noise \sigma_{mem}','Lapse \theta'};

 
corr_param2 = 'PC1'; % Define which paramter to correlate with
corr_param2_titles = {'PC1 score'}; % give this a title 
corr_type = 'Pearson'; % Choose 'Pearson' or 'Spearman'

% Load IDs 
hc=behav_hc_final;
mci=behav_mci_final;
mci=mci(~ismember(mci,'02'));  % exclude subject 2 here because one test of the CERAD test battery was missing 
cog_def=behav_cog_def_final;
subj = horzcat (hc,mci,cog_def);
n=[length(hc) length(mci) length(cog_def)];
nn = [1,length(subj); 1,length(hc);  length(hc)+1, length(hc)+length(mci)]; 
noise = []; 
lapse = []; 
lapse_noise = []; 
criterion = []; 
WMacc = []; 
PC1 = []; 

groups={'All','OHC','MCI'}; 
% Pull data
for s = 1:length(subj) 
     l = find(strcmp(cerad_data.ID,subj{s})); 
       
           noise(s,1)=cerad_data.noise(l); 
           lapse(s,1)=cerad_data.lapse(l); 
           criterion(s,1)=cerad_data.criterion(l); 
           PC1(s,1)=cerad_data.PC1_score_zsubj_all(l); 
           WMacc(s,1)=cerad_data.WMacc(l); 
      
           load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/lapse+noise/',subj{s},'_lapse_noise.mat']); 
           lapse_noise (s,1) = lapse_noise_single_subj; 
           clear lapse_noise_single_subj
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
fig_h = 25.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

mycolors(1).type = colors.jetblack; 
mycolors(2).type = colors.sky; 
mycolors(3).type = colors.rosered; 
myedgecolors(1).type = colors.black;
myedgecolors(2).type = colors.teal; 
myedgecolors(3).type = colors.lipstick; 

scattersize = 15; 

% Make plot where each row contains the correlations for a single parameter
% to correlate and each column contains the groups (all, ohc, mci) 

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
       if strcmp(corr_type,'Pearson') && pval <=0.05
           eval(['sl', num2str(uu), '=plot(', corr_param1,'(nn(uu,1):nn(uu,2)),f,lw,1,color,mycolors(uu).type),hold on'])
       end 
       eval(['xlim([min(',corr_param1,'),','max(',corr_param1,')]);'])
       eval(['ylim([min(',corr_param2,'),','max(',corr_param2,')]);'])
       ylabel('PC1 score',fs,7)
      
       if uu == 2
          xlabel(corr_param1_titles{u},fs,7)
       end 
      
       if pval<0.05
          title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,3))],fs,7); 
       else  title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,3))],fs,7,'fontweight','normal'); 
       end 
    
       if u == 1 
          legend(groups{uu})
       end 
    end 
    
    eval(['corr_diffs_',corr_param1,'=permcompare_corrs(',corr_param1,'(nn(3,1):nn(3,2)),',corr_param2,'(nn(3,1):nn(3,2)),',corr_param1,'(nn(2,1):nn(2,2)),',corr_param2,'(nn(2,1):nn(2,2)),10000);']) % Compare correlations OHC vs MCI
    
end 

saveas(figure(1),['/Users/ginamonov/Servers/mountpoint1/final_figures/Figure_4/figure4.fig']) % Figure 4, Extended Data Figure 4-1A-B
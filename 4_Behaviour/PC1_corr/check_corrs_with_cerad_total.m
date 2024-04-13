% Perform correlational analyses of behavioural measures derived from the WM task with CERAD total score 
% (Chandler et al.) and compare the correlations between OHC and MCI subgroups 
% --> Extended Data Figure 4-1C 
% Gina Monov, UKE, 2024

clear all 
close all
clc
addpath '/Users/ginamonov/Desktop/new_codes_before_submission/'
%load colors
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat']);

load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')
addpath('/Users/ginamonov/Servers/mountpoint1/functions/')

load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);


% Specify parameters to correlate with PC1
corr_param1_all = {'WMacc','lapse_noise','criterion','noise','lapse'}; % lapse_noise noise lapse WMacc PC1 criterion 
corr_param1_titles = {'WM accuracy','Noise \sigma_{mem} + Lapse \theta','Threshold \delta','Noise \sigma_{mem}','Lapse \theta'};

 
corr_param2 = 'cerad_total'; % Define which paramter to correlate with
corr_param2_titles = {'Cerad total score'}; % give this a title 
corr_type = 'Pearson'; % Choose 'Pearson' or 'Spearman'

% Load IDs 
hc=behav_hc_final;
mci=behav_mci_final;
mci=mci(~ismember(mci,'02'));  % exclude subject 2 here because one test of the cerad was missing 
cog_def=behav_cog_def_final;
subj = horzcat (hc,mci,cog_def);
n=[length(hc) length(mci) length(cog_def)];
nn = [1,length(subj); 1,length(hc);  length(hc)+1, length(hc)+length(mci)]; 
noise = []; 
lapse = []; 
lapse_noise = []; 
criterion = []; 
WMacc = []; 
cerad_total = []; 

groups={'All','OHC','MCI'}; 
% Pull data
for s = 1:length(subj) 
     l = find(strcmp(cerad_data.ID,subj{s})); 
       
           noise(s,1)=cerad_data.noise(l); 
           lapse(s,1)=cerad_data.lapse(l); 
           criterion(s,1)=cerad_data.criterion(l); 
           cerad_total(s,1)=cerad_data.totalscore(l); 
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
fig_w = 18; % figure width
fig_h = 7; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

mycolors(1).type = colors.jetblack; 
mycolors(2).type = colors.sky; 
mycolors(3).type = colors.rosered; 
myedgecolors(1).type = colors.black;
myedgecolors(2).type = colors.teal; 
myedgecolors(3).type = colors.lipstick; 

scattersize = 15; 

% Make plot where each subplot contains correlations for a single parameter
% to correlate with and each bar contains the correlation coefficient for
% each group (all, ohc, mci) with corresponding p-values on top of each bar

counts = 0; 
%% 
for u = 1:length(corr_param1_all) % loop over paramaters to correlate with cerad total score 
     corr_param1 = corr_param1_all{u}; 
     counts = counts +1; 
     eval(['s',num2str(counts),'=subplot(1,length(corr_param1_all),counts), hold on']) % make subplot
    
     for uu = 1:length(groups) 
        eval(['[r_', groups{uu},', pval_',groups{uu},'] = corr(',corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2,'(nn(uu,1):nn(uu,2)),rows,pw,tt,corr_type)']);  
     end 
    
     eval(['corr_diffs_',corr_param1,'=permcompare_corrs(',corr_param1,'(nn(3,1):nn(3,2)),',corr_param2,'(nn(3,1):nn(3,2)),',corr_param1,'(nn(2,1):nn(2,2)),',corr_param2,'(nn(2,1):nn(2,2)),10000);']) %Compare correlations
     eval(['p_diff = corr_diffs_', corr_param1, '.Delta_p;'])
     eval(['r_diff = corr_diffs_', corr_param1, '.Delta_r;'])
     p_group = ones(length(groups),length(groups)); 
     p_group(2,3) = p_diff; 
     p_group(3,2) = p_diff; 
   
     superbar([r_All, r_OHC, r_MCI],'BarFaceColor', [myedgecolors(1).type;myedgecolors(2).type;myedgecolors(3).type],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
     set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'All','OHC','MCI'}, 'Fontsize',7), xlim([0.3 3.7]), ylim([-1 1]), hold on 
     title(corr_param1_titles{u}, 'Fontsize',8), 
     
     if u == 1, ylabel('r'), end 
     
     ps = [pval_All, pval_OHC, pval_MCI]; 
     rs = [r_All, r_OHC, r_MCI]; 
     xc = 0.08; 
     
     for ll = 1:length(groups)
         if rs(ll)>0
            p_pos(ll)=rs(ll)+xc; 
         elseif rs(ll)<0
             p_pos(ll)=rs(ll)-xc; 
         end 
     end 
     for w = 1:length(groups) 
        if ps(w) <0.05
             fweight = 'bold';
        else fweight = 'normal';
        end 
        text(w-0.6, p_pos(w), ['p=',num2str(round(ps(w),3))],'Fontsize',7,'fontweight',fweight)%, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
     end 
     if p_diff < 0.05
         fweight = 'bold';
     else fweight = 'normal';
     end 
     text(0.9, -0.92,  ['\Deltar=',num2str(round(r_diff,2)), ', p=',num2str(round(p_diff,2))],'Fontsize',7,'fontweight',fweight)
     clear p_diff r_diff p_pos
end 

saveas(figure(1),['/Users/ginamonov/Servers/mountpoint1/final_figures/Figure_4/check_cerad_total.fig']) % Extended Data Figure 4-1C
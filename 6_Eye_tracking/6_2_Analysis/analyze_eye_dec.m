% Analyze eye-tracking decoding results 
% Plot time courses of decoding precision for all subjects, YHC, OHC and
% MCI separately 
% Compute the mean decoding precision during delay and compare these
% between the groups 
% Correlate the decoding precision with other behavioural measures (WM
% accuracy, fitted model parameters)
% Permutation test: https://github.com/rudyvdbrink/Tools/tree/master/statistics
% cluster-based permutation test: Edden M. Gerber (2023). permutest, MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/71737-permutest
% --> Figure 9B-E
% Gina Monov, UKE, 2024 

%% Load subject IDs, define paths etc. 

clear all 
close all
load('/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat')
load('/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat');
load('/Users/ginamonov/Servers/mountpoint1/meg_analysis/datasets_overview.mat') % trial exclusion
addpath '/Users/ginamonov/Servers/mountpoint1/functions/'
addpath '/Users/ginamonov/Desktop/new_codes_before_submission/'
load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')

hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final;
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']); %YHC IDs
mci_subj = horzcat (hc,pat,cog_def);
allsubj = horzcat (mci_subj,yhc); 

allsubj = allsubj(~ismember(allsubj,{'15','37'})); % Exclude subjects with technical problems in the recording  

sampling_rate = '250'; % Choose sampling rate to analyze (250 Hz or 1000 Hz (only older subjects))

% Find subjects with <60 trials after preprocessing & exclude them 
few_trials = {}; 
n_t = []; 
n_t_older = []; 
for s = 1:length(allsubj) 
   load(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/data4decoding/trials4decoding_x_y/', allsubj{s}, '_trial_data_x_y_',sampling_rate,'.mat']);
   n_trials(s) = size(trial_data_x_y,1); 
   if n_trials(s) < 60
      few_trials = horzcat(few_trials, allsubj{s}); 
   else n_t(end+1,1) = size(trial_data,1); 
       if ismember(allsubj{s},mci_subj) 
           n_t_older(end+1,1) = size(trial_data,1); 
       end 
   end 
  
end 
save(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/data4decoding/few_trials.mat'],'few_trials') %save these IDs for later 
allsubj = allsubj(~ismember(allsubj,few_trials)); 

% Find the subgroup IDs of included subjects 
older_subj = {}; 
yhc_ids_new = {}; 
mci_ids_new = {}; 
ohc_ids_new = {}; 
unc_ids_new = {}; 

for s = 1:length(allsubj) 
    if ismember(allsubj{s},mci_subj) 
       older_subj = horzcat(older_subj,allsubj{s}); 
    end 
    if ismember(allsubj{s},yhc) 
       yhc_ids_new = horzcat(yhc_ids_new,allsubj{s}); 
    elseif ismember(allsubj{s},pat) 
       mci_ids_new = horzcat(mci_ids_new,allsubj{s}); 
    elseif ismember(allsubj{s},hc) 
       ohc_ids_new = horzcat(ohc_ids_new,allsubj{s}); 
    elseif ismember(allsubj{s},cog_def) 
       unc_ids_new = horzcat(unc_ids_new,allsubj{s});  
    end 
end 
save(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/eye_subj_ids.mat'],'yhc_ids_new','ohc_ids_new','mci_ids_new','unc_ids_new') %save these IDs for later 

allsubj = {}; 
allsubj = horzcat(ohc_ids_new,mci_ids_new,unc_ids_new,yhc_ids_new); 
n_mci = length(mci_ids_new); 
n_ohc = length(ohc_ids_new); 
n_unc = length(unc_ids_new); 
n_yhc = length(yhc_ids_new); 

nn = []; 
nn(1,:) = [1, length(allsubj)]; % ALL 
nn(2,:) = [n_ohc+n_mci+n_unc+1 , n_ohc+n_mci+n_unc+n_yhc]; % YHC 
nn(3,:) = [1, n_ohc]; % OHC 
nn(4,:) = [n_ohc+1, n_ohc+n_mci];  % MCI
nn(5,:) = [n_ohc+n_mci+1, n_ohc+n_mci+n_unc]; % UNC 
clear s

%% Load fitted model parameters and WM accuracy 
PC1 = []; % Initialize variables
noise = []; 
lapse_noise = [];
lapse = []; 
criterion = []; 
WMacc = []; 

for s = 1:length(allsubj) 
      if ismember(allsubj{s},yhc)
          load(['/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/sM+bound+theta/', allsubj{s}, '.mat']); 
      else 
          load(['/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/sM+bound+theta/', allsubj{s}, '.mat']); 
          ii = find(strcmp(cerad_data.ID,allsubj{s})); 
          if strcmp(allsubj{s},cerad_data.ID(ii)) ~= 1 % just a check
              disp error 
          end 
          PC1(end+1,1) = cerad_data.PC1_score_zsubj_all(ii); 
      end 
           noise(s,1)=pm(1); 
           criterion(s,1)=pm(2);
           lapse(s,1)=pm(3); 
            
           WMacc(s,1)=mean(all_behav(:,6)); 

      load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/lapse+noise/',allsubj{s},'_lapse_noise.mat']); 
      lapse_noise (s,1) = lapse_noise_single_subj; 
      clear lapse_noise_single_subj all_behav pm 
    
end 

clear s
dec_prec_older = []; 

% Initialize time courses 
tc = []; 
tc_yhc = []; 
tc_mci = []; 
tc_ohc = []; 
tc_older = [];
for s = 1:length(allsubj) 
    dec_res = readtable(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/decoding_results/decode_', allsubj{s}, '_', sampling_rate,'_ndr_x_y.csv']);
    times = dec_res.latency; 
    mean_dec_prec(s,1) = mean(dec_res.test_correlation(times>=0.5 & times<1.5)); 
    eye_dec_prec = mean_dec_prec(s,1); 
    save(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/decoding_results/eye_mean_dec_prec/',allsubj{s},'_eye_mean_dec_prec_',sampling_rate,'.mat'],'eye_dec_prec')
    clear eye_dec_prec
    tc(s,:) = dec_res.test_correlation; 

    clear dec_res
end 
clear s

% Extract time courses for each subgroup 
tc_yhc = tc(nn(2,1):nn(2,2),:); 
tc_ohc = tc(nn(3,1):nn(3,2),:); 
tc_mci = tc(nn(4,1):nn(4,2),:); 
tc_unc = tc(nn(5,1):nn(5,2),:); 

%% Compare decoding precision during the delay 
% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 5; % figure width
fig_h = 5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


group_alloc = zeros(length(allsubj),1); 

for s = 1:length(allsubj)
    if ismember(allsubj{s},yhc) 
        group_alloc(s) = 1; 
    elseif ismember(allsubj{s},hc) 
        group_alloc(s) = 2; 
    elseif ismember(allsubj{s},pat) 
        group_alloc(s) = 3; 
    end 
end 


P = zeros(3,3); 
[h, p, diff, diff_null]=permtest2(mean_dec_prec(group_alloc == 1),mean_dec_prec(group_alloc == 2), 10000); 
P(2,1) = p; 
P(1,2) = p; 
[h, p, diff, diff_null]=permtest2(mean_dec_prec(group_alloc == 1),mean_dec_prec(group_alloc == 3), 10000); 
P(3,1) = p; 
P(1,3) = p; 
[h, p, diff, diff_null]=permtest2(mean_dec_prec(group_alloc == 2),mean_dec_prec(group_alloc == 3), 10000); 
P(3,2) = p; 
P(2,3) = p; 

for r = 1:3
errorbars(1,r) = std(mean_dec_prec(group_alloc == r),1)./sqrt(size(mean_dec_prec(group_alloc == r),1)); 
end

superbar([mean(mean_dec_prec(group_alloc == 1),1),mean(mean_dec_prec(group_alloc == 2),1),mean(mean_dec_prec(group_alloc == 3),1)],'E',errorbars,'P',P,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',true, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC' , 'OHC','MCI'}, 'Fontsize',8), xlim([0.2 3.8]), ylim([-0.1,0.6])
ax = gca; ax.XColor=[0 0 0];
ax = gca; ax.YColor=[0 0 0];
ylabel('Decoding precision during delay','Fontsize',8,'Color',colors.black)
% test decoding precision against zero 

[h, p, diff, diff_null]=permtest(mean_dec_prec(group_alloc == 1),0,10000,0.05,'right');
p_yhc = p; 
[h, p, diff, diff_null]=permtest(mean_dec_prec(group_alloc == 2),0,10000,0.05,'right');
p_ohc = p; 
[h, p, diff, diff_null]=permtest(mean_dec_prec(group_alloc == 3),0,10000,0.05,'right');
p_mci = p; 

if p_yhc < 10^-4
  h = text(0.5,0.01,['p<10^{-4}'],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8) 
elseif p_yhc < 0.05
  h = text(0.5,0.01,['p=',num2str(round(p_yhc,4))],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8)
else
  h = text(0.5,0.01,['p=',num2str(round(p_yhc,4))],'HorizontalAlignment','left','Fontweight','normal','Fontsize',8)
end 
set(h,'Rotation',90)
clear h 

if p_ohc < 10^-4
  h = text(1.5,0.01,['p<10^{-4}'],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8) 
elseif p_ohc < 0.05
  h = text(1.5,0.01,['p=',num2str(round(p_ohc,4))],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8)
else
  h = text(1.5,0.01,['p=',num2str(round(p_ohc,4))],'HorizontalAlignment','left','Fontweight','normal','Fontsize',8)
end 
set(h,'Rotation',90)
clear h 

if p_mci < 10^-4
  h = text(2.5,0.01,['p<10^{-4}'],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8) 
elseif p_mci < 0.05
  h = text(2.5,0.01,['p=',num2str(round(p_mci,4))],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8)
else
  h = text(2.5,0.01,['p=',num2str(round(p_mci,4))],'HorizontalAlignment','left','Fontweight','normal','Fontsize',8)
end 
set(h,'Rotation',90)
clear h 


 savefig(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/eye_figures/group_comparison',sampling_rate,'.fig']) % Figure 9C
 close all
    
%% Plot decoding precision time courses 
% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16; % figure width
fig_h = 5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
%% Subplot 1
% Perform cluster-based permutation test against zero for decoding time course 

data4stats = []; 
data4stats = tc'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);

    % Plotting 
    
subplot(1,2,1), hold on 
    s1=plot(times,mean(tc,1),'Color',colors.black), hold on 
    shadedErrorBar(times,mean(tc,1),std(tc,[],1)./sqrt(size(tc,1)),{'Color',colors.black},1)
         
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
             plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.025,'Color',colors.black,'LineWidth',1.5), hold on 
             end
         end 
         
clear clusters                 
onsets = [0 0.5]; %mem on and mem off
xline(onsets(1),'--','Color',colors.black,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.black,'LineWidth',1), hold on 
yline(0,'-','Color',colors.black), hold on % Reference line around zero
ylim([-0.1  0.4])
xlim([-0.1,1.5])
ylabel('Correlation coefficient','Fontsize',8,'color',colors.black)         
xlabel('Time from stimulus onset (s)','Fontsize',8,'color',colors.black)   
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.1 0.2 0.3 0.4])
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{ '0'  '0.1' '0.2'  '0.3' '0.4'}) 
ax = gca; ax.XColor=[0 0 0];
ax = gca; ax.YColor=[0 0 0];
legend([s1],['All subjects, N=', num2str(length(allsubj))])

 %% Subplot 2        
subplot(1,2,2), hold on          
% Perform cluster-based permutation test against zero for decoding time course 

data4stats = []; 
data4stats = tc_yhc'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);

    % Plotting 
     
    s2=plot(times,mean(tc_yhc,1),'Color',colors.yellow), hold on 
    shadedErrorBar(times,mean(tc_yhc,1),std(tc_yhc,[],1)./sqrt(size(tc_yhc,1)),{'Color',colors.yellow},1)
         
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
             plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.1,'Color',colors.yellow,'LineWidth',1.5), hold on 
             end
         end 
clear clusters

         data4stats = []; 
data4stats = tc_mci'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);

    % Plotting 
     
    s3=plot(times,mean(tc_mci,1),'Color',colors.rosered), hold on 
    shadedErrorBar(times,mean(tc_mci,1),std(tc_mci,[],1)./sqrt(size(tc_mci,1)),{'Color',colors.rosered},1)
         
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
             plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.13,'Color',colors.rosered,'LineWidth',1.5), hold on 
             end
         end 
clear clusters
         
         data4stats = []; 
data4stats = tc_ohc'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);

    % Plotting 
     
    s4=plot(times,mean(tc_ohc,1),'Color',colors.sky), hold on 
    shadedErrorBar(times,mean(tc_ohc,1),std(tc_ohc,[],1)./sqrt(size(tc_ohc,1)),{'Color',colors.sky},1)
         
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
             plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.16,'Color',colors.sky,'LineWidth',1.5), hold on 
             end
         end 
clear clusters

set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7])
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{ '0'  '0.1' '0.2'  '0.3' '0.4' '0.5' '0.6' '0.7'})         
ax = gca; ax.XColor=[0 0 0];
ax = gca; ax.YColor=[0 0 0];
onsets = [0 0.5]; %mem on and mem off
xline(onsets(1),'--','Color',colors.black,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.black,'LineWidth',1), hold on 
yline(0,'-','Color',colors.black), hold on % Reference line around zero
ylim([-0.2  0.7])
ylabel('Correlation coefficient','Fontsize',8)         
xlabel('Time from stimulus onset (s)','Fontsize',8)      
legend([s2,s3,s4],['YHC, N=', num2str(n_yhc)],['MCI, N=', num2str(n_mci)],['OHC, N=', num2str(n_ohc)])         
xlim([-0.1,1.5])

savefig(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/eye_figures/eye_dec_time_courses',sampling_rate,'.fig']) % Figure 9B and Figure 9D 
close all
    
  
  
  %% Run correlation with behavioral parameters%%%%%%%%%%%% 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify parameters to correlate with WM accuracy and fitted model
% paramaters 
corr_param1_all = {'WMacc','lapse_noise','criterion'}; %,'noise','lapse'}; % lapse_noise noise lapse WMacc PC1 criterion 
corr_param1_titles = {'WM accuracy','Noise \sigma_{mem} + Lapse \theta','Threshold \delta'};%,'Noise \sigma_{mem}','Lapse \theta'};

corr_param2 = 'mean_dec_prec'; % Define which paramter to correlate with
corr_param2_titles = {'Decoding precision during delay'}; % give this a title 


% Specifiy colors 
mycolors(1).type = colors.jetblack; %all
mycolors(2).type = colors.yellow; %yhc
mycolors(3).type = colors.sky; %ohc
mycolors(4).type = colors.rosered; %mci
myedgecolors(1).type = colors.black;
myedgecolors(2).type = colors.honey; 
myedgecolors(3).type = colors.teal; 
myedgecolors(4).type = colors.lipstick; 

rows = 'Rows'; 
pw = 'Pairwise'; 
tt = 'type'; 
mec = 'MarkerEdgeColor'; 
mfc = 'MarkerFaceColor'; 
color = 'color';
lw = 'linewidth'; 
fs = 'Fontsize'; 
corr_type = 'Pearson'; % Choose 'Pearson' or 'Spearman'

% Make plot where each row contains the correlations for a single parameter
% to correlate and each column contains the groups (all, ohc, mci) 
% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16; % figure width
fig_h = 5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
counts = 0; 
groups={'All','YHC','OHC','MCI'}; 

for u = 1:length(corr_param1_all) % loop over paramaters to correlate with cerad total score 
    corr_param1 = corr_param1_all{u}; 
    counts = counts + 1; 
    eval(['s',num2str(counts),'=subplot(1,length(corr_param1_all),counts), hold on']) % make subplot
    
            for uu = 1:length(groups) 
               eval(['[r_', groups{uu},', pval_',groups{uu},'] = corr(',corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2,'(nn(uu,1):nn(uu,2)),rows,pw,tt,corr_type)']);  
            end 
    
     eval(['corr_diffs_',corr_param1,'=permcompare_corrs(',corr_param1,'(nn(3,1):nn(3,2)),',corr_param2,'(nn(3,1):nn(3,2)),',corr_param1,'(nn(4,1):nn(4,2)),',corr_param2,'(nn(4,1):nn(4,2)),10000);'])
     eval(['p_diff = corr_diffs_', corr_param1, '.Delta_p;'])
     eval(['r_diff = corr_diffs_', corr_param1, '.Delta_r;'])
     p_group = ones(length(groups),length(groups)); 
     p_group(2,3) = p_diff; 
     p_group(3,2) = p_diff; 
     %'P',p_group,
     superbar([r_All, r_YHC, r_OHC, r_MCI],'BarFaceColor', [myedgecolors(1).type;myedgecolors(2).type;myedgecolors(3).type;myedgecolors(4).type],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
     set(gca,'TickDir','out','XTick',1:4,'XTickLabel',{'All','YHC','OHC','MCI'}, 'Fontsize',7), xlim([0.3 4.7]), ylim([-0.45 0.45]), hold on 
     ax = gca; ax.XColor=[0 0 0];
     ax = gca; ax.YColor=[0 0 0];
     title(corr_param1_titles{u}, 'Fontsize',8), 
     if u ==1, ylabel('r'), end 
     
     ps = [pval_All, pval_YHC, pval_OHC, pval_MCI]; 
     rs = [r_All,r_YHC, r_OHC, r_MCI]; 
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

savefig(['/Users/ginamonov/Servers/mountpoint1/eye_tracking_data_analysis/eye_figures/eye_dec_bar_corrs_',num2str(sampling_rate),'.fig']) %Figure 9E 
close all


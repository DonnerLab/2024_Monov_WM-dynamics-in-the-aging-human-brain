% Compute R-squared for behavioral model, inlcuding the data an model
% predictions of all subjects included in the MCI study (YHC,OHC, MCI, UNC)
% This is performed separately for the two types of errors, i.e. false
% alarms (non-match trials) and misses (match trials) 
% --> Extended Data Figure 3-1E
% Gina Monov, UKE, 2023

clear all 
close all

% Load colors 
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')

% Specify paths 
modelpath_mci = '/mnt/homes/home028/gmonov/Modelling/Dec_Rule/sM+bound+theta';
modelpath_scz = '/mnt/homes/home028/gmonov/SCZ/WM_modeling/Dec_Rule/sM+bound+theta';

% Load subject IDs 
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/behav_subj.mat']);
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/yhc.mat']);
hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final; 
n=[length(pat) length(hc) length(yhc)];
subj = horzcat (hc,pat,cog_def,yhc);

delays = [1,3,9]; 

for s = 1:length(subj)    
    if sum(ismember(subj{s},yhc))==1
    fullpath = [modelpath_scz,filesep, subj{s}];
    else fullpath = [modelpath_mci,filesep, subj{s}];
    end
    load([fullpath,'.mat']);
    n_trials(s) = length(all_behav(:,6)); 
    % Match trials 
    mt_pred_1(s,1) = mean(CPs(all_behav(:,3)==1 & all_behav(:,4)==1)); 
    mt_pred_3(s,1) = mean(CPs(all_behav(:,3)==3 & all_behav(:,4)==1)); 
    mt_pred_9(s,1) = mean(CPs(all_behav(:,3)==9 & all_behav(:,4)==1)); 
    % model predicts p(different), this transformation of responses
    % required for match trials 
    acc = all_behav(:,6); 
    mt_behav_1(s,1) = 1-(mean(acc(all_behav(:,3)==1 & all_behav(:,4)==1)));
    mt_behav_3(s,1) = 1-(mean(acc(all_behav(:,3)==3 & all_behav(:,4)==1)));
    mt_behav_9(s,1) = 1-(mean(acc(all_behav(:,3)==9 & all_behav(:,4)==1)));
    
    % Far non-match trials 
    fnmt_pred_1(s,1) = mean(CPs(all_behav(:,3)==1 & all_behav(:,4)==3)); 
    fnmt_pred_3(s,1) = mean(CPs(all_behav(:,3)==3 & all_behav(:,4)==3)); 
    fnmt_pred_9(s,1) = mean(CPs(all_behav(:,3)==9 & all_behav(:,4)==3)); 
    % model predicts p(different)
    fnmt_behav_1(s,1) = mean(acc(all_behav(:,3)==1 & all_behav(:,4)==3));
    fnmt_behav_3(s,1) = mean(acc(all_behav(:,3)==3 & all_behav(:,4)==3));
    fnmt_behav_9(s,1) = mean(acc(all_behav(:,3)==9 & all_behav(:,4)==3));
    
    % Near non-match trials 
    nnmt_pred_1(s,1) = mean(CPs(all_behav(:,3)==1 & all_behav(:,4)==2)); 
    nnmt_pred_3(s,1) = mean(CPs(all_behav(:,3)==3 & all_behav(:,4)==2)); 
    nnmt_pred_9(s,1) = mean(CPs(all_behav(:,3)==9 & all_behav(:,4)==2)); 
    % model predicts p(different)
    nnmt_behav_1(s,1) = mean(acc(all_behav(:,3)==1 & all_behav(:,4)==2));
    nnmt_behav_3(s,1) = mean(acc(all_behav(:,3)==3 & all_behav(:,4)==2));
    nnmt_behav_9(s,1) = mean(acc(all_behav(:,3)==9 & all_behav(:,4)==2));     
end 

% Fit linear regression model for match trials and non-match trials
% separately 
mdl_mt = fitlm([mt_behav_1;mt_behav_3;mt_behav_9],[mt_pred_1;mt_pred_3;mt_pred_9]); 
mdl_nmt = fitlm([fnmt_behav_1;fnmt_behav_3;fnmt_behav_9;nnmt_behav_1;nnmt_behav_3;nnmt_behav_9],[fnmt_pred_1;fnmt_pred_3;fnmt_pred_9;nnmt_pred_1;nnmt_pred_3;nnmt_pred_9]); 

% Extract R-sqaured for both 
Model_statistics.MT_Rsquared_adj = mdl_mt.Rsquared.Adjusted;
Model_statistics.MT_Rsquared_ord = mdl_mt.Rsquared.Ordinary;
Model_statistics.MT_Rsquared_adj = mdl_mt.Rsquared.Adjusted;
Model_statistics.MT_Fstat = mdl_mt.ModelFitVsNullModel.Fstat;
Model_statistics.MT_p_val_Fstat = mdl_mt.ModelFitVsNullModel.Pvalue;

Model_statistics.nmt_Rsquared_adj = mdl_nmt.Rsquared.Adjusted;
Model_statistics.nmt_Rsquared_ord = mdl_nmt.Rsquared.Ordinary;
Model_statistics.nmt_Rsquared_adj = mdl_nmt.Rsquared.Adjusted;
Model_statistics.nmt_Fstat = mdl_nmt.ModelFitVsNullModel.Fstat;
Model_statistics.nmt_p_val_Fstat = mdl_nmt.ModelFitVsNullModel.Pvalue;

% Correlate 
[r_mt,p_val_mt]=corr([mt_behav_1;mt_behav_3;mt_behav_9],[mt_pred_1;mt_pred_3;mt_pred_9],'Type','Pearson'); 
[r_nmt,p_val_nmt]=corr([fnmt_behav_1;fnmt_behav_3;fnmt_behav_9;nnmt_behav_1;nnmt_behav_3;nnmt_behav_9],[fnmt_pred_1;fnmt_pred_3;fnmt_pred_9;nnmt_pred_1;nnmt_pred_3;nnmt_pred_9],'Type','Pearson'); 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 10; % figure width
fig_h = 4.66; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

s1 = subplot(1,2,1), hold on 
axis square
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')
scatter([mt_behav_1;mt_behav_3;mt_behav_9],[mt_pred_1;mt_pred_3;mt_pred_9],'MarkerFaceColor',colors.black,'MarkerEdgeColor',colors.black), hold on 

% Add stats into the figure 
if Model_statistics.MT_p_val_Fstat <10^-4
text(0.43,0.16,['R^2(adjusted)=',num2str(round(Model_statistics.MT_Rsquared_adj,2)); 'p(F-test)<10^{-4}'],'Fontsize',7,'fontname','Helvetica'), hold on 
end 
alpha(s1,.5)
ylabel('Model predictions')
xlabel('Data')
title('Match trials')
h = refline(1,0)
h.Color = colors.silver; 

s2 = subplot(1,2,2), hold on 
axis square
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')
scatter([fnmt_behav_1;fnmt_behav_3;fnmt_behav_9;nnmt_behav_1;nnmt_behav_3;nnmt_behav_9],[fnmt_pred_1;fnmt_pred_3;fnmt_pred_9;nnmt_pred_1;nnmt_pred_3;nnmt_pred_9],'MarkerFaceColor',colors.black,'MarkerEdgeColor',colors.black)
% Add stats
if Model_statistics.nmt_p_val_Fstat <10^-4
text(0.56,0.2,['R^2(adjusted)=',num2str(round(Model_statistics.nmt_Rsquared_adj,2)); 'p(F-test)<10^{-4} '],'Fontsize',7,'fontname','Helvetica')
end 
alpha(s2,.5) 
title('Non-match trials')
h = refline(1,0); 
h.Color = colors.silver

w = 0.3; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
h = 0.3; hgap = 0.08; hoff = 0.15;

savefig(['/mnt/homes/home028/gmonov/final_figures/R2_behavioral_model.fig']) % Extended Data Figure 3-1E
% Script for producing behvaioral plots + model predictions + comparison between model paramaters for the paper &
% perform mixed anova analysis: 
% --> Figure 2A-B, Figure 3B, Extended Data Figure 2-1A-B
% Additionally script computes lapse+noise parameters and saves this for
% each subject
% Mixed ANOVA analysis: Laurent Caplette. Simple RM/Mixed ANOVA for any design (https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design), MATLAB Central File Exchange. Retrieved November 20, 2019.
% Permutation test: https://github.com/rudyvdbrink/Tools/tree/master/statistics
% Bar plots: Scott Lowe (2022). superbar (https://github.com/scottclowe/superbar), GitHub.
% Shaded error bars: Rob Campbell (2019). raacampbell/shadedErrorBar (https://github.com/raacampbell/shadedErrorBar), GitHub.
% Gina Monov, UKE, 2022

clear all 
close all

%Load colors
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat']);
 
loadpath1 = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'; 
loadpath2 = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_analysis/clean_WM_data/'; 
addpath('/Users/ginamonov/Servers/mountpoint1/functions/')

modelpath_mci = '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/sM+bound+theta';
modelpath_scz = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/sM+bound+theta';

% Load subject IDs
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
hc=behav_hc_final;
pat=behav_mci_final;
cog_def = behav_cog_def_final; 
n=[length(pat) length(hc) length(yhc)];
subj = horzcat (hc,pat,cog_def,yhc);

delays = [1,3,9]; 
 
errors=zeros(length(subj),length(delays));
CP_all=zeros(length(subj),3);

% loop through Subjects
for s = 1:length(subj)
    if s <=length(hc)+length(pat)+length(cog_def)
    fullpath = [loadpath1, subj{s},filesep,'S1',filesep,subj{s}];
    
    load([fullpath,'_1_clean_allbehav.mat']);
    
    else    fullpath = [loadpath2, subj{s}];
    
    load([fullpath,'clean_allbehav.mat']);
    end 
    
    % exclude block 1 due to performance
    if strcmp(subj{s},'02')
       blockcount = allbehav(:,11);
       allbehav(blockcount == 1,:)=[]; 
    end 
    
    n_trials (s,1) = length(allbehav(:,1)); % pull number of trials
    
    % Pull behaviour
    delays = unique(allbehav(:,3));
    D = allbehav(:,3); T = allbehav(:,4); A = allbehav(:,6);
    
    % Pull error rates for each subjects across all trial types
    errors(s,1) = 1-(nansum(A(D==1)))./length(A(D==1));
    errors(s,2) = 1-(nansum(A(D==3)))./length(A(D==3));
    errors(s,3) = 1-(nansum(A(D==9)))./length(A(D==9));
    
    % Pull error rates specifically for the near distance trials (nfa) 
    nfa(s,1) = 1-(nansum(A(T==2 & D==1)))./length(A(T==2 & D==1));
    nfa(s,2) = 1-(nansum(A(T==2 & D==3)))./length(A(T==2 & D==3));
    nfa(s,3) = 1-(nansum(A(T==2 & D==9)))./length(A(T==2 & D==9));
    
        % Pull error rates specifically for the match trials
    misses(s,1) = 1-(nansum(A(T==1 & D==1)))./length(A(T==1 & D==1));
    misses(s,2) = 1-(nansum(A(T==1 & D==3)))./length(A(T==1 & D==3));
    misses(s,3) = 1-(nansum(A(T==1 & D==9)))./length(A(T==1 & D==9));
    
        % Pull error rates specifically for the far distance trials (ffa) 
    ffa(s,1) = 1-(nansum(A(T==3 & D==1)))./length(A(T==3 & D==1));
    ffa(s,2) = 1-(nansum(A(T==3 & D==3)))./length(A(T==3 & D==3));
    ffa(s,3) = 1-(nansum(A(T==3 & D==9)))./length(A(T==3 & D==9));
    
%     if length(A(isnan(A)))>0 % Check whether there are unexpected nans
%         disp yes
%     end 

    % Pull model parameters & model predictions for each subjects
    pm = []; 
    CPs = []; 
    if sum(ismember(subj{s},yhc))==1
    fullpath = [modelpath_scz,filesep, subj{s}];
    else fullpath = [modelpath_mci,filesep, subj{s}];
    end
    load([fullpath,'.mat']);
    
    delay = allbehav(:,3); trials = allbehav(:,4); ttypes = unique(trials); dtypes = unique(delay);
    
    % Pull model predictions 
    % Convert CPs to error rates 
    for o = 1:length(CPs) 
        if trials(o) > 1
            CPs(o) = 1-CPs(o); 
        end 
    end 
    
    for d = 1:length(dtypes)
           CP_av = []; 
           CP_av = mean(CPs(delay==dtypes(d)));
           CPs_all_ttypes(s,d)=CP_av;
    end
     
    for t = 1:length(ttypes)
        for d = 1:length(dtypes)
             CP_av = []; 
             CP_av = mean(CPs(trials==ttypes(t) & delay==dtypes(d)));
             CP_all(s,d,t)=CP_av;
        end
    end
    
    % Pull fitted paramaters 
    noise(s,1) = pm(1); 
    criterion(s,1) = pm(2); 
    lapse(s,1) = pm(3); 
end 

z_noise = zscore(noise); 
z_lapse = zscore(lapse); 

lapse_noise = z_noise + z_lapse; 

% Save combined lapse and noise parameter for all analyses

for ww = 1:length(subj) 
    lapse_noise_single_subj = []; 
    lapse_noise_single_subj = lapse_noise(ww); 
    save(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/lapse+noise/' subj{ww}, '_lapse_noise.mat'],'lapse_noise_single_subj')
end 

% Get all behavuoir and model predictions per trial type for each subgroup individually 

% NFA 
nfa_hc = nfa(1:length(hc),:); 
nfa_mci = nfa(length(hc)+1:length(hc)+length(pat),:); 
nfa_unclear = nfa(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def),:); 
nfa_yhc = nfa(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc),:); 

% FFA
ffa_hc = ffa(1:length(hc),:); 
ffa_mci = ffa(length(hc)+1:length(hc)+length(pat),:); 
ffa_unclear = ffa(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def),:); 
ffa_yhc = ffa(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc),:); 

% Misses
misses_hc = misses(1:length(hc),:); 
misses_mci = misses(length(hc)+1:length(hc)+length(pat),:); 
misses_unclear = misses(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def),:); 
misses_yhc = misses(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc),:); 


% Predictions for the different trial types 
nfa_preds = CP_all(:,:,2); 
ffa_preds = CP_all(:,:,3); 
misses_preds = CP_all(:,:,1); 

pred_nfa_hc = nfa_preds(1:length(hc),:); 
pred_nfa_mci = nfa_preds(length(hc)+1:length(hc)+length(pat),:); 
pred_nfa_unclear = nfa_preds(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def),:); 
pred_nfa_yhc = nfa_preds(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc),:); 

pred_ffa_hc = ffa_preds(1:length(hc),:); 
pred_ffa_mci = ffa_preds(length(hc)+1:length(hc)+length(pat),:); 
pred_ffa_unclear = ffa_preds(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def),:); 
pred_ffa_yhc = ffa_preds(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc),:); 

pred_misses_hc = misses_preds(1:length(hc),:); 
pred_misses_mci = misses_preds(length(hc)+1:length(hc)+length(pat),:); 
pred_misses_unclear = misses_preds(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def),:); 
pred_misses_yhc = misses_preds(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc),:); 

% Noise 
noise_hc = noise(1:length(hc)); 
noise_mci = noise(length(hc)+1:length(hc)+length(pat)); 
noise_unclear = noise(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def)); 
noise_yhc = noise(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc)); 

% Criterion 
criterion_hc = criterion(1:length(hc)); 
criterion_mci = criterion(length(hc)+1:length(hc)+length(pat)); 
criterion_unclear = criterion(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def)); 
criterion_yhc = criterion(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc)); 


% Lapse 
lapse_hc = lapse(1:length(hc)); 
lapse_mci = lapse(length(hc)+1:length(hc)+length(pat)); 
lapse_unclear = lapse(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def)); 
lapse_yhc = lapse(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc)); 

% Lapse + Noise
lapse_noise_hc = lapse_noise(1:length(hc)); 
lapse_noise_mci = lapse_noise(length(hc)+1:length(hc)+length(pat)); 
lapse_noise_unclear = lapse_noise(length(hc)+length(pat)+1:length(hc)+length(pat)+length(cog_def)); 
lapse_noise_yhc = lapse_noise(length(hc)+length(pat)+length(cog_def)+1:length(hc)+length(pat)+length(cog_def)+length(yhc)); 


% Run stats on group differences for near-fas and fitted model parameters
% Concatenate only the three groups that should be compared (mci, hc, yhc) 
group_alloc = [ones(length(yhc),1);ones(length(behav_hc_final),1)+1;(ones(length(behav_mci_final),1)+2)]; % Create info abut group allocation 
nfa3 = [nfa_yhc; nfa_hc; nfa_mci]; 
noise3 = [noise_yhc; noise_hc; noise_mci]; 
criterion3 = [criterion_yhc; criterion_hc; criterion_mci]; 
lapse3 = [lapse_yhc; lapse_hc; lapse_mci]; 
lapse_noise3 = [lapse_noise_yhc; lapse_noise_hc; lapse_noise_mci]; 
pred_nfa3 = [pred_nfa_yhc; pred_nfa_hc; pred_nfa_mci]; 

   tst = combntns([1:length(unique(group_alloc))],2);
   
   % Initialize matrix that stores p_values 
   P_nfa = zeros (3,3);
   P_noise = zeros (3,3);
   P_criterion = zeros (3,3);
   P_lapse = zeros (3,3);
   P_lapse_noise = zeros (3,3);
   
for l = 1:3
   P_nfa(l,l)=nan;  
   P_noise(l,l)=nan;  
   P_criterion(l,l)=nan;  
   P_lapse(l,l)=nan;  
   P_lapse_noise(l,l)=nan; 
   P_pred_nfa(l,l)=nan; 
end 

% Run permutation tests for different combinations

for l = 1:length(tst(:,1)) 
    % NFA
      p = []; 
   [h, p, diff, diff_null]=permtest2(nfa3(group_alloc == tst(l,1),1),nfa3(group_alloc == tst(l,2),1),10000);
   P_nfa_1(tst(l,1),tst(l,2)) = p; 
   P_nfa_1(tst(l,2),tst(l,1)) = p; 
   
      p = []; 
   [h, p, diff, diff_null]=permtest2(nfa3(group_alloc == tst(l,1),2),nfa3(group_alloc == tst(l,2),2),10000);
   P_nfa_3(tst(l,1),tst(l,2)) = p; 
   P_nfa_3(tst(l,2),tst(l,1)) = p; 
   
      p = []; 
   [h, p, diff, diff_null]=permtest2(nfa3(group_alloc == tst(l,1),3),nfa3(group_alloc == tst(l,2),3),10000);
   P_nfa_9(tst(l,1),tst(l,2)) = p; 
   P_nfa_9(tst(l,2),tst(l,1)) = p; 

   % Noise
      p = []; 
   [h, p, diff, diff_null]=permtest2(noise3(group_alloc == tst(l,1)),noise3(group_alloc == tst(l,2)),10000);
   P_noise(tst(l,1),tst(l,2)) = p; 
   P_noise(tst(l,2),tst(l,1)) = p; 
   
   % Criterion
      p = []; 
   [h, p, diff, diff_null]=permtest2(criterion3(group_alloc == tst(l,1)),criterion3(group_alloc == tst(l,2)),10000);
   P_criterion(tst(l,1),tst(l,2)) = p; 
   P_criterion(tst(l,2),tst(l,1)) = p; 
   
   % Lapse 
      p = []; 
   [h, p, diff, diff_null]=permtest2(lapse3(group_alloc == tst(l,1)),lapse3(group_alloc == tst(l,2)),10000);
   P_lapse(tst(l,1),tst(l,2)) = p; 
   P_lapse(tst(l,2),tst(l,1)) = p; 
   
   % Lapse + Noise
      p = []; 
   [h, p, diff, diff_null]=permtest2(lapse_noise3(group_alloc == tst(l,1)),lapse_noise3(group_alloc == tst(l,2)),10000);
   P_lapse_noise(tst(l,1),tst(l,2)) = p; 
   P_lapse_noise(tst(l,2),tst(l,1)) = p; 
   
end

% Create errorbars for plotting 

for r = 1:3 
    % Data
    err_nfa_1(1,r) = std(nfa3(group_alloc == r,1),1)./sqrt(size(nfa3(group_alloc == r,1),1)); % 1s Delay
    err_nfa_3(1,r) = std(nfa3(group_alloc == r,2),1)./sqrt(size(nfa3(group_alloc == r,2),1)); % 3s Delay
    err_nfa_9(1,r) = std(nfa3(group_alloc == r,3),1)./sqrt(size(nfa3(group_alloc == r,3),1)); % 9s Delay 
    % Model predictions
    err_pred_nfa_1(1,r) = std(pred_nfa3(group_alloc == r,1),1)./sqrt(size(pred_nfa3(group_alloc == r,1),1)); % 1s Delay
    err_pred_nfa_3(1,r) = std(pred_nfa3(group_alloc == r,2),1)./sqrt(size(pred_nfa3(group_alloc == r,2),1)); % 3s Delay
    err_pred_nfa_9(1,r) = std(pred_nfa3(group_alloc == r,3),1)./sqrt(size(pred_nfa3(group_alloc == r,3),1)); % 9s Delay 
    % Fitted parameters
    err_noise(1,r) = std(noise3(group_alloc == r),1)./sqrt(size(noise3(group_alloc == r),1)); % Noise
    err_criterion(1,r) = std(criterion3(group_alloc == r),1)./sqrt(size(criterion3(group_alloc == r),1)); % Criterion
    err_lapse(1,r) = std(lapse3(group_alloc == r),1)./sqrt(size(lapse3(group_alloc == r),1)); % Lapse 
    err_lapse_noise(1,r) = std(lapse_noise3(group_alloc == r),1)./sqrt(size(lapse_noise3(group_alloc == r),1)); % Lapse + Noise
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Mixed Anova Analysis for effect of delay duration%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errors4sma = errors; 
errors4sma(length(pat)+length(hc)+1:length(pat)+length(hc)+length(cog_def),:)= []; %Delete the cognitive deficit group
datamat (:,:,1) = errors4sma;
subject = [1:length(subj)-length(cog_def)]; 
between_factors = zeros(length(subj)-length(cog_def),1);
between_factors(subject>length(hc)) = 1;
between_factors(subject>(length(hc)+length(pat))) = 2;

%___ Naming the factors
within_factor_names = {'Delay'};
between_factor_names = {'Group'};

%___ Mixed Anova

tbl = simple_mixed_anova(datamat, between_factors,{'Delay'},{'Group'});

% Run analyses only for misses and ffa (no effect here) and nfa alone 

miss4sma = misses;
miss4sma(length(pat)+length(hc)+1:length(pat)+length(hc)+length(cog_def),:)= []; %Delete the cognitive deficit group
tbl_miss = simple_mixed_anova(miss4sma, between_factors,{'Delay'},{'Group'});

ffa4sma = ffa;
ffa4sma(length(pat)+length(hc)+1:length(pat)+length(hc)+length(cog_def),:)= []; %Delete the cognitive deficit group
tbl_ffa = simple_mixed_anova(ffa4sma, between_factors,{'Delay'},{'Group'});

nfa4sma = nfa;
nfa4sma(length(pat)+length(hc)+1:length(pat)+length(hc)+length(cog_def),:)= []; %Delete the cognitive deficit group
tbl_nfa = simple_mixed_anova(nfa4sma, between_factors,{'Delay'},{'Group'});


% Run ANOVA to find effect for delay and distance (use all subjects here) 
% Delay 
new_data_mat(:,:,1) = errors; % all subjects included 
between_factor2=[ones(length(pat),1);ones(length(hc),1)*2; ones(length(cog_def),1)*3; ones(length(yhc),1)*4]; 
new_tbl_del = simple_mixed_anova(new_data_mat,between_factor2,{'Delay'},{'Group'});
clear new_data_mat

% Distance
new_data_mat(:,:,1) = ffa;
new_data_mat(:,:,2) = nfa;
new_tbl_dist = simple_mixed_anova(new_data_mat,between_factor2,{'Delay','Distance'},{'Group'});


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%PLOTTING%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot averaged error rates across all groups 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
l1=subplot(1,1,1), hold on 

% Data 
errorbar(1:3,mean(errors,1),std(errors,1)./sqrt(size(errors,1)),'MarkerSize',0.1,'Color',colors.black,'LineWidth',1), hold on 
s1=scatter(1:3,mean(errors,1)); set(s1,'MarkerFaceColor',num2str(colors.black),'MarkerEdgeColor','k');

% Model predictions
shadedErrorBar(1:3,mean(CPs_all_ttypes,1),std(CPs_all_ttypes,[],1)./sqrt(size(CPs_all_ttypes,1)),{'Color',colors.black},1)
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',delays, 'Fontsize',7), ylim([0 0.3]), xlim([0.8 3.2])
xlabel('Delay (s)', 'Fontsize',7), ylabel('Error rate', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',1,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];

w = 0.2; wgap = 0.019; woff = 0.05;
h = .65; hgap = 0.08; hoff = 0.12;
set(l1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
legend
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/error_functions_all.fig']) % Figure 2A
close all

% Plot bar graphs for near false alarms only + model predictions 
figure(1), hold on 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];

l1 = subplot(1,3,1), hold on 
s1 = superbar([mean(nfa_yhc(:,1),1),mean(nfa_hc(:,1),1),mean(nfa_mci(:,1),1)],'E',err_nfa_1,'P',P_nfa_1,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on

set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.5]), ylim([0,0.6])
title('1s Delay', 'Fontsize',7), ylabel('FA rate')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
l2 = subplot(1,3,2), hold on 
superbar([mean(nfa_yhc(:,2),1),mean(nfa_hc(:,2),1),mean(nfa_mci(:,2),1)],'E',err_nfa_3,'P',P_nfa_3,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.5]), ylim([0,0.6])
title('3s Delay', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
l3 = subplot(1,3,3), hold on 
superbar([mean(nfa_yhc(:,3),1),mean(nfa_hc(:,3),1),mean(nfa_mci(:,3),1)],'E',err_nfa_9,'P',P_nfa_9,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.5]), ylim([0,0.6])
title('9s Delay', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];

% Model predictions
subplot(1,3,1)
superbar([1+0.35./2,2+0.35./2,3+0.35./2]+0.35./2./2,[mean(pred_nfa_yhc(:,1),1),mean(pred_nfa_hc(:,1),1),mean(pred_nfa_mci(:,1),1)],'E',err_pred_nfa_1,'BarFaceColor', [colors.gold;colors.teal;colors.lipstick],'PStarFontSize',20,'BarWidth',0.35./2,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.5]), ylim([0,0.6])
title('1s Delay', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
subplot(1,3,2)
superbar([1+0.35./2,2+0.35./2,3+0.35./2]+0.35./2./2,[mean(pred_nfa_yhc(:,2),1),mean(pred_nfa_hc(:,2),1),mean(pred_nfa_mci(:,2),1)],'E',err_pred_nfa_3,'BarFaceColor', [colors.gold;colors.teal;colors.lipstick],'PStarFontSize',20,'BarWidth',0.35./2,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.5]), ylim([0,0.6])
title('3s Delay', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
subplot(1,3,3)
superbar([1+0.35./2,2+0.35./2,3+0.35./2]+0.35./2./2,[mean(pred_nfa_yhc(:,3),1),mean(pred_nfa_hc(:,3),1),mean(pred_nfa_mci(:,3),1)],'E',err_pred_nfa_9,'BarFaceColor', [colors.gold;colors.teal;colors.lipstick],'PStarFontSize',20,'BarWidth',0.35./2,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.5]), ylim([0,0.6])
title('9s Delay', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
legend

w = 0.2; wgap = 0.019; woff = 0.05;
h = .65; hgap = 0.08; hoff = 0.12;
set(l1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
set(l2, 'Position', [woff+wgap+w, hoff, w, h])   % [left bottom width height]
set(l3, 'Position', [woff+(wgap+w)*2, hoff, w, h])   % [left bottom width height]

savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/nfa_comparisons.fig']) %Figure 2B
close all

% Plot differences in fitted model parameters 

% Plot averaged error rates across all groups 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

l1 = subplot(1,4,1)
superbar([mean(noise_yhc(:,1),1),mean(noise_hc(:,1),1),mean(noise_mci(:,1),1)],'E',err_noise,'P',P_noise,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.2]), ylim([0,6])
title('Noise', 'Fontsize',7), ylabel('\sigma_{mem}')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
l2 = subplot(1,4,2)
superbar([mean(lapse_yhc(:,1),1),mean(lapse_hc(:,1),1),mean(lapse_mci(:,1),1)],'E',err_lapse,'P',P_lapse,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.2]), ylim([0,0.065])
title('Lapse Probability', 'Fontsize',7), ylabel('\theta')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
l3 = subplot(1,4,3)
superbar([mean(lapse_noise_yhc(:,1),1),mean(lapse_noise_hc(:,1),1),mean(lapse_noise_mci(:,1),1)],'E',err_lapse_noise,'P',P_lapse_noise,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.2]), ylim([-1,1])
title('Noise + Lapse', 'Fontsize',7), ylabel('\sigma_{mem}+\theta')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
l4 = subplot(1,4,4)
superbar([mean(criterion_yhc(:,1),1),mean(criterion_hc(:,1),1),mean(criterion_mci(:,1),1)],'E',err_criterion,'P',P_criterion,'BarFaceColor', [colors.yellow;colors.sky;colors.rosered],'PStarFontSize',15,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'YHC','OHC','MCI'}, 'Fontsize',7), xlim([-0.2 4.2]), ylim([0,17])
title('Criterion', 'Fontsize',7), ylabel('\delta')

set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];

w = 0.2; wgap = 0.03; woff = 0.05;
h = .65; hgap = 0.08; hoff = 0.12;
set(l1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
set(l2, 'Position', [woff+wgap+w, hoff, w, h])   % [left bottom width height]
set(l3, 'Position', [woff+(wgap+w)*2, hoff, w, h])   % [left bottom width height]
set(l4, 'Position', [woff+(wgap+w)*3, hoff, w, h])   % [left bottom width height]
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/model_param_group_diff.fig']) % Figure 3B
close all
save(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/behavior_w_yhc/stats_behav4paper/behav_stat.mat'],'P_*') % Save all p-values for group-wise comparisons

%% Make error plots for each trial type for Extended Data Figure 2

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

s1 = subplot(1,2,1), hold on % Match trials 

%Data 
errorbar(1:3,mean(misses_hc,1),std(misses_hc,1)./sqrt(size(misses_hc,1)),'MarkerSize',0.1,'Color',colors.sky,'LineWidth',1), hold on 
errorbar(1:3,mean(misses_mci,1),std(misses_mci,1)./sqrt(size(misses_mci,1)),'MarkerSize',0.1,'Color',colors.rosered,'LineWidth',1), hold on 
errorbar(1:3,mean(misses_yhc,1),std(misses_yhc,1)./sqrt(size(misses_yhc,1)),'MarkerSize',0.1,'Color',colors.yellow,'LineWidth',1), hold on 

%Model predictions
shadedErrorBar(1:3,mean(pred_misses_hc,1),std(pred_misses_hc,[],1)./sqrt(size(pred_misses_hc,1)),{'Color',colors.sky},1)
shadedErrorBar(1:3,mean(pred_misses_mci,1),std(pred_misses_mci,[],1)./sqrt(size(pred_misses_mci,1)),{'Color',colors.rosered},1)
shadedErrorBar(1:3,mean(pred_misses_yhc,1),std(pred_misses_yhc,[],1)./sqrt(size(pred_misses_yhc,1)),{'Color',colors.yellow},1)

set(gca,'TickDir','out','XTick',1:3,'XTickLabel',delays, 'Fontsize',7), ylim([0 0.5]), xlim([0.8 3.2])
xlabel('Delay (s)', 'Fontsize',7), ylabel('Miss rate', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];


s2 = subplot(1,2,2), hold on % Non-match trials

%Data 
errorbar(1:3,mean(nfa_hc,1),std(nfa_hc,1)./sqrt(size(nfa_hc,1)),'MarkerSize',0.1,'Color',colors.sky,'LineWidth',1), hold on 
errorbar(1:3,mean(nfa_mci,1),std(nfa_mci,1)./sqrt(size(nfa_mci,1)),'MarkerSize',0.1,'Color',colors.rosered,'LineWidth',1), hold on 
errorbar(1:3,mean(nfa_yhc,1),std(nfa_yhc,1)./sqrt(size(nfa_yhc,1)),'MarkerSize',0.1,'Color',colors.yellow,'LineWidth',1), hold on 

%Model predictions
shadedErrorBar(1:3,mean(pred_nfa_hc,1),std(pred_nfa_hc,[],1)./sqrt(size(pred_nfa_hc,1)),{'Color',colors.sky},1)
shadedErrorBar(1:3,mean(pred_nfa_mci,1),std(pred_nfa_mci,[],1)./sqrt(size(pred_nfa_mci,1)),{'Color',colors.rosered},1)
shadedErrorBar(1:3,mean(pred_nfa_yhc,1),std(pred_nfa_yhc,[],1)./sqrt(size(pred_nfa_yhc,1)),{'Color',colors.yellow},1)

%Data 
errorbar(1:3,mean(ffa_hc,1),std(ffa_hc,1)./sqrt(size(ffa_hc,1)),'MarkerSize',0.1,'Color',colors.teal,'LineWidth',1), hold on 
errorbar(1:3,mean(ffa_mci,1),std(ffa_mci,1)./sqrt(size(ffa_mci,1)),'MarkerSize',0.1,'Color',colors.lipstick,'LineWidth',1), hold on 
errorbar(1:3,mean(ffa_yhc,1),std(ffa_yhc,1)./sqrt(size(ffa_yhc,1)),'MarkerSize',0.1,'Color',colors.gold,'LineWidth',1), hold on 

%Model predictions
shadedErrorBar(1:3,mean(pred_ffa_hc,1),std(pred_ffa_hc,[],1)./sqrt(size(pred_ffa_hc,1)),{'Color',colors.teal},1)
shadedErrorBar(1:3,mean(pred_ffa_mci,1),std(pred_ffa_mci,[],1)./sqrt(size(pred_ffa_mci,1)),{'Color',colors.lipstick},1)
shadedErrorBar(1:3,mean(pred_ffa_yhc,1),std(pred_ffa_yhc,[],1)./sqrt(size(pred_ffa_yhc,1)),{'Color',colors.gold},1)


set(gca,'TickDir','out','XTick',1:3,'XTickLabel',delays, 'Fontsize',7), ylim([0 0.5]), xlim([0.8 3.2])
xlabel('Delay (s)', 'Fontsize',7), ylabel('FA rate', 'Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];

legend
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/supp_error_functions.fig']) % Extended Data Figure 2A-B
close all
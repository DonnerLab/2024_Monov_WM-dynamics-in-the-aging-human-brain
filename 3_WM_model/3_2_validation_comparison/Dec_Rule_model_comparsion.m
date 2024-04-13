% Comparison of BIC scores in algorithmic model of WM including decision
% rule (sigmoid function) with different sets of free parameters 
% --> Extended Data Figure 3-1C
% Gina Monov, UKE, 2022

clear all 
close all
clc
addpath('/Users/ginamonov/Servers/mountpoint1/functions/')

load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
allsubj = horzcat(behav_hc_final,behav_mci_final,behav_cog_def_final); 
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
group2analyze = 'all'; 
if strcmp(group2analyze,'all')
subj = horzcat(allsubj,yhc); % both dataset 1 and dataset 2 
end 

loadpath_mci = '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/';
loadpath_scz = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/';

% Load colors
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat']);

% Set which models to compare 
model_types = {'sM+bound','sM+bound+theta','sM+bound+lambda','sM+sD+bound','sM+bound+theta+lambda','sM+sD+bound+theta+lambda'}; 

% Initialize BIC scores matrix (nsubj*nmodels)
bic_scores = zeros(length(subj),length(model_types));

% Intialization for checking 
ces = zeros(length(subj),length(model_types));
max_cp = zeros(length(subj),length(model_types));
min_cp = zeros(length(subj),length(model_types));

for t = 1:length(model_types) 
   type = model_types(t); %current type/set of free parameters to be tested 
  
   for f = 1:length(subj)
    
    % Load modelling results 
    if sum(ismember(subj{f},allsubj))==1
       fullpath = [loadpath_mci,type{1,1}];
    elseif sum(ismember(subj{f},yhc))==1
       fullpath = [loadpath_scz,type{1,1}];
    end
    
    load([fullpath,filesep,subj{f},'.mat']); 

    bic_scores(f,t) = 2.*ce + length(pm).*log(length(all_behav(:,1)));
    
    % Checking whether there are unexpected modeling results for any subject 
    if ce < 0
       ces(f,t) = 1;
    end 
   
    if max(CPs)>1
       max_cp(f,t) = 1;  
    end 
    
    if min(CPs)<0
       min_cp(f,t) = 1; 
    end 
    
   end
end 

% How many combinations to test?
tst = combntns([1:length(model_types)],2); 

% Finding the best fitting model per participant
winning_model = zeros(size(bic_scores)); 

for f= 1:length(subj)
    winning_model(f,bic_scores(f,:) == min(bic_scores(f,:))) = 1; 
end

% Test significance in difference of BIC score (wilcoxon signed rank test) 
P1 = zeros (length(model_types),length(model_types));
for l = 1:length(model_types) 
    P1(l,l)=nan;  
end 

for l = 1:length(tst(:,1)) 
   p = []; 
   p = signrank(bic_scores(:,tst(l,1)), bic_scores(:,tst(l,2)));
   P1(tst(l,1),tst(l,2)) = p; 
   P1(tst(l,2),tst(l,1)) = p; 
end

% Make plot for comparison: BIC scores relative to Model 2

rel_bic_scores = bic_scores-bic_scores(:,2); 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 7; % figure width
fig_h = 6; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

s1 = subplot(1,1,1), hold on 

yyaxis left
ax = gca; 
ax.YColor = colors.black; 
ylabel('BIC(x) - BIC(2)','FontSize',7)
ylim([min(min(rel_bic_scores))-2,max(max(rel_bic_scores))+2])

hold on

for jj = 1:length(rel_bic_scores(:,1))
    if ismember(subj{jj},yhc)
        use_color = colors.black; 
    elseif ismember(subj{jj},behav_mci_final)
        use_color = colors.rosered; 
    elseif ismember(subj{jj},behav_hc_final)
        use_color = colors.black; 
    elseif ismember(subj{jj},behav_cog_def_final)
        use_color = colors.black; 
    end 
    scatter(1:length(model_types),rel_bic_scores(jj,:),1,'MarkerFaceColor',use_color,'MarkerEdgeColor',use_color), hold on 
    plot(1:length(model_types),rel_bic_scores(jj,:), '-','Color',use_color,'LineWidth',0.2), hold on 

end 
yyaxis right 

ax = gca;
ax.YColor = colors.black; 
bic_MCI = bic_scores(length(behav_hc_final)+1:length(behav_hc_final)+length(behav_mci_final),:);
for jj = 1:length(bic_scores(1,:))
    % add the mean BIC scores for whole group and for MCI patients only
    plot([jj-0.25,jj-0.15],[mean(bic_scores(:,jj)),mean(bic_scores(:,jj))],'-','Color',colors.black,'linewidth',3);
    plot([jj-0.15,jj-0.05],[mean(bic_MCI(:,jj)),mean(bic_MCI(:,jj))],'-','Color',colors.rosered,'linewidth',3);
end

ylim([100,150])
ylabel('BIC score')

xlim([0.5,6.5])
yline(0,'--k'), hold on 
set(gca,'XTick',[1:length(model_types)]);
set(gca,'XTickLabel',{'1','2','3','4','5','6'},'FontSize',7,'TickDir','out','box','off')
xlabel('Model No.','FontSize',8)
cd '/Users/ginamonov/Servers/mountpoint1/final_figures'

savefig(figure(1),['relative_model_comparison_', group2analyze, '.fig']) % Extended Data Figure 3-1C
close all
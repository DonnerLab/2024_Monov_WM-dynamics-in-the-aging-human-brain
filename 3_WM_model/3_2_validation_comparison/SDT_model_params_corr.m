% Script to relate SDT metrics to fitted model parameters 
% --> Extended Data Figure 3-1F
% Gina Monov, UKE, 2022

clear all 
close all

%Load colors
load(['/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat']);

loadpath1 = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/miss-nfa/'; 
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

% loop through Subjects
for s = 1:length(subj)
    miss_rate = []; 
    fa_rate = []; 
    nfa_rate = []; 
    ffa_rate = []; 
    H = []; 
    
    load([loadpath1, subj{s},'.mat'])
    miss_rate_all(s,1) = miss_rate; 
    fa_rate_all(s,1) = fa_rate; 
    H_all(s,1) = H;
    nfa_ffa(s,1) = nfa_rate-ffa_rate;
    
    % Pull model parameters & model predictions for each subjects
    pm = []; 
    CPs = []; 
    if sum(ismember(subj{s},yhc))==1
    fullpath = [modelpath_scz,filesep, subj{s}];
    else fullpath = [modelpath_mci,filesep, subj{s}];
    end
    load([fullpath,'.mat']);
    
    % Pull fitted paramaters 
    noise(s,1) = pm(1); 
    criterion(s,1) = pm(2); 
    lapse(s,1) = pm(3); 
    
    
    
end 

H = []; 
H = zscore(H_all); 
fas = zscore(fa_rate_all); 


d_prime = H-fas; 
c = H+fas; 

tt = 'type'; 
corr_type = 'Pearson';
lw = 'Linewidth';
color = 'color';

% Plot averaged error rates across all groups 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 17.5; % figure width
fig_h = 10.5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

subplot(1,4,1), hold on 
axis square 
scatter(criterion,c,1,'MarkerFaceColor',colors.ink,'MarkerEdgeColor',colors.black), hold on 
[r,pval] = corr(criterion,c,tt,corr_type); 
p = polyfit(criterion,c,1); 
f = polyval(p,criterion); 
plot(criterion,f,lw,1,color,colors.black), hold on 
if pval <10^-4
    title(['r=',num2str(round(r,2)),', p<10^{-4}'],'Fontsize',7,'Fontweight','bold'), hold on 
elseif  pval <0.05
    title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],'Fontsize',7,'Fontweight','bold'), hold on 
else  title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],'Fontsize',7,'Fontweight','regular'), hold on 
end 
set(gca,'fontname','Helvetica','LineWidth',0.5, 'Fontsize',7,'TickDir','out','box','off')
ylabel('c', 'Fontsize',7)
xlabel('Threshold \delta','Fontsize',7)

subplot(1,4,2), hold on 
axis square 
scatter(noise,d_prime,1,'MarkerFaceColor',colors.ink,'MarkerEdgeColor',colors.black), hold on 
[r,pval] = corr(noise,d_prime,tt,corr_type); 
p = polyfit(noise,d_prime,1); 
f = polyval(p,noise); 
plot(noise,f,lw,1,color,colors.black), hold on 
if pval <10^-4
    title(['r=',num2str(round(r,2)),', p<10^{-4}'],'Fontsize',7,'Fontweight','bold'), hold on 
elseif  pval <0.05
    title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],'Fontsize',7,'Fontweight','bold'), hold on 
else  title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],'Fontsize',7,'Fontweight','regular'), hold on 
end 
set(gca,'fontname','Helvetica','LineWidth',0.5, 'Fontsize',7,'TickDir','out','box','off')
ylabel('d prime', 'Fontsize',7)
xlabel('Noise \sigma_{mem}','Fontsize',7)

subplot(1,4,3), hold on 
axis square 
scatter(lapse,d_prime,1,'MarkerFaceColor',colors.ink,'MarkerEdgeColor',colors.black), hold on 
[r,pval] = corr(lapse,d_prime,tt,corr_type); 
p = polyfit(lapse,d_prime,1); 
f = polyval(p,lapse); 
plot(lapse,f,lw,1,color,colors.black), hold on 
if pval <10^-4
    title(['r=',num2str(round(r,2)),', p<10^{-4}'],'Fontsize',7,'Fontweight','bold'), hold on 
elseif  pval <0.05
    title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],'Fontsize',7,'Fontweight','bold'), hold on 
else  title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,4))],'Fontsize',7,'Fontweight','regular'), hold on 
end 
set(gca,'fontname','Helvetica','LineWidth',0.5, 'Fontsize',7,'TickDir','out','box','off')
ylabel('d prime', 'Fontsize',7)
xlabel('Lapse \theta','Fontsize',7)


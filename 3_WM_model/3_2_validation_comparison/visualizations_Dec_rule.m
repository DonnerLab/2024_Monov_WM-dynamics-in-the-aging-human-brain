% Visualizations of algorithmic model using sigmoid decision function: 
% --> Figure 3A 
% Correlations/histograms of fitted parameters (model 2):
% --> Extended Data Figure 3-1D
% Gina Monov, UKE, 2022

clear all
close all

% Parameters for visualization plots 

delays = [1 3 9];       % delay durations
deltas = [0,13.8,27.7,41.5,55.4,69.2,83.1,96.9,110.8,124.6,138.5];  % memorandum-target differences (in degrees; 0 corresponds to match trials)
dt = 0.05; 
nsims = 250000; 

modeltype = {'sM+bound','sM+bound+theta','sM+bound+lambda','sM+sD+bound','sM+bound+theta+lambda','sM+sD+bound+theta+lambda'}; %add all parameter settings

model2plot = modeltype{2}; % which model no. is the winning model according to BIC model comparison
load('/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat'); 
% Test simulations for a single subject and compare to behavior 

%Define subjects/paths 
addpath('/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule_Modelling/DecRule_model_wHR/')
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
loadpath_mci = '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/'; 
loadpath_scz = '/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/'; 

load('/Users/ginamonov/Servers/mountpoint1/meg_analysis/datasets_overview.mat')


load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
allsubj = horzcat(behav_hc_final,behav_mci_final,behav_cog_def_final); 

subj = horzcat(allsubj,yhc); 


% create a plot for visualization purposes 
sigma = 3; 
bound = 11.023; 
% Run simulations
ts = 0:dt:max(delays);
rnds = normrnd(0,sigma*sqrt(dt),nsims,length(ts));  % draw random diffusion increments
rnds = cumsum(rnds,2);  % calculate cumulative sum (i.e. Wiener diffusion trajectories) over time steps (overwriting to save memory)
      
%Probability density functions
          
 for d = 1:length(delays) 
        time_idx = find(ts ==(delays(d)));
        times(d) = time_idx;
        mem_trace(d) = rnds(1,time_idx); 
        pdf_degrees = fitdist(rnds(:,time_idx),'Normal');
        x_values = -30:1:30;
        if d == 1
        y1 = pdf(pdf_degrees,x_values);
        elseif d == 2
        y3 = pdf(pdf_degrees,x_values);  
        elseif d == 3
        y9 = pdf(pdf_degrees,x_values);
        end 
 end 
 
 %extract correct and incorrect trials 
 correct = []; 
 error = []; 
 correct_trace = []; 
 error_trace = []; 
 for iu = 1:length(delays)
     if rnds(1,times(iu))>= -bound && rnds(1,times(iu))<= +bound
         correct = [correct,delays(iu)];
         correct_trace = [correct_trace,mem_trace(iu)]; 
     else error = [error,delays(iu)];
          error_trace = [error_trace,mem_trace(iu)]; 
     end 
 end 
     
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
ls1 = subplot(1,1,1), hold on 

z1=plot(ts+0.5,rnds(1,:),'Color',colors.silver,'LineWidth',1.5), hold on
s0=xline(0.5,'-k', {'Memorandum Off'},'LineWidth',1,'Fontsize',8),hold on
s1=xline(0,'-k', {'Memorandum On'},'LineWidth',1,'Fontsize',8),hold on
s1=xline(1+0.5,'-k', {'1s Delay'},'LineWidth',2,'Fontsize',8,'Color',colors.lavender),hold on
s2=xline(3+0.5,'-k', {'3s Delay'},'LineWidth',2,'Fontsize',8,'Color',colors.magentapurple),hold on
s3=xline(9+0.5,'-k', {'9s Delay'},'LineWidth',2,'Fontsize',8,'Color',colors.plum),hold on

s4=yline(-bound,'--','Color',colors.tiger,'LineWidth',2),hold on
s5=yline(bound,'--','Color',colors.tiger,'LineWidth',2),hold on
s6=yline(0,'--','Color',colors.black,'LineWidth',1),hold on
z2=scatter([0.25],zeros(1,1),80,'MarkerFaceColor',colors.black,'MarkerEdgeColor',colors.black,'Linewidth',1) %add presentation of stimuli

z3=scatter(correct+0.5,correct_trace,80,'MarkerFaceColor',colors.white,'MarkerEdgeColor',colors.silver,'Linewidth',1)
z4=scatter(error+0.5,error_trace,80,'MarkerFaceColor',colors.white,'MarkerEdgeColor',colors.silver,'Linewidth',1)

z5=scatter([1,3,9]+0.5,zeros(1,length(delays)),80,'MarkerFaceColor',colors.white,'MarkerEdgeColor',colors.black,'Linewidth',1) 
uistack(z1,'top')
uistack(z2,'top')
uistack(z3,'top')
uistack(z4,'top')
uistack(z5,'top')
xlabel('time(s)','FontSize',7),ylabel('delta','FontSize',7)


% Plot probability density functions for all simulations 

y1 = y1*10; 
y3 = y3*10; 
y9 = y9*10; 
   b = plot(y1+delays(1)+0.51,x_values,':','Color',colors.silver,'LineWidth',1.5)
   
   uistack(b,'top')
    b = plot(y3+delays(2)+0.51,x_values,':','Color',colors.silver,'LineWidth',1.5)
   
   uistack(b,'top')
    b = plot(y9+delays(3)+0.51,x_values,':','Color',colors.silver,'LineWidth',1.5)
   
   uistack(b,'top')
   
   lgd = legend([z2,z5,z1,b,s4],'Memorandum','Test stimulus','Memory trace','Probability density across simulations','Criterion \delta')
   set(lgd,'Fontsize', 7); 
   ylim([-20,20])
   % Arrange plots 
w = 0.35; wgap = 0.008; ngap = 0.012; cw = 0.008; cl = 0.2; woff = 0.05;
h = 0.35; hgap = 0.08; hoff = 0.12;

set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1.5 3.5 9.5])
set(gca,'XTickLabel',{'0' '0.5' '1.5' '3.5' '9.5'})
set(ls1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
cd '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule_Modelling/Plots/'   
savefig(['match_trial_visualization.fig']) % Figure 3, left 
% --> note that this only served as a rough template for a subsequently heavily edited illustration, which is only for visualization purposes 
close all

% Make plot of STEP FUNCTION and p("different") by delay for the visualizations in the paper
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
ls1 = subplot(1,1,1), hold on 

% Loop through delay durations/mem-target deltas and calculate "same" response probabilities
% Choose values
sigmaM = 3; 
sigmaD = 0;
bound = 11.023;
theta = 0.1;
lambda = [];

deltas_hr = linspace(0,max(deltas)); % deltas in higher resolution for visualization purposes
% Run p(different)
p = run_wd_DecRule(sigmaM,sigmaD,bound,theta,delays,deltas_hr); 

if ~isempty(lambda) % only perform this if lambda is a free parameter
        for d = 1:length(delays)
           p(d,:) = p(d,:).*(1-(1-exp(-lambda*delays(d)))) + (ones(size(p(d,:))).*0.5).*(1-exp(-lambda*delays(d)));
        end
end 


delay_colors(1).type = colors.lavender;
delay_colors(2).type = colors.magentapurple;
delay_colors(3).type = colors.plum;

figure(1), hold on 
c = 'Color';
lw = 'LineWidth';
for d = 1:length(delays)
 eval(['s',num2str(d),'=plot(deltas_hr,p(d,:),c,delay_colors(d).type,lw,1), hold on'])
end 

legend([s1,s2,s3],'1s Delay','3s Delay','9s Delay')
xline(bound,'--','Color',colors.tiger,'LineWidth',1),hold on
ylim([0,1])
xlim([0,2*bound])
xlabel('Memorandum-Test Stimulus Distance','FontSize',7)
ylabel('p(different)','FontSize',7)
  % Arrange plots 
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')
w = 0.09; wgap = 0.008; ngap = 0.012; cw = 0.008; cl = 0.2; woff = 0.05;
h = 0.35; hgap = 0.08; hoff = 0.12;

set(ls1, 'Position', [woff, hoff, w, h])   % [left bottom width height]

cd '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule_Modelling/Plots/'   
savefig(['CPs4paper.fig']) %Figure 3, right 
close all

% Step function 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
ls1 = subplot(1,1,1), hold on 

deltas_hr = linspace(0,max(deltas),100000);
Dec_f = theta + (1-2*theta)./(1+exp(-(deltas_hr-bound)./sigmaD)); 


idx = find(Dec_f == 1-theta); 
idx = idx(1); 
plot(deltas_hr,Dec_f,'Color', colors.black, 'linewidth',1), hold on 
xline(bound,'--','Color',colors.tiger,'LineWidth',1),hold on
text(bound/2,theta+0.05,'\theta','Fontsize',10); 
text(bound+bound/4,1-theta+0.05,'1-\theta','Fontsize',10); 
xlabel('Memorandum-Test Stimulus Distance','FontSize',7)
ylabel('DF','FontSize',7)
ylim([0,1])
xlim([0,2*bound])
  % Arrange plots 
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')
w = 0.09; wgap = 0.008; ngap = 0.012; cw = 0.008; cl = 0.2; woff = 0.05;
h = 0.35; hgap = 0.08; hoff = 0.12;

set(ls1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
cd '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule_Modelling/Plots/'   
savefig(['DF4paper.fig']) % Figure 3, middle 
close all


%Plot correlations of model parameters for winning model  
sigmaM_y = []; %Initialize all possible free parameters
sigmaD_y = [];
bound_y = [];
theta_y = [];
lambda_y = [];
for f = 1:length(subj)
    
    %load parameters
    if sum(ismember(subj{f},allsubj))==1
    fullpath = ['/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/' model2plot];
    elseif sum(ismember(subj{f},yhc))==1
    fullpath = ['/Users/ginamonov/Servers/mountpoint1/SCZ/WM_modeling/Dec_Rule/' model2plot];
    end
    
    load([fullpath,filesep,subj{f},'.mat']);
 
      sigmaM = []; %Initialize all possible free parameters
      sigmaD = [];
      bound = [];
      theta = [];
      lambda = [];
      for r=1:length(pm)
  
         if isempty(sigmaM) & ismember('sM',model2plot)
         sigmaM_y(f) = pm(r);
         sigmaM = pm(r); 
         elseif isempty(sigmaD) & ismember('sD',model2plot)
         sigmaD_y(f) = pm(r);
         sigmaD = pm(r); 
         elseif isempty(bound) & ismember('bound',model2plot)
         bound_y(f) = pm(r);
         bound = pm(r); 
         elseif isempty(theta) & ismember('theta',model2plot)
         theta_y(f) = pm(r);
         theta = pm(r); 
         elseif isempty(lambda) & ismember('lambda',model2plot)
         lambda_y(f) = pm(r);
         lambda = pm(r); 
        end 
     
    end 
    
    allbehav = []; 
  if sum(ismember(subj{f},allsubj))==1
     for u = 1:length(datasets_overview)
        if strcmp(subj{f}, datasets_overview(u).subj) && datasets_overview(u).sess == 1
           allbehav = datasets_overview(u).clean_behavior;
           if strcmp(subj{f},'02')
               allbehav(allbehav(:,11)==1,:) = [];
           end 
        end 
     end
  else 
  load(['/Users/ginamonov/Servers/mountpoint1/SCZ/WM_analysis/clean_WM_data/',subj{f}, 'clean_allbehav.mat']);
  end 
  resp = allbehav(:,5); 
  resp_same(f,1) = length(allbehav(resp==1))/length(allbehav);%Exctract how many same responses
    
 
    
end 


% Get logical which free parameters exist in this model
model_log = zeros(1,5); 
if ismember('sM',model2plot)
model_log(1) = 1; 
end 
if ismember('sD',model2plot)
model_log(2) = 1; 
end
if ismember('bound',model2plot)
model_log(3) = 1; 
end
if ismember('theta',model2plot)
model_log(4) = 1; 
end
if ismember('lambda',model2plot)
model_log(5) = 1; 
end

mycolors(1).type=colors.grey; 
mycolors(2).type=colors.cherry;  
mycolors(3).type=colors.sky; 
mycolors(4).type=colors.magentapink;
mycolors(5).type=colors.orange;
mycolors(6).type=colors.emerald;
for uz=1:length(modeltype)
    if strcmp(model2plot,modeltype{uz})
        color4plot = mycolors(uz).type; 
    end 
end 

  params2corr = [];
  sb_title = {};
% Get parameters
   if model_log(1) == 1
    params2corr(:,end+1) = sigmaM_y;
    sb_title = horzcat(sb_title,'Noise \sigma_{mem}')
   end 
   if model_log(2) == 1
    params2corr(:,end+1) = sigmaD_y;
    sb_title = horzcat(sb_title,'Decision Noise \sigma_{dec}')
   end 
   if model_log(3) == 1
    params2corr(:,end+1) = bound_y;
    sb_title = horzcat(sb_title,'Criterion \delta')
   end 
   if model_log(4) == 1
    params2corr(:,end+1) = theta_y;
    sb_title = horzcat(sb_title,'Lapse \theta')
   end 
   if model_log(5) == 1
    params2corr(:,end+1) = lambda_y;
    sb_title = horzcat(sb_title,'Memory Lapse \lambda')
   end 
   
   num_param = length(model_log(model_log==1)); 
   
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 11.5; % figure width
fig_h = 10.5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
for l = 1:num_param*num_param

 for o = 1:num_param
 if l <= num_param*o
    y_ind = o; 
    break
 end
 end 
    x_ind = mod(l-1,num_param)+1; 
    
    if (l-(num_param*y_ind-num_param+y_ind))<1 %don't plot redundant correlations   
    eval(['s',num2str(l),'=subplot(num_param,num_param,l); hold on'])
 
if l >= num_param*num_param-num_param+1
   xlabel(sb_title{num_param*num_param-l-num_param+x_ind*2},'fontsize',7); 
end 
if mod(l-1,num_param) == 0 
   ylabel(sb_title{y_ind},'fontsize',7); 
end 

if x_ind == y_ind  % histograms 
   histogram(params2corr(:,x_ind),20,'FaceColor',color4plot)
else

    [rho,pval]=corr(params2corr(:,x_ind),params2corr(:,y_ind),'type','Pearson');
    p = polyfit(params2corr(:,x_ind),params2corr(:,y_ind),1); 
    f = polyval(p,params2corr(:,x_ind)); 
    plot(params2corr(:,x_ind),params2corr(:,y_ind),'o','MarkerFaceColor',color4plot,'MarkerEdgeColor',colors.black,'Color',colors.black,'MarkerSize',4), hold on 
    txt = {['r=' num2str(round(rho,2))]; ['p=' num2str(round(pval,4))]};
    if pval <= 0.05
        plot(params2corr(:,x_ind),f,'LineWidth',1,'Color',colors.black)
    text(0.01,max(params2corr(:,y_ind)-max(params2corr(:,y_ind)/10)),txt,'Color',colors.black,'fontsize',8,'fontweight','bold'), hold on 
    else
    text(0.01,max(params2corr(:,y_ind)-max(params2corr(:,y_ind)/10)),txt,'Color',colors.black,'fontsize',8), hold on 
    end %correlations
    end
    end 
 set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')
end 

sgtitle(['Correlations of fitted model parameters, N=',num2str(length(subj))])
cd '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule_Modelling/Plots/'
savefig(figure(1),['correlations',model2plot,'.fig']) %Extended Data Figure 3-1D

close all

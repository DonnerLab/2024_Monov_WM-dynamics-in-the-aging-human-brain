% Script for analyzing result of parameter recovery for algorithmic model of
% WM with decision rule 
% Plots histograms of fitted model parameters (Extended Data Figure 3-1A)
% and the width of the probability distribution for memory noise and
% threshold parameters 
% --> Extended Data Figure 3-1B
% Gina Monov, UKE, 2022

clear all 
close all

load('/mnt/homes/home028/gmonov/functions/colors/colors.mat') % load colors 

load(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/model_instances/median_fits.mat']) %load model instances as used for generating the datsets 

modeltype = {'sM+bound','sM+bound+theta','sM+bound+lambda','sM+sD+bound','sM+bound+theta+lambda','sM+sD+bound+theta+lambda'}; % specify model types 
model2corr = 'sM+bound+theta'; %select for which model to correlate free parameters
% specify colors for each model 
mycolors(1).type=colors.grey; 
mycolors(2).type=colors.cherry;  
mycolors(3).type=colors.sky; 
mycolors(4).type=colors.magentapink;
mycolors(5).type=colors.orange;
mycolors(6).type=colors.emerald;

% cell array of all data sets 
    allsubj={};
    for z = 1:100
        allsubj = horzcat(allsubj,num2str(z)); %create cell with dataset names 1 to 100 
    end 
    
check = zeros(length(allsubj),length(modeltype));     
width = zeros(length(modeltype),length(gen_pm));   
width(width==0) = nan; 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 11.5; % figure width
fig_h = 10.5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


 for model = 1:length(modeltype) % loop over all models 
     
     modelspec = modeltype{model}; 

     model_log = zeros(1,5); % logical of free parameters in specified model


   % Get logical which free parameters exist in this model
   if ismember('sM',modelspec)
     model_log(1) = 1; 
   end 
   if ismember('sD',modelspec)
     model_log(2) = 1; 
   end
   if ismember('bound',modelspec)
   model_log(3) = 1; 
   end
   if ismember('theta',modelspec)
   model_log(4) = 1; 
   end
   if ismember('lambda',modelspec)
   model_log(5) = 1; 
   end
   
   if strcmp(model2corr,modelspec) 
       model_log_sav = model_log; 
   end 
   sigmaM_y = []; %Initialize all possible free parameters
   sigmaD_y = [];
   bound_y = [];
   theta_y = [];
   lambda_y = [];
   loadpath = ['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/',modelspec, '/']; 

    for s = 1:length(allsubj) % loop over all gen data sets 
        
          load([loadpath,filesep,allsubj{s},'.mat']);
          input_data = dat; 
          load([loadpath,filesep,'Fits',filesep,allsubj{s},'.mat']);
          %Check whether input data and model data are identical 
          if sum(input_data(:,2)==dat(:,2)) == length(input_data(:,2))
             check(s,model) = 0; 
          else check(s,model) = 1; 
          end 
 
      sigmaM = []; %Initialize all possible free parameters
      sigmaD = [];
      bound = [];
      theta = [];
      lambda = [];
    
       for r=1:length(pm)
  
         if isempty(sigmaM) & ismember('sM',modeltype{model})
         sigmaM_y(s) = pm(r);
         sigmaM = pm(r);
         elseif isempty(sigmaD) & ismember('sD',modeltype{model})
         sigmaD_y(s) = pm(r);
          sigmaD = pm(r);
         elseif isempty(bound) & ismember('bound',modeltype{model})
         bound_y(s) = pm(r);
         bound = pm(r);
         elseif isempty(theta) & ismember('theta',modeltype{model})
         theta_y(s) = pm(r);
         theta = pm(r);
         elseif isempty(lambda) & ismember('lambda',modeltype{model})
         lambda_y(s) = pm(r);
         lambda = pm(r);
         end 
     
       end 
      
 
    end 
    
    %Concatenate parameters 
    
    params2corr = []; 
    sb_title={};
    if model_log(1) == 1
    params2corr(:,end+1) = sigmaM_y;
    sb_title = horzcat(sb_title,'Memory Noise')
   end 
   if model_log(2) == 1
    params2corr(:,end+1) = sigmaD_y;
    sb_title = horzcat(sb_title,'Decision Noise')
   end 
   if model_log(3) == 1
    params2corr(:,end+1) = bound_y;
    sb_title = horzcat(sb_title,'Criterion')
   end 
   if model_log(4) == 1
    params2corr(:,end+1) = theta_y;
    sb_title = horzcat(sb_title,'Lapse Probability')
   end 
   if model_log(5) == 1
    params2corr(:,end+1) = lambda_y;
    sb_title = horzcat(sb_title,'Memory Lapse Probability')
   end 
   
   if strcmp(model2corr,modelspec) 
       sb_title_sav = sb_title; 
       params2corr_sav = params2corr; 
   end 
   
   
    set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')
    for p = 1:length(model_log(model_log==1)) 
       
           ind = find(model_log==1); 
       
      
           ind = ind(p); 
           eval(['s',num2str(ind+length(model_log).*(model-1)),'=subplot(length(modeltype),length(gen_pm),ind+length(model_log).*(model-1)); hold on'])
          

             xline(gen_pm(ind),'--','LineWidth',2,'Color',colors.black)
             h=histfit(params2corr(:,p),10,'normal'), hold on 
          if strcmp(sb_title{p},'Memory Noise')
               xlim([0,10]); 
          elseif  strcmp(sb_title{p},'Decision Noise')
             xlim([-5,15]);
          elseif  strcmp(sb_title{p},'Criterion')
             xlim([5,15]);
          elseif strcmp(sb_title{p},'Lapse Probability')
              xlim([-0.04,0.06]);
          elseif strcmp(sb_title{p},'Memory Lapse Probability')
             xlim([-0.02,0.04]);
          end 
          
          h(1).FaceColor = mycolors(model).type; 
          h(2).Color = colors.black;
          uistack(h(2),'top')
          hold on
          
          %Get width of the distribution
          
          x_values = min(params2corr(:,p))-min(params2corr(:,p)):0.1:max(params2corr(:,p)).*2;
          pd1 = fitdist(params2corr(:,p),'Normal');
          y=pdf(pd1,x_values);
          halfHeight = (min(y) + max(y)) / 2;
          % Find left edge
          index1 = find(y >= halfHeight, 1, 'first');
          x1 = x_values(index1);
          % Find right edge
          index2 = find(y >= halfHeight, 1, 'last');
          x2 = x_values(index2);
         
          width(model,ind) = x2-x1; 
          
          set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')

        end
       
        
    end 
     cd '/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/Plots/'
     savefig(figure(1),['parameter_recovers_',modelspec,'.fig']) % Extended Data Figure 3-1A
     close all
     
  % Renaming again 
  
  params2corr = params2corr_sav; 
  model_log = model_log_sav; 
  sb_title = sb_title_sav; 
  num_param = length(model_log(model_log==1)); 
   
figure(1), hold on 
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
   histogram(params2corr(:,x_ind),20,'FaceColor',colors.cherry)
else

    [rho,pval]=corr(params2corr(:,x_ind),params2corr(:,y_ind),'type','Spearman');
    p = polyfit(params2corr(:,x_ind),params2corr(:,y_ind),1); 
    f = polyval(p,params2corr(:,x_ind)); 
    plot(params2corr(:,x_ind),params2corr(:,y_ind),'o','MarkerFaceColor',colors.cherry,'MarkerEdgeColor',colors.black,'Color',colors.black,'MarkerSize',5), hold on 
    txt = {['Spearman''s rho=' num2str(round(rho,3))]; ['p=' num2str(round(pval,5))]};
    if pval <= 0.05
    text(min(params2corr(:,x_ind)+min(params2corr(:,x_ind)/20)),max(params2corr(:,y_ind)-max(params2corr(:,y_ind)/10)),txt,'Color',colors.black,'fontsize',7,'fontweight','bold'), hold on 
    else
    text(min(params2corr(:,x_ind)+min(params2corr(:,x_ind)/20)),max(params2corr(:,y_ind)-max(params2corr(:,y_ind)/10)),txt,'Color',colors.grey,'fontsize',7), hold on 
    end %correlations
    end
    end 
    set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')

end 
  cd '/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/Plots/'
  savefig(figure(1),['corr_parameter_recovers_',modelspec,'.fig'])
  close all
              
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 6; % figure width
fig_h = 4; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


              
colors_con = [];
for r = 1:length(modeltype) 
    
    colors_con = [colors_con;mycolors(r).type]; 
end 

width_n = [width(:,1),width(:,3)]; 
 
   for p = 1:2 %only plotting memory noise and criterion  
       
       subplot(1,2,p), hold on 
    
     superbar(width_n(:,p),'BarFaceColor', colors_con,'BarEdgeColor',colors.black,'BarWidth',0.75), hold on
     set(gca,'TickDir','out','XTick',1:length(modeltype),'XTickLabel',[1:length(modeltype)], 'Fontsize',7), xlim([0.5 6.5])
        
      
    if p == 1
        title ('Memory Noise \sigma_{mem}','fontsize',7) 
        ylabel('Width of probability distribution','fontsize',7)
        xlabel('Model No.','fontsize',7)
    elseif p == 2
        title ('Criterion \delta','fontsize',7)
    end 
    set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')

   end 
    
    
   cd '/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/Plots/'
   savefig(figure(1),['width.fig']) % Extended Data Figure 3-1B 
   close all
    

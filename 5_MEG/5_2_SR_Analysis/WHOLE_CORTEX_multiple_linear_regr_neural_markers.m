% Perform multiple linear regression of fitted model parameters, WM task performance
% and cogntive status on neural markers (within- and across-trial variability, decoding precision, power modulation)
% across the entrie cortex as predictors 
% Same analysis perfomed across all subjects and OHC, MCI groups separately 
% Bar plots: Scott Lowe (2022). superbar (https://github.com/scottclowe/superbar), GitHub.
% --> Figure 7
% Gina Monov, UKE, 2023 

clear all
close all
addpath('/mnt/homes/home028/gmonov/functions/')
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
mode = 'whole_cortex'; % Perform regression analysis on 'whole_cortex' or 'sig_dec_prec' (only signifcant decoding precision parcels) 
appl_fdr = 'yes'; % correct p-values for multiple comparisons? 

if strcmp(mode,'whole_cortex') 
Areas = readtable('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/Glasser_labels.csv'); Areas_labels=Areas.Label;

ROIs={}; 
for l = 1:length(Areas_labels)./2
    roi= Areas_labels(l); 
    ROIs{l} = roi{1,1}(3:end); 
end 
% Naming from pymeg is slightly different 
new_ROIs={};
for g = 1:length(ROIs) 
    new_ROIs{g} = replace(ROIs{g},'_','-'); 
end 
end 

if strcmp(mode,'sig_dec_prec')
    load('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/sign_dec_rois.mat') %load previously saved rois that exhibited significant decoding presicion after fdr-correction 
    new_ROIs = dec_prec_rois; 
end

delay = 1; %only run this for delay activity in the first second 

groups = {'all','hc','mci'}; %'mci', 'hc', 'all'

neural_m_titles = {'Decoding precision','Power modulation','Across-trial variability','Within-trial variability'}; 

combine_lapse_noise = 'yes'; % Want to combine lapse+noise into one parameter reflecting behavioural stochasticity? 

%load table with model parameters, cognitive tests, WM performance
load('/mnt/homes/home028/gmonov/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')

% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 12.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

counter = 0; 
for group2plot= 1:length(groups) 
    group = groups{group2plot};
    allsubj ={}; 
    decoding = []; 
    neural_markers = [];
    
        
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/meg_subj.mat']);

if strcmp(group,'mci')
   subj2test = meg_mci_final; 
elseif strcmp(group,'hc')
    subj2test = meg_hc_final; 
elseif strcmp(group,'cog_def')
    subj2test = meg_cog_def_final;
elseif strcmp(group,'all')
    subj2test = meg_all_final;
end 

allsubj = subj2test; 
WMacc = []; 
noise = []; 
lapse = []; 
lapse_noise = []; 
criterion = []; 
PC1 = []; 
trialcount = []; 

% Create matrix with 4 neural markers for all subjects: neural_markers

 loadpath = '/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/'; 
 
for s = 1:length(allsubj) 
    
    % find out trial count for each subject 
    for iii = 1:length(datasets_overview) 
        if datasets_overview(iii).sess == 1 && strcmp(datasets_overview(iii).subj,allsubj{s})
            trialcount(s) = length(datasets_overview(iii).meg_trial_ids_1s(:,1)); 
            
        end 
    end 
    
    decoding = []; %decoding precision
    wt = []; %within trial
    at = []; %across trial
    dpower = []; %delay power
   
       for r = 1:length(new_ROIs) % Loop over all rois instead and load all neural markers for each subj 
              dec_prec  = []; 
          
              load([loadpath, 'dec_prec/' allsubj{s},'_', new_ROIs{r}]);
       
                     decoding(r) = dec_prec;
                
              load([loadpath, 'across_trial/' allsubj{s},'_', new_ROIs{r}]); 
           
                     at(r) = across_trial; 
                  
              load([loadpath, 'within_trial/' allsubj{s},'_', new_ROIs{r}]);
              
                     wt(r) = within_trial; 
                 
              load([loadpath, 'mean_TFR/' allsubj{s},'_', new_ROIs{r}]);
            
                     dpower(r) = mean_TFR; 
                  
         
       end 
                   
           % average results over the glasser parcels 
           neural_markers(s,1) = mean(decoding); 
           neural_markers(s,2) = mean(dpower); 
           neural_markers(s,3) = mean(at); 
           neural_markers(s,4) = mean(wt);
 
     
    % Pull to be predicted variables 
     ii = find(strcmp(cerad_data.ID,allsubj{s})); 
         strcmp(allsubj{s},cerad_data.ID(ii))
         PC1(s,1) = cerad_data.PC1_score_zsubj_all(ii); 
         noise(s,1) = cerad_data.noise(ii); 
         criterion(s,1) = cerad_data.criterion(ii);  
         lapse(s,1) = cerad_data.lapse(ii); 
         WMacc(s,1) = cerad_data.WMacc(ii); 
         
    
    
    load(['/mnt/homes/home028/gmonov/behavioral_analysis/lapse+noise/',allsubj{s},'_lapse_noise.mat']); 
    lapse_noise (s,1) = lapse_noise_single_subj; 
    clear lapse_noise_single_subj
     
     
end 

neural_markers(:,end+1) = trialcount; 

% Normalize the scales of neural markers 
for uii = 1:length(neural_markers(1,:)) 
    neural_markers(:,uii) = zscore(neural_markers(:,uii));
end 


% Check correlation between neural markers 

if strcmp(combine_lapse_noise,'no') 
to_predict = {'WMacc','noise','lapse','criterion','PC1'}; % PC1
to_predict_titles = {'WM task accuracy','Noise \sigma_{mem}','Lapse \theta','Criterion \delta','Cognitive score'};
elseif strcmp(combine_lapse_noise,'yes')
to_predict = {'WMacc','lapse_noise','criterion','PC1'}; % PC1
to_predict_titles = {'WM task accuracy','Noise \sigma_{mem} + Lapse \theta','Criterion \delta','Cognitive score'};
end 
if strcmp(group,'mci')
mycolors(1).type = colors.rosered; 
mycolors(2).type = colors.lipstick; 
elseif strcmp(group,'hc')
mycolors(1).type = colors.sky; 
mycolors(2).type = colors.teal;
elseif strcmp(group,'cog_def')
mycolors(1).type = colors.purple; 
elseif strcmp(group,'all')
mycolors(1).type = colors.ink; 
mycolors(2).type = colors.black; 
end

for p = 1:length(to_predict) 
    counter = counter +1; 
    topred = to_predict{p}; 
    eval(['mdl = fitlm(neural_markers,',topred,');']) 
    
    eval(['s',num2str(counter),'=subplot(4,7,counter), hold on']) 
 
    coefficients = mdl.Coefficients.Estimate(2:end-1); % exclude intercept & trial count 
    ebars = mdl.Coefficients.SE(2:end-1); % exclude intercept & trial count 
    p_values = mdl.Coefficients.pValue(2:end-1); % exclude intercept & trial count 
    R_squared = mdl.Rsquared.Ordinary; % ordinary R-squared
    R_squared_adj = mdl.Rsquared.Adjusted; % adjusted R-squared 
    p_ftest = mdl.ModelFitVsNullModel.Pvalue; % P-value of F-test 
    if strcmp(appl_fdr,'yes') %apply fdr-correction if defined at the top of script
        fdr_corr = fdr(p_values,0.05); 
        for tr = 1:length(fdr_corr)
            if fdr_corr(tr)==1 
                p_values(tr) = p_values(tr);
            elseif fdr_corr(tr)==0
                p_values(tr) = 1; % set p_value to 1 to make it n.s. if it does not survive correction
            end 
        end 
    end 

     b = superbar(1:4,coefficients,'E',ebars,'P',p_values,'BarFaceColor', mycolors(2).type,'PStarFontSize',10,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
    
     if  group2plot == length(groups)
        
     set(gca,'TickDir','out','XTick',1:4,'XTickLabel',neural_m_titles, 'Fontsize',7),
     else
      set(gca,'TickDir','out','XTick',1:4,'XTickLabel',[{}], 'Fontsize',7),
     end 
     if group2plot == round(length(groups)./2,0) && p == 1
        ylabel('Coefficient estimates')
     end 
    
    if strcmp(to_predict{p},'WMacc')
    ylim([-0.12 0.1])
    elseif strcmp(to_predict{p},'noise')
        ylim([-2.3 2.3])
    elseif strcmp(to_predict{p},'lapse')
        ylim([-0.05 0.05])   
    elseif strcmp(to_predict{p},'lapse_noise')
        ylim([-1.6 2])   
    elseif strcmp(to_predict{p},'criterion')
        ylim([-6 8])
    elseif strcmp(to_predict{p},'PC1')
        ylim([-1.7 2])     
    end 

       % Add R-squared values to plots 
       text(3.5,max(ylim)-max(ylim)/4.5,{['R^{2}=' num2str(round(R_squared,2))];['p(F-test)=' num2str(round(p_ftest,2))]},'Fontsize',5)


    if group2plot == 1
    title (to_predict_titles{p})
    end 


set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off')
end 
legend
hold on 
end 

% Arrange plots 
w = 0.1; wgap = 0.019; woff = 0.05;
h = 0.16; hgap = 0.05; hoff = 0.01;
if strcmp(combine_lapse_noise,'yes')
set(s1, 'Position', [woff, hoff+(3*(hgap+h)), w, h])   % [left bottom width height]
set(s2, 'Position', [woff+(w+wgap),hoff+(3*(hgap+h)), w, h])
set(s3, 'Position', [woff+(w+wgap)*2, hoff+(3*(hgap+h)), w, h])
set(s4, 'Position', [woff+(w+wgap)*3, hoff+(3*(hgap+h)), w, h])
set(s5, 'Position', [woff, hoff+(2*(hgap+h)), w, h])
set(s6, 'Position', [woff+(w+wgap), hoff+(2*(hgap+h)), w, h])
set(s7, 'Position', [woff+(w+wgap)*2, hoff+(2*(hgap+h)), w, h])
set(s8, 'Position', [woff+(w+wgap)*3, hoff+(2*(hgap+h)), w, h])
set(s9, 'Position', [woff, hoff+(1*(hgap+h)), w, h])
set(s10, 'Position', [woff+(w+wgap), hoff+(1*(hgap+h)), w, h])
set(s11, 'Position', [woff+(w+wgap)*2, hoff+(1*(hgap+h)), w, h])
set(s12, 'Position', [woff+(w+wgap)*3, hoff+(1*(hgap+h)), w, h])

elseif strcmp(combine_lapse_noise,'no')
set(s1, 'Position', [woff, hoff+(3*(hgap+h)), w, h])   % [left bottom width height]
set(s2, 'Position', [woff+(w+wgap),hoff+(3*(hgap+h)), w, h])
set(s3, 'Position', [woff+(w+wgap)*2, hoff+(3*(hgap+h)), w, h])
set(s4, 'Position', [woff+(w+wgap)*3, hoff+(3*(hgap+h)), w, h])
set(s5, 'Position', [woff+(w+wgap)*4, hoff+(3*(hgap+h)), w, h])
set(s6, 'Position', [woff, hoff+(2*(hgap+h)), w, h])
set(s7, 'Position', [woff+(w+wgap), hoff+(2*(hgap+h)), w, h])
set(s8, 'Position', [woff+(w+wgap)*2, hoff+(2*(hgap+h)), w, h])
set(s9, 'Position', [woff+(w+wgap)*3, hoff+(2*(hgap+h)), w, h])
set(s10, 'Position', [woff+(w+wgap)*4, hoff+(2*(hgap+h)), w, h])
set(s11, 'Position', [woff, hoff+(1*(hgap+h)), w, h])
set(s12, 'Position', [woff+(w+wgap), hoff+(1*(hgap+h)), w, h])
set(s13, 'Position', [woff+(w+wgap)*2, hoff+(1*(hgap+h)), w, h])
set(s14, 'Position', [woff+(w+wgap)*3, hoff+(1*(hgap+h)), w, h])
set(s15, 'Position', [woff+(w+wgap)*4, hoff+(1*(hgap+h)), w, h])
end 

cd '/mnt/homes/home028/gmonov/final_figures/Figure_7/'  


savefig([mode,'_regression_combine_lapse_noise_',combine_lapse_noise,'_.fig']) % Figure 7
close all

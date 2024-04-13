% Script performs PCA on CERAD test battery data, saves the result which will be pulled for later
% analyses 
% Produces the following panels from the manuscript: 
% --> Figure 1B, Extended Data Figure 1-1A-B
% Gina Monov, UKE, 2022

clear all
close all

% Load colors 
colorpath = '/Users/ginamonov/Servers/mountpoint1/functions/colors/'
load([colorpath,'colors']);
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
behav_mci_final = behav_mci_final(~ismember(behav_mci_final,'02'));  

addpath '/Users/ginamonov/Servers/mountpoint1/functions'
addpath '/Users/ginamonov/Servers/mountpoint1/psych_data/'
WM_path = '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'
model_path = '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/sM+bound+theta/'
save_path = '/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/'
cd '/Users/ginamonov/Servers/mountpoint1/psych_data'
cerad_data = readtable('Cerad_T1.xlsx');

% % Neuropsychological data in EXCEL file 
% % subtests / performance scores: 
% cerad_data.VF % Semantic Fluency 
% cerad_data.BNT % Boston Naming Test
% cerad_data.MMST % Mini Mental State Examination
% cerad_data.WL % Word List Learning 
% cerad_data.WLAB % Word List Recall 
% cerad_data.WLIN % Intrusions
% cerad_data.Sav_WL % Word list savings 
% cerad_data.discrim % Word List Recognition (calculation below) --> .WLWE_Ja & .WLWE_Nein
% cerad_data.VK_1 % Constructional Praxis
% cerad_data.VK_2 % Constructional Praxis Recall 
% cerad_data.Sav_Fig % Savings of drawn figures (constructional praxis) 
% cerad_data.Uhr % Clock drawing test
% cerad_data.TMTA % Trail Making Test A
% cerad_data.TMTB % Trail Making Test B 
% cerad_data.TMT_B_A % Trail Making Test B/A
% cerad_data.PF % Phonemic Fluency

% fixing problem with ID cells
cerad_data.ID{1} = '01';
cerad_data.ID{2} = '02';
cerad_data.ID{3} = '03';
cerad_data.ID{4} = '04';
cerad_data.ID{5} = '05';
cerad_data.ID{6} = '06';
cerad_data.ID{7} = '07';
cerad_data.ID{8} = '08';
cerad_data.ID{9} = '09';
cerad_data.ID{10} = '10';
cerad_data.ID{11} = '11';
cerad_data.ID{12} = '12';
cerad_data.ID{13} = '13';
cerad_data.ID{14} = '14';
cerad_data.ID{15} = '15';
cerad_data.ID{16} = '16';
cerad_data.ID{17} = '17';
cerad_data.ID{18} = '18';
cerad_data.ID{19} = '19';
cerad_data.ID{20} = '20';
cerad_data.ID{21} = '21';
cerad_data.ID{22} = '22';
cerad_data.ID{23} = '23';
cerad_data.ID{24} = '24';
cerad_data.ID{25} = '25';
cerad_data.ID{26} = '26';
cerad_data.ID{27} = '27';
cerad_data.ID{28} = '28';
cerad_data.ID{29} = '29';
cerad_data.ID{30} = '30';
cerad_data.ID{31} = '31';
cerad_data.ID{32} = '32';
cerad_data.ID{33} = '33';
cerad_data.ID{34} = '34';
cerad_data.ID{35} = '35';
cerad_data.ID{36} = '36';
cerad_data.ID{37} = '37';
cerad_data.ID{38} = '38';
cerad_data.ID{39} = '39';
cerad_data.ID{40} = '40';
cerad_data.ID{41} = '41';
cerad_data.ID{42} = '42';
cerad_data.ID{43} = '43';
cerad_data.ID{44} = '44';
cerad_data.ID{45} = '45';
cerad_data.ID{46} = '46';
cerad_data.ID{47} = '47';
cerad_data.ID{48} = '48';
cerad_data.ID{49} = '49';
cerad_data.ID{50} = '50';
cerad_data.ID{51} = '51';
cerad_data.ID{52} = '52';
cerad_data.ID{53} = '53';
cerad_data.ID{54} = '54';

idx = zeros(height(cerad_data),1); 

for r = 1:height(cerad_data) 
     if ~ismember(cerad_data.ID(r),behav_all_final) || strcmp(cerad_data.ID(r),'02')
         idx(r) = 1; 
     end 
end 
idx = find(idx == 1); 

 % calculation of BNT total
 for t=1:height(cerad_data)
    cerad_data.BNT(t) = cerad_data.BNT_1(t)+cerad_data.BNT_2(t)+cerad_data.BNT_3(t);
 end 
 
% Remove all subjects not included in final analysis 
cerad_data(idx,:)=[]; 

% Calculation of discriminability (i.e. word list recogntion, which should go to the PCA instead of WLWE_Ja and WLWE_Nein) (and the total
% score+total z-score) 
for t=1:height(cerad_data)
    % adding WM performance and fitted model parameters 
    
    ID=cerad_data.ID(t);
    cd '/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/clean_allbehav_data/'
    if isfolder(ID{1,1}) == 1
       load([ID{1,1},filesep,'S1',filesep,ID{1,1},'_1_clean_allbehav.mat']);
        % exclude one block due to performance issues
       if strcmp(ID{1,1},'02') % 02 first block with performance <50%
         blockcount = allbehav(:,11);
         allbehav(blockcount == 1,:)=[]; 
       end 

       A = allbehav(:,6); % accuracy 
       cerad_data.WMacc(t)=nansum(A)./length(A);
    else cerad_data.WMacc(t) = nan; 
    end 
    
    cd '/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule/sM+bound+theta/'
    if isfile([ID{1,1},'.mat']) == 1
        load([ID{1,1},'.mat']);
        cerad_data.noise(t) = pm(1); 
        cerad_data.criterion(t) = pm(2);
        cerad_data.lapse(t) = pm(3); 
    else
        cerad_data.noise(t) = nan;
        cerad_data.criterion(t) = nan;
        cerad_data.lapse(t) = nan;
    end 
    % calculatimg discriminability for each participant
    cerad_data.discrim(t)=(1-(((10-cerad_data.WLWE_Ja(t)+(10-cerad_data.WLWE_Nein(t)))/20)))*100;
    
    % discriminability as in d' like in chandler total score
    cerad_data.chandler_discrim(t) = cerad_data.WLWE_Ja(t)-(10-cerad_data.WLWE_Nein(t)); 
    
    % calculating Chandler et al. total score 
    if cerad_data.VF(t)>24 % ceiling for Semantic Fluency
        cerad_data.totalscore(t)=24+cerad_data.BNT_1(t)+cerad_data.BNT_2(t)+cerad_data.BNT_3(t)+cerad_data.WL(t)+cerad_data.VK_1(t)+cerad_data.WLAB(t)+cerad_data.chandler_discrim(t); % ceiling for Semantic Fluency
    else cerad_data.totalscore(t)=cerad_data.VF(t)+cerad_data.BNT_1(t)+cerad_data.BNT_2(t)+cerad_data.BNT_3(t)+cerad_data.WL(t)+cerad_data.VK_1(t)+cerad_data.WLAB(t)+cerad_data.chandler_discrim(t);
    end
   
end

   
   % Z-scoring the data of subjects in the dataset 
    cerad_data.VFz = zscore(cerad_data.VF); 
    cerad_data.PFz = zscore(cerad_data.PF); 
    cerad_data.BNTz = zscore(cerad_data.BNT);
    cerad_data.WLz = zscore(cerad_data.WL); 
    cerad_data.WLABz = zscore(cerad_data.WLAB); 
    cerad_data.discrimz = zscore(cerad_data.discrim); 
    cerad_data.chandler_discrimz = zscore(cerad_data.chandler_discrim); 
    cerad_data.VK_1z = zscore(cerad_data.VK_1); 
    cerad_data.VK_2z = zscore(cerad_data.VK_2);
    cerad_data.MMSEz = zscore(cerad_data.MMST); 
    cerad_data.TMTAz = -zscore(cerad_data.TMTA); %Flip sign here so that all postive values indicate high performance
    cerad_data.TMTBz = -zscore(cerad_data.TMTB); %Flip sign here so that all postive values indicate high performance
    
%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCA ON ALL TESTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%


% concatenate metrics for PCA
forpca = [cerad_data.VFz cerad_data.BNTz cerad_data.WLz cerad_data.WLABz cerad_data.discrimz cerad_data.VK_1z ...  
          cerad_data.VK_2z cerad_data.TMTAz cerad_data.TMTBz cerad_data.PFz];

pcalabels = {'Semantic Fluency','BNT','WL Learning','WL Recall','WL Recognition','Constructional Praxis',...
             'CP Recall','TMTA','TMTB', 'Phonemic Fluency'};
      
% set any missing values to nan 
forpca(forpca==-99) = nan;

% run PCA
[coeff,score,latent,~,explained] = pca(forpca);

%% Run further analyses and make plots for paper 

% ROC curve vor diagnostic properties of cognitive integrity scores 
cerad_data.PC1_score_zsubj_all = score(:,1); 
cerad_data.PC2_score_zsubj_all = score(:,2);
cerad_data.PC3_score_zsubj_all = score(:,3);

for s = 1:length(behav_hc_final)
    for l = 1:height(cerad_data) 
        if strcmp(cerad_data.ID(l),behav_hc_final{s})
           totalscore_HC(s) = cerad_data.totalscore(l); 
           PC1_scores_zsubj_all_HC(s) = cerad_data.PC1_score_zsubj_all(l); 
        end 
    end 
    
end 


for s = 1:length(behav_mci_final)
    for l = 1:height(cerad_data) 
        if strcmp(cerad_data.ID(l),behav_mci_final{s})
           totalscore_pat(s) = cerad_data.totalscore(l); 
           PC1_scores_zsubj_all_pat(s) = cerad_data.PC1_score_zsubj_all(l); 
        end 
    end 
    
end 


% ROC curves for diagnosis 

[xtotal,ytotal,~,auctotal] = perfcurve([ones(length(totalscore_HC),1); zeros(length(totalscore_pat),1)],[totalscore_HC'; totalscore_pat'],1); 
[xPC1zsubjall,yPC1zsubjall,~,aucPC1zsubjall] =  perfcurve([ones(length(PC1_scores_zsubj_all_HC),1); zeros(length(PC1_scores_zsubj_all_pat),1)],[PC1_scores_zsubj_all_HC'; PC1_scores_zsubj_all_pat'],1); 


h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

s1 = subplot(1,1,1), hold on 
axis square
subplot(1,1,1)
plot(xPC1zsubjall,yPC1zsubjall,'Color', colors.grey,'LineWidth',1), hold on 
x1 = linspace (0,1); 
y1 = linspace (0,1); 
plot(x1,y1,'--','color',colors.black,'LineWidth',1)
legend('PC1 score','Fontsize',7)
txt = ['AUROC = ' num2str(round(aucPC1zsubjall,4))]; 
text(0.6,0.3-0.1,txt,'Color',colors.grey,'Fontsize',7); 
xlabel('1-Specificity','Fontsize',7) 
ylabel('Sensitivity','Fontsize',7)
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out')
w = 0.31; wgap = 0.019; woff = 0.05;
h = .4; hgap = 0.08; hoff = 0.12;
set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]

savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/PCA_ROC_paper.fig']) % Figure 1B right 
close all


% Initialize figure 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

s1 = subplot(1,1,1), hold on 
x = categorical ({'1','2','3','4','5','6','7','8','9','10'})
x=reordercats(x,{'1','2','3','4','5','6','7','8','9','10'})
b=bar(x,explained,'FaceColor','flat','EdgeColor',colors.black), hold on
ylabel('Variance explained (%)','Fontsize',15), hold on
xlabel('Principal Components','Fontsize',15), hold on

ylim([0,(max(explained)+0.7)])
b.CData(1,:) = colors.grey;

for o=2:10
b.CData(o,:) = colors.white;

end 
for i1=1:numel(explained)
    if i1 == 1
    text(x(i1),explained(i1),num2str(explained(i1),'%0.2f'),'Fontsize',7,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    end 
end
set(gca,'FontName','Helvetica','LineWidth',1,'FontSize',7,'TickDir','out')
w = 0.31; wgap = 0.019; woff = 0.05;
h = 0.42; hgap = 0.08; hoff = 0.12;
set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/PCA_var_exp_paper.fig']) % Figure 1B left 
close all


h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

s1 = subplot(1,1,1), hold off

heatmap('PC1 Loadings',pcalabels, round(coeff(:,1),2),'Fontsize',7), 
colormap(flipud(gray))

set(gca,'FontName','Helvetica','FontSize',7)
w = 0.19; wgap = 0.019; woff = 0.05;
h = .17; hgap = 0.08; hoff = 0.12;

savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/PCA_loadings_heatmap_paper.fig']) % Figure 1B left, inset 
close all

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
s1 = subplot(1,1,1), hold on 
axis square
scatter(cerad_data.totalscore,score(:,1),15,'MarkerFaceColor',colors.black, 'MarkerEdgeColor',colors.black); 
[r,pval]=corr(cerad_data.totalscore, score(:,1),'Rows','pairwise');
xlim([57,100])
p = polyfit(cerad_data.totalscore,score(:,1),1); 
f = polyval(p,cerad_data.totalscore); 
l1 = plot(cerad_data.totalscore,f,'Linewidth',1,'Color',colors.black)
if pval >= 10^-4
txt = {['r=' num2str(round(r,4))], ['p=' num2str(round(pval,4))]}; 
text(80,-3,txt,'Fontsize',7), hold on
else txt = {['r=' num2str(round(r,4))], ['p<10^{-4}']}; 
text(80,-3,txt,'Fontsize',7), hold on
end 
ylabel('PC1 Score','Color',colors.black,'Fontsize',7); 
xlabel('CERAD Total Score','Color',colors.black,'Fontsize',7); 

set(gca,'FontName','Helvetica','FontSize',7,'TickDir','out')
w = 0.31; wgap = 0.019; woff = 0.05;
h = 0.42; hgap = 0.08; hoff = 0.12;
set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]


savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/PCA_correlations_paper.fig']) % Figure 1B, middle 
close all

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 16.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
s1 = subplot(1,1,1), hold on 
axis square
scatter(cerad_data.MMST,score(:,1),15,'MarkerFaceColor',colors.black, 'MarkerEdgeColor',colors.black); 
[r,pval]=corr(cerad_data.MMST, score(:,1),'Rows','pairwise');
xlim([25,30])
p = polyfit(cerad_data.MMST,score(:,1),1); 
f = polyval(p,cerad_data.MMST); 
l1 = plot(cerad_data.MMST,f,'Linewidth',1,'Color',colors.black)
if pval >= 10^-4
txt = {['r=' num2str(round(r,4))], ['p=' num2str(round(pval,4))]}; 
text(24.1,4,txt,'Fontsize',7), hold on
else txt = {['r=' num2str(round(r,4))], ['p<10^{-4}']}; 
text(24.1,4,txt,'Fontsize',7), hold on
end 
ylabel('PC1 Score','Color',colors.black,'Fontsize',7); 
xlabel('MMSE','Color',colors.black,'Fontsize',7); 

set(gca,'FontName','Helvetica','FontSize',7,'TickDir','out')
w = 0.31; wgap = 0.019; woff = 0.05;
h = 0.42; hgap = 0.08; hoff = 0.12;
set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]


savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/PCA_correlations_MMSE_paper.fig']) % Extended Data Figure 1-1B 
close all

 

save([save_path,'cerad_table_PCA.mat'],'cerad_data') % Save results

% Subtest performance for each group (OHC, MCI, UNC) 
for s = 1:length(behav_mci_final)
    for t = 1:height(cerad_data) 
        if strcmp(cerad_data.ID(t),behav_mci_final{s})
           
          VF_z_mci(s,1) = cerad_data.VFz(t);
          BNT_z_mci(s,1) = cerad_data.BNTz(t);
          TMTA_z_mci(s,1) = cerad_data.TMTAz(t);
          TMTB_z_mci(s,1) = cerad_data.TMTBz(t);
          MMSE_z_mci (s,1)= cerad_data.MMSEz(t);
          WL_z_mci(s,1) = cerad_data.WLz(t);
          WLAB_z_mci(s,1) = cerad_data.WLABz(t);
          VK1_z_mci(s,1) = cerad_data.VK_1z(t); 
          Discr_z_mci(s,1) =  cerad_data.discrimz(t); 
          PF_z_mci(s,1) = cerad_data.PFz(t);
          VK2_z_mci(s,1) = cerad_data.VK_2z(t);

        end 
    end 
    
end 

for s = 1:length(behav_hc_final)
    for t = 1:height(cerad_data) 
        if strcmp(cerad_data.ID(t),behav_hc_final{s})
           
          VF_z_hc(s,1) = cerad_data.VFz(t);
          BNT_z_hc(s,1) = cerad_data.BNTz(t);
          TMTA_z_hc(s,1) = cerad_data.TMTAz(t);
          TMTB_z_hc(s,1) = cerad_data.TMTBz(t);
          MMSE_z_hc (s,1)= cerad_data.MMSEz(t);
          WL_z_hc(s,1) = cerad_data.WLz(t);
          WLAB_z_hc(s,1) = cerad_data.WLABz(t);
          VK1_z_hc(s,1) = cerad_data.VK_1z(t); 
          Discr_z_hc(s,1) =  cerad_data.discrimz(t); 
          PF_z_hc(s,1) = cerad_data.PFz(t);
          VK2_z_hc(s,1) = cerad_data.VK_2z(t);

        end 
    end 
    
end 

for s = 1:length(behav_cog_def_final)
    for t = 1:height(cerad_data) 
        if strcmp(cerad_data.ID(t),behav_cog_def_final{s})
           
          VF_z_cog_def(s,1) = cerad_data.VFz(t);
          BNT_z_cog_def(s,1) = cerad_data.BNTz(t);
          TMTA_z_cog_def(s,1) = cerad_data.TMTAz(t);
          TMTB_z_cog_def(s,1) = cerad_data.TMTBz(t);
          MMSE_z_cog_def (s,1)= cerad_data.MMSEz(t);
          WL_z_cog_def(s,1) = cerad_data.WLz(t);
          WLAB_z_cog_def(s,1) = cerad_data.WLABz(t);
          VK1_z_cog_def(s,1) = cerad_data.VK_1z(t); 
          Discr_z_cog_def(s,1) =  cerad_data.discrimz(t); 
          PF_z_cog_def(s,1) = cerad_data.PFz(t);
          VK2_z_cog_def(s,1) = cerad_data.VK_2z(t);

        end 
    end 
    
end 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 12.5; % figure width
fig_h = 4.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

titles = {'Semantic Fluency','BNT','WL Learning','WL Recall','WL Recognition','Constructional Praxis',...
              'CP Recall','TMTA','TMTB', 'Phonemic Fluency','MMSE'};
          
          
subplot(1,3,1), hold on           
y = [mean(VF_z_hc),mean(BNT_z_hc),mean(WL_z_hc),mean(WLAB_z_hc),mean(Discr_z_hc),...
    mean(VK1_z_hc),mean(VK2_z_hc),mean(TMTA_z_hc),mean(TMTB_z_hc),...
    mean(PF_z_hc),mean(MMSE_z_hc)];


err = [std(VF_z_hc,1)./sqrt(size(VF_z_hc,1));std(BNT_z_hc,1)./sqrt(size(BNT_z_hc,1));std(WL_z_hc,1)./sqrt(size(WL_z_hc,1));std(WLAB_z_hc,1)./sqrt(size(WLAB_z_hc,1));std(Discr_z_hc,1)./sqrt(size(Discr_z_hc,1));...
    std(VK1_z_hc,1)./sqrt(size(VK1_z_hc,1));std(VK2_z_hc,1)./sqrt(size(VK2_z_hc,1));std(TMTA_z_hc,1)./sqrt(size(TMTA_z_hc,1));std(TMTB_z_hc,1)./sqrt(size(TMTB_z_hc,1));...
    std(PF_z_hc,1)./sqrt(size(PF_z_hc,1));std(MMSE_z_hc,1)./sqrt(size(MMSE_z_hc,1))];


s1 = superbar(y,'E',err,'BarFaceColor', [colors.sky;colors.sky;colors.sky],'BarWidth',0.7, 'ErrorbarLineWidth',0.3), hold on


set(gca,'TickDir','out','XTick',1:length(titles),'XTickLabel',titles, 'Fontsize',7), xlim([-0.2 length(titles)+1.2]), ylim([-1.5,1.5])
title('OHC', 'Fontsize',7,'fontweight','normal'), ylabel('performance z-score')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0]; 

subplot(1,3,2), hold on 

y = [mean(VF_z_mci),mean(BNT_z_mci),mean(WL_z_mci),mean(WLAB_z_mci),mean(Discr_z_mci),...
    mean(VK1_z_mci),mean(VK2_z_mci),mean(TMTA_z_mci),mean(TMTB_z_mci),...
    mean(PF_z_mci),mean(MMSE_z_mci)];


err = [std(VF_z_mci,1)./sqrt(size(VF_z_mci,1));std(BNT_z_mci,1)./sqrt(size(BNT_z_mci,1));std(WL_z_mci,1)./sqrt(size(WL_z_mci,1));std(WLAB_z_mci,1)./sqrt(size(WLAB_z_mci,1));std(Discr_z_mci,1)./sqrt(size(Discr_z_mci,1));...
    std(VK1_z_mci,1)./sqrt(size(VK1_z_mci,1));std(VK2_z_mci,1)./sqrt(size(VK2_z_mci,1));std(TMTA_z_mci,1)./sqrt(size(TMTA_z_mci,1));std(TMTB_z_mci,1)./sqrt(size(TMTB_z_mci,1));...
    std(PF_z_mci,1)./sqrt(size(PF_z_mci,1));std(MMSE_z_mci,1)./sqrt(size(MMSE_z_mci,1))];



s2 = superbar(y,'E',err,'BarFaceColor', [colors.rosered;colors.rosered;colors.rosered],'BarWidth',0.7, 'ErrorbarLineWidth',0.3), hold on


set(gca,'TickDir','out','XTick',1:length(titles),'XTickLabel',titles, 'Fontsize',7), xlim([-0.2 length(titles)+1.2]), ylim([-1.5,1.5])
title('MCI', 'Fontsize',7,'fontweight','normal')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];


         
subplot(1,3,3), hold on   

y = [mean(VF_z_cog_def),mean(BNT_z_cog_def),mean(WL_z_cog_def),mean(WLAB_z_cog_def),mean(Discr_z_cog_def),...
    mean(VK1_z_cog_def),mean(VK2_z_cog_def),mean(TMTA_z_cog_def),mean(TMTB_z_cog_def),...
    mean(PF_z_cog_def),mean(MMSE_z_cog_def)];


err = [std(VF_z_cog_def,1)./sqrt(size(VF_z_cog_def,1));std(BNT_z_cog_def,1)./sqrt(size(BNT_z_cog_def,1));std(WL_z_cog_def,1)./sqrt(size(WL_z_cog_def,1));std(WLAB_z_cog_def,1)./sqrt(size(WLAB_z_cog_def,1));std(Discr_z_cog_def,1)./sqrt(size(Discr_z_cog_def,1));...
    std(VK1_z_cog_def,1)./sqrt(size(VK1_z_cog_def,1));std(VK2_z_cog_def,1)./sqrt(size(VK2_z_cog_def,1));std(TMTA_z_cog_def,1)./sqrt(size(TMTA_z_cog_def,1));std(TMTB_z_cog_def,1)./sqrt(size(TMTB_z_cog_def,1));...
    std(PF_z_cog_def,1)./sqrt(size(PF_z_cog_def,1));std(MMSE_z_cog_def,1)./sqrt(size(MMSE_z_cog_def,1))];

s3 = superbar(y,'E',err,'BarFaceColor', [colors.tortilla;colors.tortilla;colors.tortilla],'BarWidth',0.7, 'ErrorbarLineWidth',0.3), hold on


set(gca,'TickDir','out','XTick',1:length(titles),'XTickLabel',titles, 'Fontsize',7), xlim([-0.2 length(titles)+1.2]), ylim([-1.5,1.5])
title('UNC', 'Fontsize',7,'fontweight','normal')
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7)
ax = gca; ax.XColor=[0 0 0];
    
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/CERAD_scores_per_group_paper.fig']) % Extended Data Figure 1-1A
close all
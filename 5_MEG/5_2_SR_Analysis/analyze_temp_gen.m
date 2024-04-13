% Run analyses on temporal generalization results 
% Depending on inputs (area and shift between train and test) 
% script produces Figure 8C-E & Extended Data Figure 8-1A-B 
% Permutation test: https://github.com/rudyvdbrink/Tools/tree/master/statistics
% cluster-based permutation test: Edden M. Gerber (2023). permutest, MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/71737-permutest
% Gina Monov, UKE, 2024 

clear all
close all
addpath '/Users/ginamonov/Servers/mountpoint1/functions/'
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/meg_subj.mat']);
load('/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat')
load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')
run_dynamic_coding = 'yes'; % Analyze dynamic WM code? 'yes' or 'no'

% Load IDs
mci = meg_mci_final; 
hc = meg_hc_final; 
cog_def = meg_cog_def_final;
subj = horzcat(hc,mci,cog_def);

groups={'All','OHC','MCI'}; 
nn = [1,length(subj); 1,length(hc);  length(hc)+1, length(hc)+length(mci)]; 


rows = 'Rows'; 
pw = 'Pairwise'; 
tt = 'type'; 
mec = 'MarkerEdgeColor'; 
mfc = 'MarkerFaceColor'; 
color = 'color';
lw = 'linewidth'; 
fs = 'Fontsize'; 
corr_type = 'Pearson'; % Choose 'Pearson' or 'Spearman'


noise = []; 
lapse = []; 
lapse_noise = []; 
criterion = []; 
WMacc = []; 
PC1 = []; 


diag_offset = 0; % define shifts for sum of off-diagonal elements, zero ~ all off-diagonal elements; 9 ~ |Test-train time| >= 0.5 s (20 consecutive time points)

% Define which area to analyze ('whole_cortex' or 'Dorsal_visual')

area='Dorsal_visual';
area_title = 'Dorsal visual cortex';
 
% area='whole_cortex';
% area_title = 'Whole cortex';

onsets = [0 0.5]; %mem on and mem off


% Load behavioral data for subjects 
for s = 1:length(subj) 
         l = find(strcmp(cerad_data.ID,subj{s})); 

               noise(s,1)=cerad_data.noise(l); 
               lapse(s,1)=cerad_data.lapse(l); 
               criterion(s,1)=cerad_data.criterion(l); 
               PC1(s,1)=cerad_data.PC1_score_zsubj_all(l); 
               WMacc(s,1)=cerad_data.WMacc(l); 

               load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/lapse+noise/',subj{s},'_lapse_noise.mat']); 
               lapse_noise (s,1) = lapse_noise_single_subj; 
               clear lapse_noise_single_subj
    
end 

%% Pull data 
% Create one logical that contains which elements fall inside a significant
% cluster 

load(['/Users/ginamonov/Desktop/new_codes_before_submission/temp_gen_stats/temp_gen_stats_all_', area, '.mat'])
cluster_stats_all_c = clusters; 
cluster_stats_all_p = p_values; 
new_cluster_stats = []; 
cou = 0; 
for ee = 1:length(p_values) 
    if p_values(ee) < 0.05
        cou = cou +1; 
       new_cluster_stats(:,:,cou) = cluster_stats_all_c{1,ee}; 
    end 
end 
new_cluster_stats = sum(new_cluster_stats,3); 
new_cluster_stats = logical(new_cluster_stats>0); 

% loop over subjects 
for s = 1:length(subj) 
   if strcmp(area,'Dorsal_visual')   
      load(['/Users/ginamonov/Servers/mountpoint1/meg_analysis/Decoding/temp_gen_decode_5_35_coarse/temp_gen_' subj{s}, '_',area,'_', '1_5_35.mat']);  
   elseif strcmp(area,'whole_cortex') 
      load(['/Users/ginamonov/Servers/mountpoint1/meg_analysis/Decoding/temp_gen_decode_5_35_whole_cortex/temp_gen_' subj{s}, '_',area,'_', '1_5_35.mat']); 
   end 

        % Initialize generalisation matrix
        gen_matrix = zeros(length(times),length(times)); 
        gen_matrix(gen_matrix == 0) = nan; 
        
        % Add indices to the times
        for t = 1:length(corrs) 
           idx_train_time(t) = find(times==train_time(t)); 
           idx_gen_time (t) = find(times==gen_time(t)); 

        end 
        
        for c = 1:length(corrs) 
            gen_matrix(idx_train_time(c),idx_gen_time(c)) = corrs(c); 
        end

        gen_matrices_all (:,:,s) = gen_matrix;
        
        % add a summed value only for elements that were significant across
        % all subjects 
        significant_gen_matrix = gen_matrix;
        significant_gen_matrix(new_cluster_stats == 0) = 0; 
        gen_matrices_sig (:,:,s) = significant_gen_matrix; 
        
        
        
        % Compute the sum of off-diagnonal correlation coefficients during the delay
        % period 
        % Pull delay period elemens 
        summed_gen_matrix = gen_matrix(times>= 0.5 & times <1.5,times>= 0.5 & times <1.5);
        summed_gen_matrix_sig = significant_gen_matrix(times>= 0.5 & times <1.5,times>= 0.5 & times <1.5);
        
        % set diagonal to zero according to the diagonal offset (increase
        % for wider shifts between train and test time)
          clear significant_gen_matrix
        for uu = 1:length(summed_gen_matrix(:,1)) 
            for zz = 1:length(summed_gen_matrix(1,:))
                if zz == uu
                   summed_gen_matrix(zz,uu) = 0;

                 if diag_offset ~= 0 
                  for oo = 1:diag_offset
                   if uu-oo > 0 
                     summed_gen_matrix(zz,uu-oo) = 0;
                   end 
                   if uu+oo <= length(summed_gen_matrix)
                     summed_gen_matrix(zz,uu+oo) = 0;
                   end 
                  end 
                 end 
                end 
            end  
        end 
        
        % Do the same for only significant elements
       for uu = 1:length(summed_gen_matrix_sig(:,1)) 
            for zz = 1:length(summed_gen_matrix_sig(1,:))
                if zz == uu
                   summed_gen_matrix_sig(zz,uu) = 0;

                 if diag_offset ~= 0 
                  for oo = 1:diag_offset
                   if uu-oo > 0 
                     summed_gen_matrix_sig(zz,uu-oo) = 0;
                   end 
                   if uu+oo <= length(summed_gen_matrix_sig)
                     summed_gen_matrix_sig(zz,uu+oo) = 0;
                   end 
                  end 
                 end 
                end 
            end  
       end 
        
        % sum over all the values in the remaining matrix 
        off_diagonal_sum(s,1) = sum(summed_gen_matrix,'all'); 
        off_diagonal_sum_sig(s,1) = sum(summed_gen_matrix_sig,'all');     

        diagonal(s,:) = diag(gen_matrix); % Pull the diagonal (train == test) 
        summed_diagonal(s) = sum(diagonal(s,times>= 0.5 & times <1.5)); 

 
        clear gen_matrix corrs gen_time train_time summed_gen_matrix
  
end 

%% Make plot for wide shift method visualization 

if diag_offset == 9 
   dummy_matrix = zeros(size(gen_matrices_all(:,:,1))); 
          for uu = 1:length(dummy_matrix(:,1)) 
            for zz = 1:length(dummy_matrix(1,:))
                if zz == uu
                   dummy_matrix(zz,uu) = 1;

                 if diag_offset ~= 0 
                  for oo = 1:diag_offset
                   if uu-oo > 0 
                     dummy_matrix(zz,uu-oo) = 1;
                   end 
                   if uu+oo <= length(dummy_matrix)
                     dummy_matrix(zz,uu+oo) = 1;
                   end 
                  end 
                 end 
                end 
            end  
          end 
        dummy_matrix(times<0.5,times<0.5) = 1; 
        dummy_matrix(times<0.5,:) = 1; 
        dummy_matrix(:,times<0.5) = 1; 
        % Plot averaged error rates across all groups 
        h =  findobj('type','figure');
        cfig = length(h)+1;
        fig_w = 5; % figure width
        fig_h = 5; % figure height
        f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
        set(gcf,'renderer','Painters'),


        axis square 
        s1 = subplot(1,1,1)
        imagesc(times,times,dummy_matrix,[0,0.1])
        colormap gray
        onsets= [0 0.5]; 
        xline(onsets(1),'--','Color',colors.silver,'LineWidth',1), hold on 
        xline(onsets(2),'--','Color',colors.silver,'LineWidth',1), hold on 
        yline(onsets(1),'--','Color',colors.silver,'LineWidth',1), hold on 
        yline(onsets(2),'--','Color',colors.silver,'LineWidth',1), hold on 
        line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
        line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
        set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
        set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})
        set(gca,'YDir','normal')
         
   
            xlabel('Test time','FontSize',8)
            ylabel('Train time','FontSize',8)
      
            % Arrange plots 
            w = 0.2; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
            h = 0.2; hgap = 0.08; hoff = 0.15;

        set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
        savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_wide_shift_vis_',num2str(diag_offset),'.fig']) %Figure 8E, left inset
      close all

end 

%% Plot group differences in summed temp gen score for all off-diagonal elements during delay and for significant-only  

% Plot averaged error rates across all groups 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 18; % figure width
fig_h = 7; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


for er = 1:1 % set to 1:2 if significant only elements should be considered in a second plot
eval(['s',num2str(er),'=subplot(1,2,', num2str(er), '), hold on']) 
axis square 
off_sum2test = []; 
        if er == 1
           off_sum2test = off_diagonal_sum; 
        elseif er == 2
           off_sum2test = off_diagonal_sum_sig; 
        end 
p = []; 
P_diag = []; 
[h, p, diff, diff_null]=permtest2(off_sum2test(1:length(hc)),off_sum2test(length(hc)+1:length(hc)+length(mci)),10000); 
P_diag(1:2,1:2) = p; 
[h, P_zero_hc, diff, diff_null] = permtest(off_sum2test(1:length(hc)),0,10000); 
[h, P_zero_mci, diff, diff_null]= permtest(off_sum2test(length(hc)+1:length(hc)+length(mci)),0,10000); 

group_alloc = zeros(length(subj),1); 
group_alloc(1:length(hc)) = 1; 
group_alloc(length(hc)+1:length(hc)+length(mci)) = 2; 
for r = 1:2
errorbars(1,r) = std(off_sum2test(group_alloc == r),1)./sqrt(size(off_sum2test(group_alloc == r),1)); 
end
superbar([mean(off_sum2test(1:length(hc)),1),mean(off_sum2test(length(hc)+1:length(hc)+length(mci)),1)],'E',errorbars,'P',P_diag,'BarFaceColor', [colors.sky;colors.rosered],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',true, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
set(gca,'TickDir','out','XTick',1:2,'XTickLabel',{'OHC','MCI'}, 'Fontsize',7), xlim([-0.2 3.2]), ylim([0,30])

% Add p-values to plot 
if P_zero_hc < 10^-4
  h = text(.5,0.01,['p<10^{-4}'],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8) 
elseif P_zero_hc < 0.05
  h = text(.5,0.01,['p=',num2str(round(P_zero_hc,4))],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8)
else
  h = text(.5,0.01,['p=',num2str(round(P_zero_hc,4))],'HorizontalAlignment','left','Fontweight','normal','Fontsize',8)
end 
set(h,'Rotation',90)
clear h 

if P_zero_mci < 10^-4
  h = text(1.5,0.01,['p<10^{-4}'],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8) 
elseif P_zero_mci < 0.05
  h = text(1.5,0.01,['p=',num2str(round(P_zero_mci,4))],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8)
else
  h = text(1.5,0.01,['p=',num2str(round(P_zero_mci,4))],'HorizontalAlignment','left','Fontweight','normal','Fontsize',8)
end 
set(h,'Rotation',90)
clear h 



if er == 1
   title('WM code generalization during delay', 'Fontsize',8), ylabel('\Sigma off-diagonal')
else
    title('WM code generalization during delay', 'Fontsize',8), ylabel('\Sigma off-diagonal (only significant)')
end 
sgtitle(area_title,'fontsize',10)

end 

        % Arrange plots 
            w = 0.25; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
            h = 0.45; hgap = 0.08; hoff = 0.15;

            set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
        %   set(s2, 'Position', [woff+w+wgap,hoff, w, h])
            
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_group_diffs_',area,'_',num2str(diag_offset),'.fig']) % Figure 8C
close all

%% Is sum of wide shift off-diagonal elements >0 and correlate sum with behavioral measures WM accuracy and Noise + Lapse 

if diag_offset == 9 
   
        h =  findobj('type','figure');
        cfig = length(h)+1;
        fig_w = 8; % figure width
        fig_h = 10; % figure height
        f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
        set(gcf,'renderer','Painters'),
        off_sum2test = off_diagonal_sum; % use all off diagonal elements within the wide shifted elements 
      
        s1 = subplot(1,1,1), hold on 
        [h, P_zero_wide_shifts_all, diff, diff_null]= permtest(off_sum2test,0,10000); 
        
        if P_zero_wide_shifts_all < 10^-4
           h = text(.6,0.01,['p<10^{-4}'],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8) 
        elseif P_zero_wide_shifts_all < 0.05
           h = text(.6,0.01,['p=',num2str(round(P_zero_wide_shifts_all,4))],'HorizontalAlignment','left','Fontweight','bold','Fontsize',8)
        else
           h = text(.6,0.01,['p=',num2str(round(P_zero_wide_shifts_all,4))],'HorizontalAlignment','left','Fontweight','normal','Fontsize',8)
        end 
        set(h,'Rotation',90)
        clear h 
        
        
        superbar([mean(off_sum2test,1)],'E',std(off_sum2test,1)./sqrt(size(off_sum2test,1)) ,'P',P_zero_wide_shifts_all,'BarFaceColor', [colors.black],'PStarFontSize',15,'BarWidth',0.5,'PStarShowNS',true, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
        set(gca,'TickDir','out','XTick',1:1,'XTickLabel',{['All, N =',num2str(length(subj))]}, 'Fontsize',8), xlim([0.4 1.6]), ylim([0,4])
         ylabel('\Sigmaoff-diagonal (wide)')
            w = 0.2; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
            h = 0.3; hgap = 0.08; hoff = 0.15;

            set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
            
            savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_all_ag_zero',area,'_',num2str(diag_offset),'.fig']) %Figure 8E, right inset 
            close all
            
      % Plot correlations of wide shift off-diagonal sum and WM accuracy and Noise +Lapse parameters 
      
              % Initialize figure 
        h =  findobj('type','figure');
        cfig = length(h)+1;
        fig_w = 13; % figure width
        fig_h = 7; % figure height
        f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
        set(gcf,'renderer','Painters'),
        
        
        mycolors(1).type = colors.jetblack; 
        myedgecolors(1).type = colors.black; 

                scattersize = 15; 
                corr_param1_all = {'WMacc','lapse_noise'}; 
                corr_param1_titles = {'WM accuracy','Noise \sigma_{mem} + Lapse \theta'}; 
                corr_param2 = 'off_diagonal_sum'; 
                corr_param2_title = '\Sigma off diagonal (wide)';

                counts = 0; 
                for u = 1:length(corr_param1_all) 
                       corr_param1 = corr_param1_all{u}; 
                       for uu = 1:1 % plot only all subjects (set to 1:(length(groups) if all groups are desired 
                       counts = counts+1;  
                       eval(['s',num2str(counts),'=subplot(length(corr_param1_all),1,counts), hold on']) 
                       axis square
                       set(gca,'fontname','Helvetica','LineWidth',0.5, 'Fontsize',7,'TickDir','out','box','off')
                       eval(['[r,pval] = corr(',corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2,'(nn(uu,1):nn(uu,2)),rows,pw,tt,corr_type);']); 
                       eval(['p = polyfit(', corr_param1,'(nn(uu,1):nn(uu,2)),', corr_param2,'(nn(uu,1):nn(uu,2)),1);']); 
                       eval(['f = polyval(p,',corr_param1,'(nn(uu,1):nn(uu,2)));']); 
                       eval(['l', num2str(counts), '=scatter(', corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2, '(nn(uu,1):nn(uu,2)),scattersize,mfc,mycolors(uu).type,mec,myedgecolors(uu).type),hold on'])
                       if strcmp(corr_type,'Pearson') && pval < 0.05
                           eval(['sl', num2str(uu), '=plot(', corr_param1,'(nn(uu,1):nn(uu,2)),f,lw,1,color,mycolors(uu).type),hold on'])
                       end 

                             eval(['xlim([min(',corr_param1,'),','max(',corr_param1,')]);'])
                      eval(['ylim([min(',corr_param2,'),','max(',corr_param2,')]);'])
                      ylabel(corr_param2_title,fs,7)

                      if uu == 1
                      xlabel(corr_param1_titles{u},fs,7)
                      end 

                    if pval<0.05
                        title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,3))],fs,7); 
                    else  title(['r=',num2str(round(r,2)),', p=',num2str(round(pval,3))],fs,7,'fontweight','normal'); 
                    end 


                    eval(['corr_diffs_',corr_param1,'=permcompare_corrs(',corr_param1,'(nn(3,1):nn(3,2)),',corr_param2,'(nn(3,1):nn(3,2)),',corr_param1,'(nn(2,1):nn(2,2)),',corr_param2,'(nn(2,1):nn(2,2)),10000);'])
                       end  
                end 


         % Arrange plots 
        w = 0.15; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
        h = 0.25; hgap = 0.08; hoff = 0.15;

        set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
        set(s2, 'Position', [woff+w*1.3+wgap,hoff, w, h])

        savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_wide_shift_corrs_',area,'.fig'])

      
end 

%% Run correlation with behavioral parameters 
% Specify parameters to correlate with PC1
corr_param1_all = {'WMacc','lapse_noise','criterion','PC1'}; %,'noise','lapse'}; % lapse_noise noise lapse WMacc PC1 criterion 
corr_param1_titles = {'WM accuracy','Noise \sigma_{mem} + Lapse \theta','Threshold \delta','PC1 score'};%,'Noise \sigma_{mem}','Lapse \theta'};

off_sum2test = []; 
off_sum2test = off_diagonal_sum; % Choose whether all off diagonal elements or only significant off diagonal elememts should be considered 
        
corr_param2 = 'off_sum2test'; % Define which paramter to correlate with
corr_param2_titles = {'\Sigma off-diagonal'}; % give this a title 

% Plot averaged error rates across all groups 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 18; % figure width
fig_h = 5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

% Specifiy colors 
mycolors(1).type = colors.jetblack; 
mycolors(2).type = colors.sky; 
mycolors(3).type = colors.rosered; 
myedgecolors(1).type = colors.black;
myedgecolors(2).type = colors.teal; 
myedgecolors(3).type = colors.lipstick; 


% Make plot where each row contains the correlations for a single parameter
% to correlate and each column contains the groups (all, ohc, mci) 

counts = 0; 

for u = 1:length(corr_param1_all)  
     corr_param1 = corr_param1_all{u}; 
    counts = counts +1; 
    eval(['s',num2str(counts),'=subplot(1,length(corr_param1_all),counts), hold on']) % make subplot
    
    for uu = 1:length(groups) 
       eval(['[r_', groups{uu},', pval_',groups{uu},'] = corr(',corr_param1,'(nn(uu,1):nn(uu,2)),',corr_param2,'(nn(uu,1):nn(uu,2)),rows,pw,tt,corr_type)']);  
    end 
    
     eval(['corr_diffs_',corr_param1,'=permcompare_corrs(',corr_param1,'(nn(3,1):nn(3,2)),',corr_param2,'(nn(3,1):nn(3,2)),',corr_param1,'(nn(2,1):nn(2,2)),',corr_param2,'(nn(2,1):nn(2,2)),10000);'])
     eval(['p_diff = corr_diffs_', corr_param1, '.Delta_p;'])
     eval(['r_diff = corr_diffs_', corr_param1, '.Delta_r;'])
     p_group = ones(length(groups),length(groups)); 
     p_group(2,3) = p_diff; 
     p_group(3,2) = p_diff; 
     %'P',p_group,
     superbar([r_All, r_OHC, r_MCI],'BarFaceColor', [myedgecolors(1).type;myedgecolors(2).type;myedgecolors(3).type],'PStarFontSize',15,'BarWidth',0.7,'PStarShowNS',false, 'ErrorbarLineWidth',0.3,'PLineWidth',0.5), hold on
     set(gca,'TickDir','out','XTick',1:3,'XTickLabel',{'All','OHC','MCI'}, 'Fontsize',7), xlim([0.3 3.7]), ylim([-1 1]), hold on 
     title(corr_param1_titles{u}, 'Fontsize',8), 
     if u ==1, ylabel('r'), end 
     
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

savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_bar_corrs_',area,'_',num2str(diag_offset),'.fig']) %Figure 8D-E
close all


%% Run analyses for testing dynamic WM coding 

if strcmp(run_dynamic_coding,'yes')

        delay_matrices = gen_matrices_all; 
         ps_tp1 = [];
         ps_tp2 = [];
         dyn = [];
         dyn_mask = []; 
         indx = []; 
        for tp1 = 1:length(times) 
           for tp2 = 1:length(times) 
               if times(tp1) ~= times(tp2) %omit diagonal for dynamic code test
                   to_test = squeeze(delay_matrices(tp1,tp2,:)); 
                   comp_tp1 = squeeze(delay_matrices(tp1,tp1,:)); 
                   comp_tp2 = squeeze(delay_matrices(tp2,tp2,:)); 
                   [h, p1, diff, diff_null]=permtest2(to_test,comp_tp1,10000); 
                   [h, p2, diff, diff_null]=permtest2(to_test,comp_tp2,10000); 
                   ps_tp1(end+1,1) = p1;
                   indx(end+1,1:2) = [tp1,tp2]; %save indices that belong to p-values 
                   ps_tp2(end+1,1) = p2;
                   if p1 < 0.05 && p2 <0.05
                      dyn_mask(tp1,tp2) = 1; % save logical of signifcant elements 
                   else dyn_mask(tp1,tp2) = 0; 
                   end 
                   clear p1 p2 to_test comp_tp1 comp_tp2
               else dyn_mask(tp1,tp2)=0;
               end 
           end 
        end
        fdr_ps_tp1 = fdr(ps_tp1,0.05); 
        fdr_ps_tp2 = fdr(ps_tp2,0.05); 

        fdr_corrs = fdr_ps_tp1 + fdr_ps_tp2; 
        fdr_corrs(fdr_corrs~=2) = 0; % set to zero if not both are significant after fdr-correctioson 
        fdr_corrs(fdr_corrs == 2) = 1; % set to 1 only when both ps were signifcant afte fdr-correction 

        for uu = 1:length(fdr_corrs) 
            if fdr_corrs(uu) == 1
               dyn_mask_fdr(indx(uu,1),indx(uu,2)) = 1; 
            else dyn_mask_fdr(indx(uu,1),indx(uu,2)) = 0; 
            end 
        end 
        
        % Initialize figure 
        h =  findobj('type','figure');
        cfig = length(h)+1;
        fig_w = 13; % figure width
        fig_h = 10; % figure height
        f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
        set(gcf,'renderer','Painters'),
        subj_counter = 0;


        for re = 1:2
        fdr_corrected_choices = {'no','yes'}; 
        fdr_corrected = fdr_corrected_choices{re}; 
        eval(['s',num2str(re),'=subplot(1,2,',num2str(re),'), hold on'])
        axis square
        contourf(times, times, mean(gen_matrices_all,3),50,'linestyle','none'), hold on 
        if strcmp(fdr_corrected,'no')
        contour(times,times,dyn_mask,1,'Color',colors.white), hold on
        else contour(times,times,dyn_mask_fdr,1,'Color',colors.white), hold on
        end 
        clims = [-0.01, 0.3];
        caxis(clims);
        cmap=brewermap(256, '*Spectral'); % take red-blue for all others 
        colormap(cmap)

        xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
        xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
        yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
        yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
        line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
        line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
        set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
        set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})

        if strcmp(fdr_corrected,'yes')
            title('FDR-corrected')
        else title('unthresholded')
        end 
        
        if re == 2
           
         
             cb = colorbar('east'), hold on 
        
             pos1 = get(gca,'position');
             pos = get(cb,'position'); 
             pos(3) = 0.01; 
             pos(4) = 0.2; 
             set(cb, 'Position', [0.95,0.3, 0.01, 0.5])
         
        end 
   
            xlabel('Test time','FontSize',8)
            ylabel('Train time','FontSize',8)
      
        
        end 
        
        %% Arrange plots 
w = 0.1638; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
h = 0.3393; hgap = 0.08; hoff = 0.15;

set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
set(s2, 'Position', [woff+w*1.7+wgap,hoff, w, h])

set(cb, 'Position', [woff+(w+wgap)*1.9+w+0.01, hoff, cw, h])

savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_dynamic_',area,'.fig']) %Extended Data Figure 8-1
end 


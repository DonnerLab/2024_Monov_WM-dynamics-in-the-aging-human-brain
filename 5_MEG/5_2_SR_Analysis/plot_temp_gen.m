% Plot temporal generalization matrices 
% Options include single subjects, all subjects, OHC, MCI and difference
% maps 
% cluster-based permutation test: Edden M. Gerber (2023). permutest, MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/71737-permutest
% --> Figure 8A-B
% Gina Monov, UKE, 2024 

clear all
close all
addpath /mnt/homes/home030/aarazi/fieldtrip-20201009/
addpath '/Users/ginamonov/Servers/mountpoint1/functions/'
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/meg_subj.mat']);
load('/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat')
load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')

mci = meg_mci_final; 
hc = meg_hc_final; 
cog_def = meg_cog_def_final;
subj = horzcat(hc,mci,cog_def);

groups={'All','OHC','MCI'}; 
nn = [1,length(subj); 1,length(hc);  length(hc)+1, length(hc)+length(mci)]; 

plot_difference_matrix = 'no'; % indicate whether to matrix of group difference
plot_single_subj = 'no'; % want to plot single subjects in a separate figure? 
compare_diag_ini = 'yes'; % want to compare diagonal to time-resolved decoding? 
plot_examples = 'yes'; % Want to plot a few single subjects as examples? 

noise = []; 
lapse = []; 
lapse_noise = []; 
criterion = []; 
WMacc = []; 
PC1 = []; 
diag_offset = 0; % define shifts for sum of off-diagonal elements, zero ~ all off-diagonal elements; 9 ~ |Test-train time| >= 0.5 s (20 consecutive time points)

% Define are: 'Dorsal_visual' or 'whole_cortex'
area='Dorsal_visual';
area_title = 'Dorsal visual cortex';

% area='whole_cortex';
% area_title = 'Whole cortex';

onsets = [0 0.5]; %mem on and mem off

% Initialize figure 
if strcmp(plot_single_subj,'yes')
    h =  findobj('type','figure');
    cfig = length(h)+1;
    fig_w = 18; % figure width
    fig_h = 25; % figure height
    f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
    set(gcf,'renderer','Painters'),
elseif strcmp(plot_examples,'yes')
    h =  findobj('type','figure');
    cfig = length(h)+1;
    fig_w = 18; % figure width
    fig_h = 7; % figure height
    f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
    set(gcf,'renderer','Painters'),
end 
    
subj_counter = 0;
s_subj_counter = 0;

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

% Compute the sum of off-diagnonal correlaion coefficients during the delay
% period 

summed_gen_matrix = gen_matrix(times>= 0.5 & times <1.5,times>= 0.5 & times <1.5);
% set diagonal to zero

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

 diagonal(s,:) = diag(gen_matrix); 
 
 summed_diagonal(s) = sum(diagonal(s,times>= 0.5 & times <1.5)); 


% sum over all the values in the remaining matrix 
off_diagonal_sum(s,1) = sum(summed_gen_matrix,'all'); 
% Save matrices of all subjects in one variable  
gen_matrices_all (:,:,s) = gen_matrix;
subj_counter = subj_counter+1; 
if strcmp(plot_single_subj,'yes')
        subplot(6,7,s)
        contourf(times, times, gen_matrix,50,'linestyle','none') 
        axis square
        clims = [-0.01, 0.3];
        caxis(clims);
      
        if ismember(subj{s},mci)
        title(['#',num2str(subj_counter),' (MCI)'])
        elseif ismember(subj{s},hc)
        title(['#',num2str(subj_counter),' (OHC)'])
        elseif ismember(subj{s},cog_def)
        title(['#',num2str(subj_counter),' (UNC)'])
        end 
        cmap=brewermap(256, '*Spectral'); 
        colormap(cmap)
        if s==1 
        xlabel('Test time','FontSize',6)
        ylabel('Train time','FontSize',6)
        end 
        xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
        xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
        yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
        yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
        line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
        line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
        set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
        set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})
elseif  strcmp(plot_examples,'yes')
    
    if subj_counter == 15 || subj_counter == 28 || subj_counter == 23 || subj_counter == 36 
        s_subj_counter = s_subj_counter+1; 
        eval(['s',num2str(s_subj_counter),'=subplot(1,4,',num2str(s_subj_counter), '),hold on']) 
        contourf(times, times, gen_matrix,50,'linestyle','none') 
        axis square
        clims = [-0.01, 0.3];
        caxis(clims);
        if ismember(subj{s},mci)
        title(['Example #',num2str(s_subj_counter),' (MCI)'])
        elseif ismember(subj{s},hc)
        title(['Example #',num2str(s_subj_counter),' (OHC)'])
        elseif ismember(subj{s},cog_def)
        title(['Example #',num2str(s_subj_counter),' (UNC)'])
        end 
        cmap=brewermap(256, '*Spectral'); 
        colormap(cmap)
        if s==1 
        xlabel('Test time','FontSize',6)
        ylabel('Train time','FontSize',6)
        end 
        xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
        xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
        yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
        yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
        line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
        line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
        set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
        set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})
         xlabel('Test time (s)','Fontsize',6)
         ylabel('Train time (s)','Fontsize',6) 
         if s_subj_counter == 4
            cb = colorbar('east'), hold on 
            pos1 = get(gca,'position');
            pos = get(cb,'position')
            pos(3) = 0.01; 
            pos(4) = 0.2; 
            set(cb,'position',pos)
            set(gca,'position',pos1); 
             
         
         % Arrange plots 
            w = 0.1638; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
            h = 0.3393; hgap = 0.08; hoff = 0.15;

            set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
            set(s2, 'Position', [woff+w+wgap,hoff, w, h])
            set(s3, 'Position', [woff+(w+wgap)*2, hoff, w, h])
            set(s4, 'Position', [woff+(w+wgap)*3, hoff, w, h])
            set(cb, 'Position', [woff+(w+wgap)*3+w+0.02, hoff, cw, h])
         end

             
    end 
    
end 
clear gen_matrix corrs gen_time train_time summed_gen_matrix
  
end 

if strcmp(plot_single_subj,'yes')
% Save single subjects figure 
sgtitle(area_title,'fontsize',10)
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_all_',area,'.fig'])
close all
end

if strcmp(plot_examples,'yes')
          % Save single subjects figure 
        sgtitle(area_title,'fontsize',10)
        savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_all_example_subjs',area,'.fig'])
        close all
end


if strcmp(compare_diag_ini,'yes')
for t = 1:length(subj)
    subplot(6,7,t)
    axis square
    plot(times,diagonal(t,:),'color',colors.green), hold on
    exact_loc = readtable(['/Users/ginamonov/Servers/mountpoint1/meg_analysis/Decoding/decode_5_35_coarse/' subj{t}, '_Dorsal_visual_1_5_35.csv']);
    plot(times,diagonal(t,:),'color',colors.green), hold on
    plot(exact_loc.latency,exact_loc.test_correlation,'color',colors.blue)
end
sgtitle(area_title,'fontsize',10)
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/temp_gen_compare2dec_',area,'.fig'])
close all
end 

%% Plot all, MCI, OHC, Difference (only if selected)
% construct train, test indices for cluster-based permutation test 
index = [idx_gen_time;idx_train_time;1:(length(times)*length(times))]';

% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 18; % figure width
fig_h = 7; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

gen_matrices_ohc = gen_matrices_all(:,:,1:length(hc)); 
gen_matrices_mci = gen_matrices_all(:,:,length(hc)+1:length(hc)+length(mci)); 


mean_gen_matrices_all = mean(gen_matrices_all,3); 
all_max = max(max(mean_gen_matrices_all)); 
all_min = min(min(mean_gen_matrices_all)); 

%% ALL

cb_count = 0; %counter for colorbars 
east = 'East';
south = 'south';
position = 'position';
s1 = subplot(1,3,1), hold on 

contourf(times, times, mean_gen_matrices_all,50,'linestyle','none') 
axis square
clims = [all_min, all_max];
clims = [-0.01, 0.3];
caxis(clims);
cmap=brewermap(256, '*Spectral'); % take red-blue for all others 
colormap(cmap)
% Run stats
null4stats = zeros(size(gen_matrices_all)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(gen_matrices_all,null4stats,true,0.05,10000,true);

%translaste clusters into indices for plotting contour
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             
             if p_values(ccc) < 0.05
                stat_mask = zeros(length(times),length(times)); 
                clust_idx = clusters{1,ccc};
                for rr=1:length(clust_idx) 
                   idx = find(clust_idx(rr) == index(:,3));
                   new_clust_idx(1,rr) = index(idx,1);
                   new_clust_idx(2,rr) = index(idx,2);
                end 
                new_clust_idx = new_clust_idx';
                for i = 1:length(new_clust_idx(:,1))
                    stat_mask(new_clust_idx(i,1),new_clust_idx(i,2))=1;
                    
                end 
               stat_mask = logical(stat_mask); 
               clusters{1,ccc} = stat_mask;
               contour(times,times,clusters{1,ccc},1,'Color','k'), hold on 
               clear stat_mask new_clust_idx
             end
         end 
         save(['/Users/ginamonov/Desktop/new_codes_before_submission/temp_gen_stats/temp_gen_stats_all_', area, '.mat'],'clusters','p_values')
         title(['All subjects, N=',num2str(length(subj))])
         xlabel('Test time (s)','Fontsize',8)
         ylabel('Train time (s)','Fontsize',8)
         line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
         line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
         xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})

%% OHC only 
s2=subplot(1,3,2), hold on 

contourf(times, times, mean(gen_matrices_ohc,3),50,'linestyle','none') 
axis square
clims = [all_min, all_max];
clims = [-0.01, 0.3];
caxis(clims);
cmap=brewermap(256, '*Spectral'); % take red-blue for all others 
colormap(cmap)

% Run stats
null4stats = zeros(size(gen_matrices_ohc)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(gen_matrices_ohc,null4stats,true,0.05,10000,true);

%translaste clusters into indices for plotting contour
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             stat_mask = zeros(length(times),length(times)); 
             if p_values(ccc) < 0.05
                clust_idx = clusters{1,ccc};
                for rr=1:length(clust_idx) 
                   idx = find(clust_idx(rr) == index(:,3));
                   new_clust_idx(1,rr) = index(idx,1);
                   new_clust_idx(2,rr) = index(idx,2);
                end 
                new_clust_idx = new_clust_idx';
                for i = 1:length(new_clust_idx(:,1))
                    stat_mask(new_clust_idx(i,1),new_clust_idx(i,2))=1;
                    
                end 
               stat_mask = logical(stat_mask); 
               clusters{1,ccc} = stat_mask;
               contour(times,times,clusters{1,ccc},1,'Color','k'), hold on 
               clear stat_mask new_clust_idx
             end
         end 
         set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
         line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
         line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on
         
         set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})
         title(['OHC, N=',num2str(length(hc))],'Fontsize',8)
         

%% MCI only 
s3 = subplot(1,3,3), hold on 

contourf(times, times, mean(gen_matrices_mci,3),50,'linestyle','none') 
axis square
clims = [all_min, all_max];
clims = [-0.01, 0.3];
caxis(clims);
cmap=brewermap(256, '*Spectral'); % take red-blue for all others 
colormap(cmap)

% Run stats
null4stats = zeros(size(gen_matrices_mci)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(gen_matrices_mci,null4stats,true,0.05,10000,true);

%translaste clusters into indices for plotting contour
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             stat_mask = zeros(length(times),length(times)); 
             if p_values(ccc) < 0.05
                clust_idx = clusters{1,ccc};
                for rr=1:length(clust_idx) 
                   idx = find(clust_idx(rr) == index(:,3));
                   new_clust_idx(1,rr) = index(idx,1);
                   new_clust_idx(2,rr) = index(idx,2);
                end 
                new_clust_idx = new_clust_idx';
                for i = 1:length(new_clust_idx(:,1))
                    stat_mask(new_clust_idx(i,1),new_clust_idx(i,2))=1;
                    
                end 
               stat_mask = logical(stat_mask); 
               clusters{1,ccc} = stat_mask;
               contour(times,times,clusters{1,ccc},1,'Color','k')
                
               clear stat_mask new_clust_idx
             end
         end 
         
             cb_count = cb_count+1; 
             eval(['cb',num2str(cb_count),'=colorbar(east), hold on'])
             pos1 = get(gca,'position');
             eval(['pos = get(cb',num2str(cb_count),',position)'])
             pos(3) = 0.01; 
             pos(4) = 0.2; 
             eval(['set(cb',num2str(cb_count),',position,pos)']) 
             set(gca,'position',pos1);          
         
             set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
           set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})
title(['MCI, N=',num2str(length(mci))],'Fontsize',8)

  line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
         line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on
         
         
         


%% OHC-MCI
if strcmp(plot_difference_matrix,'yes')
     
s4 = subplot(1,4,4), hold on 

contourf(times, times, mean(gen_matrices_ohc,3)-mean(gen_matrices_mci,3),50,'linestyle','none') 
axis square
clims = [all_min, all_max];
clims = [-0.1, 0.1];
caxis(clims);
cmap=brewermap(256, '*Spectral'); % take red-blue for all others 
colormap(cmap)

% Run stats

[clusters,p_values,t_sums,permutation_distribution] = permutest(gen_matrices_ohc,gen_matrices_mci,false,0.05,10000,true);

%translaste clusters into indices for plotting contour
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             stat_mask = zeros(length(times),length(times)); 
             if p_values(ccc) < 0.05
                clust_idx = clusters{1,ccc};
                for rr=1:length(clust_idx) 
                   idx = find(clust_idx(rr) == index(:,3));
                   new_clust_idx(1,rr) = index(idx,1);
                   new_clust_idx(2,rr) = index(idx,2);
                end 
                new_clust_idx = new_clust_idx';
                for i = 1:length(new_clust_idx(:,1))
                    stat_mask(new_clust_idx(i,1),new_clust_idx(i,2))=1;
                    
                end 
               stat_mask = logical(stat_mask); 
               clusters{1,ccc} = stat_mask;
               contour(times,times,clusters{1,ccc},1,'Color','k'), hold on 
               
               clear stat_mask new_clust_idx
             end
         end 
         
             cb_count = cb_count+1; 
             eval(['cb',num2str(cb_count),'=colorbar(east), hold on'])
             pos1 = get(gca,'position');
             eval(['pos = get(cb',num2str(cb_count),',position)'])
             pos(3) = 0.01; 
             pos(4) = 0.2; 
             eval(['set(cb',num2str(cb_count),',position,pos)']) 
             set(gca,'position',pos1); 
         
 title('OHC-MCI','Fontsize',8)    

 set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.5 1 1.5])
 set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{'0' '0.5' '1' '1.5'})
         
 sgtitle(area_title,'Fontsize',10,'Fontweight','bold');       
 end 
%% Arrange plots 
w = 0.1638; wgap = 0.025; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.08;
h = 0.3393; hgap = 0.08; hoff = 0.15;

set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
set(s2, 'Position', [woff+w+wgap,hoff, w, h])
set(s3, 'Position', [woff+(w+wgap)*2, hoff, w, h])

set(cb1, 'Position', [woff+(w+wgap)*2+w+0.01, hoff, cw, h])

if strcmp(plot_diff_matrix,'yes')
   set(s4, 'Position', [woff+(w+wgap)*3+1.5*wgap, hoff, w, h])
   set(cb2, 'Position', [woff+(w+wgap)*3+1.5*wgap+w+0.01, hoff, cw, h])
end 
  line(zeros(size(times(times>=0 & times<=0.5)))+min(times),times(times>=0 & times<=0.5),'LineWidth',1, 'Color',colors.silver)
         line(times(times>=0 & times<=0.5),zeros(size(times(times>=0 & times<=0.5)))+min(times),'LineWidth',1, 'Color',colors.silver)
xline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(1),'--','Color',colors.white,'LineWidth',1), hold on 
yline(onsets(2),'--','Color',colors.white,'LineWidth',1), hold on

%% Save groups figures 
savefig(['/Users/ginamonov/Servers/mountpoint1/final_figures/',area,'_temp_gen_subgroups.fig'])



% Make the brain map plots for neural markers across cortex for 22 Glasser
% groups 
% Specify which neural marker to plot, whether a significance threshold
% should be applied (FDR-correction) 
% Permutation test: https://github.com/rudyvdbrink/Tools/tree/master/statistics
% --> Figure 5B, Figure 6B 
% Gina Monov, UKE, 2023

clear all
close all
% loading support files and templates 
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
addpath /mnt/homes/home030/aarazi/fieldtrip-20201009/
addpath /mnt/homes/home028/gmonov/functions/
ft_defaults
atlas = ft_read_cifti('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/support_files/Glasser_atlas.dlabel.nii');
Areas = readtable('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/Glasser_labels.csv'); Areas_labels=Areas.Label;
addpath /mnt/homes/home030/aarazi/Confirmation_bias/code/ColorMap_2
%load table with model parameters, cognitive tests, WM performace
load('/mnt/homes/home028/gmonov/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')
%Subject IDs
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/meg_subj.mat']);
mci = meg_mci_final; 
hc = meg_hc_final; 
cog_def = meg_cog_def_final;
allsubj = horzcat(mci,hc,cog_def); 

% Specify which neural marker to plot: dec_prec, mean_TFR, across_trial, within_trial (thresholded or unthresholded?)
plot_significant_only = 'yes'; % set to 'yes' to only plot parcels that are significant after fdr correction, 'no' to plot unthresholded results 
plot_diff = 'no'; % choose 'yes' if group difference maps are desired
if strcmp(plot_diff,'yes')
    n_maps = 4; 
else n_maps = 3; 
end 
% Select which neural markers to plot 
%neural_markers = {'across_trial'}; % Across-trial variability 
%neural_markers = {'mean_TFR'}; % Power modulations 
%neural_markers = {'within_trial'}; %Within-trial variability 
neural_markers = {'dec_prec'}; %Decoding precision 

% Load 22 groups defined in supplementary material of glasser et al. 2016 
load(['/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/glasser_22_groups.mat']);

% get the lables for finegrained ROIs for plotting later 
ROIs180={}; 
for l = 1:length(Areas_labels)./2
    roi= Areas_labels(l); 
    ROIs180{l} = roi{1,1}(3:end); 
end 
% Naming from pymeg is slightly different 
new_ROIs180={};
for g = 1:length(ROIs180) 
    new_ROIs180{g} = replace(ROIs180{g},'_','-'); 
end 


n_vertices=64984; % n of vertices in the template 
counts = 0; % initialize counter for brain maps
cb_count = 0; % initialize counter for colorbars 

% Initialize figure 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 20.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),



for markers = 1:length(neural_markers) % loop over neural markers: just define 1 above, since the plot arrangement deals with only 1 row 
Vertices_values=NaN(n_vertices, 1);
vals2plot = nan(length(allsubj),length(ROIs)); % Intitialite vals2plot

% Pull values to plot

vals2plot=nan(length(allsubj),22); %initialize values only for 1 hemi set of rois 
 if strcmp(plot_significant_only,'no')
     pvals_rois = nan(length(ROIs),1); % Initialize p-values vector if only significant parcels should be plotted 
 end 
for o = 1:length(ROIs(1,:)) 
    for s = 1:length(allsubj)
       
       vals_in_parcel = []; % Initialize values for a certain cluster that will be averaged 
       for tui = 1:length(glasser_22(:,1))
           if ~isempty(glasser_22{tui,o})
              val=load(['/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/',neural_markers{markers},'/',allsubj{s},'_',glasser_22{tui,o},'.mat']);  %Load value of neural marker that is of interest 
              eval(['vals_in_parcel(tui) = val.',neural_markers{markers},'(1,1);']) 
              clear val
           else break
           
           end   
       end 

       vals2plot(s,o)=mean(vals_in_parcel); % average values for the subparcels in this ROI and add to Subject X ROI matrix 
        
      
    end
end  

     % Test significance 
    if strcmp(plot_significant_only,'yes')
        p_thres = 0.05; %threshold for permutation test 
        if strcmp(neural_markers,'dec_prec') % specify tail for decoding precision 
            tail = 'right';
        else
            tail = 'both';
        end 
        for stat = 1:length(ROIs) %loop over rois to make stats for each instance 
        p = []; 
        [h, p, diff, diff_null] = permtest(vals2plot(:,stat),0,10000,p_thres,tail); 
        pvals_rois_all(stat) = p; 
        p = []; 
        [h, p, diff, diff_null] = permtest(vals2plot(1:length(mci),stat),0,10000,p_thres,tail); 
        pvals_rois_mci(stat) = p; 
        [h, p, diff, diff_null] = permtest(vals2plot(length(mci)+1:length(mci)+length(hc),stat),0,10000,p_thres,tail); 
        pvals_rois_hc(stat) = p; 
        [h, p, diff, diff_null] = permtest2(vals2plot(1:length(mci),stat),vals2plot(length(mci)+1:length(mci)+length(hc),stat),10000); 
        pvals_rois_group_diff(stat) = p;         
        end 
     
        fdr_all = fdr(pvals_rois_all,0.05); 
        fdr_mci = fdr(pvals_rois_mci,0.05); 
        fdr_hc = fdr(pvals_rois_hc,0.05); 
        fdr_group_diff = fdr(pvals_rois_group_diff,0.05);

        fdr_all = [fdr_all,fdr_all]; 
        fdr_mci = [fdr_mci,fdr_mci]; 
        fdr_hc = [fdr_hc,fdr_hc];  
        fdr_group_diff = [fdr_group_diff,fdr_group_diff]; 
    end 
    


    %Simply double to get the same data on both hemispheres since the data
    %has been averaged over both hemispheres already 

 vals2plot = [vals2plot,vals2plot]; 

 % Pull the averaged data across all subjects, mci, ohc, and ohc-mci 
 mean_all = mean(vals2plot,1); 
 mean_mci = mean(vals2plot(1:length(mci),:),1);
 mean_hc = mean(vals2plot(length(mci)+1:length(mci)+length(hc),:),1);
 diff_mci_hc = mean_hc-mean_mci; 



 % If only significant should be plotted, set the data to nan according to
 % results of fdr-correction 
 if strcmp(plot_significant_only,'yes')
    for ee = 1:length(vals2plot(1,:))
        if fdr_all(ee) == 1 
            mean_all(ee) = mean_all(ee); 
        else mean_all(ee) = nan; 
        end
        if fdr_mci(ee) == 1 
            mean_mci(ee) = mean_mci(ee); 
        else mean_mci(ee) = nan; 
        end
        if fdr_hc(ee) == 1 
            mean_hc(ee) = mean_hc(ee); 
        else mean_hc(ee) = nan; 
        end 
        if fdr_group_diff(ee) == 1 
            diff_mci_hc(ee) = diff_mci_hc(ee); 
        else diff_mci_hc(ee) = nan; 
        end        
                
    end 

    % Save which of the clusters are significant and compare the
    % mean of these between mci and ohc 

      dec_significant_parcels = zeros(length(mean_all),1); 
      dec_significant_parcels(~isnan(mean_all)) = 1;
      dec_significant_parcels(isnan(mean_all)) = 0; 
      dec_significant_parcels=dec_significant_parcels(1:22);
      dec_prec_all_sig = vals2plot(:,1:22);
      dec_prec_all_sig = dec_prec_all_sig(:,dec_significant_parcels==1); 
      dec_prec_all_sig = mean(dec_prec_all_sig,2);
      dec_prec_rois = {}; 
      for iu = 1:length(ROIs)
          if dec_significant_parcels(iu) == 1
            dec_prec_rois = horzcat(dec_prec_rois,ROIs{iu});
          end 
      end 
      [h, p_sig_parcel_group_diff, diff, diff_null] = permtest2(dec_prec_all_sig(1:length(mci)),dec_prec_all_sig(length(mci)+1:length(mci)+length(hc)),10000); 
 
 end 
 
 for u = 1:n_maps   % loop for plotting all (1), hc (2), mci (3) and difference maps (4) next to each other 
     
    % redefine vals2plot as the group 2 be plotted 
    if u == 1, vals2plot = mean_all; 
    elseif u == 2, vals2plot = mean_hc; 
    elseif u == 3, vals2plot = mean_mci; 
    elseif u == 4, vals2plot = diff_mci_hc; 
    end 

    %Tranform vals2plot from 22 parcels into 180 ROIs 
     vals2plot_new = zeros(1,length(ROIs180)); 
     vals2plot_new(vals2plot_new==0)=nan; 

     for puz = 1:length(new_ROIs180) 
         for ruz = 1:length(ROIs)
             for quz = 1:length(glasser_22(:,1))
               if strcmp(new_ROIs180{puz},glasser_22{quz,ruz})
                   % disp (glasser_22{quz,ruz}) % Check if it works
                   vals2plot_new(puz)= vals2plot(ruz);
                   break
               end 
               
             end 
         end 
     end 
vals2plot = [vals2plot_new,vals2plot_new]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i=1:size(Areas_labels,1)
    ROI=cell2mat(Areas_labels(i));
    id_Area=find(strcmp(Areas_labels, ROI));
    if ~isempty(id_Area)
        Atlasid=atlas.indexmax;
        idAr=find(Atlasid==id_Area);
        Vertices_values(idAr)=vals2plot(i);
    end
 end

 %Change clims appropriately for each neural marker
 if u < 4 % all limits should be the same for the first three maps 
     
   find_clims = [abs(min(mean_all)),abs(max(mean_all)),abs(min(mean_mci)),abs(max(mean_mci)),abs(min(mean_hc)),abs(max(mean_hc))]; % find appropriate values
   clim = [-max(find_clims)-0.05 max(find_clims)+0.05]; 
     
 elseif u == 4 % new limits for the difference maps 
   
   find_clims = [abs(min(diff_mci_hc)),abs(max(diff_mci_hc))]; 
   clim = [-max(find_clims)-0.01 max(find_clims)+0.01];

 end 
 
  % Get rid of negative limits for trial variability measures & specify
  % color map 
 if strcmp(neural_markers{markers},'across_trial') || strcmp(neural_markers{markers},'within_trial') 
     if u < 4 % don't do that for the difference maps 
         new_clim = min([min(mean_all),min(mean_mci),min(mean_hc)]); 
         clim = [new_clim,max(clim)];
         cmap=brewermap(256, '*Reds');
         cmap = flipud(cmap); 
     end 
 else cmap=brewermap(256, '*RdBu'); % take red-blue for all others 
 end 

 Vertices_values(isnan(Vertices_values))=100; % missing valus (NaN) gets very high values, to be plotted as gray 


 cmap=[cmap; colors.silver]; % defining color map; last row is set to be gray
if isnan(clim(1))
   clim = [-0.1,0.1]; % if no data is plotted just put any limits for the code to run 
end 
  

 % plot Left hemisphere, by taking the first half of "Vertices_values"
for x = 1:2 % Plot twice to rotate one of them later to display medial brain 
     counts = counts+1; 

     eval(['s',num2str(counts),'=subplot(length(neural_markers),4.*2,counts), hold on'])
     

     Plot_scalp_map('L', cmap, clim, Vertices_values(1:n_vertices/2)); 
     
     
     if x == 1
     view([270 0]); 
     camlight infinite
     lighting flat
     material dull 
     end 
     if x == 2
         
      view([90 0]); 
      camlight infinite
      lighting flat
      material dull
     end
     xlim([-100 50]); zlim([-70 80]); ylim([-105 70])
     if x == 1
     % add titles in the first row     
         if markers == 1
           if u == 1
            title(['All subjects, N=',num2str(length(allsubj))],'Fontsize',8) 
           elseif u == 2
            title(['OHC, N=',num2str(length(hc))],'Fontsize',8) 
           elseif u == 3
            title(['MCI, N=',num2str(length(mci))],'Fontsize',8) 
           elseif u == 4
            title('OHC-MCI','Fontsize',8) 
           end 
         end 

     end 
     east = 'East'; 
     position = 'Position';
         if u >= 3 && x == 2
             cb_count = cb_count+1; 
             eval(['cb',num2str(cb_count),'=colorbar(east), hold on'])
             pos1 = get(gca,'position');
             eval(['pos = get(cb',num2str(cb_count),',position)'])
             pos(3) = 0.01; 
             pos(4) = 0.2; 
             eval(['set(cb',num2str(cb_count),',position,pos)']) 
             set(gca,'position',pos1); 

        end 

     
end

    

 end 
end 

% Arrange plots 
w = 0.1638; wgap = -0.08; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.05;
h = 0.3393; hgap = 0.08; hoff = 0.12;

set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
set(s2, 'Position', [woff+w-ngap,hoff, w, h])
set(s3, 'Position', [woff+w*2+wgap*1, hoff, w, h])
set(s4, 'Position', [woff+w*3+wgap*1-ngap, hoff, w, h])
set(s5, 'Position', [woff+w*4+wgap*2, hoff, w, h])
set(s6, 'Position', [woff+w*5+wgap*2-ngap, hoff, w, h])

if strcmp(plot_diff,'yes')
set(s7, 'Position', [woff+w*6.2+wgap*3, hoff, w, h])
set(s8, 'Position', [woff+w*7.2+wgap*3-ngap, hoff, w, h])
set(cb2, 'Position', [woff+w*8.25+wgap*3-ngap,hoff,cw,h])
end 

set(cb1, 'Position', [woff+w*5.95+wgap*2-ngap,hoff,cw,h])

% Save Figure 
savefig(['/mnt/homes/home028/gmonov/final_figures/brain_plots22_',plot_significant_only,'_',neural_markers{markers},'.fig']) % Figure 5B, Figure 6B
close all

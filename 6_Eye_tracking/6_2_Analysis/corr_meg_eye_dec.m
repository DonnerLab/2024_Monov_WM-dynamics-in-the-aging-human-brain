% Analyze relationship of MEG decoding precision to decoding precision
% from gaze positions for each of the 22 Glasser groups


clear all
close all

% loading support files and templates 
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
addpath /mnt/homes/home030/aarazi/fieldtrip-20201009/
addpath /mnt/homes/home028/gmonov/functions/
addpath /mnt/homes/home028/gmonov/meg_analysis/Brain_plot/
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

%Exclude subjects that had to be excluded for analysis of eyetracking data
load('/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_subj_ids.mat')
eye_ids = horzcat(mci_ids_new,ohc_ids_new,unc_ids_new); 

allsubj = allsubj(ismember(allsubj,eye_ids)); 

plot_significant_only = 'no'; % set to 'yes' to only plot parcels that are significant after fdr correction, 'no' to plot unthresholded results 

neural_markers = {'dec_prec'}; 

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
%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%

Vertices_values=NaN(n_vertices, 1);
vals2plot = nan(length(allsubj),length(ROIs)); % Intitialite vals2plot

% Pull values to plot

vals2plot=nan(length(allsubj),22); %initialize values only for 1 hemi set of rois 
 if strcmp(plot_significant_only,'no')
     pvals_rois = nan(length(ROIs),1); % Initialize p-values vector if only significant parcels should be plotted 
 end 
for o = 1:length(ROIs(1,:)) 
    for s = 1:length(allsubj)
       load(['/mnt/homes/home028/gmonov/eye_tracking_data_analysis/decoding_results/eye_mean_dec_prec/',allsubj{s},'_eye_mean_dec_prec_x_y_1000.mat'])
       eye_dec(s,1) = eye_dec_prec; 
       clear eye_dec_prec
       vals_in_parcel = []; % Initialize values for a certain cluster that will be averaged 
       for tui = 1:length(glasser_22(:,1))
           if ~isempty(glasser_22{tui,o})
              val=load(['/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/dec_prec/',allsubj{s},'_',glasser_22{tui,o},'.mat']);  %Load value of neural marker that is of interest 
              eval(['vals_in_parcel(tui) = val.dec_prec(1,1);']) 
              clear val
           else break
           
           end   
       end 

       vals2plot(s,o)=mean(vals_in_parcel); % average values for the subparcels in this ROI and add to Subject X ROI matrix 
        
      
    end
end  

for g = 1:length(glasser_22)
    [r,p] = corr(vals2plot(:,g),eye_dec)
    if strcmp(plot_significant_only,'yes')
        if p <0.05
          corrs2plot(g) = r;    
        else corrs2plot(g) = nan;  
        end 
    else
    corrs2plot(g) = r; 
    end 
    corr_p(g) = p; 
    clear r p 
end 



 %Tranform vals2plot from 22 parcels into 180 ROIs 
     corrs2plot_new = zeros(1,length(ROIs180)); 
     corrs2plot_new(corrs2plot_new==0)=nan; 

     for puz = 1:length(new_ROIs180) 
         for ruz = 1:length(ROIs)
             for quz = 1:length(glasser_22(:,1))
               if strcmp(new_ROIs180{puz},glasser_22{quz,ruz})
                   corrs2plot_new(puz)= corrs2plot(ruz);
                   break
               end 
               
             end 
         end 
     end 

corrs2plot = [corrs2plot_new,corrs2plot_new]; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i=1:size(Areas_labels,1)
    ROI=cell2mat(Areas_labels(i));
    id_Area=find(strcmp(Areas_labels, ROI));
    if ~isempty(id_Area)
        Atlasid=atlas.indexmax;
        idAr=find(Atlasid==id_Area);
        Vertices_values(idAr)=corrs2plot(i);
    end
 end


 %Change clims 

 clim = [-0.5,0.5]; 
  % Get rid of negative limits for trial variability measures & specify
  % color map 

 cmap=brewermap(256, '*RdBu'); % take red-blue for all others 
 

 Vertices_values(isnan(Vertices_values))=100; % missing valus (NaN) gets very high values, to be plotted as gray 


 cmap=[cmap; colors.silver]; % defining color map; last row is set to be gray

  counts = 0; 

 % plot Left hemisphere, by taking the first half of "Vertices_values"
for x = 1:2 % Plot twice to rotate one of them later to display medial brain 
     counts = counts+1; 

     eval(['s',num2str(counts),'=subplot(1,2,counts), hold on'])
     

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
    
     if x == 2
     cb = colorbar('east'), hold on ; 
     end       
    
end

% Arrange plots 
w = 0.1638; wgap = -0.08; ngap = 0.05; cw = 0.008; cl = 0.2; woff = 0.05;
h = 0.3393; hgap = 0.08; hoff = 0.12;

set(s1, 'Position', [woff, hoff, w, h])   % [left bottom width height]
set(s2, 'Position', [woff+w-ngap,hoff, w, h])
set(cb, 'Position', [woff+w*2.95+wgap*2-ngap,hoff,cw,h])


savefig(['/mnt/homes/home028/gmonov/eye_tracking_data_analysis/eye_figures/corr_meg_eye_dec.fig'])
close all
% Script to plot cortical areas on brain maps in defined colors 
% Produces inset brain maps of Figure 5 
% Gina Monov, UKE, 2022 

clear all
close all
% loading support files and templates
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')
addpath /mnt/homes/home032/aarazi/fieldtrip-20201009
ft_defaults
atlas = ft_read_cifti('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/support_files/Glasser_atlas.dlabel.nii');
Areas = readtable('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/Glasser_labels.csv'); Areas_labels=Areas.Label;
addpath /mnt/homes/home032/aarazi/Confirmation_bias/code/ColorMap_2
load('/mnt/homes/home028/gmonov/meg_analysis/Brain_plot/m1ind.mat')

% Specify surface to be shown
surf2show = 'medial'; % lateral, medial, occipital

parcels={{'V1'},...
    {'V2', 'V3', 'V4'},...
    {'V6', 'V3A', 'V7', 'IPS1', 'V3B', 'V6A'},...
    {'MST', 'LO1', 'LO2', 'MT','V4t', 'FST', 'LO3', 'V3CD', 'PH'},...
    {'V8', 'FFC', 'PIT', 'VMV1', 'VMV3', 'VMV2', 'VVC'},...
    {'6a', '6d'},...
    {'FEF', 'PEF'},...
    {'6v', '6r'},...
    {'4'}};
names={'V1', 'V2-V4', 'Dorsal visual', 'MT+', 'Ventral visual','PMd', 'Eye fields', 'PMv', 'M1'};
m1idx = [m1ind;m1ind]; %Indices of JW's M1 hand area

% Specify colors 
tcols = [colors.blueberry;...
         colors.cobalt; ...
         colors.sky;...
         colors.teal;...
         colors.purple;...
         colors.seafoam;...
         colors.juniper;...
         colors.green;...
         colors.lime];
     
for i=1:length(names)
    clusters(i).name=names{i};
    clusters(i).roi=parcels{i};
    clusters(i).colors = tcols(i,:); 
  
end
n_vertices=64984; % n of vertices in the template 


     
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


% replace primary motor area (4) with JW's hand area 
idx = find(strcmp(ROIs,'4')); 
atlas.indexmax(atlas.indexmax==idx) = nan; 
atlas.indexmax(atlas.indexmax==idx+length(ROIs)) = nan; 
for m1h = 1:length(atlas.indexmax) 
    if m1idx(m1h) == 1
       atlas.indexmax(m1h) = idx; 
    end
    
end 

vals2plot=NaN(length(ROIs), 1);
Vertices_values=NaN(n_vertices, 1);  %initialize values only for 1 hemi set of rois 

for r = 1:length(ROIs)
    for c = 1:length(clusters) 
        for p = 1:length(clusters(c).roi)
        if strcmp(ROIs{r},clusters(c).roi{p})
            vals2plot(r) = c; 
        end 
        end 
    end 
end 
vals2plot = [vals2plot;vals2plot];

% add 
 for i=1:size(Areas_labels,1)
    ROI=cell2mat(Areas_labels(i));
    id_Area=find(strcmp(Areas_labels, ROI));
    if ~isempty(id_Area)
        Atlasid=atlas.indexmax;
        idAr=find(Atlasid==id_Area);
        Vertices_values(idAr)=vals2plot(i);
    end
 end
 
 Vertices_values(isnan(Vertices_values))=length(clusters)+1; % missing valus (NaN) gets very high values, to be plotted as white 


cmap=tcols; cmap=[cmap; colors.white];  % defining color map; last row is set to be white
clim=[1,length(clusters)+1]; 

% Initialize figure 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 6.9; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),


 x = 1; 
 eval(['s',num2str(1),'=subplot(3,1,x),hold on'])

 Plot_scalp_map('L', cmap, clim, Vertices_values(1:n_vertices/2)); 
    

     view([270 0]); 
     lighting gouraud
     material dull
     camlight infinite
     light('Style','Infinite')
     if strcmp(surf2show,'occipital')
 
         view([360 0]); 
       
     elseif strcmp(surf2show,'medial')

         view([90 0]);
     
     end
      xlim([-100 50]); zlim([-70 80]); ylim([-105 70])



%Arrange plots 
w = 0.25; wgap = 0.019; woff = 0.01;
h = 0.3; hgap = 0.01; hoff = 0.01;

set(s1, 'Position', [woff, hoff+(2*(hgap+h)), w, h])   % [left bottom width height]


savefig(['/mnt/homes/home028/gmonov/final_figures/roi_depiction_',surf2show,'.fig'])
close all

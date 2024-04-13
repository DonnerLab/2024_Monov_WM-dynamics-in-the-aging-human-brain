% Plot decoding results for each clustered brain region and perform
% cluster-based permutation test against zero: Edden M. Gerber (2023). permutest, MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/71737-permutest
% --> Figure 6A
% Gina Monov, UKE, 2023

clear all 
close all

% Load colors 
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')
% Load subject IDs
delay = '1'; %Specify delay duration 
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/meg_subj.mat']);
mci = meg_mci_final; 
hc = meg_hc_final; 
cog_def = meg_cog_def_final;
allsubj = horzcat(mci,hc,cog_def);

no_results = {}; % to check if any data are missing 

names={'V1', 'V2-V4', 'Dorsal_visual', 'MT', 'Ventral_visual','PMd', 'Eye_fields', 'PMv', 'JWG_M1'};
titlenames={'V1', 'Early visual', 'Dorsal visual', 'MT+ & neighbours', 'Ventral visual','PMd', 'Eye fields', 'PMv', 'M1'};

% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 14; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

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
    clusters_all_subj(i).name=names{i};
    clusters_all_subj(i).titlenames=titlenames{i};
    clusters_all_subj(i).roi='1'; % Do not specify rois since the decoding was performed in the entire cluster 
    clusters_all_subj(i).colors = tcols(i,:); 
  
end

for cc = 1:length(clusters_all_subj) 
    for s = 1:length(allsubj) 

         decoding_results = []; 
         decoding_results = readtable(['/mnt/homes/home028/gmonov/meg_analysis/Decoding/decode_5_35_coarse/' allsubj{s} '_' clusters_all_subj(cc).name '_' delay '_5_35.csv']); 


     if ~isempty(decoding_results)
        if sum(isnan(decoding_results.test_correlation(:,1:end)))>0 % Check whether there are unexpected nans
           no_results = horzcat(no_results,allsubj{s});
           clusters_all_subj(cc).corr(s,:) = nan; 
        else 

             times = decoding_results.latency;
             v = 1:length(times); 
             oddidx = @(v) v(1:2:end);           % Addressing Odd-Indexed Element (results accidentally save twice for each time point in the decoding)
             y1 = oddidx(v); 

        clusters_all_subj(cc).corr(s,:) = decoding_results.test_correlation(y1); 
        times = times(y1); 
        end 
     end 


     % Extract number of principal components for one example cluster (V1) 
     if cc == 1
         n_pcs(s) = mean(decoding_results.n_comps_80);  
     end 
             
    end
    
    
              
% Perform cluster-based permutation test against zero for decoding time course 

data4stats = []; 
data4stats = clusters_all_subj(cc).corr'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);
 
    % Plotting 
    eval(['s',num2str(cc),'=subplot(3,10,cc); hold on']) 
    plot(times,mean(clusters_all_subj(cc).corr,1),'Color',tcols(cc,:)), hold on 
    shadedErrorBar(times,mean(clusters_all_subj(cc).corr,1),std(clusters_all_subj(cc).corr,[],1)./sqrt(size(clusters_all_subj(cc).corr,1)),{'Color',tcols(cc,:)},1)
         
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
                plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.01,'Color',tcols(cc,:),'LineWidth',1.5), hold on 
             end
         end 
                
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.1 0.2 0.3 0.4])
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{ '0'  '0.1' '0.2'  '0.3' '0.4'})
onsets = [0 0.5]; % mem on and mem off
xline(onsets(1),'--','Color',colors.black,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.black,'LineWidth',1), hold on 
yline(0,'-','Color',colors.black), hold on % Reference line around zero
title(clusters_all_subj(cc).titlenames,'Fontsize',8,'Color', clusters_all_subj(cc).colors)

if cc == 1
   xlabel('Time from memorandum (s)','Fontsize',8)
   ylabel('Decoding precision','Fontsize',8)
end 

xlim([-0.1 1.5]); ylim([-0.05  0.4])
   
end 

% Arrange plots 
w = 0.0933; wgap = 0.019; woff = 0.05;
h = 0.14; hgap = 0.06; hoff = 0.12;

set(s1, 'Position', [woff, hoff+(1*(hgap+h)), w, h])   % [left bottom width height]
set(s2, 'Position', [woff+(w+wgap),hoff+(1*(hgap+h)), w, h])

set(s3, 'Position', [woff+(w+wgap)*2, hoff+(2*(hgap+h)), w, h])
set(s4, 'Position', [woff+(w+wgap)*2, hoff+(1*(hgap+h)), w, h])
set(s5, 'Position', [woff+(w+wgap)*2, hoff, w, h])

set(s6, 'Position', [woff+(w+wgap)*3, hoff+(2*(hgap+h)), w, h])
set(s7, 'Position', [woff+(w+wgap)*3, hoff+(1*(hgap+h)), w, h])
set(s8, 'Position', [woff+(w+wgap)*3, hoff, w, h])

set(s9, 'Position', [woff+(w+wgap)*4, hoff+(1*(hgap+h)), w, h])


savefig(['/mnt/homes/home028/gmonov/final_figures/decoding.fig']) % Figure 6A
close all
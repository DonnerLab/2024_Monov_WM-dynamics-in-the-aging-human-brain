% 3 and 9s TFRs of power modulations for all subjects included in MEG analysis 
% Dorsal visual cortex only
% --> Extended Data Figure 5-1A
% Gina Monov, UKE, 2022 

clear all

% Load colors 
load('/mnt/homes/home028/gmonov/functions/colors/colors.mat')

% Load subject IDs
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/Final/meg_subj.mat']);
mci = meg_mci_final; 
hc = meg_hc_final; 
cog_def = meg_cog_def_final;

allsubj = horzcat(mci,hc,cog_def);

parcels={{'V6', 'V3A', 'V7', 'IPS1', 'V3B', 'V6A'}};
names={'Dorsal visual'};

% Specify colors 
tcols = [colors.sky];
     
for i=1:length(names)
    clusters_all_subj(i).name=names{i};
    clusters_all_subj(i).roi=parcels{i};
    clusters_all_subj(i).colors = tcols(i,:);   
end

% Fieldtrip 
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home032/aarazi/fieldtrip-20201009'
ft_defaults

% Define path to load subject data from 
loadpath = ['/mnt/homes/home028/gmonov/meg_analysis/task_analysis/TFR/output_ave_TFR/'];

% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 22.5; % figure width
fig_h = 14; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),
   delays = [3,9];
   
for loop = 1:2 % loop over the two delay durations 
    
    delay = delays(loop); %define delay duration to plot 
    foi = [2 120];  % bounds on frequencies for plotting/cluster-corrected stats
    bl = [-0.4 -0.2]; %baseline times
    clustalpha = 0.05;  % intial threhsold for forming clusters
    trlwin = [bl(1) delay+0.5];
for cc = 1:length(clusters_all_subj) % Loop over the clusters to plot     
for s = 1:length(allsubj)
    freqs = [];
    times = [];
    all_parcels(s).tfrs = [];
    tfrs = [];
    clusters = [];
    load( ['/mnt/homes/home028/gmonov/meg_analysis/task_analysis/TFR/output_ave_TFR_all_glasser/', allsubj{s},'_', num2str(delay) '_F_AVE_TFR_all_glasser.mat'])

        freq = freqs(1,:);
        trltimes = times(1,:);  
        new_freqs = freqs(freqs>=foi(1) & freqs <=foi(2));
        parcel_tfrs = [];
            for pp = 1:length(clusters_all_subj(cc).roi) %Loop over the parcels in the cluster 
                      idx = find(strcmp(clusters.roi,clusters_all_subj(cc).roi{pp}));

                      parcel_tfrs = [parcel_tfrs;clusters.parcel_data{1,idx}];

            end 
            parcel_tfrs = mean(parcel_tfrs,1); 

            tfrs = [tfrs;parcel_tfrs];

        
        tfrs = mean(tfrs,1);
        

    tfrs_full = tfrs; 

    tfrs = tfrs(:,freq>=foi(1) & freq<=foi(2),trltimes>=trlwin(1) & trltimes<=trlwin(2));

    % Save full length data of all subjects
    all_subj_data(:,:,s) = squeeze(mean(tfrs,1)); 
    
    % Create structure for this subject to mimic output from FT_FREQANALYSIS
    cstruct = struct;
    cstruct.label = unique(clusters_all_subj(cc).name);  % arbitrary label for motor lateralization 'channel'
    cstruct.fsample = 20;  
    cstruct.freq = freq(freq>=foi(1) & freq<=foi(2));
    cstruct.time = trltimes(trltimes>=trlwin(1) & trltimes<=trlwin(2));
    cstruct.dimord = 'chan_freq_time';
    trl_avg = mean(tfrs_full,1); % Pull data for this subject only again 
    % Add structures to group-level arrays
    if ~isnan(tfrs)
    cstruct.powspctrm = squeeze(reshape(trl_avg(1,freq>=foi(1) & freq<=foi(2),trltimes>=trlwin(1) & trltimes<=trlwin(2)),[1,size(trl_avg(1,freq>=foi(1) & freq<=foi(2),trltimes>=trlwin(1) & trltimes<=trlwin(2)))]));
    cstruct.powspctrm = reshape(cstruct.powspctrm,[size(trl_avg(1,freq>=foi(1) & freq<=foi(2),trltimes>=trlwin(1) & trltimes<=trlwin(2)))]);
    else cstruct.powspctrm = zeros(1,56,39);
         cstruct.powspctrm(cstruct.powspctrm==0) = nan; 
    end 
    all_avF{s} = cstruct;
    cstruct.powspctrm = zeros(size(cstruct.powspctrm)); allBnull{s} = cstruct;
 end 

try load(['/mnt/homes/home028/gmonov/meg_analysis/task_analysis/TFR/stats4FigureS5/',num2str(delay),'_',clusters_all_subj(cc).name,'_','stats.mat'])
    
catch 
    
% Set up cfg stucture for cluster-based permutation tests
cfg_stat = [];
cfg_stat.latency     = 'all';
cfg_stat.frequency   = 'all';
cfg_stat.avgovertime = 'no';
cfg_stat.avgoverfreq = 'no';
cfg_stat.parameter   = 'powspctrm';
cfg_stat.method      = 'montecarlo';
cfg_stat.statistic   = 'depsamplesT';  
cfg_stat.alpha       = 0.05;
cfg_stat.correctm    = 'cluster';
cfg_stat.clusteralpha = clustalpha;  % intial threhsold for forming clusters
cfg_stat.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics
cfg_stat.correcttail = 'prob';
cfg_stat.numrandomization = 10000;
cfg_stat.neighbours = struct([]); % empty neighbour structure
cfg_stat.channel = 1; 
Nsub = length(allsubj);
cfg_stat.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg_stat.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg_stat.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg_stat.uvar                = 2; % the 2nd row in cfg.design contains the subject number

% Calculate cluster-corrected stats
stats = ft_freqstatistics(cfg_stat,all_avF{:},allBnull{:});

% Save cluster correctes stats for this cluster to save time 
save(['/mnt/homes/home028/gmonov/meg_analysis/task_analysis/TFR/stats4FigureS5/',num2str(delay),'_',clusters_all_subj(cc).name,'_','stats.mat'],'stats')
end

% Plotting 
if cc == 1 && loop == 1
find_clims = [abs(min(mean(all_subj_data,3))),max(mean(all_subj_data,3))];
clims = [-max(find_clims) max(find_clims)];
cmap = colormapsmoothPM(500);
end 
 caxis(clims)
eval(['s',num2str(loop),'=subplot(1,2,loop); hold on']) 
S = contourf(times(times>=trlwin(1)&times<=trlwin(2)),freqs(freqs>=2 & freqs <=120),squeeze(nanmean(all_subj_data,3)),300,'linestyle','none'); hold on
colormap(cmap), caxis(clims)
contour(stats.time(stats.time>=trlwin(1)&stats.time<=trlwin(2)),stats.freq,squeeze(stats.mask),1,'Color','k'), hold on  
if delay == 3
 set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[2 7 20 50 120],'yscale','log')
 set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2' '2.5' '3' '3.5'},'YTickLabel',{'2' '7' '20' '50' '120'})

 
elseif delay == 9 
  set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',6,'TickDir','out','box','off','XTick',[0 0.5 1.5 3.5 9.5],'YTick',[2 7 20 50 120],'yscale','log')
 set(gca,'XTickLabel',{'0' '0.5' '1.5' '3.5' '9.5'},'YTickLabel',{'2' '7' '20' '50' '120'})
  
end 
onsets = [0 0.5]; %mem on and mem off
xline(onsets(1),'--','Color',colors.black,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.black,'LineWidth',1), hold on 

title(clusters_all_subj(cc).name,'Fontsize',8,'Color', clusters_all_subj(cc).colors)

if cc == 1

if loop == 1
xlabel('Time from memorandum (s)','Fontsize',8)
ylabel('Frequency (Hz)','Fontsize',8)
end 
end 

if loop == 2
cb=colorbar('South'); ylabel(cb,'Power modulation (dB)','FontSize',8)
end 


 xlim([-0.1 delay+0.5]);

end 
clear tfrs_full tfrs cstruct all_subj_data all_avF all_Bnull parcels_tfrs stats trl_avg
end 

% Arrange plots 
w = 0.0933; wgap = 0.019; woff = 0.05;
h = 0.14; hgap = 0.06; hoff = 0.12;


set(s1, 'Position', [woff+(w+wgap)*2, hoff+(2*(hgap+h)), w, h])
set(s2, 'Position', [woff+(w+wgap)*3, hoff+(2*(hgap+h)), w*2.5, h])


set(cb, 'Position', [woff+(w+wgap)*4,hoff,w,0.025])

savefig(['/mnt/homes/home028/gmonov/final_figures/FS5_tfrs.fig']) % Extended Data Figure 5-1A
close all
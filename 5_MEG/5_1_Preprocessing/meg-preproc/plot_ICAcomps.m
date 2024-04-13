function [] = plot_ICAcomps(filein, artifact_type)
% PLOT_ICACOMPS plots ICA components
%   ICA components and the correpsonding timelocked data, frequency
%   decomposotion and spectrum. Components are ordered by the coherence value
%   The 25 components with the highest values are plotted (5 in each block)
%   and saved in the folder comp_ICA

% Current ID?
loadpath = '/mnt/homes/home028/gmonov/meg_analysis/Dataset_info/final_masterfile/';
load([loadpath,'Master_file.mat']);
    info_EL_blocks = zeros(length(file_info),7);
    files = {zeros(length(file_info),1)};
for f = 1:length(file_info)
        info_EL_blocks(f,1)=file_info(f).blocknums(1);
        info_EL_blocks(f,2)=file_info(f).blocknums(end);
        info_EL_blocks(f,3)=length(file_info(f).blocknums);
        info_EL_blocks(f,4)=file_info(f).sess;
        info_EL_blocks(f,5)=file_info(f).ELex(1);
        info_EL_blocks(f,6)=file_info(f).ELex(2);
        if length(file_info(f).blocknums)==4
            info_EL_blocks(f,7:8)=file_info(f).ELex(3:4);
        elseif length(file_info(f).blocknums)==3
            info_EL_blocks(f,7)=file_info(f).ELex(3);
            info_EL_blocks(f,8)=nan;
        else info_EL_blocks(f,7:8)=nan;
        end 
       
        files{f,1}=file_info(f).name; 
    end 
        
%idx_file to process?
idx_file = find(strcmp(files, filein));
ID = [file_info(idx_file).subj,'_',num2str(file_info(idx_file).sess)];
name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/comp_ICA/' ID '/'];
cd(name)
load([name 'comp_' ID '.mat'])

name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/ICA/comp_ICA/' ID '/Matrices/'];
cd(name)
        
load([name artifact_type '_fdcomp_' ID '.mat']);
load([name artifact_type '_timelock_' ID '.mat']);
load([name artifact_type '_freq_' ID '.mat']);

top_num = 25;

[val_cor,idx_coh] = sort(nanmean(abs(fdcomp.cohspctrm),2),'descend');
top_components = idx_coh(1:top_num);
%art_stop_components = sort(top_components);
%man_components = [artifact_type '_' num2str(top_num) 'TopCoh_' ID '.mat'];
%save(man_components,'art_stop_components')

fprintf(['\n ---------- Plot components (' artifact_type ') ----------\n'])

% Plot timelocked data with trial avergaed artifactual channel and topography

% Set some plotting configurations
cfg = [];
cfg.layout = 'CTF275.lay';
cfg.style = 'straight';
cfg.shading = 'interp';
cfg.marker = 'off';
cfg.comment = 'no';

%         correct_labels = ft_channelselection('meg', data.label);
%         comp_labels = comp.label;
%         comp.label = correct_labels(1:length(comp_labels));

pCount = [1 6 11 16 21];

comp_count = 1;

for f = 1:5
    myfig = figure; clf
    p = uipanel('Parent',myfig);
    p.Title = [artifact_type ' - Top #' num2str(1+(f-1)*5) '-' num2str(f*5) ' components'];
    p.FontWeight = 'bold'; p.FontSize = 12; p.TitlePosition = 'centertop';
     
    for i = 1:5
        fprintf(['\n... Plotting component ' num2str(comp_count) ' of 25\n'])
        subplot(5,5,pCount(i),'Parent',p); plot(timelock.time, timelock.avg(top_components(comp_count)+1,:), 'linewidth',2) % Component
        hold on
        plot(timelock.time, timelock.avg(1,:)./max(timelock.avg(1,:))*max(timelock.avg(top_components(comp_count)+1,:))) % adjusted trial-averaged artifactual channel
          
        title({['Time series top #' num2str(comp_count) ' , c=' num2str(round(val_cor(comp_count),3))];''})
        ax = gca; ax.FontSize = 8; xlim([min(timelock.time) max(timelock.time)]);
        if i>4, xlabel('Time in seconds'); end
        %legend('component', 'artifact','Location','southwest','FontSize',8); legend('boxoff')
        subplot(5,5,pCount(i)+1,'Parent',p); plot(fdcomp.freq, abs(fdcomp.cohspctrm(top_components(comp_count),:))); % coherence spectrum
        title({'Coherence spectrum';''})
        ax = gca; ax.FontSize = 8;
        if i>4, xlabel('Frequency in Hz'); end
        subplot(5,5,pCount(i)+2,'Parent',p); plot(freq.freq, squeeze(nanmean(abs(freq.fourierspctrm(:,top_components(comp_count)+1,:))))); % frequency spectrum
        title({'Frequency spectrum';''})
        ax = gca; ax.FontSize = 8; set(gca,'xscale','log');
        if i>4, xlabel('Frequency in Hz'); end
        subplot(5,5,pCount(i)+3,'Parent',p); scatter(1:size(timelock.trialvar,1), timelock.trialvar(:,top_components(comp_count)+1),'.k') % Variance per trial
        title({'Variance';''})
        if i>4, xlabel('Trials'); end
        ax = gca; ax.FontSize = 8;
        axis tight
        subplot(5,5,pCount(i)+4,'Parent',p); cfg.component = top_components(comp_count); ft_topoplotIC(cfg,comp) % topography
        ax = gca; ax.FontSize = 8;
        comp_count = comp_count + 1;
    end

        fig_name = ['/mnt/homes/home028/gmonov/meg_analysis/meg_preprocessing/identify_ICA_comps/' ID '/'];
      
        if 7==exist(fig_name,'dir')
            cd(fig_name)
        else
            mkdir(fig_name)
            cd(fig_name)
        end
    
    coh_top = [artifact_type '_' num2str(f) 'TopCoh_' ID '.fig'];
    saveas(gca, coh_top, 'fig')
    close all
end
end

